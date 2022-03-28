#include "RandomFieldGenerator.h"
#include "Region.h"
#include "Exception.h"
#include "ErrorHandler.h"

#include <chrono>

using namespace std;

namespace csmp {

template<uint32_t dim>
RandomFieldGenerator<dim>::RandomFieldGenerator()
 : m1_(0,0), m2_(0,0), m3_(0,0), pi_(csmp::PI), k_(0), fname("random_permeability_field.txt")
{

} 

template<uint32_t dim>
Matrix& RandomFieldGenerator<dim>::UniformRandomMatrix(size_t m, size_t n)
{
   // ascertain that we get a different distribution everytime we call this
    std::random_device rd;
    std::mt19937::result_type seed = rd() ^ (
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::system_clock::now().time_since_epoch()
                ).count() +
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now().time_since_epoch()
                ).count() );

  std::mt19937 gen(seed);
  std::uniform_real_distribution<> rndist(0,1);
  
  m1_.Resize(m,n);
  for(size_t i{0U}; i<m; i++) {
      for(size_t j{0U}; j<n; j++)
          m1_(i,j) = rndist(gen);
    } 
  return m1_;
}


template<uint32_t dim>
Matrix& RandomFieldGenerator<dim>::NormalRandomMatrix(size_t m, size_t n)
{    
    std::random_device rd;
    std::mt19937::result_type seed = rd() ^ (
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::system_clock::now().time_since_epoch()
                ).count() +
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now().time_since_epoch()
                ).count() );

  std::mt19937 gen(seed);
  std::normal_distribution<> rndist(0,1);
  
  m3_.Resize(m,n);
  for (size_t i{0U}; i<m; i++) {
      for (size_t j{0U}; j<n; j++) {
          m3_(i,j) = rndist(gen);
      }
  }
  return m3_;
}

template<uint32_t dim>
double RandomFieldGenerator<dim>::TheoreticalStandardDeviation(size_t m, double xl, double yl, double Lx, double Ly)
{    
  const size_t cnt(m+2);
  const double t1(pi_*xl/Lx), t2(pi_*yl/Ly), pi2(2.0*pi_);
  Matrix M1(cnt, cnt), M2(cnt, cnt);
  
  for ( size_t i{0U}; i<cnt; i++ )
      for (size_t j{0U}; j<cnt; j++ ) {
          M1(i,j) = static_cast<double>(j)*t1;
          M2(i,j) = static_cast<double>(i)*t2;
        }

  M1 *= M1;
  M2 *= M2;

  double sigma(0.0);
  for ( size_t i{0U}; i<cnt; i++ )
      for (size_t j{0U}; j<cnt; j++ )
        sigma += std::exp(-(M1(i,j)+M2(i,j))/pi2);
  
  return std::sqrt(sigma);

}


/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
void RandomFieldGenerator<dim>::RandomElementField2D( Model<dim>& mdl, const char* variable, double mean, double sigma, double xlength, double ylength, 
                                                      bool logarithmic, const char* region, size_t iterations )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Generates a Gaussian random field for a scalar variable placed on the elements
of the entire Model using the Karhuem Love decomposition. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->
A reference to the Model, the name of the variable for which the random
field shall be generated, its mean value, the standard variation, the
correlation length in x and y direction. The flag bool determines if the input
mean is logarithmic (e.g., for the permeability) in which case the exponent
of the computed random field is stored back to the nodes. The second to last 
argument is the name of the Model subregion ("Model") by default. The last argument 
sets the number of iterations to compute the initial Gaussian random field.
The default is 50 and cannot be set to a smaller number as, otherwise, the
iteration will become inaccurate.

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->


<H4>Implementation:</H4><!------------------------------------------------>


<H4>Application:</H4><!--------------------------------------------------->


<H4>Messages:</H4><!------------------------------------------------------>

<!------------------------------------------------------------------------>
tested: O.K. */
template<uint32_t dim>
void RandomFieldGenerator<dim>::RandomElementField2D( Model<dim>& mdl, const char* variable, double mean, double sigma, double xlength, double ylength, 
													  bool logarithmic, const char* region, size_t iterations )
{
  
  Index         key(mdl.Database().StorageKey(variable));
  
  if ( key.place != ELEMENT or key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomElementField2D",
                     variable,    "Property must be a scalar variable placed on the elements");
  if ( dim != 2 )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomElementField2D",
                     "Methods works only in two dimensions");
  if ( xlength <= 0.0 or ylength <= 0.0 )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomElementField2D",
                     "Correlation must be larger than zero");
  if ( sigma <= 0.0 )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomElementField2D",
                     "Standard deviation must be larger than zero");
  if ( iterations < 50 ) {
     throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::RandomElementField2D",
                     "Number of iterations less than 50, Monte Carlo iteration possibly incorrect, resetting iterations to 50");
     iterations = 50;
   }

  // get normal random matrix
  const size_t m(iterations);
  Matrix R(m+1,m+1);
  R = NormalRandomMatrix(m+1,m+1);
  
  // define compute model dimensions and some other constant variables
  Point<dim> xy_min, xy_max;
  mdl.MinMaxCoordinates( xy_min, xy_max );
  const double Lx(xy_max[0]-xy_min[0]);
  const double Ly(xy_max[1]-xy_min[1]);
  const double sqrt2(std::sqrt(2.0));
  const double sqrtLxLy(std::sqrt(Lx*Ly));
  const double sqrtLxLy2(2.0/sqrtLxLy);
  const double term(sqrt2/sqrtLxLy);
  const double pi2(pi_*2.0);
  const double sigma2(TheoreticalStandardDeviation(iterations, xlength, xlength, Lx, Ly));
  const double min(1.0e-25);
  
  
  double   sc, sc2;
  Point<dim> bc;
  double   v1, v2, v3, t1, t2, t3, t4;
  
  // resize the storage vector for the random permeability field
  Region<dim>& mref = mdl.Region(region);
  k_.resize(mdl.Region("Model").Elements());

  cout << "\nRandomFieldGenerator<dim>::RandomElementField2D: Generating random field for '" << variable << "' in region " << region;
  cout << "\nUsing mean: " << mean << ", standard deviation: " << sigma << ", correlation length x: " << xlength << ", correlation length y: " << ylength << endl;
    
  // loop over elements
  for ( auto it = mref.ElementsBegin(); it != mref.ElementsEnd(); it++ ) {

      sc = (1.0/sqrtLxLy*R(0,0));
      bc = (*it)->BaryCenter(); // xy coordinates of bary centre
      v1 = v2 = v3 = 1.0;
      for (size_t i{0U};i<m; i++) {
          t1 = std::cos(pi_ * v1 * bc[0]/Lx);
          t2 = std::exp(-std::pow((pi_ * v1 * xlength/Lx),2.0)/pi2);
          t3 = std::cos(pi_ * v1 * bc[1]/Ly);
          t4 = std::exp(-std::pow((pi_ * v1 * ylength/Ly),2.0)/pi2);
          sc += term * t1 * t2 * R(0,i+1);
          sc += term * t3 * t4 * R(i+1,0);
          v1 += 1.0;
        }
      for (size_t j{0U}; j<m; j++ ) {
          v3 = 1.0;
          for (auto k=0; k<m; k++ ) {
              t1 = std::cos(pi_ * v2 * bc[0]/Lx);
              t2 = std::cos(pi_ * v3 * bc[1]/Ly);
              t3 = std::pow((pi_ * v2 * xlength/Lx),2.0) + std::pow((pi_ * v3 * ylength/Ly),2.0);
              t4 = std::exp(-t3/pi2);
              sc += sqrtLxLy2 * t1 * t2 * t4 * R(j+1,k+1); 
              v3 += 1.0;
            }
          v2 += 1.0;
        } 
      // scale with standard deviation  
      sc *= sigma/sigma2;
      sc += mean;
      if ( std::fabs(sc) < min ) sc = min;
      
      if ( logarithmic ) sc2 = std::pow(10.0,sc);
      else               sc2 = sc;
      (*it)->Store( key, makeScalar(PLAIN,sc2) );
      k_[(*it)->Idx()] = sc2;
                                                   
    }                                                    


}


template<uint32_t dim>
/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
void RandomFieldGenerator<dim>::RandomNodeField2D( Model<dim>& mdl, const char* variable, double mean, double sigma, double xlength, double ylength, 
                                                  bool logarithmic, const char* region, size_t iterations )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Generates a Gaussian random field for a scalar variable placed on the nodes
of the entire Model using the Karhuem Love decomposition. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->
A reference to the Model, the name of the variable for which the random
field shall be generated, its mean value, the standard variation, the
correlation length in x and y direction. The flag bool determines if the input
mean is logarithmic (e.g., for the permeability) in which case the exponent
of the computed random field is stored back to the nodes. The second to last 
argument is the name of the Model subregion ("Model") by default. The last argument 
sets the number of iterations to compute the initial Gaussian random field.
The default is 50 and cannot be set to a smaller number as, otherwise, the
iteration will become inaccurate.

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->


<H4>Implementation:</H4><!------------------------------------------------>


<H4>Application:</H4><!--------------------------------------------------->


<H4>Messages:</H4><!------------------------------------------------------>

<!------------------------------------------------------------------------>
tested: O.K. */
void RandomFieldGenerator<dim>::RandomNodeField2D( Model<dim>& mdl, const char* variable, double mean, double sigma, double xlength, double ylength, 
												  bool logarithmic, const char* region, size_t iterations )
{
  

  Index         key(mdl.Database().StorageKey(variable));
  
  if ( key.place != NODE or key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomNodeField2D",
                     variable,    "Property must be a scalar variable placed on the nodes");
  if ( dim != 2 )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomNodeField2D",
                     "Methods works only in two dimensions");
  if ( xlength <= 0.0 or ylength <= 0.0 )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomNodeField2D",
                     "Correlation must be larger than zero");
  if ( sigma <= 0.0 )
     throw csmp::Exception( FATAL_ERROR, "RandomFieldGenerator<double,2>::RandomNodeField2D",
                     "Standard deviation must be larger than zero");
  if ( iterations < 50 ) {
     throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::RandomNodeField2D",
                     "Number of iterations less than 50, Monte Carlo iteration possibly incorrect, resetting iterations to 50");
     iterations = 50;
   }

  // get normal random matrix
  const size_t m(iterations);
  Matrix R(m+1,m+1);
  R = NormalRandomMatrix(m+1,m+1);
  
  // define compute model dimensions and some other constant variables
  Point<dim> xy_min, xy_max;
  mdl.MinMaxCoordinates( xy_min, xy_max );
  const double Lx(xy_max[0]-xy_min[0]);
  const double Ly(xy_max[1]-xy_min[1]);
  const double sqrt2(std::sqrt(2.0));
  const double sqrtLxLy(std::sqrt(Lx*Ly));
  const double sqrtLxLy2(2.0/sqrtLxLy);
  const double term(sqrt2/sqrtLxLy);
  const double pi2(pi_*2.0);
  const double sigma2(TheoreticalStandardDeviation(iterations, xlength, xlength, Lx, Ly));
  const double min(1.0e-25);
  
  
  double   sc, sc2;
  Point<dim> bc;
  double   v1, v2, v3, t1, t2, t3, t4;

  // resize the storage vector for the random permeability field
  Region<dim>& mref = mdl.Region(region);
  k_.resize( mdl.Region("Model").Nodes());

  cout << "\nRandomFieldGenerator<dim>::RandomNodeField2D: Generating random field for '" << variable << "' in region " << region;
  cout << "\nUsing mean: " << mean << ", standard deviation: " << sigma << ", correlation length x: " << xlength << ", correlation length y: " << ylength << endl;
    
  // loop over nodes
  for ( auto it = mref.NodesBegin(); it != mref.NodesEnd(); it++ ) {

      sc = (1.0/sqrtLxLy*R(0,0)); 
      bc[0] = (*it)->x();
      bc[1] = (*it)->y();
      v1 = v2 = v3 = 1.0;
      for (size_t i{0U};i<m; i++) {
          t1 = std::cos(pi_ * v1 * bc[0]/Lx);
          t2 = std::exp(-std::pow((pi_ * v1 * xlength/Lx),2.0)/pi2);
          t3 = std::cos(pi_ * v1 * bc[1]/Ly);
          t4 = std::exp(-std::pow((pi_ * v1 * ylength/Ly),2.0)/pi2);
          sc += term * t1 * t2 * R(0,i+1);
          sc += term * t3 * t4 * R(i+1,0);
          v1 += 1.0;
        }
      for (size_t j{0U}; j<m; j++ ) {
          v3 = 1.0;
          for (auto k=0; k<m; k++ ) {
              t1 = std::cos(pi_ * v2 * bc[0]/Lx);
              t2 = std::cos(pi_ * v3 * bc[1]/Ly);
              t3 = std::pow((pi_ * v2 * xlength/Lx),2.0) + std::pow((pi_ * v3 * ylength/Ly),2.0);
              t4 = std::exp(-t3/pi2);
              sc += sqrtLxLy2 * t1 * t2 * t4 * R(j+1,k+1); 
              v3 += 1.0;
            }
          v2 += 1.0;
        } 
      // scale with standard deviation  
      sc *= sigma/sigma2;
      sc += mean;
      if ( std::fabs(sc) < min ) sc = min;
      
      if ( logarithmic ) sc2 = std::pow(10.0,sc);
      else               sc2 = sc;
      (*it)->Store( key, makeScalar(PLAIN,sc2) );
      
      k_[(*it)->Idx()] = sc2;
    }                                                    

}


template<uint32_t dim>
void RandomFieldGenerator<dim>::OutputRandomElementField( Model<dim>& mdl )
{

  if ( k_.size() == 0 ) {
      throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::OutputRandomElementField2D",
                      "Random field was not generated, no field is saved");
      return;
    }
 
  const Region<dim>& mref = mdl.Region("Model");
  if ( k_.size() != mref.Elements() ) {
      throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::OutputRandomElementField2D",
                      "Size of the random field does not correspond to number of finite elements in Model");
      return;
    }
    
  std::ofstream ofs;
  ofs.open( fname.c_str(), ios::out|ios::trunc );  
  
  for ( auto it = mref.ElementsBegin(); it != mref.ElementsEnd(); it++ )
    ofs << k_[(*it)->Idx()] << endl;

  ofs.close();
  
  cout << "\nRandomFieldGenerator::OutputRandomElementField2D: Random element field saved to file '" << fname << "'" << endl;
  
}

template<uint32_t dim>
void RandomFieldGenerator<dim>::OutputRandomNodeField( Model<dim>& mdl )
{

  if ( k_.size() == 0 ) {
      throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::OutputRandomNodeField2D",
                      "Random field was not generated, no field is saved");
      return;
    }
  
  const Region<dim>& mref = mdl.Region("Model");
  if ( k_.size() != mref.Nodes() ) {
      throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::OutputRandomNodeField2D",
                      "Size of the random field does not correspond to number of nodes in Model");
      return;
    }
   
  std::ofstream ofs;
  ofs.open( fname.c_str(), ios::out|ios::trunc );  

  for ( auto it = mref.NodesBegin(); it != mref.NodesEnd(); it++ )
    ofs << k_[(*it)->Idx()] << endl;

  ofs.close();
  
  cout << "\nRandomFieldGenerator::OutputRandomNodeField2D: Random element field saved to file '" << fname << "'" << endl;
  
}



template<uint32_t dim>
void RandomFieldGenerator<dim>::InputRandomElementField( Model<dim>& mdl, const char* variable )
{
  
  std::ifstream ifs;
  std::string text_line;
  std::vector<double> k_temp;
  ifs.open( fname.c_str() );
  
  if ( !ifs.is_open() ) { 
       throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::InputRandomElementField2D",
                      "Input file for random element field could not be located");
      return;
    }
 
  while ( getline( ifs, text_line ) ) k_temp.push_back(std::atof(text_line.c_str()));
   
  ifs.close();

  Region<dim>& mref = mdl.Region("Model");
  if ( k_temp.size() != mref.Elements() ) {
      throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::InputRandomElementField2D",
                      "Size of the random field does not correspond to number of finite elements in Model");
      return;
    }
    
  const Index    key(mdl.Database().StorageKey(variable));
  ScalarVariable sc;

  for ( auto it = mref.ElementsBegin(); it != mref.ElementsEnd(); it++ ) {
      sc() = k_temp[(*it)->Idx()];
      (*it)->Store( key, sc );
    }
                                                   

  cout << "\nRandomFieldGenerator::InputRandomElementField2D: Random element field successfull input from file '" << fname << "'" << endl;
  
}

template<uint32_t dim>
void RandomFieldGenerator<dim>::InputRandomNodeField( Model<dim>& mdl, const char* variable )
{

  std::ifstream ifs;
  std::string text_line;
  std::vector<double> k_temp;
  ifs.open( fname.c_str() );
  
  if ( !ifs.is_open() ) { 
       throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::InputRandomNodeField2D",
                      "Input file for random element field could not be located");
      return;
    }
 
  while ( getline( ifs, text_line ) ) k_temp.push_back(std::atof(text_line.c_str()));
   
  ifs.close();

  Region<dim>& mref = mdl.Region("Model");
  if ( k_temp.size() != mref.Nodes() ) {
      throw csmp::Exception( ERROR, "RandomFieldGenerator<double,2>::InputRandomNodeField2D",
                      "Size of the random field does not correspond to number of nodes in Model");
      return;
    }
    
  const Index    key(mdl.Database().StorageKey(variable));
  ScalarVariable sc;
  
  for ( auto it = mref.NodesBegin(); it != mref.NodesEnd(); it++ ) {
      sc() = k_temp[(*it)->Idx()];
      (*it)->Store( key, sc );
    }


  cout << "\nRandomFieldGenerator::InputRandomNodeField2D: Random element field successfull input from file '" << fname << "'" << endl;
  
}


template class RandomFieldGenerator<2U>;

} // csp
