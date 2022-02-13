#include "CSMP_mathUtilities.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

#include <chrono>

using namespace std;

namespace csmp {

/** rsqrt -- computes reciprocal square root */
double
rsqrt( double val )
{
    return 1.0 / sqrt(val);
}


/** _v_norm1 -- computes (scaled) 1-norms of vectors */
double vector_norm1( vector<double>& x, vector<double>& scale )
{
	size_t	 i, dim = x.size();
	double s, sum(0.0);

    if ( scale.empty() ) {
         scale.resize(dim);
         for ( i=0; i<dim; i++ ) sum += fabs(x[i]);
      }
      
    for ( i=0; i<dim; i++ ) {	
         s = scale[i];
		 sum += ( s == 0.0 ) ? fabs(x[i]) : fabs(x[i]/s);
	  }

	return sum;
}


double vector_norm2( vector<double>& x, vector<double>& scale )
{
	size_t	 i, dim = x.size();
	double s, sum(0.0);

    if ( scale.empty() ) {
         scale.resize(dim);
         for ( i=0; i<dim; i++ ) sum += x[i] * x[i];
      }

    for ( i=0; i<dim; i++ ) {	
         s = scale[i];
		 sum += ( s == 0.0 ) ? x[i] * x[i] : (x[i]/s) * (x[i]/s);
	  }
	  
	return sqrt(sum);
}





/**
 
Converts (global) interpolation function derivative matrix into form where it
can be multiplied, for instance, with a 2DOF material property matrix or
for computation of 2D stiffness matrix.  

Original B: 
 
    | dN0dx dN1dx ... dN(nodes)dx |
B = | dN0dy dN1dy ... dN(nodes)dy |


Transformed B (see Zienkiewicz Vol.1 p. 48, eq. 3.10b): 
 
    | dN0dx   0   dN1dx  ... dN(nodes)dx   0   |
B = |   0   dN0dy   0    ...   0   dN(nodes)dy |
    | dN0dy dN0dx dN1dy  ... dN5dy dN(nodes)dx |
  
 
@section arguments Input Arguments 

The shape function derivative matrix DN for 1 degree of freedom per node. 

@return The expanded matrix is returned into the input matrix.

@section application Application

To pre-process DN matrices in multi-DOF computations.  

    @test OK SKM refactored and retested 2/12/2014
*/
void dN_To2DOF( uint32_t nodes, DenseMatrix<DM_MIN>& B )
 {
    assert( B.Cols() == nodes );
    const auto dof(B.Cols() * 2U);

    // farming old matrix out into entries of new one
    // starting with last row
    B.Resize(3,B.Cols() * 2U);
    for ( auto i=0; i<nodes; i++ )
      {
         B(2,i*2)   = B(1,i);
         B(2,i*2+1) = B(0,i);
      }
    
    // spreading out first row
    for ( int i=nodes; i>0; i-- ) B(0,(i-1)*2) = B(0,(i-1));
    // zeroing intermediate positions
    for ( auto i=1; i<dof; i+=2 ) B(0,i) = 0.;
      
    // spreading out second row
    for ( int i=nodes; i>0; i-- ) B(1,(i-1)*2+1) = B(1,(i-1));
    // zeroing intermediate positions
    for ( auto i=0; i<dof; i+=2 ) B(1,i) = 0.;
      
 } // end DNto2DOF




/**
    Transformation of interpolation function derivative matrix to 
    3 vectorial degrees of freedom per node as shown in Zienkewicz,
    vol. 1, p.133.
    
    @test OK SKM retested 2/12/2014
*/
void  dN_To3DOF( uint32_t nodes, DenseMatrix<DM_MIN>& B )
 {
    assert( B.Cols() == nodes );
    const auto dof = B.Cols() * 3U;
    const auto old_nodes = B.Cols();

    B.Resize(6,dof);
    
    // inserting values into the extra three bottom rows
    for ( auto i=0; i<old_nodes; i++ )
      {
         // fourth row
         B(3,i*3)   = B(1,i); // d/dy
         B(3,i*3+1) = B(0,i); // d/dx
         B(3,i*3+2) = 0.0;    // 0
         
         // fifth row
         B(4,i*3)   = 0.0;    // 0
         B(4,i*3+1) = B(2,i); // d/dz
         B(4,i*3+2) = B(1,i); // d/dy

         // sixth row
         B(5,i*3)   = B(2,i); // d/dz
         B(5,i*3+1) = 0.0;    // 0
         B(5,i*3+2) = B(0,i); // d/dx
      }
      
    // spreading out the values in the first three rows
    // (going backward in order not to overwrite values in old
    //  position, while these are still needed)
    for ( auto i=(old_nodes-1U); (i+1U)>0U; i-- )
      {
         // first row
         B(0,i*3)   = B(0,i);
         B(0,i*3+1) = 0.0;
         B(0,i*3+2) = 0.0; // zeros in the first row
         
         // second row
         B(1,i*3+1) = B(1,i);
         B(1,i*3)   = 0.0;
         B(1,i*3+2) = 0.0; // zeros in the first row
         
         // third row
         B(2,i*3+2) = B(2,i);
         B(2,i*3)   = 0.0;
         B(2,i*3+1) = 0.0; // zeros in the first row
      }
    
 } // end dN_To3DOF




/** _v_norm_inf -- computes (scaled) infinity-norm (supremum norm) of vectors */
double vector_norm_inf(  vector<double>& x, vector<double>& scale )
{
	size_t	 i, dim = x.size();
	double sum(0.0);

    if ( scale.empty() ) {
         scale.resize(dim);
         for ( i=0; i<dim; i++ ) sum += fabs(x[i]);
      }

cout <<"\nvector_norm_inf (in CSMP_math): function not implemented yet."<< endl;
return 0.;
/*
	double s, tmp, maxval = 0.0;
	if ( scale == (VEC *)NULL )
		for ( i = 0; i < dim; i++ )
		{	tmp = fabs(x->ve[i]);
			maxval = ms_max(maxval,tmp);
		}
	else if ( scale->dim < dim )
		error(E_SIZES,"_v_norm_inf");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			tmp = ( s== 0.0 ) ? fabs(x->ve[i]) : fabs(x->ve[i]/s);
			maxval = ms_max(maxval,tmp);
		}

	return maxval;
*/
}


void  vector_randomize( vector<double>& x, double scale_fac )
 {
   // get a different seed every time 
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
    std::uniform_real_distribution<> rndist(0,scale_fac);
 
    for ( vector<double>::iterator it=x.begin();
          it!=x.end(); it++ )
      *it = rndist(gen);

 } // end vector_randomize




template<typename Var>
void out( const Var& obj )
 {
   cout <<"\ncontainer: "<< obj.size() << endl;
   
   int row_break(1), split_after(10); 
   for ( auto i=0; i<obj.size(); i++, row_break++ )
     {
        if ( obj[i] > 0 ) cout <<" ";
        cout << obj[i] <<" ";
        if ( row_break == split_after )
          {
             cout << endl;
             row_break = 0;
          }
      }
   cout << endl;
 }




/// approximation of the complementary error function using a Chebyshev polynomial
double  erfc_Chebyshev( double x )
{
	double z=fabs(x);
	double t=1.0/(1.0+0.5*z);
	double ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	
	return std::max( x >= 0.0 ? ans : 2.0-ans, 0. );
}


double  erf_Chebyshev( double x )
{
	double z=fabs(x);
	double t=1.0/(1.0+0.5*z);
	double ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	
	// ascertain that the result is positive
	return std::max( 1. - (x >= 0. ? ans : 2. - ans), 0. );
}




// template explicit instantiations

template void out( const vector<ONE_BYTE_NUMBER>& );
template void out( const vector<uint32_t>& );
template void out( const vector<int32_t>& );
template void out( const vector<double>& );



/// calculates arithmetic average of variables
template<typename Var> void average( const vector<Var>& var_vec, Var& var )
 {
    typename vector<Var>::const_iterator  it(var_vec.begin());
    var = *it++;
    while ( it != var_vec.end() ) var += *it++;
    var /= static_cast<double>(var_vec.size());
 }
   
template void average( const vector<ScalarVariable >&, ScalarVariable& );

template void average( const vector<VectorVariable<1U> >&, VectorVariable<1U>& );
template void average( const vector<VectorVariable<2U> >&, VectorVariable<2U>& );
template void average( const vector<VectorVariable<3U> >&, VectorVariable<3U>& );

template void average( const vector<TensorVariable<1U> >&, TensorVariable<1U>& );
template void average( const vector<TensorVariable<2U> >&, TensorVariable<2U>& );
template void average( const vector<TensorVariable<3U> >&, TensorVariable<3U>& );





/// spline interpolation of values
double splineValue( double x, double x1, double x2, double y1, double y2, double k1, double k2)
{
    const double a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double t = ( x - x1) / ( x2 - x1 );

    return (1. - t)*y1 + t*y2 + t*(1.-t)*( a*(1.-t) + b*t);
}



double splineDerivative( double x, double x1, double x2, double y1, double y2, double k1, double k2)
{
    const double a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double t = ( x - x1) / ( x2 - x1 );

    return (y2-y1)/( x2-x1 ) + (1.-2.*t)*( a*(1.-t)+b*t)/(x2-x1) + t*(1.-t)*(b-a)/(x2-x1);
}



double splineSecondDerivative( double x, double x1, double x2, double y1, double y2, double k1, double k2)
{
    const double a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double t = ( x - x1) / ( x2 - x1 );

    return 2.*( b-2.*a +(a-b)*3.*t)/(x2-x1)/(x2-x1);
}




/** 
    Finding the root of the function by the Secant method

*/
double secant_method( double xmin, double xmax, double (*function)( double ), double tolerance)
{
  double xm{std::numeric_limits<double>::quiet_NaN()}, x0{std::numeric_limits<double>::quiet_NaN()};
  double c;
  
  if (function(xmin) * function(xmax) < 0) {  // check the range for finding root
    do {
      
      x0 = (xmin * function(xmax) - xmax * function(xmin)) / (function(xmax) - function(xmin));
      c = function(xmin) * function(x0);

      xmin = xmax;
      xmax = x0;

      if (c == 0) break;
      
      xm = (xmin * function(xmax) - xmax * function(xmin)) / (function(xmax) - function(xmin));
      
    } while (fabs(xm - x0) >= tolerance); // repeat the loop
    
  } else
  {
    cout << "Error at CSMP_mathUtilities:the_secant_method Can not find a root in the given interval";
    xm =  std::numeric_limits<double>::quiet_NaN() ;
  }
  
  return xm;
}




/**
 
 finding maximum of a function using golden-section search, see:
 
 https://en.wikipedia.org/wiki/Golden-section_search
 
 */
double maximum_of_function( double xmin, double xmax, double (*function)( double ), double tolerance )
{
  
  double gr((sqrt(5.)+1.0)/2.0);  // golden ratio
  
  double xup(xmax - (xmax - xmin) / gr);
  double xlow(xmin + (xmax - xmin) / gr);
  
  do {
    if (function(xup) > function(xlow))   // check the function is maximized
    {
      xmax = xlow;
    } else
    {
      xmin = xup;
    }
    
    xup  = xmax - (xmax - xmin) / gr;
    xlow = xmin + (xmax - xmin) / gr;
    
  } while(fabs(xup - xlow) >= tolerance) ;
  
  return (xup + xlow)/2. ;
  
}
 
  
  
  
  
/**
 
 finding minimum of a function using golden-section search, see:
 
 https://en.wikipedia.org/wiki/Golden-section_search
 
 */
double minimum_of_function( double xmin, double xmax, double (*function)( double ), double tolerance )
{
  
  double gr((sqrt(5.)+1.0)/2.0); // golden ratio
  
  double xup(xmax - (xmax - xmin) / gr);
  double xlow(xmin + (xmax - xmin) / gr);
  
  do {
    if (function(xup) < function(xlow))  // check the function is minimized
    {
      xmax = xlow;
    } else
    {
      xmin = xup;
    }
    
    xup  = xmax - (xmax - xmin) / gr;
    xlow = xmin + (xmax - xmin) / gr;
    
  } while(fabs(xup - xlow) >= tolerance) ;
  
  return (xup + xlow)/2. ;
  
}

  
  
  
/**
 
 finding the x value, in the range [xmin,xmax],  where function the df(x)/dx-(f(x)-f(x1))/(x-x1).. this line is the secant and tangent to f(x) at the shock point
 Mahyar: this is a basic function for two phase shock velocity analysis.

*/
double g( double x, double x1, double (*function)( double ), double (*dfunction)( double ))   // "g(x) = df(x)/dx-(f(x)-f(x1))/(x-x1)"
  {
    return dfunction(x)-(function(x)-function(x1))/(x-x1) ;
  }
  
  
  
double secant_line( double x1, double xmin, double xmax, double (*function)( double ), double (*dfunction)( double ), double tolerance)
  {
    
    double xm, x0;
    double c;
    
     if (g(xmin,x1,function ,dfunction) * g(xmax,x1,function ,dfunction) < 0) {  // check the range for finding root
       do {
        
        x0 = (xmin * g(xmax,x1,function ,dfunction) - xmax * g(xmin,x1,function ,dfunction)) / (g(xmax,x1,function ,dfunction) - g(xmin,x1,function ,dfunction));
        c  = g(xmin,x1,function ,dfunction) * g(x0,x1,function ,dfunction);
        
        xmin = xmax;
        xmax = x0;
        
        if (c == 0.0) break;
        
        xm = (xmin * g(xmax,x1,function ,dfunction) - xmax * g(xmin,x1,function ,dfunction)) / (g(xmax,x1,function ,dfunction) - g(xmin,x1,function ,dfunction));
        
       } while (fabs(xm - x0) >= tolerance); // repeat the loop
      
     }
    else
    {
      cout << "Error at CSMP_mathUtilities:the_secant_method Can not find a root in the given interval";
      xm =  std::numeric_limits<double>::quiet_NaN() ;
    }
    
    return xm;
  }





} // csmp









