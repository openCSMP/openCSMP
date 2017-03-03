#include <sstream>
#include "TensorVariable.h"
#include "StressRotate.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp
{

/**
   constructs class from the orientations (vectors) of sigma1 and 3, 
   and the magnitudes of the principal stresses 1,2 and 3
   
   @attention the orientations of sigma1 and sigma must be perpendicular to eachother.
*/
StressRotate::StressRotate( const Point<3U>& n1, 
                            const Point<3U>& n3, 
                            double64 sigma1, double64 sigma2, double64 sigma3 )
  : original_n1_(n1), 
    original_n3_(n3),
    n1_(n1), n3_(n3), 
    sigma1_(sigma1), sigma2_(sigma2), sigma3_(sigma3), 
    x_(0.), y_(0.), z_(0.)
  {
  }


/// copy constructor
StressRotate::StressRotate( const StressRotate& sr )
  : original_n1_(sr.original_n1_), original_n3_(sr.original_n3_),
    n1_(sr.n1_), n3_(sr.n3_),
    sigma1_(sr.sigma1_), sigma2_(sr.sigma2_), sigma3_(sr.sigma3_),
    x_(sr.x_), y_(sr.y_), z_(sr.z_)
 {
 }



StressRotate::~StressRotate() {}

  
/// converts to state initialized by default constructor
void StressRotate::Reset() 
  {
    n1_ = original_n1_;
    n3_ = original_n3_;
    x_ = 0.;
    y_ = 0.;
    z_ = 0.;
  }



// inline functions

// maximum (compressive) stress (Pa/m2), corresponding to max Eigenvalue of stress tensor
double64 StressRotate::Sigma1() const
  { return sigma1_; }
  
// intermediate principal stress (Pa/m2)
double64 StressRotate::Sigma2() const
  { return sigma2_; }
  
// minimum compressive stress (Pa/m2) = mimimum Eigenvalue of stress tensor
double64 StressRotate::Sigma3() const
  { return sigma3_; }
  





/**
    Finds the the Cartesion stresses from the principal stresses
    using the formulae (6.86, 6.87, and 6.88) from Pollard & Fletcher's
   "Fundamentals of Structural Geology", see p. 225
*/
void StressRotate::CartesianStressTensor( TensorVariable<3U>& stress ) const 
  {
    Point<3U> n2 (crossProduct(n3_,n1_));
    n2.NormalizeLengthTo(1.);
    
    const double64 mxx ( n1_[0] ), mxy ( n1_[1] ), mxz ( n1_[2] ) ; 
    const double64 myx ( n2[0] ), myy ( n2[1] ), myz ( n2[2] ) ; 
    const double64 mzx ( n3_[0] ), mzy ( n3_[1] ), mzz ( n3_[2] ) ; 
    
    const double64 sxx (sigma1_*mxx*mxx + sigma2_*myx*myx + sigma3_*mzx*mzx);    
    const double64 sxy (sigma1_*mxx*mxy + sigma2_*myx*myy + sigma3_*mzx*mzy);   
    const double64 syy (sigma1_*mxy*mxy + sigma2_*myy*myy + sigma3_*mzy*mzy);   
    const double64 szz (sigma1_*mxz*mxz + sigma2_*myz*myz + sigma3_*mzz*mzz);   
    const double64 syz (sigma1_*mxy*mxz + sigma2_*myy*myz + sigma3_*mzy*mzz);   
    const double64 szx (sigma1_*mxz*mxx + sigma2_*myz*myx + sigma3_*mzz*mzx);
    
    stress(0,0) = sxx;
    stress(0,1) = sxy;
    stress(0,2) = szx; //sxz
    
    stress(1,0) = sxy; //syx
    stress(1,1) = syy;
    stress(1,2) = syz;
    
    stress(2,0) = szx;
    stress(2,1) = syz; //szy
    stress(2,2) = szz;
    
  } // end CartesianStressTensor
 
  
  
 /**   CartesianStressTensor
 
Same as previous method but with principal stresses augmented by isostatic stress
in order to take into account offset of current location from where the 
stress was measured.
      
 */
 void StressRotate::CartesianStressTensor( TensorVariable<3U>& stress, double64 iso_stress_offset ) const
  {
    Point<3U> n2 (crossProduct(n3_,n1_));
    n2.NormalizeLengthTo(1.);
    
    const double64 mxx ( n1_[0] ), mxy ( n1_[1] ), mxz ( n1_[2] ) ; 
    const double64 myx ( n2[0] ), myy ( n2[1] ), myz ( n2[2] ) ; 
    const double64 mzx ( n3_[0] ), mzy ( n3_[1] ), mzz ( n3_[2] ) ; 
    
    const double64 s1_offs = sigma1_ + iso_stress_offset;
    const double64 s2_offs = sigma2_ + iso_stress_offset;
    const double64 s3_offs = sigma3_ + iso_stress_offset;
    
    const double64 sxx (s1_offs * mxx*mxx + s2_offs * myx*myx + s3_offs * mzx*mzx);
    const double64 sxy (s1_offs * mxx*mxy + s2_offs * myx*myy + s3_offs * mzx*mzy);
    const double64 syy (s1_offs * mxy*mxy + s2_offs * myy*myy + s3_offs * mzy*mzy);
    const double64 szz (s1_offs * mxz*mxz + s2_offs * myz*myz + s3_offs * mzz*mzz);
    const double64 syz (s1_offs * mxy*mxz + s2_offs * myy*myz + s3_offs * mzy*mzz);
    const double64 szx (s1_offs * mxz*mxx + s2_offs * myz*myx + s3_offs * mzz*mzx);
    
    stress(0,0) = sxx;
    stress(0,1) = sxy;
    stress(0,2) = szx; //sxz
    
    stress(1,0) = sxy; //syx
    stress(1,1) = syy;
    stress(1,2) = syz;
    
    stress(2,0) = szx;
    stress(2,1) = syz; //szy
    stress(2,2) = szz;
    
  } // end CartesianStressTensor (with isostatic offset)
 
  
  
  
/**
    Rotates the internally stored stress state around one of the principal coordinate
    axes, x, y or z by the parameter angle (in degrees).
     
    Method applies the Standard Rotation Matrices to each principal stress vector individually 
    Weisstein, Eric W. "Rotation Matrix." From MathWorld--A Wolfram Web Resource.
    http://mathworld.wolfram.com/RotationMatrix.html 
    
    @param axes   is either of the x, y or z coordinate axes in the Kartesian system
    @param angle  is the angle in degress (0..360)
*/
void StressRotate::Rotate( char axis, double64 angle ) 
  {
    if ( axis != 'x' and axis != 'y' and axis != 'z' )
      throw Exception( CSMP_ERROR, "StressRotate::Rotate", "axis parameter (1) could not be identified; should be x,y or z");
  
    if ( angle < 0. or angle > 360. )
      throw Exception( CSMP_ERROR, "StressRotate::Rotate", "rotation angle is in degrees and must be between 0 and 360");

    //conversion from degrees to radians
    //conversion from yaw-pitch-roll to right hand side rule
    const double64 angle_r( -angle * PI / 180. );
    
    TensorVariable<3U> rotation;
    
    if( axis == 'x' )
    {
       rotation(0,0)=1.;
       rotation(0,1)=0.;
       rotation(0,2)=0.;
       
       rotation(1,0)=0.;
       rotation(1,1)=cos(angle_r);
       rotation(1,2)=sin(angle_r);
       
       rotation(2,0)=0.;
       rotation(2,1)=-sin(angle_r);
       rotation(2,2)=cos(angle_r);
       
       x_ += angle;
    }
    else if( axis == 'y' )
    {
       rotation(0,0)=cos(angle_r);
       rotation(0,1)=0.;
       rotation(0,2)=-sin(angle_r);
       
       rotation(1,0)=0.;
       rotation(1,1)=1.;
       rotation(1,2)=0.;
       
       rotation(2,0)=sin(angle_r);
       rotation(2,1)=0.;
       rotation(2,2)=cos(angle_r);
       
       y_ += angle;
    }
    else if( axis == 'z' )
    {
       rotation(0,0)=cos(angle_r);
       rotation(0,1)=sin(angle_r);
       rotation(0,2)=0.;
       
       rotation(1,0)=-sin(angle_r);;
       rotation(1,1)=cos(angle_r);;
       rotation(1,2)=0.;
       
       rotation(2,0)=0.;
       rotation(2,1)=0.;
       rotation(2,2)=1.;
       
       z_ += angle;
    }
    else
    throw csmp::Exception( CSMP_ERROR, "StressRotate<dim>::Rotate",
                          "axis could not be recognised: possibilities are x,y,and z");      
 
    n1_ = (rotation * VectorVariable<3U>(n1_.Coordinates())).P();
    n3_ = (rotation * VectorVariable<3U>(n3_.Coordinates())).P();
    
  } // end Rotate
  


/// outputs state of object
void StressRotate::Out(std::ostream& os) const
  {
    os << "\nStressRotate: ";
    os << "\nSigma1: " << sigma1_ << " Sigma2: " << sigma2_ << " Sigma3: " << sigma3_;
    os << "\nOriginal Principal Vectors:";
    os << "\nN1: " << original_n1_[0] << ", " << original_n1_[1] << ", " << original_n1_[2]; 
    os << "\nN3: " << original_n3_[0] << ", " << original_n3_[1] << ", " << original_n3_[2]; 
    os << "\nPrincipal Vectors:";
    os << "\nN1: " << n1_[0] << ", " << n1_[1] << ", " << n1_[2]; 
    os << "\nN3: " << n3_[0] << ", " << n3_[1] << ", " << n3_[2]; 
    os << "\nX Axis: " << x_ << " Y Axis: " << y_ << " Z Axis: " << z_;
    os << "\nCartesian Stress Tensor:"<< endl;
    TensorVariable<3U>  stress;
    CartesianStressTensor(stress);
    os << stress;
    os.flush();
  }
 
  
/// by default n1, n2 and n3 are normalized to sigma1, sigma2, sigma3 respectively
void StressRotate::OutputToVTK( const string& filename, 
                                const Point<3U>& location,
                                bool normalize_by_sigma3 ) const
  { 
    //calculate n2
    Point<3U> n1 (n1_);
    Point<3U> n2 (crossProduct(n3_,n1_));
    Point<3U> n3 (n3_);
    
    //normalize to sigma1, sigma2, sigma3 respectively
    n1.NormalizeLengthTo( fabs(sigma1_) );
    n2.NormalizeLengthTo( fabs(sigma2_) );
    n3.NormalizeLengthTo( fabs(sigma3_) );

    if ( normalize_by_sigma3 ) {
         //normalize all to sigma3
         n1 /= fabs(sigma3_);
         n2 /= fabs(sigma3_);
         n3 /= fabs(sigma3_);
      }
    
    ofstream ofs((filename+".vtk").c_str(), ios::out);
    
    const char* red   ("1. 0. 0.");
    const char* green ("0. 1. 0.");
    const char* blue  ("0. 0. 1.");
    const char* black ("0. 0. 0.");
    
    ofs << "# vtk DataFile Version 2.0 \
            \nPrincipal Stresses       \
            \nASCII\
            \nDATASET UNSTRUCTURED_GRID\
            \nPOINTS 3 float";
    ofs << endl << location[0] << " " << location[1] << " " << location[2];
    ofs << endl << location[0] << " " << location[1] << " " << location[2];
    ofs << endl << location[0] << " " << location[1] << " " << location[2];
    ofs << "\nPOINT_DATA 3\
            \nVECTORS vectors float";
    ofs << endl << n1[0] << " " << n1[1] << " " << n1[2];//sigma1
    ofs << endl << n2[0] << " " << n2[1] << " " << n2[2];//sigma2
    ofs << endl << n3[0] << " " << n3[1] << " " << n3[2];//sigma3
    ofs << "\nCOLOR_SCALARS colors 3";
    ofs << endl << red << endl << green << endl << (sigma3_<0.?blue:black);
    
    cout <<"\nStressRotate::OutputToVTK: file '"<< filename <<"' written successfully.\n";
  }



/**
    Normalizes and scales the vectors by the scale factor at the same time.
*/
void StressRotate::OutputToVTK( const string& filename, 
                                const Point<3U>& location,
                                double scale_factor ) const
  { 
    //calculate n2
    Point<3U> n1 (n1_);
    Point<3U> n2 (crossProduct(n3_,n1_));
    Point<3U> n3 (n3_);
    
    // normalization to sigma1, sigma2, sigma3 respectively
    n1.NormalizeLengthTo( fabs(sigma1_) );
    n2.NormalizeLengthTo( fabs(sigma2_) );
    n3.NormalizeLengthTo( fabs(sigma3_) );
    n1 /= (1./scale_factor) * fabs(sigma3_);
    n2 /= (1./scale_factor) * fabs(sigma3_);
    n3 /= (1./scale_factor) * fabs(sigma3_);
    
    ofstream ofs((filename+".vtk").c_str(), ios::out);
    
    const char* red   ("1. 0. 0.");
    const char* green ("0. 1. 0.");
    const char* blue  ("0. 0. 1.");
    const char* black ("0. 0. 0.");
    
    ofs << "# vtk DataFile Version 2.0 \
            \nPrincipal Stresses       \
            \nASCII\
            \nDATASET UNSTRUCTURED_GRID\
            \nPOINTS 3 float";
    ofs << endl << location[0] << " " << location[1] << " " << location[2];
    ofs << endl << location[0] << " " << location[1] << " " << location[2];
    ofs << endl << location[0] << " " << location[1] << " " << location[2];
    ofs << "\nPOINT_DATA 3\
            \nVECTORS vectors float";
    ofs << endl << n1[0] << " " << n1[1] << " " << n1[2];//sigma1
    ofs << endl << n2[0] << " " << n2[1] << " " << n2[2];//sigma2
    ofs << endl << n3[0] << " " << n3[1] << " " << n3[2];//sigma3
    ofs << "\nCOLOR_SCALARS colors 3";
    ofs << endl << red << endl << green << endl << (sigma3_<0.?blue:black);
    
    cout <<"\nStressRotate::OutputToVTK: file '"<< filename <<"' written successfully.\n";
  }




} // end csmp

