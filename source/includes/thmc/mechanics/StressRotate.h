#ifndef CSMP_STRESS_ROTATE_H
#define CSMP_STRESS_ROTATE_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class TensorVariable;


class StressRotate {
  public:
    StressRotate( const Point<3U>& n1, 
                  const Point<3U>& n3, 
                  double64 sigma1, double64 sigma2, double64 sigma3 );
  
    StressRotate( const StressRotate& );
                  
    ~StressRotate();
  
    /// returns the full Cartesian stress tensor rotated to align principal axes according to internally stored trend
    void CartesianStressTensor( TensorVariable<3U>& stress ) const;
  
    /// as previous, but augmentation modification of principal stresses by user-supplied isostatic stress
    void CartesianStressTensor( TensorVariable<3U>& stress, double64 iso_stress_offset ) const;
    
    /// counter-clockwise rotation of stress vectors (n1,n3) looking down on the rotation axis (use 180-angle to get cw.)
    void Rotate( char axis, double64 angle ); 
  
    /// loss free restoration of stress to state before rotation
    void Reset();
   
    /// maximum (compressive) stress (Pa/m2), corresponding to max Eigenvalue of stress tensor
    double64 Sigma1() const;
    /// intermediate principal stress (Pa/m2)
    double64 Sigma2() const;
    /// minimum compressive stress (Pa/m2) = mimimum Eigenvalue of stress tensor
    double64 Sigma3() const;
    
    void Out() const { Out(std::cout); }
    void Out(std::ostream& os) const;
  
    /// output 3 vectors representing the principal axes of the stress tensor
    void OutputToVTK( const std::string& filename, 
                      const Point<3U>& location, 
                      bool normalize_by_sigma3=true ) const;

    /// output 3 vectors representing the principal axes of the stress tensor and scale them by factor
    void OutputToVTK( const std::string& filename,
                      const Point<3U>& location, 
                      double scale_factor ) const;
  
  protected:
    Point<3U> original_n1_, original_n3_;
    Point<3U> n1_, n3_;
    double64 sigma1_, sigma2_, sigma3_;
    
    //rotations x,y,z
    double64 x_;
    double64 y_;
    double64 z_;
};


// inline functions

// maximum (compressive) stress (Pa/m2), corresponding to max Eigenvalue of stress tensor
inline double64 StressRotate::Sigma1() const
  { return sigma1_; }
  
// intermediate principal stress (Pa/m2)
inline double64 StressRotate::Sigma2() const
  { return sigma2_; }
  
// minimum compressive stress (Pa/m2) = mimimum Eigenvalue of stress tensor
inline double64 StressRotate::Sigma3() const
  { return sigma3_; }
  




} // end csmp

#endif
