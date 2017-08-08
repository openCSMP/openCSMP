#ifndef CSMP_MATH_UTILITIES_H
#define CSMP_MATH_UTILITIES_H

#include "CSMP_definitions.h"
#include "CSMP_random.h"
#include "DenseMatrix.h"

// minor convienient additions to standard C++ functionality

namespace csmp {

/*! \file CSMP_mathUtilities.h */

/**
@defgroup CSMPglobalVariables Global CSMP Variables
@{
*/

const double64 PI( 3.14159265358979324 );

/// all trigonomic functions in C++ take arguments in radians
inline double64 degreesToRadians( double64 deg ) { return deg * PI/180.; }
inline double64 radiansToDegrees( double64 rad ) { return rad * 180./PI; }

/**
@}
*/

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// sign function to determine a positive or negative multiplier= -1 0 1 extracting the sign of a number
template <typename T>
T sign( T val ) {
    return static_cast<T>((T(0) < val) - (val < T(0)));
}

/// Square function
template<typename T>
T square( T val ) {
    return val*val;
}

/// Reciprocal square root function
double64
rsqrt( double64 val );

/// for conversion of numbers to strings use to_string() function

/// erf() and erfc() approximated with Chebyshev polynomials
double64  erf_Chebyshev( double64 );
double64  erfc_Chebyshev( double64 );

/// L1 norm
double64 vector_norm1( std::vector<double64>& x, std::vector<double64>& scale );

/// L2 norm
double64 vector_norm2( std::vector<double64>& x, std::vector<double64>& scale );

/// L-infinity norm
double64 vector_norm_inf( std::vector<double64>& x, std::vector<double64>& scale );

/// randomly perturbs the values stored in the supplied floating-point vector
void     vector_randomize( random_generator& rng, std::vector<double64>& x, double64 scale_fac=1. );

/// averaging the floating-point values of scalar, vector and tensor variables in CSMP
template<typename Var> void average( const std::vector<Var>&, Var& );

/// FEM matrix transformations, for the case where the solution variable is a vector in 2D, Zienkewicz, volume 1, p. 22
void  dN_To2DOF( size_t nodes, DenseMatrix<DM_MIN>& DN );

/// FEM matrix transformations, for the case where the solution variable is a vector in 3D, Zienkewicz, volume 1, p. 133
void  dN_To3DOF( size_t nodes, DenseMatrix<DM_MIN>& DN );

/**
@}
*/

/// console output of arguments like STL vectors etc.
template<typename Var> void out( const Var& obj );

/// scalar product of any two STL vectors (loop is unrolled)
template<unsigned int dim, typename T>
class vectorDotProduct {
   public:
      /// inline function is defined within class declaration
     static T Result( typename std::vector<T>::const_iterator& ait, 
                      typename std::vector<T>::const_iterator& bit ) {  
          // look at the template specifier: it gets decremented here                    
          return (*ait) * (*bit) + vectorDotProduct<dim-1U,T>::Result( ++ait, ++bit );
       }
 };

/// partial specialization of previous template needed to stop recursion
template<typename T>
class vectorDotProduct<1U,T> {
   public:
     static T Result( typename std::vector<T>::const_iterator& ait, 
                      typename std::vector<T>::const_iterator& bit ) { 
          return (*ait) * (*bit);
       }
 };

/// inlined convenience function to facilitate multiplication of STL vectors
template<unsigned int dim, typename T>
T vector_product( const std::vector<T>& a, 
                  const std::vector<T>& b ) {
     typename std::vector<T>::const_iterator ait(a.begin()); 
     typename std::vector<T>::const_iterator bit(b.begin());
     return vectorDotProduct<dim,T>::Result( ait, bit );
  }


} // csmp

#endif
