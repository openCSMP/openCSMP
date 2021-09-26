#ifndef CSMP_MATH_UTILITIES_H
#define CSMP_MATH_UTILITIES_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"

// minor convienient additions to standard C++ functionality

namespace csmp {

/*! \file CSMP_mathUtilities.h */

/**
@defgroup CSMPglobalVariables Global CSMP Variables
@{
*/

const double64 PI( 3.14159265358979324 );

/**
@}
*/

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// all trigonomic functions in C++ take arguments in radians
inline double64 degreesToRadians( double64 deg ) { return deg * PI/180.; }
inline double64 radiansToDegrees( double64 rad ) { return rad * 180./PI; }


/**
    Shave off decimal places from floating point values so that numbers that contain
    noise match for a given target precision.
    
    @param scale 10^decimal places, i.e. 1e-5 to get 5 decimal places
    
*/
inline double64 quantiseToScale(double64 x, double64 scale)
{
    double ipart;
    const double frac = std::modf(x, &ipart);
    return ipart + scale * (int64_t)(frac / scale + 0.5);
} 


/// sign function to determine a positive or negative multiplier= -1 0 1 extracting the sign of a number
template <typename T>
inline T sign( T val ) {
    return static_cast<T>((T(0) < val) - (val < T(0)));
}

/// Square function
template<typename T>
inline T square( T val ) {
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
void     vector_randomize( std::vector<double64>& x, double64 scale_fac=1. );

/// averaging the floating-point values of scalar, vector and tensor variables in CSMP
template<typename Var> void average( const std::vector<Var>&, Var& );

/// FEM matrix transformations, for the case where the solution variable is a vector in 2D, Zienkewicz, volume 1, p. 22
void  dN_To2DOF( size_t nodes, DenseMatrix<DM_MIN>& DN );

/// FEM matrix transformations, for the case where the solution variable is a vector in 3D, Zienkewicz, volume 1, p. 133
void  dN_To3DOF( size_t nodes, DenseMatrix<DM_MIN>& DN );


/// auxiliary functions for spline interpolation
double64 splineValue( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2);
double64 splineDerivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2);
double64 splineSecondDerivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2);
  

// ROOT FINDING

/// finding root of function f(x) by Secant method within the x range [xmin,xmax]
double64 secant_method( double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 tolerance);

/// in the range [xmin,xmax], finds the intersection of the nonlinear function (flinear) with a linear function (fnonlinear) with the given x-axis intercept
double64 secant_line( double64 x_intercept, double64 xmin, double64 xmax, double64 (*flinear)( double64 ), double64 (*fnonlinear)( double64 ), double64 tolerance );
  
/// finding the x value, in the range [xmin,xmax],  where function f(x) is maximum by using golden-section search
double64 maximum_of_function( double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 tolerance);

/// finding the x value, in the range [xmin,xmax],  where function f(x) is minimum by using golden-section search
double64 minimum_of_function( double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 tolerance);

/// finding the x value, in the range [xmin,xmax],  where function the "df(x)/dx=(f(x)-f(x1))/(x-x1)".. this line is the secant and tangent to f(x) at that point
double64 g( double64 x, double64 x1, double64 (*function)( double64 ), double64 (*dfunction)( double64 ) );
  

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


/// linearly interpolate between two values
template<typename T>
T lerp(const double64 t, T x0, T x1) {
    return (1.0-t) * x0 + t * x1;
}


} // csmp

#endif
