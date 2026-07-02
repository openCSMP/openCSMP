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

constexpr double PI = 3.14159265358979323846;

/**
@}
*/

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// all trigonomic functions in C++ take arguments in radians
inline double degreesToRadians( double deg ) { return deg * PI/180.; }
inline double radiansToDegrees( double rad ) { return rad * 180./PI; }


/**
    Shave off decimal places from floating point values so that numbers that contain
    noise match at  the given target precision.
    
    @param scale 10^decimal places, i.e. 1e-5 to get 5 decimal places
*/
inline double quantiseToScale(double x, double scale)
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
double
rsqrt( double val );

/// for conversion of numbers to strings use to_string() function

/// erf() and erfc() approximated with Chebyshev polynomials
double  erf_Chebyshev( double );
double  erfc_Chebyshev( double );

/// L1 norm
double vector_norm1( std::vector<double>& x, std::vector<double>& scale );

/// L2 norm
double vector_norm2( std::vector<double>& x, std::vector<double>& scale );

/// L-infinity norm
double vector_norm_inf( std::vector<double>& x, std::vector<double>& scale );

/// randomly perturbs the values stored in the supplied floating-point vector
void     vector_randomize( std::vector<double>& x, double scale_fac=1. );

/// averaging the floating-point values of scalar, vector and tensor variables in CSMP
template<typename Var> void average( const std::vector<Var>&, Var& );

/// FEM matrix transformations, for the case where the solution variable is a vector in 2D, Zienkewicz, volume 1, p. 22
void  dN_To2DOF( uint32_t nodes, DenseMatrix<DM_MIN>& DN );

/// FEM matrix transformations, for the case where the solution variable is a vector in 3D, Zienkewicz, volume 1, p. 133
void  dN_To3DOF( uint32_t nodes, DenseMatrix<DM_MIN>& DN );


/// auxiliary functions for spline interpolation
double splineValue( double x, double x1, double x2, double y1, double y2, double k1, double k2);
double splineDerivative( double x, double x1, double x2, double y1, double y2, double k1, double k2);
double splineSecondDerivative( double x, double x1, double x2, double y1, double y2, double k1, double k2);
  

// ROOT FINDING 

/// finding root of function f(x) by Secant method within the x range [xmin,xmax]
double secant_method( double xmin, double xmax, double (*function)( double ), double tolerance);

/// in the range [xmin,xmax], finds the intersection of the nonlinear function (flinear) with a linear function (fnonlinear) with the given x-axis intercept
double secant_line( double x_intercept, double xmin, double xmax, double (*flinear)( double ), double (*fnonlinear)( double ), double tolerance );
  
/// finding the x value, in the range [xmin,xmax],  where function f(x) is maximum by using golden-section search
double maximum_of_function( double xmin, double xmax, double (*function)( double ), double tolerance);

/// finding the x value, in the range [xmin,xmax],  where function f(x) is minimum by using golden-section search
double minimum_of_function( double xmin, double xmax, double (*function)( double ), double tolerance);

/// finding the x value, in the range [xmin,xmax],  where function the "df(x)/dx=(f(x)-f(x1))/(x-x1)".. this line is the secant and tangent to f(x) at that point
double g( double x, double x1, double (*function)( double ), double (*dfunction)( double ) );
  

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
T lerp(const double t, T x0, T x1) {
    return (1.0-t) * x0 + t * x1;
}

/// helper for matrix indices
void print( std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&  v );


} // csmp

#endif
