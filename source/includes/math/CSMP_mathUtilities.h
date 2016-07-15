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

const double64 PI( 3.14159265358979324F );

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
template <typename T> int sgn(T val) {
    return static_cast<T>((T(0) < val) - (val < T(0)));
}

/// conversion of numbers to strings (in >=C++11, use to_string() function)
void         uintToString( unsigned int i, std::string& );
std::string  uintToString( unsigned int i );
void         ulongToString( unsigned long i, std::string& );


/// rounds a double to an integer
int32  rint( double64 ); 

/// rounds to an unsigned long
long64  lrint( double64 );

/// 
template<typename var> void average( const std::vector<var>&, var& );

/// erf() and erfc() approximated with Chebyshev polynomials
double64  erf_Chebyshev( double64 );
double64  erfc_Chebyshev( double64 );

double64 vector_norm1( std::vector<double64>& x, std::vector<double64>& scale );
double64 vector_norm2( std::vector<double64>& x, std::vector<double64>& scale );
double64 vector_norm_inf( std::vector<double64>& x, std::vector<double64>& scale );
void     vector_randomize( std::vector<double64>& x, double64 scale_fac=1. );


/// min for 3 argument values
template <typename T> inline const T& min(const T& a, const T& b, const T& c)
{
    double64 cmin(std::min(a,b));
    return std::min( cmin, c );
}

/// max for 3 argument values
template <typename T> inline const T& max(const T& a, const T& b, const T& c)
{
    double64 cmax(std::max(a,b));
    return std::max( cmax, c );
}


// interpolation function matrix conversions for when the result variable is not a scalar, see Zienkewicz

/// Zienkewicz, volume 1, p. 22
void  dN_To2DOF( size_t nodes, DenseMatrix<DM_MIN>& DN );

/// Zienkewicz, volume 1, p. 133
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
