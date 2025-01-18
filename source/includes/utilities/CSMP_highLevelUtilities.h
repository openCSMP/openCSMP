#ifndef CSMP_HIGH_LEVEL_UTILITIES_H
#define CSMP_HIGH_LEVEL_UTILITIES_H

#include "CSMP_definitions.h"

namespace csmp {

/**
@file CSMP_highLevelUtilities.h
@brief Leftover global functions for which no suitable place has yet been found in the code.
@author S.K. Matthai
*/

class Standard_IO_Handler;
template<uint32_t> class Point;
template<uint32_t> class Element;

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// utility that tokenises string into substrings using the supplied delimiter(s).
std::vector<std::string> splitString( std::string str, char delimiter );

/// to fill whitespace in strings with character of choice, for instance '_'
void replaceWhiteSpaceBy( std::string&, char ascii_char );

/// converts types to strings without loss of precision (std::to_string limits to 6 significant digits)
template <typename T>
inline std::string number_to_string( const T& value ) {
  std::stringstream sstr;
  sstr << std::setprecision( 16 ) << std::ios::scientific << value;
  return sstr.str();
}


/// scientific-format string with appropriate digits (without loss of precision (std::to_string limits to 6 significant digits)
template <typename T>
inline std::string numberToString( const T& value, bool verbose=false ) {
  if constexpr( std::is_same<T,double>::value ) {
      // int snprintf ( char * s, size_t n, const char * format, ... );
      //                sign + digit +  dp +       digits          + e + sign + expo + \0 (15 significant digits)
      std::string dstring( 1 + 1 +      1  + (DBL_DECIMAL_DIG - 1) + 1 + 1    + 5, '\0');
      auto slen = snprintf( &dstring[0], dstring.size(), "%1.14e", value );
      dstring.resize(slen);
      if ( verbose ) std::cout <<"\nmain: floating point value: "<< dstring <<" vs. "<< value << std::endl;
      return dstring;
    }
  else if constexpr ( std::is_same<T,float>::value )
    {
      //                sign + digit +  dp +       digits          + e + sign + expo + \0 (8 significant digits)
      std::string fstring( 1 + 1 +      1  + (FLT_DECIMAL_DIG - 1) + 1 + 1    + 3, '\0');
      auto slen = snprintf( &fstring[0], fstring.size(), "%1.7e", value );
      fstring.resize(slen);
      if ( verbose ) std::cout <<"\nmain: floating point value: "<< fstring <<" vs. "<< value << std::endl;
      return fstring;
    }
  return number_to_string( value );
}


/// evaluates the distance between 2 points
bool areFartherApartThan( const double* pn, const double* pw, double distance );


/// renumber node index values for a range of elements in
/// a continuous sequence from (0..n-1); the number of unique nodes is returned
size_t  renumberElementNodes( std::vector<Element<1U>*>::iterator first,
                              std::vector<Element<1U>*>::iterator last );

size_t  renumberElementNodes( std::vector<Element<2U>*>::iterator first,
                              std::vector<Element<2U>*>::iterator last );

size_t  renumberElementNodes( std::vector<Element<3U>*>::iterator first,
                              std::vector<Element<3U>*>::iterator last );


/**
@brief Vector element subtraction.

@author  P. Lang
@date    06/11/2012

@param [in,out]  toRemoveFrom  vector from which elements are deleted.
@param [in,out]  toRemove      vector which contains elements to be deleted from toRemoveFrom

@return  new size of vector.

*/
template<typename T>
size_t removeVectorElements( std::vector<T>& toRemoveFrom, std::vector<T>& toRemove )
{
  std::sort( toRemoveFrom.begin(), toRemoveFrom.end() );
  std::sort( toRemove.begin(), toRemove.end() );
  std::vector<T> toRemoveFromNew;
  set_difference( toRemoveFrom.begin(), toRemoveFrom.end(), toRemove.begin(), toRemove.end(), std::back_inserter( toRemoveFromNew ) );
  const size_t deleted( toRemoveFrom.size() - toRemoveFromNew.size() );
  toRemoveFrom.swap( toRemoveFromNew );
  std::vector<T> cache( toRemoveFrom );
  toRemoveFrom.swap( cache );
  return deleted;
}

/// find all possible combinations of the numbers provided in 'sequence', where 'samples' specified how many numbers shall be combined; combinations are returned into deque.
template<typename intType>
size_t createUniqueCombinations( std::vector<intType>& sequence, intType samples,
                                 std::deque<std::vector<intType> >& combinations );


/// Returns true if file on ifstream is empty.
bool isInputFileEmpty( std::ifstream& );


#if defined _MSC_VER || defined __MINGW32__
char * strptime( const char *s, const char *format, struct tm *tm );
#endif


} // end csmp


#endif
