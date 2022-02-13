#ifndef CSMP_HIGH_LEVEL_UTILITIES_H
#define CSMP_HIGH_LEVEL_UTILITIES_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"

namespace csmp {

/**
@file CSMP_highLevelUtilities.h
@brief global csmp funtions
@author S.K. Matthai
*/

class Standard_IO_Handler;
template<uint32_t> class Point;
template<uint32_t> class Visitor;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Model;
template<uint32_t> class VSet;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Region;
template<uint32_t> class MeshManager;

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// utility that tokenises string into substrings using the supplied delimiter(s).
std::vector<std::string> splitString( std::string str, char delimiter );

/// to fill whitespace in strings with character of choice, for instance '_'
void replaceWhiteSpaceBy( std::string&, char ascii_char );

/// converts integral types to strings without loss of precision (std::to_string limits to 6 significant digits)
template <typename T>
inline std::string number_to_string( const T& value ) {
  std::stringstream sstr;
  sstr << std::setprecision( 16 ) << value;
  return sstr.str();
}

/// evaluates the distance between 2 points
bool areFartherApartThan( const double* pn, const double* pw, double distance );

// finds node by point coordinate; returns -1 if not found; @attention tolerance needs to account for single-precision of CAD tools
template<uint32_t dim>
long  findNode( const Model<dim>&, const Point<dim>& pxyz, double tolerance, bool verbose = false );

/// find node by its position as identified from its coordinates: tolerance should take into account single-precision of CAD tools
size_t  findNode( const Model<1U>&,
                  double nx, double tolerance );
/// 2D version
size_t  findNode( const Model<2U>&,
                  double nx, double ny, double tolerance );
/// 3D version
size_t  findNode( const Model<3U>&,
                  double nx, double ny, double nz, double tolerance );

/// prints sorted global element node numbers in a compact way
template<uint32_t dim, template<uint32_t> class CELL>
void printNodes( const CELL<dim>& );

/// linear interpolation to point between 2 points in 1D
double linearInterpolate( const std::pair<Point<1U>, double>&,
                            const std::pair<Point<1U>, double>&,
                            const Point<1U>& );

/// bilinear interpolation to point between 2 points in 2D
double linearInterpolate( const std::pair<Point<2U>, double>&,
                            const std::pair<Point<2U>, double>&,
                            const Point<2U>& );

/// 3D trilinear interpolation
double linearInterpolate( const std::pair<Point<3U>, double>&,
                            const std::pair<Point<3U>, double>&,
                            const Point<3U>& );

/// just a stub (ignores the bilinear)
double bilinearInterpolate( size_t idx_x, size_t idx_y,
                              const Point<1U>& xy1,
                              const Point<1U>& xy2,
                              const Point<1U>& coord,
                              double p1, double p2, double, double );

double bilinearInterpolate( size_t idx_x, size_t idx_y,
                              const Point<2U>& xy1,
                              const Point<2U>& xy2,
                              const Point<2U>& coord,
                              //                       val@x0,y0  val@x1,y0  val@x1,y1  val@x0,y1
                              double p1, double p2, double p3, double p4 );

double bilinearInterpolate( size_t idx_x, size_t idx_y,
                              const Point<3U>& xy1,
                              const Point<3U>& xy2,
                              const Point<3U>& coord,
                              //                       val@x0,y0  val@x1,y0  val@x1,y1  val@x0,y1
                              double p1, double p2, double p3, double p4 );

/// extrapolate element property values stored in VSet to the vertices
template<uint32_t dim, class VarType>
void extrapolateElementToNodeProperty( csmp::VSet<dim>&             vset,
                                       const std::vector<VarType>&  elmnt_values,
                                       std::vector<VarType>&        nodal_values );

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
size_t createUniqueCombinations( std::vector<int64_t>& sequence, size_t samples,
                                 std::deque<std::vector<int64_t> >& combinations );




/// Returns true if file on ifstream is empty.
bool isInputFileEmpty( std::ifstream& );


#if defined _MSC_VER || defined __MINGW32__
char * strptime( const char *s, const char *format, struct tm *tm );
#endif


} // end csmp


#endif
