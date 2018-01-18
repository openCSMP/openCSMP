#ifndef CSMP_HIGH_LEVEL_UTILITIES_H
#define CSMP_HIGH_LEVEL_UTILITIES_H

#include "CSMP_definitions.h"
#include "CSMP_random.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "Exception.h"
#include "FiniteElement.h"
#include "binaryReadWrite.h"

namespace csmp {

/**
@file CSMP_highLevelUtilities.h
@brief global csmp funtions
@author S.K. Matthai
*/

class Standard_IO_Handler;
template<size_t> class Point;
template<size_t> class Visitor;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Model;
template<size_t> class VSet;
template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Region;
template<size_t> class MeshManager;

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// to fill whitespace in strings with character of choice, for instance '_'
void replaceWhiteSpaceBy( std::string&, char ascii_char );

/// converts integral types to strings without loss of precision (std::to_string limits to 6 significant digits)
template <typename T>
inline std::string number_to_string( const T& value ) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

/// returns intermediate (true) or maximum (false) model dimensions
template<size_t  dim>
double64  printModelDimensions( const Model<dim>&, bool intermed_or_max=false ); 

/// calculates the center of gravity of the model
template<size_t  dim>
Point<dim>  centerOfGravity( const Model<dim>& );

/// prints range to screen; returns either min(arg=false) or maximum variable value (default)
template<size_t  dim>
double64  printRangeOfVariable( const Model<dim>&, 
                                const char* var, bool print_maximum=true );

/// prints range of target variable in model to screen and logs it to IO handler
template<size_t  dim>
double64  printRangeOfVariable( const Model<dim>&, 
                                Standard_IO_Handler& io, const char* var,
                                bool max_instead_of_min=true );

/// prints range of target variable within specific model subdomain
template<size_t  dim>
double64  printRangeOfVariable( const Model<dim>&, 
                                const char* group, const char* var, bool max_or_min=true );

/// prints range of target variable within specific model subdomain and logs it to IO handler
template<size_t  dim>
double64  printRangeOfVariable( const Model<dim>&, 
                                Standard_IO_Handler&,
                                const char* group, const char* var,
                                bool max_instead_of_min=true );

/// prints min/max values stored in supplied vector
void printRangeOf( const std::vector<std::pair<double64,double64> >& );

/// prints min/max values stored in supplied vector of vectors
void printRangeOfVectorOfVectors( const std::vector<std::vector<double64> >&  );

/// convert the flag(s) of a variable into integer values stored in its number part
template<size_t dim> 
void flagToNumber( Model<dim>&, const char* variable );
 
/// convert the flag(s) of first variable into integer values stored in the second variable
template<size_t dim> 
void flagToNumber( Model<dim>&, const char* flag_variable, const char* number_variable );

/// using random number generator, adds percentage of Gaussian noise to variable values
template<size_t dim> 
void randomPerturb( random_generator& gen, Model<dim>&, const char* prop, double64 by_percent_of_max_value );

/// compares the mesh connectivity in the model with that of the input vset; returns true if both have the same
template<size_t dim>
bool compareConnectivity( const Model<dim>&, const VSet<dim>& );

/// evaluates the distance between 2 points
bool areFartherApartThan( const double64* pn, const double64* pw, double64 distance );

/// Returns true if file on ifstream is empty.
bool isInputFileEmpty( std::ifstream& file );
  
/// container of element pointers and local face ids of elements contacting each other across a split boundary
typedef std::pair<std::pair<Element<3U>*,size_t>,std::pair<Element<3U>*,size_t> > OppositeElements;
/// find all elements in a model that contact eachother across split interfaces and are node-matched
bool findSplitInterfaceElements( const Region<3U>&, std::set<OppositeElements>& );

/// finds node by point coordinate; returns -1 if not found; @attention tolerance needs to account for single-precision of CAD tools
template<size_t dim>
long  findNode( const Model<dim>&, const Point<dim>& pxyz, double64 tolerance, bool verbose=false );

/// find node by its position as identified from its coordinates: tolerance should take into account single-precision of CAD tools
size_t  findNode( const Model<1U>&, 
                  double64 nx, double64 tolerance );
/// 2D version
size_t  findNode( const Model<2U>&, 
                  double64 nx, double64 ny, double64 tolerance );
/// 3D version
size_t  findNode( const Model<3U>&, 
                  double64 nx, double64 ny, double64 nz, double64 tolerance );

/// retrieves and returns the first contiguous element patch that can be reached by mesh traversal from the starting element
template<size_t dim> 
void floodFill( Element<dim>* const eptr, std::set<Element<dim>*>& output_contiguous_subset );

/// retrieves and returns the element ids of the first contiguous element patch that can be reached by mesh traversal from the starting element
template<size_t dim> 
void floodFillViaIndexes( const Region<dim>&, size_t starting_idx,
                          std::set<size_t>& output_contiguous_subset );
  
/// determines whether mesh in model is built from finite elements with a local coordinate system
template<size_t dim>
bool isoparametricElementMesh( const Model<dim>& );

/// checks region for whether it contains elements of the same dimensionality
template<size_t dim>
bool containsElementsOfTtype( const Region<dim>&, ELEMENT_DIMENSION );

/// recreates neighbor connectivity among all equidimensional elements (volumetric-, surfacic- and line elements) 
template<size_t dim>
void  establishNeighborConnectivity( std::vector<Element<dim>*>& simplexVector, bool unassign_neighbors_outside = false );

template<size_t dim>
void  establishNeighborConnectivity( std::vector<Face<dim>*>& simplexVector, bool unassign_neighbors_outside = false );

template<size_t dim>
void  establishNeighborConnectivity( std::vector<InterFace<dim>*>& simplexVector, bool unassign_neighbors_outside = false );

/// checks all elements of the surface region for whether their neighbor elements have normals that deviate less than 90o from their normals
template<size_t dim>
bool checkNeighborNormalsForConsistentOrientation( const Region<dim>& );

/// calculates the absolute value of the maximum difference between distributed variable values with the same type and placement
template<size_t dim>
double64 maximumResidual( const Model<dim>&,
                          const char* new_property, const char* old_property, 
                          bool absolute );

/// in target region, element variable is extrapolated to node and back as many times as indicated by n_smoothing_cycles
template<size_t dim>  
void smoothElementVariable( Model<dim>&, const char*, const char* element_var, const char* temp_node_var, size_t n_smoothing_cycles );

/// replaces no-data values of target variable with nearest-neighbor values until there are none left, by default NAN's are no-data values
template<size_t dim>  
void nearestNeighborFill( Model<dim>&, const char* target_region, const char* variable, double64 no_data_value );

/// imposes either an upper- or a lower limit on the variable in the region of interest
template<size_t dim>
void imposeLimitOn( Model<dim>& model, const char* region, const char* variable, bool upper_limit, double64 limit_value );

/// barycentre-to-node distances for parent elements returned into vector [e1,e2...e_n,e_sum] with a length of parent elements+1
template<size_t dim>
void distanceWeights( typename std::vector<Node<dim>*>::const_iterator nodes_begin,
                      typename std::vector<Node<dim>*>::const_iterator nodes_end,
                      std::vector<std::vector<double64> >& distances_and_weight );


/// linear interpolation to point between 2 points in 1D
double64 linearInterpolate( const std::pair<Point<1U>,double64>&, 
                            const std::pair<Point<1U>,double64>&,
                            const Point<1U>& );            

/// bilinear interpolation to point between 2 points in 2D
double64 linearInterpolate( const std::pair<Point<2U>,double64>&,
                            const std::pair<Point<2U>,double64>&,
                            const Point<2U>& );            

/// 3D trilinear interpolation
double64 linearInterpolate( const std::pair<Point<3U>,double64>&, 
                            const std::pair<Point<3U>,double64>&,
                            const Point<3U>& );            

/// just a stub (ignores the bilinear)
double64 bilinearInterpolate( size_t idx_x, size_t idx_y, 
                              const Point<1U>& xy1, 
                              const Point<1U>& xy2,
                              const Point<1U>& coord, 
                              double64 p1, double64 p2, double64, double64 );

double64 bilinearInterpolate( size_t idx_x, size_t idx_y, 
                              const Point<2U>& xy1, 
                              const Point<2U>& xy2,
                              const Point<2U>& coord, 
     //                       val@x0,y0  val@x1,y0  val@x1,y1  val@x0,y1
                              double64 p1, double64 p2, double64 p3, double64 p4 );

double64 bilinearInterpolate( size_t idx_x, size_t idx_y, 
                               const Point<3U>& xy1, 
                               const Point<3U>& xy2,
                               const Point<3U>& coord, 
      //                       val@x0,y0  val@x1,y0  val@x1,y1  val@x0,y1
                               double64 p1, double64 p2, double64 p3, double64 p4 );

/// extrapolate element properties from VSet to node properties
template<size_t dim,class VarType>
void extrapolateElementToNodeProperty( csmp::VSet<dim>&             vset,
                                       const std::vector<VarType>&  elmnt_data,
                                       std::vector<VarType>&        nodal_data );

/// mapping node coordinates to a node variable
template<size_t dim>
void assignNodeCoordinatesTo( Model<dim>& sg, const char coordinate, const char* nodal_variable );

/// for Triangulator meshes: for meshes created from pixels, finds outlier triangles (3 nodes at region boundaries) and flips their property values
void stripDomainEdgesFor( Model<2U>&, const char* el_prop );

/// renumber node index values for a range of elements in
/// a continuous sequence from (0..n-1); the number of unique nodes is returned
size_t  renumberElementNodes( std::vector<Element<1U>*>::iterator first,
                              std::vector<Element<1U>*>::iterator last );

size_t  renumberElementNodes( std::vector<Element<2U>*>::iterator first,
                              std::vector<Element<2U>*>::iterator last );
                                 
size_t  renumberElementNodes( std::vector<Element<3U>*>::iterator first,
                              std::vector<Element<3U>*>::iterator last );

/// additions by P. Lang (2012); SKM @todo explain what this is good for
namespace femDataOutputDispatch {

template<class Var>
inline void initVariable( csmp::Index, Var& ) {}

template<>
inline void initVariable( csmp::Index key, ArrayVariable& var )
  { var.Resize( key.dataDepth ); }

template<>
inline void initVariable( csmp::Index key, FlaggedArrayVariable& var )
  { var.Resize( key.dataDepth ); }

} // femDataOutputDispatch



/// domain (Model, Region...) variables binary IO
template<class V, class D, size_t dim>
bool variablesOut( FILE* fp, const D& domain, const PropertyDatabase<dim>& pref, VARIABLE_TYPE vtype )
  {
  size_t vcount( pref.VariableCount( domain.Placement(), vtype ) );
  fwrite( (void*)&vcount, sizeof(size_t), 1, fp );
  std::set<std::string> propList;
  pref.ListProperties( domain.Placement(), vtype, propList );
  for( std::set<std::string>::const_iterator it( propList.begin() ); it != propList.end(); ++it )
    {
    V var;
    Index key( pref.StorageKey( it->c_str() ) );
    femDataOutputDispatch::initVariable( key, var );
    domain.Read( key, var );
    skm_C_fwrite( fp, it->c_str() );
    if( !var.Out(fp) )
      return false;
    }
  return true;
  }


template<class V, class D, size_t dim>
bool variablesIn( FILE* fp, D& domain, const PropertyDatabase<dim>& pref, VARIABLE_TYPE )
  {
  size_t vcount(-1);
  fread( (void*)&vcount, sizeof(size_t), 1, fp );
  for( size_t i(0); i < vcount; ++i )
    {
    V var;
    char propName[200];
    skm_C_fread( fp, propName );
    Index key( pref.StorageKey(propName) );
    femDataOutputDispatch::initVariable( key, var );
    if( !var.In(fp) )
      return false;
    domain.Store( key, var );
    }
  return true;
  }

template<class D, size_t dim>
bool domainVariablesOut( FILE* fp, const D& domain, const PropertyDatabase<dim>& pref )
  {
    if( !variablesOut<ScalarVariable>       ( fp, domain, pref, SCALAR ) )
      return false;
    if( !variablesOut<VectorVariable<dim> > ( fp, domain, pref, VECTOR ) )
      return false;
    if( !variablesOut<TensorVariable<dim> > ( fp, domain, pref, TENSOR ) )
      return false;
    if( !variablesOut<ArrayVariable>        ( fp, domain, pref, ARRAY ) )
      return false;
    if( !variablesOut<FlaggedArrayVariable> ( fp, domain, pref, FLAGGEDARRAY ) )
      return false;
    return true;
  }

template<class D, size_t dim>
bool domainVariablesIn( FILE* fp, D& domain, const PropertyDatabase<dim>& pref )
  {
  if( !variablesIn<ScalarVariable>          ( fp, domain, pref, SCALAR ) )
    return false;
  if( !variablesIn<VectorVariable<dim> >    ( fp, domain, pref, VECTOR ) )
    return false;
  if( !variablesIn<TensorVariable<dim> >    ( fp, domain, pref, TENSOR ) )
    return false;
  if( !variablesIn<ArrayVariable>           ( fp, domain, pref, ARRAY ) )
    return false;
  if( !variablesIn<FlaggedArrayVariable>    ( fp, domain, pref, FLAGGEDARRAY ) )
    return false;
  return true;
  }





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
  std::sort(toRemoveFrom.begin(), toRemoveFrom.end());
  std::sort(toRemove.begin(), toRemove.end());
  std::vector<T> toRemoveFromNew;
  set_difference(toRemoveFrom.begin(),toRemoveFrom.end(), toRemove.begin(),toRemove.end(), std::back_inserter(toRemoveFromNew)); 
  const size_t deleted( toRemoveFrom.size() - toRemoveFromNew.size() );
  toRemoveFrom.swap(toRemoveFromNew);
  std::vector<T> cache(toRemoveFrom);
  toRemoveFrom.swap(cache);
  return deleted;
}


#if defined _MSC_VER || defined __MINGW32__
char * strptime(const char *s, const char *format, struct tm *tm);
#endif

/**
@}
*/

} // end namespace csmp

#endif
