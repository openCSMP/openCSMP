#include "CSMP_highLevelUtilities.h"
#include "compareFloats.h"
#include "Box.h"
#include "VSet.h"
#include "Node.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "ModelTime.h"
#include "TextInterface.h"
#include "Matrix.h"
#include "convertColorToPermeability.h"
#include "Triangulator.h"
#include "MeshManager.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearLineElement.h"

using namespace std;

namespace csmp {

// strings
void replaceWhiteSpaceBy( string& p, char ascii_char )
 {
    if ( !p.empty() )
      for ( size_t i=0U; i<p.size(); i++ )
        if ( p[i] == ' ' || p[i] == '\t' ||
             p[i] == '\n' || p[i] == '\r' ) p[i] = ascii_char;
 }



template<size_t dim>
bool isoparametricElementMesh( const Model<dim>& sg )
 {
   std::string  etype(parseFiniteElementType((*sg.Region("Model").ElementsBegin())->FE()->ElementType()));

     if ( etype.find("ISOPARAMETRIC") != std::string::npos ) return true;
     return false;

 }// end isoparametricElementMesh
template bool isoparametricElementMesh<1U>( const csmp::Model<1U>& );
template bool isoparametricElementMesh<2U>( const Model<2U>& );
template bool isoparametricElementMesh<3U>( const Model<3U>& );


/**

Evaluates distance between 2 points. If it is smaller that supplied
value, the method returns false, else it returns true.

@section arguments Input Arguments 

The locations of the points that shall be evaluated.

@return The minimum distance that the points should be apart from one another.
*/
bool areFartherApartThan( const double64* pn, const double64* pw, double64 distance )
 {
    double64 dx = pn[0] - pw[0];
    double64 dy = pn[1] - pw[1];
    double64 dz = pn[2] - pw[2];
    
    double64 separation = sqrt( dx*dx + dy*dy + dz*dz );
    
    if ( distance < separation ) return true;
    
    return false;
    
 } // end

/**
 * Returns true if file on ifstream is empty.(Aug 2014)
 * @author Julian E. Mindel
 */
bool isInputFileEmpty( std::ifstream& pFile )
{
   return pFile.peek() == std::ifstream::traits_type::eof();
}

/**

Returns the ID (0..n-1) of any node which is located within the tolerance of the
target coordinates. If the node cannot be found it returns max() of index type. 

@section arguments Input Arguments 

The current model that shall be searched for the node, the node 
coordinates, and the tolerance which shall be applied in comparing the
supplied coordinates with those of the actual node points.  

@return The Idx (0..n-1) of the node of interest or UINTMAX 
(if this node does not exist) will
be returned.  

@section application Application

The method is used to retrieve point locations from the mesh in order to
identify points that cannot be grouped into individual families using
the ANSYS mesher.  

*/
size_t  findNode( const Model<3U>& sg, double64 nx, double64 ny, double64 nz, 
                  double64 tolerance )
 {
    const Region<3>&  sgroup(sg.Region("Model"));
 
    for ( vector<csmp::Node<3U>*>::const_iterator 
          it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ ) {
//         if ( approximatelyEqual( nx, (*it)->x(), tolerance ) &&
//              approximatelyEqual( ny, (*it)->y(), tolerance ) &&
//              approximatelyEqual( nz, (*it)->z(), tolerance ) )
         if ( fabs(nx-(*it)->x()) <= tolerance and
              fabs(ny-(*it)->y()) <= tolerance and
              fabs(nz-(*it)->z()) <= tolerance )
           return (*it)->Idx();
      }
 
    stringstream  out("The targeted node with the coordinate (x,y,z): ");
    out << nx <<" "<< ny <<" "<< nz <<" could not be found; ";
    out <<" returning node index="<< UINT_MAX << endl;
    throw csmp::Exception( WARNING, "findNode", out.str() );
    
    return UINT_MAX;
     
} // end find_node



size_t  findNode( const Model<2U>& sg, double64 nx, double64 ny,  
                  double64 tolerance )
 {
    const Region<2>&  sgroup(sg.Region("Model"));
 
    for ( vector<Node<2U>*>::const_iterator 
          it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ ) {
//         if ( approximatelyEqual( nx, (*it)->x(), tolerance ) &&
//              approximatelyEqual( ny, (*it)->y(), tolerance ) )
           if ( fabs(nx-(*it)->x()) <= tolerance and
                fabs(ny-(*it)->y()) <= tolerance )
             return (*it)->Idx();
      }
 
    stringstream  out("The targeted node with the coordinates (x,y): ");
    out << nx <<" "<< ny <<" could not be found; ";
    out <<" returning node index="<< UINT_MAX << endl;
    throw csmp::Exception( WARNING, "findNode", out.str() );
    
    return UINT_MAX;
     
} // end find_node



size_t  findNode( const Model<1U>& sg, double64 nx, double64 tolerance )
 {
    const Region<1>&  sgroup(sg.Region("Model"));
 
    for ( vector<Node<1U>*>::const_iterator 
          it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ ) {
//         if ( approximatelyEqual( nx, (*it)->x(), tolerance ) )
         if ( fabs(nx-(*it)->x()) <= tolerance ) return (*it)->Idx();
      }
 
    stringstream  out("The targeted node with the coordinate (x): ");
    out << nx <<" could not be found; ";
    out <<" returning node index="<< UINT_MAX << endl;
    throw csmp::Exception( WARNING, "findNode", out.str() );
    
    return UINT_MAX;
     
} // end find_node


/**
    Generic version for 1-3 dimensions, using Point object to identify the node
    location.
    
    If a node is found its local Idx() number is returned 
    (care has to be taken that this index is a valid number.
    
    If the node cannot be found, -1, is returned.
    
    @attention if verbose is on and the point cannot be found this is reported.
    
    @author SKM 22/9/2014.
*/
template<size_t dim>
long  findNode( const Model<dim>& sg, const Point<dim>& pxyz, double64 tolerance, bool verbose )
 {
    const Region<dim>&  sgroup(sg.Region("Model"));
 
    for ( typename vector<Node<dim>*>::const_iterator
          it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ )
      if ( pxyz.CoincidesWithWithinTolerance( (*it)->Coordinate(), tolerance ) )
        return (*it)->Idx();

    ErrorHandler& csmp_error(ErrorHandler::Instance());
   
    if ( verbose ) {
         stringstream  out("The targeted node with the coordinate (x): ");
         out << pxyz <<" could not be found; ";
         out <<" returning node index="<< -1 << endl;
         csmp_error.notice( WARNING, "findNode:", out.str() );
      }
    return -1;
     
} // end find_node

template long findNode( const Model<1U>&, const Point<1U>&, double64, bool );
template long findNode( const Model<2U>&, const Point<2U>&, double64, bool );
template long findNode( const Model<3U>&, const Point<3U>&, double64, bool );



/// prints sorted global element node numbers in a compact way
template<size_t dim, template<size_t> class CELL>
void printNodes( const CELL<dim>& c )
 {
    set<size_t> nodes;
    for ( size_t i=0U; i<c.Nodes(); i++ ) nodes.insert( c.N(i)->Idx() );
    cout <<" "<< c.Idx() <<": ";
    for ( auto it : nodes ) cout << it <<",";
    cout <<" ";
 }

template void printNodes( const Element<1U>& );
template void printNodes( const Element<2U>& );
template void printNodes( const Element<3U>& );
template void printNodes( const Face<1U>& );
template void printNodes( const Face<2U>& );
template void printNodes( const Face<3U>& );
template void printNodes( const InterFace<1U>& );
template void printNodes( const InterFace<2U>& );
template void printNodes( const InterFace<3U>& );




size_t  renumberElementNodes( vector<Element<1U>*>::iterator first,
                              vector<Element<1U>*>::iterator last )
 {
    assert( first != last );

    set<size_t>  node_numbers;
    size_t       counts(0);
    
    while ( first != last ) {
         for ( vector<Node<1U>*>::iterator
               nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); nit++ ) {
              pair<set<size_t>::iterator,bool>
              sit=node_numbers.insert(counts);
              if ( sit.second == true ) (*nit)->Idx( counts++ );
              else assert( (*nit)->Idx() == (*sit.first) );
           }
         first++;
      }
    
    return node_numbers.size();
    
 } // end renumberElementNodes





size_t  renumberElementNodes( vector<Element<2U>*>::iterator first,
                              vector<Element<2U>*>::iterator last )
 {
    assert( first != last );

    set<size_t>  node_numbers;
    size_t       counts(0);
    
    while ( first != last ) {
         for ( vector<Node<2U>*>::iterator
               nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); nit++ ) {
              pair<set<size_t>::iterator,bool>
              sit=node_numbers.insert(counts);
              if ( sit.second == true ) (*nit)->Idx( counts++ );
              else assert( (*nit)->Idx() == (*sit.first) );
           }
         first++;
      }
    
    return node_numbers.size();
    
 } // end renumberElementNodes




size_t  renumberElementNodes( vector<Element<3U>*>::iterator first,
                              vector<Element<3U>*>::iterator last )
 {
    assert( first != last );

    set<size_t>  node_numbers;
    size_t       counts(0);
    
    while ( first != last ) {
         for ( vector<Node<3U>*>::iterator
               nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); nit++ ) {
              pair<set<size_t>::iterator,bool>
              sit=node_numbers.insert(counts);
              if ( sit.second == true ) (*nit)->Idx( counts++ );
              else assert( (*nit)->Idx() == (*sit.first) );
           }
         first++;
      }
    
    return node_numbers.size();
    
 } // end renumberElementNodes





void printRangeOfVectorOfVectors( const vector<vector<double64> >&  data )
 {
    assert( !data.empty() );
    double64 vmin(data[0][0]), 
              vmax(data[0][0]);
    
    for ( vector<vector<double64> >::const_iterator it=data.begin(); it!=data.end(); it++ )
      for ( vector<double64>::const_iterator  dit=(*it).begin(); dit!=(*it).end(); dit++ ) {
           vmin = std::min( vmin, (*dit) );
           vmax = std::max( vmax, (*dit) );
        }
    
    cout <<"\nprintRangeOfVectorOfVectors: Data range: "<< vmin <<" to "<< vmax << endl;
    
 } // end printRangeOfVectorOfVectors



 
void printRangeOf( const vector<pair<double64,double64> >&  data )
 {
    assert( !data.empty() );
    double64 vmin(data[0].first), 
              vmax(data[0].second);
    
    for ( vector<pair<double64,double64> >::const_iterator it=data.begin(); it!=data.end(); it++ ) {
           vmin = std::min( vmin, (*it).first );
           vmax = std::max( vmax, (*it).second );
        }
    
    cout <<"\nprintRangeOf: Data range: "<< vmin <<" to "<< vmax << endl;
    
 } // end printRangeOf
 



/**
    replaces no-data values of target variable with nearest-neighbor values until there are none left, by default NAN's are no-data values
*/
template<size_t dim>  
void nearestNeighborFill( Model<dim>& model, const char* target_region, const char* variable, double64 no_data_value )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( !model.Database().IsDefined(variable) ) {
        csmp_error.notice( WARNING, "nearestNeighborFill:", variable, "is not defined, nothing could be done.");
        return;
      }
    csmp::Index var_key = model.Database().StorageKey(variable);
    if ( var_key.type != SCALAR ) {
        csmp_error.notice( WARNING, "nearestNeighborFill:", "method currently only handles scalars, nothing could be done.");
        return;
      }
    if ( var_key.place != ELEMENT && var_key.place != NODE ) {
        csmp_error.notice( WARNING, "nearestNeighborFill:", "method currently only handles variables placed on element or nodes, nothing could be done.");
        return;
      }
    if ( !model.ContainsRegion(target_region) ) {
        csmp_error.notice( WARNING, "nearestNeighborFill:", target_region, "is not defined, nothing could be done.");
        return;
      }
    Region<dim>& gref(model.Region(target_region));
  
    // -------------------------------------------------------------------------------------------------------------
    // 1. processing element variables
    // -------------------------------------------------------------------------------------------------------------
    if ( var_key.place == ELEMENT ) {
         // counting no-data values, and memorizing pointers to elements with such values
         set<Element<dim>*> elementsMissingDataValues, filledValues;
         for ( typename vector<Element<dim>*>::const_iterator it=gref.ElementsBegin(); it!=gref.ElementsEnd(); it++ )
           if ( isnan((*it)->Read(var_key)) || fabs(no_data_value-(*it)->Read(var_key)) < numeric_limits<double64>::epsilon() )
             elementsMissingDataValues.insert( (*it) );
          
         if ( elementsMissingDataValues.empty() ) {
              csmp_error.notice( INFO, "nearestNeighborFill:", "all elements have valid data values, nothing was done.");
              return;
           }
         // starting nearest neighbor-fill loop
         else {
              cout <<"\nnearestNeighborFill: detected "<< elementsMissingDataValues.size() <<" elements with no-data values; filling these now.\n";
              do {
                 // looping over the element neighbors, collecting data values for later weighted averaging,
                 // using the inverse of the barycenter to barycenter distance as weighting factor
                 for ( typename set<Element<dim>*>::iterator
                       it=elementsMissingDataValues.begin(); it!=elementsMissingDataValues.end(); it++ )
                   {
                      set<pair<double64,double64> > valuesAndWeights;
                      for ( size_t i=0U; i<(*it)->Neighbors(); i++ )
                        if ( (*it)->Neighbor(i) != NULL )
                          {
                             Point<dim> bctr = (*it)->BaryCenter();
                             double64   nval = (*it)->Neighbor(i)->Read(var_key);
                             // if the neighbor element exists and has a valid variable value, the barycentric distance is determined and stored
                             if ( !isnan(nval) && fabs(no_data_value-nval) > numeric_limits<double64>::epsilon() ) {
                                  double64   distance  = bctr.DistanceTo( (*it)->Neighbor(i)->BaryCenter() );
                                  valuesAndWeights.insert( make_pair(nval,1./distance) );
                               }
                          }
                      // assigned a weighted average to the element if possible
                      if ( !valuesAndWeights.empty() ) {
                           double64  sumOfWeights(0.), sumOfWeightedVals(0.);
                           for ( set<pair<double64,double64> >::const_iterator
                                 vit=valuesAndWeights.begin(); vit!=valuesAndWeights.end(); vit++ ) {
                                sumOfWeightedVals += (*vit).first * (*vit).second;
                                sumOfWeights      += (*vit).second;
                             }
                           (*it)->Store( var_key, makeScalar(PLAIN,sumOfWeightedVals/sumOfWeights) );
                           filledValues.insert( (*it) );
                        }
                   }
                  // removing the successfully filled elements from elementsMissingDataValues
                  for ( typename set<Element<dim>*>::const_iterator it=filledValues.begin(); it!=filledValues.end(); it++ )
                    elementsMissingDataValues.erase( (*it) );
                  filledValues.clear();
                }
              while ( !elementsMissingDataValues.empty() );
          }
        return;
      } // end element variable processing


    // -------------------------------------------------------------------------------------------------------------
    // 2. processing node variables
    // -------------------------------------------------------------------------------------------------------------
    if ( var_key.place == NODE ) {
         // counting no-data values, and memorizing pointers to elements with such values
         set<Node<dim>*> nodesMissingDataValues, filledValues;
         for ( typename vector<Node<dim>*>::const_iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
           if ( isnan((*it)->Read(var_key)) || fabs(no_data_value-(*it)->Read(var_key)) < numeric_limits<double64>::epsilon() )
             nodesMissingDataValues.insert( (*it) );
          
         if ( nodesMissingDataValues.empty() ) {
              csmp_error.notice( INFO, "nearestNeighborFill:", "all nodes have valid data values, nothing was done.");
              return;
           }
         // starting nearest neighbor-fill loop
         else {
              cout <<"\nnearestNeighborFill: detected "<< nodesMissingDataValues.size() <<" nodes with no-data values; filling these now.\n";
              do {
                 // looping over the element neighbors, collecting data values for later weighted averaging,
                 // using the inverse of the barycenter to barycenter distance as weighting factor
                 for ( typename set<Node<dim>*>::iterator
                       it=nodesMissingDataValues.begin(); it!=nodesMissingDataValues.end(); it++ )
                   {
                      set<pair<double64,double64> > valuesAndWeights;
                      for ( size_t i=0U; i<(*it)->Neighbors(); i++ )
                        if ( (*it)->Neighbor(i) != NULL )
                          {
                             Point<dim> nxyz = (*it)->Coordinate();
                             double64   nval = (*it)->Neighbor(i)->Read(var_key);
                             // if the neighbor element exists and has a valid variable value, the barycentric distance is determined and stored
                             if ( !isnan(nval) && fabs(no_data_value-nval) > numeric_limits<double64>::epsilon() ) {
                                  double64   distance  = nxyz.DistanceTo( (*it)->Neighbor(i)->Coordinate() );
                                  valuesAndWeights.insert( make_pair(nval,1./distance) );
                               }
                          }
                      // assigned a weighted average to the element if possible
                      if ( !valuesAndWeights.empty() ) {
                           double64  sumOfWeights(0.), sumOfWeightedVals(0.);
                           for ( set<pair<double64,double64> >::const_iterator
                                 vit=valuesAndWeights.begin(); vit!=valuesAndWeights.end(); vit++ ) {
                                sumOfWeightedVals += (*vit).first * (*vit).second;
                                sumOfWeights      += (*vit).second;
                             }
                           (*it)->Store( var_key, makeScalar(PLAIN,sumOfWeightedVals/sumOfWeights) );
                           filledValues.insert( (*it) );
                        }
                   }
                  // removing the successfully filled elements from elementsMissingDataValues
                  for ( typename set<Node<dim>*>::const_iterator it=filledValues.begin(); it!=filledValues.end(); it++ )
                    nodesMissingDataValues.erase( (*it) );
                  filledValues.clear();
                }
              while ( !nodesMissingDataValues.empty() );
          }
        return;
      } // end node variable processing

 } // end nearestNeighborFill
 
template void nearestNeighborFill( Model<1U>&, const char*, const char*, double64 );
template void nearestNeighborFill( Model<2U>&, const char*, const char*, double64 );
template void nearestNeighborFill( Model<3U>&, const char*, const char*, double64 );



double64 bilinearInterpolate( size_t idx_x, size_t, 
                              const Point<1U>& xy1, 
                              const Point<1U>& xy2,
                              const Point<1U>& coord, 
                              double64 p1, double64 p2, double64, double64 )
{
      // If min-coords. are equivalent to max-coords. the boundary-value average 
      // is assigned.
      if ( xy1 == xy2 ) 
        {
           cout <<"\nbilinearInterpolate: min/max coordinates are identical:"<< endl;
           xy1.Out();
           xy2.Out();
           return (p1+p2) / 2.0;      
        }

      if ( p1 == p2 ) return p1;
    
      // 4. computing interpolation function. Num. Recip. p. 105
      double64 t = (coord[idx_x] - xy1[idx_x] ) / (xy2[idx_x] - xy1[idx_x] );
      
      // 5 bi-linear interpolation
      return (1. - t) * p1 + t * p2;   

} // end bilinearInterpolate





/**

interpolateXY() uses bilinear interpolation to find the value of a
property specified at the corner points of a rectangle, at the coordinates
of a point located inside of this rectangle. 

Since interpolateXY() carries out an interpolation on a planar surface 
in 3D space, it also needs the integer indices which give the axis
in the reference coordinate system. 

@section arguments Input Arguments 

The first two arguments of interpolateXY() specify the coordinate axis
indices of the following Point<dim> arguments which shall be used in 
the interpolation. For instance, for i=0, j=1, the interpolation will be 
carried out in the XY plane. 

The three following VectorVariable<2U> arguments specify the lower left and 
upper right right corners of the rectangle in which the variable value shall 
be interpolated at a point given by the third VectorVariable<2U> argument. 

The last four floating point arguments define the values of the variable
which is to be interpolated. They are the corner points of the rectangle
listed in counter-clockwise fashion (e.g., lower left, lower right, upper
right and upper left corners, respectively). 

@return The result of the interpolation is returned into a double64 type variable. 

@section application Application

Function is used by AssignBoundaryValues(). 

@section messages Messages 

The function will report an error and return the average value of the 
four cornerpoints if their coordinates are identical. 
*/
double64 bilinearInterpolate( size_t idx_x, size_t idx_y, 
                              const Point<2U>& xy1,
                              const Point<2U>& xy2,
                              const Point<2U>& coord, 
                              double64 p1, double64 p2, double64 p3, double64 p4 )
{
      // If min-coords. are equivalent to max-coords. the boundary-value average 
      // is assigned.
      if ( xy1 == xy2 ) {
           cout <<"\nbilinearInterpolate: min/max coordinates are identical:"<< endl;
           xy1.Out();
           xy2.Out();
           return (p1+p2+p3+p4) / 4.0;      
        }

      if ( p1 == p2 && p2 == p3 && p3 == p4 ) return p1;
    
      // 4. computing interpolation functions. Num. Recip. p. 105
      double64 t = (coord[idx_x]-xy1[idx_x] ) / (xy2[idx_x] - xy1[idx_x] );
      double64 u = (coord[idx_y]-xy1[idx_y] ) / (xy2[idx_y] - xy1[idx_y] );
      
      // 5 bi-linear interpolation
      return (1.-t) * (1.-u) * p1 + t * (1.-u) * p2 + t * u * p3 + (1.-t) * u * p4;   

} // end bilinearInterpolate




/**

interpolateXY() uses bilinear interpolation to find the value of a
property specified at the corner points of a rectangle, at the coordinates
of a point located inside of this rectangle. 

Since interpolateXY() carries out an interpolation on a planar surface 
in 3D space, it also needs the integer indices which give the axis
in the reference coordinate system. 

@section arguments Input Arguments 

The first two arguments of interpolateXY() specify the coordinate axis
indices of the following VectorVariable<dim> arguments which shall be used in 
the interpolation. For instance, for i=0, j=1, the interpolation will be 
carried out in the XY plane. 

The three following VectorVariable<dim> arguments specify the lower left and 
upper right right corners of the rectangle in which the variable value shall 
be interpolated at a point given by the third VectorVariable<dim> argument. 

The last four floating point arguments define the values of the variable
which is to be interpolated. They are the corner points of the rectangle
listed in counter-clockwise fashion (e.g., lower left, lower right, upper
right and upper left corners, respectively). 

@return The result of the interpolation is returned into a double64 type variable. 

@section application Application

Function is used by AssignBoundaryValues(). 

@section messages Messages 

The function will report an error and return the average value of the 
four cornerpoints if their coordinates are identical. 
*/
double64 bilinearInterpolate( size_t idx_x, size_t idx_y, 
                              const Point<3U>& xy1,
                              const Point<3U>& xy2,
                              const Point<3U>& coord, 
                              double64 p1, double64 p2, double64 p3, double64 p4 )
{
      // If min-coords. are equivalent to max-coords. the boundary-value average 
      // is assigned.
      if ( xy1 == xy2 ) {
           cout <<"\nbilinearInterpolate: min/max coordinates are identical:"<< endl;
           xy1.Out();
           xy2.Out();
           return (p1+p2+p3+p4) / 4.0;      
        }

      if ( p1 == p2 && p2 == p3 && p3 == p4 ) return p1;
    
      // 4. computing interpolation functions. Num. Recip. p. 105
      double64 t = (coord[idx_x]-xy1[idx_x] ) / (xy2[idx_x] - xy1[idx_x] );
      double64 u = (coord[idx_y]-xy1[idx_y] ) / (xy2[idx_y] - xy1[idx_y] );
      
      // 5 bi-linear interpolation
      return (1.-t) * (1.-u) * p1 + t * (1.-u) * p2 + t * u * p3 + (1.-t) * u * p4;   

} // end bilinearInterpolate


/// interpolate along boundaries of 2D rectangle-shaped model
double64 linearInterpolate( const pair<Point<1U>,double64>& p1, // endpoint1, value1
                            const pair<Point<1U>,double64>& p2, // endpoint2, value2
                            const Point<1U>& pt )                // current point x,y,z
{
    // endmember value range
    double64 dval = p2.second - p1.second;
    
    // distance between endpoints
    double64 dx = p2.first[0] - p1.first[0];
    
    // distance between current point and point 1
    double64 dist = pt[0] - p1.first[0];

    // computing the interpolated value (y = b + mx)
    //     min     normalized distance    gradient
    return p1.second + (dist * dval) / dx;

} // end



/// interpolate along boundaries of 2D rectangle-shaped model
double64 linearInterpolate( const pair<Point<2U>,double64>& p1, // endpoint1, value1
                            const pair<Point<2U>,double64>& p2, // endpoint2, value2
                            const Point<2U>& pt )                // current point x,y,z
{
    // endmember value range
    double64 dval = p2.second - p1.second;
    
    // distance between endpoints
    double64 dx   = sqrt( (p2.first[0]-p1.first[0])*(p2.first[0]-p1.first[0]) +
                           (p2.first[1]-p1.first[1])*(p2.first[1]-p1.first[1]) );
    
    // distance between current point and point 1
    double64 dist = sqrt( (pt[0]-p1.first[0])*(pt[0]-p1.first[0]) +
                           (pt[1]-p1.first[1])*(pt[1]-p1.first[1]) );

    // computing the interpolated value (y = b + mx)
    //     min     normalized distance    gradient
    return p1.second + (dist * dval) / dx;

} // end


/// interpolate along boundaries of 2D rectangle-shaped model
double64 linearInterpolate( const pair<Point<3U>,double64>& p1, // endpoint1, value1
                             const pair<Point<3U>,double64>& p2, // endpoint2, value2
                             const Point<3U>& pt )                // current point x,y,z
{
    // endmember value range
    double64 dval = p2.second - p1.second;
    
    // distance between endpoints
    double64 dx   = sqrt( (p2.first[0]-p1.first[0])*(p2.first[0]-p1.first[0]) +
                           (p2.first[1]-p1.first[1])*(p2.first[1]-p1.first[1]) +
                           (p2.first[2]-p1.first[2])*(p2.first[2]-p1.first[2]) );
    
    // distance between current point and point 1
    double64 dist = sqrt( (pt[0]-p1.first[0])*(pt[0]-p1.first[0]) +
                           (pt[1]-p1.first[1])*(pt[1]-p1.first[1]) +
                           (pt[2]-p1.first[2])*(pt[2]-p1.first[2]) );

    // computing the interpolated value (y = b + mx)
    //     min     normalized distance    gradient
    return p1.second + (dist * dval) / dx;

} // end


template<size_t dim,class VarType>
void extrapolateElementToNodeProperty( csmp::VSet<dim>&             vset,
                                       const std::vector<VarType>& elmnt_data,
                                       std::vector<VarType>&       nodal_data )
{
    /// mesh specs
    const size_t num_elements( vset.Elements() );
    const size_t num_nodes( vset.Vertices() );
    size_t eid;
    size_t nid;
    size_t enodes;

    /// fem data
    int fem_type;
    csmp::Point<dim> pt;
    double volume;
    size_t lnid;

    /// apply volume averaging to the element property
    VarType var;
    var = 0.0;
    nodal_data.resize( num_nodes, var );
    std::vector<double> total_volume( num_nodes, 0.0 );
    for ( eid = 0; eid<num_elements; ++eid )
    {
        enodes = vset.PlistSize( eid );
        if( enodes > 1 )
        {
            /// get element data
            var = elmnt_data[ eid ];
            /// calculate volume of element
            fem_type = ( vset.HybridElementTypeMesh() ? vset.ElementType( eid ) : vset.ElementType( 0 ) );
            switch( fem_type )
            {
            case csmp::ISOPARAMETRIC_LINEAR_HEXAHEDRON:
                {
                    csmp::IsoparametricLinearHexahedron     fem;
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            case csmp::ISOPARAMETRIC_LINEAR_PYRAMID:
                {
                    csmp::IsoparametricLinearPyramid        fem;
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            case csmp::ISOPARAMETRIC_LINEAR_PRISM:
                {
                    csmp::IsoparametricLinearPrism fem;
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            case csmp::ISOPARAMETRIC_LINEAR_TETRAHEDRON:
                {
                    csmp::IsoparametricLinearTetrahedron fem;
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            case csmp::ISOPARAMETRIC_LINEAR_QUADRILATERAL:
                {
                    csmp::IsoparametricLinearQuadrilateral  fem(3);
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            case csmp::ISOPARAMETRIC_LINEAR_TRIANGLE:
                {
                    csmp::IsoparametricLinearTriangle fem(3);
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            case csmp::ISOPARAMETRIC_LINEAR_BAR:
                {
                    csmp::IsoparametricLinearLineElement    fem(3);
                    fem.XY.Resize( enodes, dim );
                    lnid = 0;
                    for ( std::vector<size_t>::iterator
                          nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit, ++lnid )
                    {
                        nid = (*nit);
                        pt[0] = vset.Px( nid );
                        if( dim > 1U )  pt[1] = vset.Py( nid );
                        if( dim == 3U ) pt[2] = vset.Pz( nid );
                        fem.XY.AssignRow( lnid, pt );
                    }
                    volume = fem.Volume();
                }
                break;
            default:
                volume = 1.0;
                break;
            }
            var *= volume;
            /// summ up data
            for ( std::vector<size_t>::iterator
                  nit = vset.PlistBegin( eid ); nit != vset.PlistEnd( eid ); ++nit )
            {
                nid = (*nit);
                nodal_data[ nid ]   += var;
                total_volume[ nid ] += volume;
            }
        }
    }
    /// calculating the nodal averages
    for ( nid = 0; nid < num_nodes; ++nid )
         nodal_data[ nid ] /= total_volume[ nid ];
}

template void extrapolateElementToNodeProperty(csmp::VSet<1U>&,const std::vector<csmp::ScalarVariable>&,std::vector<csmp::ScalarVariable>&);
template void extrapolateElementToNodeProperty(csmp::VSet<1U>&,const std::vector<csmp::VectorVariable<1U> >&,std::vector<csmp::VectorVariable<1U> >&);
template void extrapolateElementToNodeProperty(csmp::VSet<1U>&,const std::vector<csmp::TensorVariable<1U> >&,std::vector<csmp::TensorVariable<1U> >&);
template void extrapolateElementToNodeProperty(csmp::VSet<1U>&,const std::vector<csmp::ArrayVariable>&,std::vector<csmp::ArrayVariable>&);
template void extrapolateElementToNodeProperty(csmp::VSet<1U>&,const std::vector<csmp::FlaggedArrayVariable>&,std::vector<csmp::FlaggedArrayVariable>&);

template void extrapolateElementToNodeProperty(csmp::VSet<2U>&,const std::vector<csmp::ScalarVariable>&,std::vector<csmp::ScalarVariable>&);
template void extrapolateElementToNodeProperty(csmp::VSet<2U>&,const std::vector<csmp::VectorVariable<2U> >&,std::vector<csmp::VectorVariable<2U> >&);
template void extrapolateElementToNodeProperty(csmp::VSet<2U>&,const std::vector<csmp::TensorVariable<2U> >&,std::vector<csmp::TensorVariable<2U> >&);
template void extrapolateElementToNodeProperty(csmp::VSet<2U>&,const std::vector<csmp::ArrayVariable>&,std::vector<csmp::ArrayVariable>&);
template void extrapolateElementToNodeProperty(csmp::VSet<2U>&,const std::vector<csmp::FlaggedArrayVariable>&,std::vector<csmp::FlaggedArrayVariable>&);

template void extrapolateElementToNodeProperty(csmp::VSet<3U>&,const std::vector<csmp::ScalarVariable>&,std::vector<csmp::ScalarVariable>&);
template void extrapolateElementToNodeProperty(csmp::VSet<3U>&,const std::vector<csmp::VectorVariable<3U> >&,std::vector<csmp::VectorVariable<3U> >&);
template void extrapolateElementToNodeProperty(csmp::VSet<3U>&,const std::vector<csmp::TensorVariable<3U> >&,std::vector<csmp::TensorVariable<3U> >&);
template void extrapolateElementToNodeProperty(csmp::VSet<3U>&,const std::vector<csmp::ArrayVariable>&,std::vector<csmp::ArrayVariable>&);
template void extrapolateElementToNodeProperty(csmp::VSet<3U>&,const std::vector<csmp::FlaggedArrayVariable>&,std::vector<csmp::FlaggedArrayVariable>&);




 
  
  
  

/**
   @todo document this function
*/
#if defined _MSC_VER || defined __MINGW32__
const char * strp_weekdays[] =
{ "sunday", "monday", "tuesday", "wednesday", "thursday", "friday", "saturday" };
const char * strp_monthnames[] =
{ "january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december" };
bool strp_atoi(const char * & s, int & result, int low, int high, int offset)
{
	bool worked = false;
	char * end;
	unsigned long num = strtoul(s, &end, 10);
	if (num >= (unsigned long)low && num <= (unsigned long)high)
	{
		result = (int)(num + offset);
		s = end;
		worked = true;
	}
	return worked;
}
char * strptime(const char *s, const char *format, struct tm *tm)
{
	bool working = true;
	while (working && *format && *s)
	{
		switch (*format)
		{
		case '%':
		{
					++format;
					switch (*format)
					{
					case 'a':
					case 'A': // weekday name
						tm->tm_wday = -1;
						working = false;
						for (size_t i = 0; i < 7; ++i)
						{
							size_t len = strlen(strp_weekdays[i]);
							if (!strnicmp(strp_weekdays[i], s, len))
							{
								tm->tm_wday = i;
								s += len;
								working = true;
								break;
							}
							else if (!strnicmp(strp_weekdays[i], s, 3))
							{
								tm->tm_wday = i;
								s += 3;
								working = true;
								break;
							}
						}
						break;
					case 'b':
					case 'B':
					case 'h': // month name
						tm->tm_mon = -1;
						working = false;
						for (size_t i = 0; i < 12; ++i)
						{
							size_t len = strlen(strp_monthnames[i]);
							if (!strnicmp(strp_monthnames[i], s, len))
							{
								tm->tm_mon = i;
								s += len;
								working = true;
								break;
							}
							else if (!strnicmp(strp_monthnames[i], s, 3))
							{
								tm->tm_mon = i;
								s += 3;
								working = true;
								break;
							}
						}
						break;
					case 'd':
					case 'e': // day of month number
						working = strp_atoi(s, tm->tm_mday, 1, 31, 0);
						break;
					case 'D': // %m/%d/%y
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_mon, 1, 12, -1);
								  if (working && *s == '/')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_mday, 1, 31, 0);
									  if (working && *s == '/')
									  {
										  ++s;
										  working = strp_atoi(s, tm->tm_year, 0, 99, 0);
										  if (working && tm->tm_year < 69)
											  tm->tm_year += 100;
									  }
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'H': // hour
						working = strp_atoi(s, tm->tm_hour, 0, 23, 0);
						break;
					case 'I': // hour 12-hour clock
						working = strp_atoi(s, tm->tm_hour, 1, 12, 0);
						break;
					case 'j': // day number of year
						working = strp_atoi(s, tm->tm_yday, 1, 366, -1);
						break;
					case 'm': // month number
						working = strp_atoi(s, tm->tm_mon, 1, 12, -1);
						break;
					case 'M': // minute
						working = strp_atoi(s, tm->tm_min, 0, 59, 0);
						break;
					case 'n': // arbitrary whitespace
					case 't':
						while (isspace((int)*s))
							++s;
						break;
					case 'p': // am / pm
						if (!strnicmp(s, "am", 2))
						{ // the hour will be 1 -> 12 maps to 12 am, 1 am .. 11 am, 12 noon 12 pm .. 11 pm
							if (tm->tm_hour == 12) // 12 am == 00 hours
								tm->tm_hour = 0;
						}
						else if (!strnicmp(s, "pm", 2))
						{
							if (tm->tm_hour < 12) // 12 pm == 12 hours
								tm->tm_hour += 12; // 1 pm -> 13 hours, 11 pm -> 23 hours
						}
						else
							working = false;
						break;
					case 'r': // 12 hour clock %I:%M:%S %p
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_hour, 1, 12, 0);
								  if (working && *s == ':')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_min, 0, 59, 0);
									  if (working && *s == ':')
									  {
										  ++s;
										  working = strp_atoi(s, tm->tm_sec, 0, 60, 0);
										  if (working && isspace((int)*s))
										  {
											  ++s;
											  while (isspace((int)*s))
												  ++s;
											  if (!strnicmp(s, "am", 2))
											  { // the hour will be 1 -> 12 maps to 12 am, 1 am .. 11 am, 12 noon 12 pm .. 11 pm
												  if (tm->tm_hour == 12) // 12 am == 00 hours
													  tm->tm_hour = 0;
											  }
											  else if (!strnicmp(s, "pm", 2))
											  {
												  if (tm->tm_hour < 12) // 12 pm == 12 hours
													  tm->tm_hour += 12; // 1 pm -> 13 hours, 11 pm -> 23 hours
											  }
											  else
												  working = false;
										  }
									  }
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'R': // %H:%M
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_hour, 0, 23, 0);
								  if (working && *s == ':')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_min, 0, 59, 0);
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'S': // seconds
						working = strp_atoi(s, tm->tm_sec, 0, 60, 0);
						break;
					case 'T': // %H:%M:%S
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_hour, 0, 23, 0);
								  if (working && *s == ':')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_min, 0, 59, 0);
									  if (working && *s == ':')
									  {
										  ++s;
										  working = strp_atoi(s, tm->tm_sec, 0, 60, 0);
									  }
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'w': // weekday number 0->6 sunday->saturday
						working = strp_atoi(s, tm->tm_wday, 0, 6, 0);
						break;
					case 'Y': // year
						working = strp_atoi(s, tm->tm_year, 1900, 65535, -1900);
						break;
					case 'y': // 2-digit year
						working = strp_atoi(s, tm->tm_year, 0, 99, 0);
						if (working && tm->tm_year < 69)
							tm->tm_year += 100;
						break;
					case '%': // escaped
						if (*s != '%')
							working = false;
						++s;
						break;
					default:
						working = false;
					}
		}
			break;
		case ' ':
		case '\t':
		case '\r':
		case '\n':
		case '\f':
		case '\v':
			// zero or more whitespaces:
			while (isspace((int)*s))
				++s;
			break;
		default:
			// match character
			if (*s != *format)
				working = false;
			else
				++s;
			break;
		}
		++format;
	}
	return (working ? (char *)s : 0);
}
#endif // _MSC_VER



/**
 @brief Find element which contains a given point.

 Note that this only searches volumetric elements. It is not recommended
 that you use this function if you need to search for many points.

 @author  A.J. Bromage
 @date    11/04/2018

 @param [in] region  region to search
 @param [in] query   query point

 @return  the element which contains the point, or NULL if no element does
 */
Element<3u>* pointInVolumeElement( const Region<3u>& region, const Point<3u>& query )
    {
      std::vector<size_t> fnids;
      fnids.reserve(4);

      auto eend = region.ElementsEnd();
      for (auto eit = region.ElementsBegin(); eit != eend; ++eit) {

        // 1. Volume elements only
 
        if (!(*eit)->IsVolumeElement()) {
          continue;
        }

        // 2. Test against axis-aligned bounding box

        double64 minx = +std::numeric_limits<double64>::max();
        double64 miny = +std::numeric_limits<double64>::max();
        double64 minz = +std::numeric_limits<double64>::max();
        double64 maxx = -std::numeric_limits<double64>::max();
        double64 maxy = -std::numeric_limits<double64>::max();
        double64 maxz = -std::numeric_limits<double64>::max();
        const size_t iNrNodes = (*eit)->Nodes();
        for ( size_t iNode = 0; iNode < iNrNodes; ++iNode ) {
          auto n = (*eit)->N(iNode);
          minx = std::min(minx, n->x());
          maxx = std::max(maxx, n->x());
          miny = std::min(miny, n->y());
          maxy = std::max(maxy, n->y());
          minz = std::min(minz, n->z());
          maxz = std::max(maxz, n->z());
        }

        if (query[0] < minx || query[0] > maxx
            || query[1] < miny || query[1] > maxy
            || query[2] < minz || query[2] > maxz) {
          continue;
        }

        // 3. Test against all faces

        bool reject = false;
        auto fe = (*eit)->FE();
        const size_t iNrFaces = (*eit)->Faces();
        for (size_t iFace = 0; iFace < iNrFaces && !reject; ++iFace) {
          fe->NodesOfFace(iFace, fnids);
          const size_t iNrFacePts = fnids.size();
          for (size_t iFacePt = 0; iFacePt < iNrFacePts; iFacePt += 2) {
            auto p0 = (*eit)->N(fnids[(iFacePt+0) % iNrFacePts])->Coordinate();
            auto p1 = (*eit)->N(fnids[(iFacePt+1) % iNrFacePts])->Coordinate();
            auto p2 = (*eit)->N(fnids[(iFacePt+2) % iNrFacePts])->Coordinate();

            auto normal = crossProduct(p2-p0, p1-p0);
            normal.NormalizeLengthTo(1.0f);
            const double64 queryDotNormal = dotProduct(query, normal);
            const double64 p0DotNormal = dotProduct(p0, normal);

            if (queryDotNormal < p0DotNormal) {
              reject = true;
              break;
            }
          }
        }
        if (!reject) {
          return *eit;
        }
        fnids.clear();
      }

      return 0;
    }



/// Utility that tokenises string into substrings using the supplied delimiter(s).
std::vector<std::string> splitString( std::string str, char delimiter )
{
  size_t pos = 0U;
  string token, s = str;
  vector<string> items;
  while ( (pos = s.find( delimiter )) != std::string::npos ) {
    token = s.substr( 0, pos );
    items.push_back( token );
    s.erase( 0, pos + 1 );
  }
  if ( !s.empty() ) items.push_back( s );
  return items;
}




/**
      reads vector<vector> from filestream where the elements of the vector are sequential 
*/   
template<typename T>   
void readVectorOfVectors( ifstream& ifs, size_t total_items, size_t entries_per_vector, deque<vector<T> >& file_records )
 {
    file_records.clear();
    // not in deque: file_records.reserve( total_items / entries_per_vector );
    // reading plist
    size_t item=0U; 
    while ( item < total_items )
      {
         // node ID's in file range 0...nodes-1
         vector<T> data;
         data.reserve( entries_per_vector );
         int32 id;
         for ( size_t i=0; i<entries_per_vector; ++i ) {
              ifs >> id;
              assert( id >= 0 && id << total_items ); // assumption that there are not more nodes that elements*nodes_per_element
              data.push_back( id );
              item++;
           }
         file_records.emplace_back( data );
      }

    if ( item != total_items )
         throw csmp::Exception( ERROR, "readVectorOfVectors", 
                               "File record of vector<vector<typename>> was not correctly read" );

 } // end read inlined vector

template void readVectorOfVectors( ifstream&, size_t, size_t, deque<vector<size_t> >& );






} // end namespace csmp
