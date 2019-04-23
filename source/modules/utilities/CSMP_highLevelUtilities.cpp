#include "CSMP_highLevelUtilities.h"
#include "compareFloats.h"
#include "Box.h"
#include "VSet.h"
#include "Node.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Boundary.h"
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



/**

Method obtains the range of values for the physical variable that it
is prompted for by the user. For vector and tensor variables, the
length and the minimum/maximum eigenvalues are returned (check whether
the latter actually happens).  

With the optional boolean argument (default=true), the user can determin
the return value. As the default, the maximum obtained value is returned;
else the minimum.  

@section arguments Input Arguments 

The model to be examined, the name of the target physical variable, and 
targeted return value (true->maximum, false->minimum of target variable).

The second version of this method also takes an I/O handler as argument
in order to log the calculated values to file etc. 

@return Either the maximum (default) or the minimum value of the target variable.
The result is printed to the screen.

*/
template<size_t  dim>
double64 printRangeOfVariable( const Model<dim>& sg, 
                               const char* var, bool max_or_min )
 {
     double64 pmin, pmax;
     sg.MinMaxOf( var, pmin, pmax );
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< sg.Database().Unit(var) <<"]: '";
     cout << var <<"': "<< pmin <<" to "<< pmax << endl;
          
     if ( !max_or_min ) return pmin;
     return pmax;
 }



template<size_t  dim>
double64 printRangeOfVariable( const Model<dim>& sg,
                               Standard_IO_Handler& io,
                               const char* var, bool max_or_min )
 {
     double64  pmin, pmax;
     sg.MinMaxOf( var, pmin, pmax );
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< sg.Database().Unit(var) <<"]: '";
     cout << var <<"': "<< pmin <<" to "<< pmax << endl;
     
     // recording the measured variable value range at given timestep
     double64& model_time( ModelTime::Instance().modelTime );
     char   info[100];
     sprintf( info, "%lf", model_time );
     string var_info(info);
     var_info += " secs, range of'";
     var_info += var;
     var_info += "' [";
     var_info += sg.Database().Unit(var);
     var_info += "]: ";
     sprintf( info, "%lf", pmin );
     var_info += info;
     var_info += " to ";
     sprintf( info, "%lf", pmax );
     var_info += info;
     
     io.RecordInformation( var_info ); 
     
     if ( !max_or_min ) return pmin;
     return pmax;
 }




/// as above, but for individual model regions
template<size_t  dim>
double64 printRangeOfVariable( const Model<dim>& sg, 
                               const char* group, const char* var, bool max_or_min )
 {
     double64  pmin, pmax;
     const PropertyDatabase<dim>& p_ref = sg.Database();
     const PLACEMENT place = sg.Database().Placement(var);
     
     if ( sg.ContainsRegion(group) && !faceVariable(place) && !interFaceVariable(place) )
       sg.Region( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsBoundary(group) ) sg.Boundary( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsSplitBoundary(group) ) sg.SplitBoundary( group ).MinMaxOf( var, pmin, pmax );
     else {
          cerr <<"\nprintRangeOfVariable: '"<< group <<"' does not exist."<< endl;
          return std::numeric_limits<double64>::quiet_NaN();
       }
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"' in subdomain of model '"<< group <<"': "<< pmin <<" to "<< pmax << endl;
          
     if ( !max_or_min ) return pmin;
     return pmax;
 }




template<size_t  dim>
double64 printRangeOfVariable( const Model<dim>& sg, 
                                Standard_IO_Handler& io, 
                                const char* group,
                                const char* var, bool max_or_min )
 {
     double64& model_time( ModelTime::Instance().modelTime );
     double64         pmin, pmax;
     const PropertyDatabase<dim>& p_ref = sg.Database();
     const PLACEMENT place = sg.Database().Placement(var);

     if ( sg.ContainsRegion(group) && !faceVariable(place) && !interFaceVariable(place) )
       sg.Region( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsBoundary(group) ) sg.Boundary( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsSplitBoundary(group) ) sg.SplitBoundary( group ).MinMaxOf( var, pmin, pmax );
     else {
          cerr <<"\nprintRangeOfVariable: '"<< group <<"' does not exist."<< endl;
          return std::numeric_limits<double64>::quiet_NaN();
       }
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"': "<< pmin <<" to "<< pmax <<" in subdomain of model '"<< group <<"'"<< endl;
     
     // recording the measured variable value range at given timestep
     char info[100];
     sprintf( info, "%lf", model_time );
     string var_info(info);
     var_info += info;
     var_info += ", region: ";
     var_info += group;
     var_info += ", secs, range of '";
     var_info += var;
     var_info += "' [";
     var_info += p_ref.Unit(var);
     var_info += "]: ";
     sprintf( info, "%lf", pmin );
     var_info += info;
     var_info += " to ";
     sprintf( info, "%lf", pmax );
     var_info += info;
     
     io.RecordInformation( var_info ); 
     
     if ( !max_or_min ) return pmin;
     return pmax;
 }






/**

Method prints the physical dimensions of the model and returns either
the maximum or the intermediate axis, depending on the value of its 
second argument.  

returns intermediate (true) or maximum (false) model dimensions.  

@section arguments Input Arguments 

The current model and a boolean flag. For 'true' the intermediate axis
is returned, if 'false' the long axis is returned. 

@return The intermediate or long axis of the current model.  
*/
template<size_t  dim>
double64  printModelDimensions( const Model<dim>& sg, bool intermed_or_max )
 {
    Point<dim> xyz_min, xyz_max; 
    sg.MinMaxCoordinates( xyz_min, xyz_max );    
    cout <<"\nprintModelDimensions: Dimensions of model (meters): "<< endl;
    cout <<"xmin, xmax (horizontal right):    "<< xyz_min[0] <<" "<< xyz_max[0] << endl;
    if ( dim != 1U ) cout <<"ymin, ymax (vertical upward):     "<< xyz_min[1] <<" "<< xyz_max[1] << endl;
    if ( dim == 3U ) cout <<"zmin, zmax (horizontal to front): "<< xyz_min[2] <<" "<< xyz_max[2] << endl << endl;

    set<double64,greater<double64> >  axis;
    axis.insert( xyz_max[0] - xyz_min[0] );
    if ( dim != 1U ) axis.insert( xyz_max[1] - xyz_min[1] );
    if ( dim == 3U ) axis.insert( xyz_max[2] - xyz_min[2] );
    
    set<double64,greater<double64> >::const_iterator  it = axis.begin();
    
    if ( !intermed_or_max ) return *it;
    
    if ( axis.size() >= 2U ) it++;
    
    return *it;

 } // end printModelDimensions




/** 
    calculates the center of gravity of the model by averaging
    the barycenter locations of all highest-dimensional elements.
*/
template<size_t  dim>
Point<dim>  centerOfGravity( const Model<dim>& model )
 {
    const Region<dim>& mref(model.Region("Model"));
    typename vector<Element<dim>*>::const_iterator it(mref.ElementsBegin());
    Point<dim>  center((*it)->BaryCenter());
    double64    counter(0.);
    it++;
   
    while( it != mref.ElementsEnd() ) {
         if ( dim == 3U ) {
               if ( (*it)->FE()->IsVolumeElement() ) {
                    center += (*it)->BaryCenter();
                    counter += 1.;
                 }
            }
         else if ( dim == 2U ) {
               if ( (*it)->FE()->IsSurfaceElement() ) {
                    center += (*it)->BaryCenter();
                    counter += 1.;
                 }
            }
         else /* 1D */ {
                    center += (*it)->BaryCenter();
                    counter += 1.;
            }
         it++;
      }
    center /= counter;
    return center;

 } // end CenterOfGravity

template Point<1U>  centerOfGravity( const Model<1U>& );
template Point<2U>  centerOfGravity( const Model<2U>& );
template Point<3U>  centerOfGravity( const Model<3U>& );



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
 



// template instantiations
template
double64  printModelDimensions( const Model<1U>& sg, bool intermed_or_max ); 

template
double64  printRangeOfVariable( const Model<1U>& sg, 
                                const char* var, bool max_or_min );
template
double64  printRangeOfVariable( const Model<1U>& sg, 
	                              Standard_IO_Handler& io, const char* var,
	                              bool max_instead_of_min );
template
double64  printRangeOfVariable( const Model<1U>& sg, 
                                const char* group, const char* var, bool max_or_min );
template
double64  printRangeOfVariable( const Model<1U>& sg, 
                                Standard_IO_Handler& io,
                                const char* group, const char* var,
                                bool max_instead_of_min );

template
double64  printModelDimensions( const Model<2U>& sg, bool intermed_or_max ); 

template
double64  printRangeOfVariable( const Model<2U>& sg, 
                                const char* var, bool max_or_min );
template
double64  printRangeOfVariable( const Model<2U>& sg, 
                                Standard_IO_Handler& io, const char* var,
                                bool max_instead_of_min );
template
double64  printRangeOfVariable( const Model<2U>& sg, 
                                const char* group, const char* var, bool max_or_min );
template
double64  printRangeOfVariable( const Model<2U>& sg, 
                                Standard_IO_Handler& io,
                                const char* group, const char* var,
                                bool max_instead_of_min );

template
double64  printModelDimensions( const Model<3U>& sg, bool intermed_or_max ); 

template
double64  printRangeOfVariable( const Model<3U>& sg, 
                                const char* var, bool max_or_min );
template
double64  printRangeOfVariable( const Model<3U>& sg, 
                                Standard_IO_Handler& io, const char* var,
                                bool max_instead_of_min );
template
double64  printRangeOfVariable( const Model<3U>& sg, 
                                const char* group, const char* var, bool max_or_min );
template
double64  printRangeOfVariable( const Model<3U>& sg, 
                                Standard_IO_Handler& io,
                                const char* group, const char* var,
                                bool max_instead_of_min );


/**
    Checks whether region of interest contains elements of the type of interest.
*/
template<size_t dim>  
bool containsElementsOfTtype( const Region<dim>& gref, ELEMENT_DIMENSION dimension )
 {
    for ( typename vector<Element<dim>*>::const_iterator
          it=gref.ElementsBegin(); it!=gref.ElementsEnd(); it++ )
      if ( parseFiniteElementDimension( (*it)->FE_Type() ) == dimension )
        return true; 
      
    return false;
      
 } // end containsElementsOfTtype

template bool containsElementsOfTtype<1U>( const Region<1>&, ELEMENT_DIMENSION );
template bool containsElementsOfTtype<2U>( const Region<2>&, ELEMENT_DIMENSION );
template bool containsElementsOfTtype<3U>( const Region<3>&, ELEMENT_DIMENSION );









/**
     Smoothes scalar element variable by extrapolating it to the nodes and back-interpolating it to barycenters
     Apart from the mname of element variable to be smoothed, the name of the temporary node variable needs to be specified.
     Uses the dummy variable 'dummy node' to store the interim result
*/   
template<size_t dim>  
void smoothElementVariable( Model<dim>& model, const char* region, const char* element_var, const char* temp_node_var, size_t n_smoothing_cycles )
 {
    csmp::Index eprop_key = model.Database().StorageKey(element_var);
    csmp::Index nprop_key = model.Database().StorageKey(temp_node_var);
    assert( nprop_key.place == NODE );

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    if ( eprop_key.place != ELEMENT ) {
         csmp_error.notice( ERROR, "smoothElementVariable", "Smoothed variable must be placed on the element" );
         return;
      }
    if ( nprop_key.place != NODE ) {
         csmp_error.notice( ERROR, "smoothElementVariable", "Temporary variable must be placed on the node" );
         return;
      }
    if ( nprop_key.type != eprop_key.type ) {
         csmp_error.notice( ERROR, "smoothElementVariable", "Smoothed and temporary variable must have the same type" );
         return;
      }
    if ( n_smoothing_cycles == 0 ) {
         csmp_error.notice( WARNING, "smoothElementVariable", "smoothing cycles=0; nothing was done" );
         return;
      }
   
    // smoothing
    Region<dim> ref = model.Region(region);
   
    for ( size_t i=0U; i<n_smoothing_cycles; i++ ) {
         ref.ExtrapolateElementToNodeProperty( element_var, temp_node_var );
         ref.InterpolateNodeToElementProperty( temp_node_var, element_var );
      }

 } // end smoothElementVariable

template void smoothElementVariable( Model<1U>&, const char*, const char*, const char*, size_t );
template void smoothElementVariable( Model<2U>&, const char*, const char*, const char*, size_t );
template void smoothElementVariable( Model<3U>&, const char*, const char*, const char*, size_t );




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

Randomly perturbs the values of a scalar target property by subtracting an
amount which varies between minus zero and the specified percentage of the 
original maximum value of the target property. 

@section arguments Input Arguments 

RandomPerturb() requires the name of the property which shall be perturbed 
and the percentage of the original maximum value of the property by which 
the property shall be perturbed, in order to operate. 

@section application Application

Processes which are critically dependent on initial conditions can profit 
from a 'noisy' input signal when one tries to simulate natural behaviour. 

@section messages Messages 

RandomPerturb() only handles scalar variables and it will therefore report 
an error and return without executing when one tries to perturb a vector or 
tensor variable. 

*/
template<size_t dim>
void randomPerturb( random_generator& rng, Model<dim>& sg, const char* prop, double64 by_percent_of_max_value )
 {
    csmp::Index prop_key = sg.Database().StorageKey(prop);
    Region<dim>&  sgroup(sg.Region("Model"));
    
    if ( prop_key.type != SCALAR )
      throw csmp::Exception( ERROR, "Model::RandomPerturb", 
                                     "Can only perturb scalar values so far" ); 		       

    double64 dmin, dmax;
    ScalarVariable  sc;
    sgroup.MinMaxOf( prop, dmin, dmax );

    std::uniform_real_distribution<double>
        rngen(0, by_percent_of_max_value * dmax * 0.01);

    switch( prop_key.place )
      {
         case NODE:
              for ( typename vector<Node<dim>*>::iterator
                    nit=sgroup.NodesBegin(); nit!=sgroup.NodesEnd(); nit++ )
                {
                   (*nit)->Read( prop_key, sc );
                   sc -= rngen(rng);
                   (*nit)->Store( prop_key, sc );
                }
           break;
         case ELEMENT_INTEGRATION_POINT:
              for ( typename vector<Element<dim>*>::iterator
                    eit=sgroup.ElementsBegin(); eit!=sgroup.ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                {
                   (*eit)->Read( i, prop_key, sc );
                   sc -= rngen(rng);
                   (*eit)->Store( i, prop_key, sc );
                }
           break;
         case ELEMENT:
              for ( typename vector<Element<dim>*>::iterator
                    eit=sgroup.ElementsBegin(); eit!=sgroup.ElementsEnd(); eit++ )
                {
                   (*eit)->Read( prop_key, sc );
                   sc -= rngen(rng);
                   (*eit)->Store( prop_key, sc );
                }
           break;
         default:
           cout <<"\nrandomPerturb: property placement not handled."<< endl;
      }
      
 } // end RandomPerturb

template void randomPerturb( random_generator& rng, Model<1U>&, const char*, double64 );
template void randomPerturb( random_generator& rng, Model<2U>&, const char*, double64 );
template void randomPerturb( random_generator& rng, Model<3U>&, const char*, double64 );






/**
    convert the flag(s) of a variable into integer values stored in its number part
    
    @attention works only for scalars and basic property placements.
    
    @author SKM 21/5/2014
*/
template<size_t dim> 
void flagToNumber( Model<dim>& model, const char* variable )
 {
    csmp::Region<dim>&  mref(model.Region("Model"));
    csmp::Index  prop_key = model.Database().StorageKey(variable);
   
    if ( prop_key.type != SCALAR )
      throw csmp::Exception( ERROR, "flagToNumber:", "method has not been implemented yet" );

    switch( prop_key.place )
      {
         case NODE:
              for ( typename vector<Node<dim>*>::iterator
                    nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ )
                {
                   // overwrites variable value with integer value of its flag enum
                   double64 value = static_cast<double64>( (*nit)->Status(prop_key) );
                   (*nit)->Store( prop_key, makeScalar( (*nit)->Status(prop_key), value ) );
                }
           break;
         case ELEMENT_INTEGRATION_POINT:
              for ( typename vector<Element<dim>*>::iterator
                    eit=mref.ElementsBegin(); eit!=mref.ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                {
                   double64 value = static_cast<double64>( (*eit)->Status(prop_key) );
                   (*eit)->Store( prop_key, makeScalar( (*eit)->Status(prop_key), value ) );
                }
           break;
         case ELEMENT:
              for ( typename vector<Element<dim>*>::iterator
                    eit=mref.ElementsBegin(); eit!=mref.ElementsEnd(); eit++ )
                {
                   double64 value = static_cast<double64>( (*eit)->Status(prop_key) );
                   (*eit)->Store( prop_key, makeScalar( (*eit)->Status(prop_key), value ) );
                }
           break;
         default:
           cout <<"\nflagToNumber: property placement not handled."<< endl;
      }
   
 } // end flagToNumber

template void flagToNumber( Model<1U>&, const char* );
template void flagToNumber( Model<2U>&, const char* );
template void flagToNumber( Model<3U>&, const char* );




/**
    convert the flag(s) of first variable into double values stored in the second variable
    
    @attention works only for node-property placement.
    
    @author SKM 8/12/2016
*/
template<size_t dim> 
void flagToNumber( Model<dim>& model, const char* flag_variable, const char* value_variable )
 {
    csmp::Region<dim>&  mref(model.Region("Model"));
    csmp::Index  flag_key = model.Database().StorageKey(flag_variable);  // input
    csmp::Index  prop_key = model.Database().StorageKey(value_variable); // output
 
    if ( flag_key.type != prop_key.type )
      throw csmp::Exception( ERROR, "flagToNumber:", "flag and value variables must be of the same type." );

    if ( flag_key.place != prop_key.place )
      throw csmp::Exception( ERROR, "flagToNumber:", "flag and value variables must have the same placement." );

    if ( flag_key.type != SCALAR and flag_key.type != VECTOR )
      throw csmp::Exception( ERROR, "flagToNumber:", "method handles only scalar and vector variables." );
  
    switch( prop_key.place )
      {
         case NODE:
              if ( flag_key.type == SCALAR ) {
                  for ( typename vector<Node<dim>*>::iterator
                        nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ )
                    {
                       // retrieves status of the flag variable
                       double64 value = static_cast<double64>( (*nit)->Status(flag_key) );
                       // overwrites value of value variable with integer value of its flag enum
                       (*nit)->Store( prop_key, makeScalar( (*nit)->Status(prop_key), value ) );
                    }
                }
              else if ( flag_key.type == VECTOR ) {
                  VectorVariable<dim> vc;
                  for ( typename vector<Node<dim>*>::iterator
                        nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ )
                    {
                       (*nit)->Read( flag_key, vc );
                       for ( size_t i=0U; i<dim; ++ i )
                         vc(i) = static_cast<double64>( vc.Flag(i) );
                       (*nit)->Store( prop_key, vc );
                    }
                }
           break;
         default:
           throw csmp::Exception( ERROR, "flagToNumber:", "property placement not handled." );
      }
   
 } // end flagToNumber

template void flagToNumber( Model<1U>&, const char*, const char* );
template void flagToNumber( Model<2U>&, const char*, const char* );
template void flagToNumber( Model<3U>&, const char*, const char* );






template<size_t dim>
double64 maximumResidual( const Model<dim>& sg,
                          const char* new_property, const char* old_property, bool absolute )
 {
   double64                  residual(0.), temp;
   ScalarVariable  new_prop, old_prop;
   
   const csmp::Region<dim>&  sgref(sg.Region("Model"));
   csmp::Index  new_key = sg.Database().StorageKey(new_property);
   csmp::Index  old_key = sg.Database().StorageKey(old_property);
   
   if ( new_key.place != old_key.place ) {
       throw csmp::Exception( ERROR, "maximumResidual", 
                             "Properties to be compared do not have the same placement, residual cannot be computed, returning 0.0...!");
     }

   if ( new_key.place == NODE ) {
       for ( typename std::vector<Node<dim>*>::const_iterator
             nit = sgref.NodesBegin(); nit != sgref.NodesEnd(); nit++ ) {
           (*nit)->Read( new_key, new_prop );
           (*nit)->Read( old_key, old_prop );
           if ( !absolute ) temp = std::fabs( ( new_prop() - old_prop() ) / new_prop() );
           else             temp = std::fabs( new_prop() - old_prop() );
           if ( temp > residual ) residual = temp;
         }
     }

   else if ( new_key.place == ELEMENT ) {
       for ( typename std::vector<Element<dim>*>::const_iterator
             eit = sgref.ElementsBegin(); eit != sgref.ElementsEnd(); eit++ ) {
           (*eit)->Read( new_key, new_prop );
           (*eit)->Read( old_key, old_prop );
           if ( !absolute ) temp = std::fabs( ( new_prop() - old_prop() ) / new_prop() );
           else             temp = std::fabs( new_prop() - old_prop() );
           if ( temp > residual ) residual = temp;
         }
     }

   else if ( new_key.place == ELEMENT_INTEGRATION_POINT ) {
       for ( typename std::vector<Element<dim>*>::const_iterator
             eit = sgref.ElementsBegin(); eit != sgref.ElementsEnd(); eit++ ) 
         for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
           {
             (*eit)->Read( i, new_key, new_prop );
             (*eit)->Read( i, old_key, old_prop );
             if ( !absolute ) temp = std::fabs( ( new_prop() - old_prop() ) / new_prop() );
             else             temp = std::fabs( new_prop() - old_prop() );
             if ( temp > residual ) residual = temp;
           }
     }
   
   else {
       throw csmp::Exception( ERROR, "maximumResidual", 
                             "Cannot identify placement of variable, residual cannot be computed, returning 0.0...!");
     }

   return residual;
} 

template double64 maximumResidual( const Model<1U>&, const char*, const char*, bool );
template double64 maximumResidual( const Model<2U>&, const char*, const char*, bool );
template double64 maximumResidual( const Model<3U>&, const char*, const char*, bool );


/**

Explores all the nodes and the elements from the mesh by mesh traversal from the starting root nodes, 
and returns the explored nodes and elements.

The method depends on correct neighbor information.

@param mesh: MeshMananger pointer
@param nodes, elmts: the method returns deques of pointers to the nodes and the elements

@section application Application

*/
template<size_t dim>
void exploreNodesAndElementsFromMesh(MeshManager<dim>* mesh, std::deque<Node<dim>*>& nodes, std::deque<Element<dim>*>& elmts)
{
	// traversal of the existing mesh nodes to find all nodes and elements
	set<Element<dim>*>		explored_elements;
	set<Node<dim>*>	discovered_nodes;
	deque<Node<dim>*>	current_nodes;
	for (size_t g = 0U; g < mesh->NodeGroups(); g++) {
		auto root_node = mesh->RootNode(g);
		// starting at the root node
		discovered_nodes.insert(root_node);
		current_nodes.push_back(root_node);
		while (!current_nodes.empty()) {
			Node<dim>*  n_ptr(*current_nodes.begin());
			// for all parent elements of the current node
			for (size_t i = 0U; i < n_ptr->Parents(); i++) {
				// for all the nodes of each parent element
				for (size_t j = 0U; j < n_ptr->Parent(i)->Nodes(); j++)
					// if this node is not the one from which we started
					if (j != n_ptr->ParentNodeNumber(i)) {
						auto new_node = discovered_nodes.insert(n_ptr->Parent(i)->N(j));
						if (new_node.second) current_nodes.push_back(n_ptr->Parent(i)->N(j));
					}
				// storing the explored element
				explored_elements.insert(n_ptr->Parent(i));
			}
			// removing the node from the discovered (but not yet explored) deque
			current_nodes.pop_front();
		}
	}

	// assigning the explored nodes and elements from the node and element pointer vectors	
	nodes.assign(discovered_nodes.begin(), discovered_nodes.end());
	elmts.assign(explored_elements.begin(), explored_elements.end());
}

template void exploreNodesAndElementsFromMesh(MeshManager<1U>*, std::deque<Node<1U>*>&, std::deque<Element<1U>*>&);
template void exploreNodesAndElementsFromMesh(MeshManager<2U>*, std::deque<Node<2U>*>&, std::deque<Element<2U>*>&);
template void exploreNodesAndElementsFromMesh(MeshManager<3U>*, std::deque<Node<3U>*>&, std::deque<Element<3U>*>&);

template<size_t dim>
void exploreNodesAndElementsFromMesh(const MeshManager<dim>* mesh, std::deque<const Node<dim>*>& nodes, std::deque<Element<dim>*>& elmts)
{
	// traversal of the existing mesh nodes to find all nodes and elements
	set<Element<dim>*>		explored_elements;
	set<const Node<dim>*>	discovered_nodes;
	deque<const Node<dim>*>	current_nodes;
	for (size_t g = 0U; g < mesh->NodeGroups(); g++) {
		auto root_node = mesh->RootNode(g);
		// starting at the root node
		discovered_nodes.insert(root_node);
		current_nodes.push_back(root_node);
		while (!current_nodes.empty()) {
			auto n_ptr(*current_nodes.begin());
			// for all parent elements of the current node
			for (size_t i = 0U; i < n_ptr->Parents(); i++) {
				// for all the nodes of each parent element
				for (size_t j = 0U; j < n_ptr->Parent(i)->Nodes(); j++)
					// if this node is not the one from which we started
					if (j != n_ptr->ParentNodeNumber(i)) {
						auto new_node = discovered_nodes.insert(n_ptr->Parent(i)->N(j));
						if (new_node.second) current_nodes.push_back(n_ptr->Parent(i)->N(j));
					}
				// storing the explored element
				explored_elements.insert(n_ptr->Parent(i));
			}
			// removing the node from the discovered (but not yet explored) deque
			current_nodes.pop_front();
		}
	}

	// assigning the explored nodes and elements from the node and element pointer vectors	
	nodes.assign(discovered_nodes.begin(), discovered_nodes.end());
	elmts.assign(explored_elements.begin(), explored_elements.end());
}

template void exploreNodesAndElementsFromMesh(const MeshManager<1U>*, std::deque<const Node<1U>*>&, std::deque<Element<1U>*>&);
template void exploreNodesAndElementsFromMesh(const MeshManager<2U>*, std::deque<const Node<2U>*>&, std::deque<Element<2U>*>&);
template void exploreNodesAndElementsFromMesh(const MeshManager<3U>*, std::deque<const Node<3U>*>&, std::deque<Element<3U>*>&);


/**

Explores all the faces from the mesh by mesh traversal from the starting root faces,
and returns the explored faces.

The method depends on correct neighbor information.

@param mesh: MeshMananger pointer
@param faces: the method returns deques of pointers to the faces

@section application Application

*/
/// retrieves and returns the faces that are explored by mesh traversal from the starting root faces
template<size_t dim>
void exploreFacesFromMesh(MeshManager<dim>* mesh, std::deque<Face<dim>*>& faces)
{
	// traversal of the existing mesh root faces to find all faces	
	set<Face<dim>*>	discovered_faces;
	deque<Face<dim>*>	current_faces;

	for (size_t g = 0U; g < mesh->FaceGroups(); g++) {
		auto root_face = mesh->RootFace(g);
		// starting at the first face
		discovered_faces.insert(root_face);
		current_faces.push_back(root_face);
		while (!current_faces.empty()) {
			Face<dim>*  n_ptr(*current_faces.begin());
			// for all neighbor faces of the current face
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;

				// if this neighbor is new one					
				auto new_face = discovered_faces.insert(n_ptr->Neighbor(i));
				if (new_face.second) current_faces.push_back(n_ptr->Neighbor(i));
			}
			// removing the face from the discovered (but not yet explored) deque
			current_faces.pop_front();
		}
	}
	
	// assigning the explored faces from the face pointer vectors	
	faces.assign(discovered_faces.begin(), discovered_faces.end());
}

template void exploreFacesFromMesh(MeshManager<1U>*, std::deque<Face<1U>*>&);
template void exploreFacesFromMesh(MeshManager<2U>*, std::deque<Face<2U>*>&);
template void exploreFacesFromMesh(MeshManager<3U>*, std::deque<Face<3U>*>&);

template<size_t dim>
void exploreFacesFromMesh(const MeshManager<dim>* mesh, std::deque<const Face<dim>*>& faces)
{
	// traversal of the existing mesh root faces to find all faces	
	set<const Face<dim>*>	discovered_faces;
	deque<const Face<dim>*>	current_faces;

	for (size_t g = 0U; g < mesh->FaceGroups(); g++) {
		auto root_face = mesh->RootFace(g);
		// starting at the first face
		discovered_faces.insert(root_face);
		current_faces.push_back(root_face);
		while (!current_faces.empty()) {
			auto n_ptr(*current_faces.begin());
			// for all neighbor faces of the current face
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;

				// if this neighbor is new one					
				auto new_face = discovered_faces.insert(n_ptr->Neighbor(i));
				if (new_face.second) current_faces.push_back(n_ptr->Neighbor(i));
			}
			// removing the face from the discovered (but not yet explored) deque
			current_faces.pop_front();
		}
	}

	// assigning the explored faces from the face pointer vectors	
	faces.assign(discovered_faces.begin(), discovered_faces.end());
}

template void exploreFacesFromMesh(const MeshManager<1U>*, std::deque<const Face<1U>*>&);
template void exploreFacesFromMesh(const MeshManager<2U>*, std::deque<const Face<2U>*>&);
template void exploreFacesFromMesh(const MeshManager<3U>*, std::deque<const Face<3U>*>&);

/**

Explores all the interfaces from the mesh by mesh traversal from the starting root interfaces,
and returns the explored interfaces.

The method depends on correct neighbor information.

@param mesh: MeshMananger pointer
@param interfaces: the method returns deques of pointers to the interfaces

@section application Application

*/
/// retrieves and returns the interfaces that are explored by mesh traversal from the starting root interfaces
template<size_t dim>
void exploreInterFacesFromMesh(MeshManager<dim>* mesh, std::deque<InterFace<dim>*>& interfaces)
{
	// traversal of the existing mesh root interfaces to find all interfaces	
	set<InterFace<dim>*>	discovered_interfaces;
	deque<InterFace<dim>*>	current_interfaces;

	for (size_t g = 0U; g < mesh->InterFaceGroups(); g++) {
		auto root_interface = mesh->RootInterFace(g);
		// starting at the first interface
		discovered_interfaces.insert(root_interface);
		current_interfaces.push_back(root_interface);
		while (!current_interfaces.empty()) {
			InterFace<dim>*  n_ptr(*current_interfaces.begin());
			// for all neighbor interfaces of the current interface
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;

				// if this neighbor is new one					
				auto new_interface = discovered_interfaces.insert(n_ptr->Neighbor(i));
				if (new_interface.second) current_interfaces.push_back(n_ptr->Neighbor(i));
			}
			// removing the interface from the discovered (but not yet explored) deque
			current_interfaces.pop_front();
		}
	}

	// assigning the explored interfaces from the interface pointer vectors	
	interfaces.assign(discovered_interfaces.begin(), discovered_interfaces.end());
}

template void exploreInterFacesFromMesh(MeshManager<1U>*, std::deque<InterFace<1U>*>&);
template void exploreInterFacesFromMesh(MeshManager<2U>*, std::deque<InterFace<2U>*>&);
template void exploreInterFacesFromMesh(MeshManager<3U>*, std::deque<InterFace<3U>*>&);

template<size_t dim>
void exploreInterFacesFromMesh(const MeshManager<dim>* mesh, std::deque<const InterFace<dim>*>& interfaces)
{
	// traversal of the existing mesh root interfaces to find all interfaces	
	set<const InterFace<dim>*>		discovered_interfaces;
	deque<const InterFace<dim>*>	current_interfaces;

	for (size_t g = 0U; g < mesh->InterFaceGroups(); g++) {
		auto root_interface = mesh->RootInterFace(g);
		// starting at the first interface
		discovered_interfaces.insert(root_interface);
		current_interfaces.push_back(root_interface);
		while (!current_interfaces.empty()) {
			auto n_ptr(*current_interfaces.begin());
			// for all neighbor interfaces of the current interface
			for (size_t i = 0U; i < n_ptr->Neighbors(); i++) {
				if (n_ptr->Neighbor(i) == NULL) continue;

				// if this neighbor is new one					
				auto new_interface = discovered_interfaces.insert(n_ptr->Neighbor(i));
				if (new_interface.second) current_interfaces.push_back(n_ptr->Neighbor(i));
			}
			// removing the interface from the discovered (but not yet explored) deque
			current_interfaces.pop_front();
		}
	}

	// assigning the explored interfaces from the interface pointer vectors	
	interfaces.assign(discovered_interfaces.begin(), discovered_interfaces.end());
}

template void exploreInterFacesFromMesh(const MeshManager<1U>*, std::deque<const InterFace<1U>*>&);
template void exploreInterFacesFromMesh(const MeshManager<2U>*, std::deque<const InterFace<2U>*>&);
template void exploreInterFacesFromMesh(const MeshManager<3U>*, std::deque<const InterFace<3U>*>&);


/**

Attempts a floodfill on the supplied set of elements, this so identified
contiguous model region is returned.   

The method depends on correct neighbor information.  

@section arguments Input Arguments 

Input element range iterator: In this region it proceeds to identify 
a contiguous domain.  

@param elements_contiguous_subset the method returns a set of pointers to
those elements forming the first contiguous domain that
it was able to reach from the supplied iterator.   

@section application Application

To break regions into contiguous subdomains.

*/
template<size_t dim>
void floodFill( Element<dim>* const eptr, set<Element<dim>*>& elements_contiguous_subset )
 {
    assert( eptr != NULL );
    // identifying the neighbors of the first element to be looked at
    deque<Element<dim>*>  neighbor_elements;
    const size_t  neighbors(eptr->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( eptr->Neighbor(i) != NULL )
        neighbor_elements.push_back( eptr->Neighbor(i) );
 
    // performing the floodfill, starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( static_cast<Element<dim>*>(eptr) );
      
    // element set for subsequent passes  
    deque<Element<dim>*>  new_neighbor_elements;
   
     while( !neighbor_elements.empty() )
       {
          // 1. loop over those neighbors that are not already part of the list
          for ( typename deque<Element<dim>*>::const_iterator
                nit=neighbor_elements.begin(); nit!=neighbor_elements.end(); ++nit )
            // if the element has not already been dealt with
            if ( elements_contiguous_subset.find( (*nit) ) == elements_contiguous_subset.end() ) {
                const size_t  neighbors((*nit)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( (*nit)->Neighbor(j) != NULL )
                        new_neighbor_elements.push_back( (*nit)->Neighbor(j) );
                elements_contiguous_subset.insert( (*nit) );
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_elements.clear();
      }
     
 } // end floodFill 


template void floodFill( Element<1U>* const, set<Element<1U>*>& );
template void floodFill( Element<2U>* const, set<Element<2U>*>& );
template void floodFill( Element<3U>* const, set<Element<3U>*>& );






/**
    Retrieves and returns the element ids of the first contiguous element patch
    that can be reached by mesh traversal from the starting element.
    
    @attention the elements in the region must have a unique numbering scheme.
*/
template<size_t dim> 
void floodFillViaIndexes( const Region<dim>& gref, size_t starting_idx,
                          std::set<size_t>& elements_contiguous_subset )
 {
    assert( starting_idx < gref.Elements() );
    assert( gref.E(starting_idx) != NULL );
    // identifying the neighbors of the first element to be looked at
    deque<size_t>  neighbor_elements;
    const size_t  neighbors(gref.E(starting_idx)->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( gref.E(starting_idx)->Neighbor(i) != NULL )
        neighbor_elements.push_back( gref.E(starting_idx)->Neighbor(i)->Idx() );
 
    // performing the floodfill, starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( starting_idx );
      
    // element set for subsequent passes  
    deque<size_t>  new_neighbor_elements;
   
     while( !neighbor_elements.empty() )
       {
          // 1. loop over those neighbors that are not already part of the list
          for ( typename deque<size_t>::const_iterator
                idx=neighbor_elements.begin(); idx!=neighbor_elements.end(); ++idx )
            // if the element has not already been dealt with
            if ( elements_contiguous_subset.find( (*idx) ) == elements_contiguous_subset.end() ) {
                const size_t  neighbors(gref.E(*idx)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( gref.E(*idx)->Neighbor(j) != NULL ) {
                        assert( gref.E(*idx)->Neighbor(j)->Idx() < gref.Elements() );
                        new_neighbor_elements.push_back( gref.E(*idx)->Neighbor(j)->Idx() );
                        elements_contiguous_subset.insert( gref.E(*idx)->Neighbor(j)->Idx() );
                    }
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_elements.clear();
      }
     
 } // end floodFillViaIndexes


template void floodFillViaIndexes( const Region<1U>&, size_t, set<size_t>& );
template void floodFillViaIndexes( const Region<2U>&, size_t, set<size_t>&  );
template void floodFillViaIndexes( const Region<3U>&, size_t, set<size_t>&  );










/**

StripDomainEdgesFor() tests the spatial distribution of a scalar input element
property for outliers and removes these. Property outliers are elements 
that constitute a property value boundary with >=two of their faces.
When such elements are detected, their property value is set to the
average of the surrounding elements. The property is changed only if 
the outside property is either smaller or greater than the element 
property. 

@section arguments Input Arguments 

The name of the property whose variations over the mesh shall be defined
by relatively smooth boundaries. 

@section implementation Implementation

The method using the connections among elements to test whether the 
element represents a property outlier. 

@section application Application

The method is used for regular meshes which were created from pixel-type
input data. In this case the method allows to capitalize on the element
splits which were created by the Triangulator meshing tool along 
property boundaries. The result are boundaries with 45o segments that
superseed the stepwise property boundaries of the original mesh. 

Numerous calls to StripDomainEdgesFor() also allow to erode regions 
defined by stepwise property variations. 

@section messages Messages 

The method will always warn the user that the model properties are 
modified. If the target property is not an element property, the 
simulation will be halted by a fatal error. 

If the target property is not a scalar variable the method will return
without modifying the target property and it will report a warning.  

 */
void stripDomainEdgesFor( Model<2U>& sg, const char* el_prop )
 {
     csmp::Index  prop_key = sg.Database().StorageKey(el_prop);

     if ( prop_key.place != ELEMENT )
       throw csmp::Exception( FATAL_ERROR, "stripDomainEdgesFor<2U>::StripDomainEdgesFor", 
                                    "The requested property is not an element variable");

     if ( prop_key.type != SCALAR ) {
          throw csmp::Exception( WARNING, "stripDomainEdgesFor<double64oat,2U>::StripDomainEdgesFor", 
                                   "only SCALAR variables are handled so far");
          return;
       }
     
    map<size_t,ScalarVariable > new_sc_data;
    ScalarVariable              sc;

    csmp::Region<2>&  super_group(sg.Region("Model"));

    for ( size_t n=0U; n<super_group.Elements(); n++ )
       {
          //  for elements that are not located at model boundary
          if ( super_group.E(n)->AtBoundary() == NOT )
            {
               // getting the scalar variable data
               super_group.E(n)->Read( prop_key, sc );
              
               // checking whether element-property should be changed
               // because the element is located at a region boundary
               // ---------------------------------------------------
               // 1. counting the surrounding values that are different from el-value
               double64     sc_sum(0U);
               unsigned int counter(0U);
               for ( size_t i=0U; i<super_group.E(n)->Neighbors(); i++ ) 
                 if ( sc() > super_group.E(n)->Neighbor(i)->Read( prop_key ) ) { 
                      sc_sum += super_group.E(n)->Neighbor(i)->Read( prop_key );
                      counter++;
                   }
               // if more than 2 neighbors have a different property value, this value
               // is assigned to the element
               if ( counter >= 2U ) sc = sc_sum / static_cast<double64>(counter);
          
               // storing the new values of only those elements that must be changed
               new_sc_data[ n ] = sc;
            }
       }
    
     // modyfying those elements that were found to be isolated
     // this implies that isolated squares are removed    
     for ( map<size_t,ScalarVariable >::iterator
           sc_it=new_sc_data.begin(); sc_it!=new_sc_data.end(); sc_it++ )
       super_group.E( (*sc_it).first )->Store( prop_key, (*sc_it).second );
       
     cout <<"\nstripDomainEdgesFor<2U>::StripDomainEdgesFor: ";
     cout <<"Boundaries of property regions have been modified..."<< endl;  
            
   } // end StripRoughDomainEdgesFor






/**
    Connects Element objects up to their same-dimensional neighbors
    in as much as is possible. Where there are no neighbors the neighbor pointers 
    will be nulled.
    
    The assumption is made that all nodes have a unique numbering.
    
    @author SKM 2012
*/

// CONNECTIVITY BETWEEN ELEMENTS

template<size_t dim>
void  establishNeighborConnectivity( vector<Element<dim>*>& simplexVector, bool unassign_neighbors_outside )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( Element )", "supplied element vector is empty; nothing was done" );
         return;
      }

    cout << "\nestablishNeighborConnectivity( Element ): Establishing CSMP FE neighbor connectivity...\n";
 
    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >  volume_neighbor_keys,
                                                           surface_neighbor_keys, 
                                                           line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<Element<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( size_t face=0U; face<(*it)->Faces(); face++ ) 
        {
           if ( (*it) == NULL ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( Element )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }
           // creating face key from idx's of face
           (*it)->FE()->NodesOfFace( face, fnids ); 
           for ( size_t j=0U; j<fnids.size(); j++ )
               key.insert( (*it)->N( fnids[j] ) );
             
           // inserting newly generated keys into multimap
           if ( (*it)->FE()->IsVolumeElement() )
               volume_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else // for all line elements
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           key.clear();

           // unassign neirghbors outside of the provided vector range
		   // JC: fix it since ElementRemeshingTrait was removed.
           //if( unassign_neighbors_outside )
               //(*it)->UnassignNeighbors();
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)

    cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(NULL);
        Element<dim>* e2Ptr(NULL);

        cout << "\n\t\tline elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                                 it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) ) 
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // line elements
    
    // 2.2 surface elements
    // --------------------
    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(NULL);
        Element<dim>* e2Ptr(NULL);

        cout << "\n\t\tsurface elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                                 it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) ) 
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                          face pointer                   nbor face idx        neighbor pointer
                    ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                    ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( dim == 3U and !volume_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(NULL);
        Element<dim>* e2Ptr(NULL);

        cout << "\n\t\tvolume elements...\n";

        //                key             n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >::iterator it1(volume_neighbor_keys.begin()),
                                                                                 it2(volume_neighbor_keys.begin());

        it2++;

        while ( it2 != volume_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ( (*it1).second.second != 0 and (*it2).second.second != 0 ) ) 
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == volume_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // dim=3
    
 } // end establishNeighborConnectivity

// explicit instantiations
template void establishNeighborConnectivity<1U>( std::vector<csmp::Element<1U>*>&, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::Element<2U>*>&, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::Element<3U>*>&, bool );






// CONNECTIVITY BETWEEN FACES

/*
template<size_t dim>
void  establishNeighborConnectivity( std::vector<csmp::Face<dim>*>& simplexVector, bool unassign_neighbors_outside )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( Face )", "supplied element vector is empty; nothing was done" );
         return;
      }

    cout << "\nestablishNeighborConnectivity( Face ): Establishing CSMP FE neighbor connectivity...\n";

    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<size_t,Face<dim>*> >  surface_neighbor_keys,
                                                        line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<Face<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( size_t face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == NULL ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( Face )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }
           // creating face key from idx's of face
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( size_t j=0U; j<fnids.size(); j++ )
               key.insert( (*it)->N( fnids[j] ) );

           // inserting newly generated key into multimap
           if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsLineElement() )
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else{
               csmp_error.notice( ERROR, "establishNeighborConnectivity( Face )", "supplied element vector contains volumetric element! nothing was done" );
               return;
           }
           key.clear();

           // unassign neighbors outside of the provided vector range
           if( unassign_neighbors_outside )
               (*it)->UnassignNeighbors();
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() ) {

        Face<dim>* e1Ptr(NULL);
        Face<dim>* e2Ptr(NULL);

        cout << "\n\t\tline elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Face<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                              it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                   (*it2).second.second->Assign( (*it2).second.first, e1Ptr );

                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=1

    // 2.2 surface elements
    // --------------------
    cout << "surface elements...";

    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        Face<dim>* e1Ptr(NULL);
        Face<dim>* e2Ptr(NULL);

        cout << "\n\t\tsurface elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Face<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                              it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                          face pointer                   nbor face idx        neighbor pointer
                    (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                    (*it2).second.second->Assign( (*it2).second.first, e1Ptr );


                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // surfaces


 } // end establishNeighborConnectivity

template void establishNeighborConnectivity<1U>( std::vector<csmp::Face<1U>*>&, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::Face<2U>*>&, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::Face<3U>*>&, bool );
*/


// CONNECTIVITY BETWEEN INTERFACES

template<size_t dim>
void  establishNeighborConnectivity( std::vector<csmp::InterFace<dim>*>& simplexVector, bool unassign_neighbors_outside )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( InterFace ):", "supplied element vector is empty; nothing was done." );
         return;
      }

    cout << "\nestablishNeighborConnectivity( InterFace ): Establishing CSMP FE neighbor connectivity...\n";

    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<size_t,InterFace<dim>*> >  surface_neighbor_keys,
                                                             line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<InterFace<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( size_t face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == NULL ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }

           // creating face key from idx's of face
           (*it)->CurrentSide( INSIDE );
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( size_t j=0U; j<fnids.size(); j++ )
           {
               key.insert( (*it)->N( fnids[j], INSIDE ) );
               key.insert( (*it)->N( fnids[j], OUTSIDE ) );
           }

           // inserting newly generated key into multimap
           if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsLineElement() )
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else{
               csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )", "supplied element vector contains volumetric element! nothing was done" );
               return;
           }
           key.clear();

           // JC: check it later! unassign neirghbors outside of the provided vector range
           //if( unassign_neighbors_outside )
               //(*it)->UnassignNeighbors();
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    cout << "line elements...";

    if ( !line_neighbor_keys.empty() ) {

        InterFace<dim>* e1Ptr(NULL);
        InterFace<dim>* e2Ptr(NULL);

        cout << "\n\t\tline elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,InterFace<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                                   it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                   (*it2).second.second->Assign( (*it2).second.first, e1Ptr );

                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=1

    // 2.2 surface elements
    // --------------------
    cout << "surface elements...";

    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        InterFace<dim>* e1Ptr(NULL);
        InterFace<dim>* e2Ptr(NULL);

        cout << "\n\t\tsurface elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,InterFace<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                                   it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                          face pointer                   nbor face idx        neighbor pointer
                    (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                    (*it2).second.second->Assign( (*it2).second.first, e1Ptr );


                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=2


 } // end establishNeighborConnectivity

template void establishNeighborConnectivity<1U>( std::vector<csmp::InterFace<1U>*>&, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::InterFace<2U>*>&, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::InterFace<3U>*>&, bool );



/**
    loops over the surface cells of the model subdomain, 
    checking whether any of the projections of the normals of the neighbor 
    cells onto the normal of the current element are negative.
    
    @return if any of the projections is negative, method returns false.
    A negative projection will mean that the normals of the neighboring cells 
    are at > to 90^o to the current cell; this should not be the case
    unless the surface has cusps with in it.
    
    @return if the region does not consist entirely of surface elements
    an error is reported.
    
    @attention even if the normals have a consistent orientation, this method
    may return false if the surface contains a cusp (>90^o kink).
    
    @test SKM 30/3/2016 for 3D only
 
*/
template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<3U>&  subdomain )
 {
    vector<double64> normal(3U), nbor_normal(3U);
   
    size_t non_surface_elements(0U);
    for ( auto it=subdomain.ElementsBegin(); it!=subdomain.ElementsEnd(); ++it )
      // this method only considers surface elements
      if ( (*it)->IsSurfaceElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( size_t i=0U; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double64 result(0.);
                  for ( size_t j=0U; j<3U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_surface_elements++;
   
    if ( non_surface_elements > 0U )
      ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation (3D):",
                                       subdomain.Name(), "region contained not only surface elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation


template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<2U>&  subdomain )
 {
    vector<double64> normal(2U), nbor_normal(2U);
   
    size_t non_line_elements(0U);
    for ( auto it=subdomain.ElementsBegin(); it!=subdomain.ElementsEnd(); ++it )
      // this method only considers line elements
      if ( (*it)->IsLineElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( size_t i=0U; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double64 result(0.);
                  for ( size_t j=0U; j<2U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_line_elements++;
   
    if ( non_line_elements > 0U )
      ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation (2D):",
                                       subdomain.Name(), "region contained not only line elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation





/**
    Stub: ID there are no boundaries so this should be a compile time assert
    
    @todo TODO: use static_assert<> here on the template argument
*/
template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<1U>&  subdomain )
 {
    ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation (1D):",
                                     subdomain.Name(), "one-dimensional models have no boundaries." );
    return false;
   
 } // end checkNeighborNormalsForConsistentOrientation






/**
     Assigns chosen node coordinate (x or y or z) to the target node variable.
*/
template<size_t dim>
void assignNodeCoordinatesTo( Model<dim>& sg, const char coordinate, const char* node_var )
 {
      csmp::Index nvar_key = sg. Database().StorageKey(node_var);
      assert( nvar_key.place == NODE );
      
      Region<dim>&  sgref(sg.Region("Model"));

      const typename vector<Node<dim>* >::iterator  nodesEnd(sgref.NodesEnd());
      
      if ( coordinate == 'x' or coordinate == 'X' )
        for ( typename vector<Node<dim>* >::iterator nit=sgref.NodesBegin(); nit!=nodesEnd; ++nit )
          (*nit)->Store( nvar_key, makeScalar( (*nit)->Status(nvar_key), (*nit)->x() ) );
        
      if ( dim > 1 and (coordinate == 'y' or coordinate == 'Y') )
        for ( typename vector<Node<dim>* >::iterator nit=sgref.NodesBegin(); nit!=nodesEnd; ++nit )
          (*nit)->Store( nvar_key, makeScalar( (*nit)->Status(nvar_key), (*nit)->y() ) );

      if ( dim > 2 and (coordinate == 'z' or coordinate == 'Z') )
        for ( typename vector<Node<dim>* >::iterator nit=sgref.NodesBegin(); nit!=nodesEnd; ++nit )
          (*nit)->Store( nvar_key, makeScalar( (*nit)->Status(nvar_key ), (*nit)->z() ) );
 
 }  // end 

template void assignNodeCoordinatesTo( Model<1U>&, const char, const char* );
template void assignNodeCoordinatesTo( Model<2U>&, const char, const char* );
template void assignNodeCoordinatesTo( Model<3U>&, const char, const char* );




/**

Function evaluates that the VSet connectivity is exactly the same
as the data in the current Model!

The Model is used as the reference case.

*/
template<size_t dim>
bool compareConnectivity( const Model<dim>& sg, const VSet<dim>& vset )
 {
    bool correct(true);
   
    const Region<dim>&  gref(sg.Region("Model")); 
    if ( gref.Elements() != vset.Elements() ) cout <<"\ncompareConnectivity: element number mismatch."<< endl;
    if ( gref.Nodes() != vset.Vertices() ) cout <<"\ncompareConnectivity: node number mismatch."<< endl;
  
    // 1. plist
    for ( uint32 i=0U; i<gref.Elements(); i++ )
      {
         for ( uint32 j=0U; j<gref.E(i)->Nodes(); j++ )
           if ( gref.E(i)->N(j)->Idx() != vset.Plist( gref.E(i)->Idx(), j ) ) {
                 cerr <<"\ncompareConnectivity: plist inconsistency: sg node id: "<< gref.E(i)->N(j)->Idx();
                 cerr <<" vs. vset nid: "<< vset.Plist( gref.E(i)->Idx(), j );
                 correct = false;
             }
      }
    
    // 2. pfverts
    for ( uint32 i=0U; i<gref.Elements(); i++ )
      {
         for ( uint32 j=0U; j<gref.E(i)->Neighbors(); j++ )
           if ( gref.E(i)->Neighbor(j) and
                static_cast<int32>(gref.E(i)->Neighbor(j)->Idx()) != vset.Pfvert( gref.E(i)->Idx(), j ) ) {
                 cerr <<"\ncompareConnectivity: plist inconsistency: sg node id: "<< gref.E(i)->Neighbor(j)->Idx();
                 cerr <<" vs. vset nid: "<< vset.Pfvert( gref.E(i)->Idx(), j );
                 correct = false;
             }
      }
   
   return correct;
    
 } // end compare

template bool compareConnectivity<2U>( const Model<2U>&, const VSet<2U>& );
template bool compareConnectivity<3U>( const Model<3U>&, const VSet<3U>& );



/**
    Imposes a user-defined upper or lower limit on the value of the variable of interest.
    
    For a vector variable, its length gets scaled to the limit value.
    For a tensor variable nothing can be done yet, so an exception is thrown.
    
    @author SKM 7/9/2014
*/
template<size_t dim>
void imposeLimitOn( Model<dim>& model, const char* region, const char* variable, bool upper_limit, double64 limit_value )
 {
    Region<dim>&  rref(model.Region(region));
    csmp::Index   prop_key(model.Database().StorageKey(variable));
    double64      min, max;
    model.Database().RangeOf( variable, min, max );
   
    if ( upper_limit && limit_value > max ) {
         cerr <<"\nIntended upper limit on variable '"<< variable <<"' exceeds that defined in database: ";
         cerr << limit_value <<" vs. "<< max << endl;
         throw csmp::Exception( ERROR, "imposeLimitOn:", "user-defined limit is out of bounds specified in variable database." );
      }
    if ( !upper_limit && limit_value < min ) {
         cerr <<"\nIntended lower limit on variable '"<< variable <<"' is lower than that defined in database: ";
         cerr << limit_value <<" vs. "<< min << endl;
         throw csmp::Exception( ERROR, "imposeLimitOn:", "user-defined limit is out of bounds specified in variable database." );
      }
   
    if ( prop_key.type == SCALAR ) {
        ScalarVariable  sc;
        if ( upper_limit )
          switch( prop_key.place )
            {
               case MODEL:
                  model.Store( prop_key, makeScalar( model.Status(prop_key), std::min(limit_value,rref.Read(prop_key)) ) );
                 break;
               case REGION:
                  rref.Store( prop_key, makeScalar( rref.Status(prop_key), std::min(limit_value,rref.Read(prop_key)) ) );
                 break;
               case ELEMENT:
                  for ( typename vector<Element<dim>*>::iterator
                        it=rref.ElementsBegin();  it!=rref.ElementsEnd(); ++it ) {
                      double64 val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::min(limit_value,val) ) );
                   }
                 break;
               case ELEMENT_INTEGRATION_POINT:
                  for ( typename vector<Element<dim>*>::iterator
                        it=rref.ElementsBegin();  it!=rref.ElementsEnd(); ++it )
                    for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ ) {
                         double64 val = (*it)->Read( i, prop_key );
                         (*it)->Store( i, prop_key, makeScalar( (*it)->Status(i,prop_key), std::min(limit_value,val) ) );
                      }
                 break;
               case NODE:
                  for ( typename vector<Node<dim>*>::iterator
                        it=rref.NodesBegin(); it!=rref.NodesEnd(); ++it ) {
                      double64 val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::min(limit_value,val) ) );
                   }
                 break;
               default:
                 throw csmp::Exception( ERROR, "imposeLimitOn:", "variable placement not recognized." );
            }
          else // if a lower limit shall be imposed
          switch( prop_key.place )
            {
               case MODEL:
                  model.Store( prop_key, makeScalar( model.Status(prop_key), std::max(limit_value,rref.Read(prop_key)) ) );
                 break;
               case REGION:
                  rref.Store( prop_key, makeScalar( rref.Status(prop_key), std::max(limit_value,rref.Read(prop_key)) ) );
                 break;
               case ELEMENT:
                  for ( typename vector<Element<dim>*>::iterator
                        it=rref.ElementsBegin();  it!=rref.ElementsEnd(); ++it ) {
                      double64 val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::max(limit_value,val) ) );
                   }
                 break;
               case ELEMENT_INTEGRATION_POINT:
                  for ( typename vector<Element<dim>*>::iterator
                        it=rref.ElementsBegin();  it!=rref.ElementsEnd(); ++it )
                    for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ ) {
                         double64 val = (*it)->Read( i, prop_key );
                         (*it)->Store( i, prop_key, makeScalar( (*it)->Status(i,prop_key), std::max(limit_value,val) ) );
                      }
                 break;
               case NODE:
                  for ( typename vector<Node<dim>*>::iterator
                        it=rref.NodesBegin(); it!=rref.NodesEnd(); ++it ) {
                      double64 val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::max(limit_value,val) ) );
                   }
                 break;
               default:
                 throw csmp::Exception( ERROR, "imposeLimitOn:", "variable placement not recognized." );
            }
      }
   
    // for a vector variable, its length gets scaled to the limit value
    else if ( prop_key.type == VECTOR ) {
        VectorVariable<dim>  vc;
        switch( prop_key.place )
          {
             case MODEL: {
                    model.Read( prop_key, vc );
                    const double64 vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    // if the vector is too long it gets scaled back
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    model.Store( prop_key, vc );
                 }
               break;
             case REGION: {
                    rref.Read( prop_key, vc );
                    const double64 vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    rref.Store( prop_key, vc );
                 }
               break;
             case ELEMENT:
                for ( typename vector<Element<dim>*>::iterator
                      it=rref.ElementsBegin();  it!=rref.ElementsEnd(); ++it ) {
                    (*it)->Read( prop_key, vc );
                    const double64 vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    (*it)->Store( prop_key, vc );
                 }
               break;
             case ELEMENT_INTEGRATION_POINT:
                for ( typename vector<Element<dim>*>::iterator
                      it=rref.ElementsBegin();  it!=rref.ElementsEnd(); ++it )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ ) {
                       (*it)->Read( i, prop_key, vc );
                        const double64 vmagnitude = vc.Length();
                        assert( vmagnitude > 0. );
                        if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                        else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                       (*it)->Store( i, prop_key, vc );
                    }
               break;
             case NODE:
                for ( typename vector<Node<dim>*>::iterator
                      it=rref.NodesBegin(); it!=rref.NodesEnd(); ++it ) {
                    (*it)->Read( prop_key, vc );
                    const double64 vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    (*it)->Store( prop_key, vc );
                 }
               break;
             default:
               throw csmp::Exception( ERROR, "imposeLimitOn:", "variable placement not recognized." );
          }
      }
    else throw csmp::Exception( ERROR, "imposeLimitOn:", "variable type not recognized." );

 } // end imposeLimitOn

template void imposeLimitOn( Model<1U>&, const char*, const char*, bool, double64 );
template void imposeLimitOn( Model<2U>&, const char*, const char*, bool, double64 );
template void imposeLimitOn( Model<3U>&, const char*, const char*, bool, double64 );






/**
    Computes barycentre-to-node distances for all parent elements of node and returns them into the supplied vector
    distances_and_weight vector [e1,e2...e_n,e_sum] with a size of parent elements+1
*/
template<size_t dim>
void distanceWeights( typename vector<Node<dim>*>::const_iterator nodes_begin,
                      typename vector<Node<dim>*>::const_iterator nodes_end,
                      vector<vector<double64> >& distances_and_weights )
 {
     distances_and_weights.resize(distance(nodes_begin,nodes_end));
     size_t  node_index(0U);
   
     while( nodes_begin != nodes_end )
       {
          double64 weight(0.);
          const Point<dim> npt((*nodes_begin)->Coordinate());
          const size_t parents((*nodes_begin)->Parents());
          distances_and_weights[node_index].reserve(parents+1U);
         
          for ( size_t i=0U; i<parents; i++ ) {
               // recording the node-to-barycentre distances
               distances_and_weights[node_index].push_back( npt.DistanceTo( (*nodes_begin)->Parent(i)->BaryCenter() ) );
               weight += distances_and_weights[node_index][i];
            }
          // storing the weights (=sum of the distances) in the last element of the vector
          distances_and_weights[node_index][parents] = weight;
          nodes_begin++;
          node_index++;
       }
   
 } // end distanceWeights

template void distanceWeights<1>( vector<Node<1>*>::const_iterator, vector<Node<1>*>::const_iterator, vector<vector<double64> >& );
template void distanceWeights<2>( vector<Node<2>*>::const_iterator, vector<Node<2>*>::const_iterator, vector<vector<double64> >& );
template void distanceWeights<3>( vector<Node<3>*>::const_iterator, vector<Node<3>*>::const_iterator, vector<vector<double64> >& );



 
  
  
  
/**
    finds those elements in a model that contact eachother across split interfaces.
    Only those elements are discovered that are node-matched.
 
    @param subdomain (non-unique) region which contains the split boundary
 
    @return returns false if none of the perimeter elements are node matched
 
    @attention method only looks at highest dimensional elements in the model;
    thus, lower dimensional elements are ignored and the neighborhood relations
    on either side of them are returned.
 
    @attention the elements inside the model that are located along the splitboundary
    do not have neighbor pointers yet.
 
    @author SKM
    @date 14/01/2018
 
    @TODO do we need to remember which side of the interface we are on?
    @todo test method on Chloe's dataset
*/
template<size_t dim>
bool findSplitInterfaceElements( const Region<dim>& subdomain,
                                 set<pair<pair<Element<dim>*,size_t>,pair<Element<dim>*,size_t> > >& interface_elmt_pairs )
 {
    // multi-element container for all elements that are located on split boundaries
    //       face search key      element      face number
    multimap<set<Point<dim> >,pair<Element<dim>*,size_t> > element_face_keys;
   
    // 1. for all elements on the perimeter of the model subdomain,
    //    generate keys from their node coordinates that are then matched with one-another
    //    in order to connect these elements
    for ( size_t n=subdomain.InteriorElements(); n<subdomain.Elements(); ++n )
        {
           // for those element faces that define the perimeter surface
           for ( size_t i=0U; i<subdomain.PerimeterFaces(n); ++i )
             {
                // get the local node numbers of the perimeter face
                vector<size_t> fnids;
                subdomain.E(n)->FE()->NodesOfFace( subdomain.PerimeterFace( n, i ), fnids );
                // add the corresponding node points to a set that will form the element face key
                pair<set<Point<dim> >,size_t> face_key;
                for ( size_t j=0U; j<fnids.size(); ++j )
                  face_key.first.insert( subdomain.E(n)->N( fnids[j] )->Coordinate() );
                // remembering the face id
                face_key.second = subdomain.PerimeterFace( n, i );
                // storing the key in the correspondance search map
                //                              node-coordinate set  element pointer   local face id
                element_face_keys.insert( make_pair( face_key.first, make_pair( subdomain.E(n), face_key.second ) ) );
             }
        }
   
    // 2. searching map for matching interface elements
    //    - a match is obtained if the element pointers are different
    //    - if there is a match, the element pair is added to the OppositeElement container
    if ( !interface_elmt_pairs.empty() ) interface_elmt_pairs.clear();
    // searching
    for ( auto it=element_face_keys.begin(); it!=element_face_keys.end(); ++it ) {
         // multimap iterator pair containing the range of shared keys
         //  face search key (point set), element(Element*), face number (size_t)
         auto result = element_face_keys.equal_range( (*it).first );
         // if more than one value was found
         if ( distance( result.first, result.second ) > 1U ) {
             // there should not be any manyfolds
             assert( distance( result.first, result.second ) == 2U );
             // advancing the result range iterator as necessary to find an Element different from (*it).second.first
             do result.first++;
             while ( (*result.first).second.first == (*it).second.first );
           // we store the matching element pair
           // set<pair<pair<Element<dim>*,size_t>,pair<Element<dim>*,size_t> > >
           interface_elmt_pairs.insert( make_pair(
                                        make_pair( (*it).second.first, (*it).second.second ),
                                        make_pair( (*result.first).second.first, (*result.first).second.second ) )
                                      );
          }
      }
  
    // 3. checking the results
    if (  interface_elmt_pairs.empty() ) return false;
 
    return true;
    
 } // end findSplitInterfaceElements

template bool findSplitInterfaceElements( const Region<3U>&, set<pair<pair<Element<3U>*,size_t>,pair<Element<3U>*,size_t> > >& );
template bool findSplitInterfaceElements( const Region<2U>&, set<pair<pair<Element<2U>*,size_t>,pair<Element<2U>*,size_t> > >& );
template bool findSplitInterfaceElements( const Region<1U>&, set<pair<pair<Element<1U>*,size_t>,pair<Element<1U>*,size_t> > >& );



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


} // end namespace csmp
