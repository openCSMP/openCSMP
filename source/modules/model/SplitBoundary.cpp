#include "SplitBoundary.h"
#include "Region.h"
#include "Boundary.h"

#include "writeVariableIf.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

#include "Model.h"
#include "MeshManagementUtilities.h"

#include "FiniteVolumeStencilManager.h"

#include "ErrorHandler.h"
//#include "CSMP_highLevelUtilities.h"
#include "PL_Utilities.h"
#include "FEM_Data.h"

#include "Visitor.h"
#include "variableOperations.h"

//#define SPLITBOUNDARY_DEBUG

using namespace std;


namespace csmp {


template<size_t dim>
SplitBoundary<dim>::SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>& pref )
  : ModelSubDomain<dim, InterFace>( splitboundaryname, pref )
{
  this->ResizePropertyStorage( this->pref_.LocalVariablesAt( Placement() ) );
}


template<size_t dim>
SplitBoundary<dim>::SplitBoundary( const SplitBoundary& ed )
  : ModelSubDomain<dim, InterFace>( ed ),
    LocalVariableStorage<dim,SplitBoundary>( ed )
{
}


template<size_t dim>
SplitBoundary<dim>::SplitBoundary( SplitBoundary&& ed )
  : ModelSubDomain<dim, InterFace>( ed ),
    LocalVariableStorage<dim,SplitBoundary>( ed )
{
}


template<size_t dim>
SplitBoundary<dim>& SplitBoundary<dim>::operator=( const SplitBoundary<dim>& ed )
{
  if ( &ed != this ) {
     ModelSubDomain<dim,InterFace>::operator=( ed ); 
     this->LVS( ed.LVS() );
  }
  return *this;
}



/** RECONSTRUCTOR of split boundary
 
    @attention the InterFace objects referred to by SubDomainInfo are expected to be allready part of the model
    @attention this also means, that these were already interconnected with one-another when the mesh was build from VSet
    
        @author SKM
        @date 21/9/21
*/
template<size_t dim>
SplitBoundary<dim>::SplitBoundary( const PropertyDatabase<dim>& pref,
                                   MeshManager<dim>& mesh,
                                   const SubDomainInfo& info )
  : ModelSubDomain<dim, InterFace>( info.name, pref )
{
  // building the interface vector (for this particular region)
  // ----------------------------------------------------------
  if ( info.interior_elmts[0] != mesh.Elements() + mesh.Faces() )
    throw csmp::Exception( ERROR, "SplitBoundary(reconstructor)",
                          "InterFace numbering is expected to start at the number of elements + faces");

  this->elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );
  for ( auto i : info.interior_elmts ) this->elmt_vec_.push_back( &(*next(mesh.InterFacesBegin(),i)) );
  for ( auto i : info.perimeter_elmts ) this->elmt_vec_.push_back( &(*next(mesh.InterFacesBegin(),i)) );

  // building the node vector
  // ------------------------
  // assigning pointers to the interior and perimeter nodes
  this->first_bd_node_ = info.interior_nodes.size();
  this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );
  for ( auto i : info.interior_nodes ) this->node_vec_.push_back( &(*next(mesh.NodesBegin(),i)) );
  for ( auto i : info.perimeter_nodes ) this->node_vec_.push_back( &(*next(mesh.NodesBegin(),i)) );

  // sorting of the pointers is necessary because the memory addresses of the new pointers will be different than in the last model
  this->SortVectors( info.interior_elmts.size(), info.interior_nodes.size() );

  // building the vector of vectors of those faces of the interfaces that lie on the subdomain perimeter
  // ---------------------------------------------------------------------------------------------------
  this->BuildPerimeterFaceVector( info.interior_elmts.size() );

  // allocating the storage for subdomain properties
  // -----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( SPLIT_BOUNDARY ) );
  
} // end Constructor

  // removing duplicates from a vector
  //sort( this->node_vec_.begin(), this->node_vec_.end() );
  //this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );






/**
Constructor of a split boundary from the set of juxtaposed elements
Asks mesh manager to construct corresponding faces.

Can be used to create Splitboundaries from node-matched interface meshes
which already contain multiplicated yet collocated nodes (can be done in ANSYS).

@author SKM 15/08/2018

*/
template<size_t dim>
SplitBoundary<dim>::SplitBoundary( std::string splitboundaryname,
                                   const PropertyDatabase<dim>& pref,
                                   const FiniteElementManager& femgr,
                                   MeshManager<dim>& mesh,
                                   const InterFaceSet<dim>& ifset )
  : ModelSubDomain<dim, InterFace>( splitboundaryname, pref )
{
  const LocalVariables&             ifvars( pref.LocalVariablesAt( INTER_FACE ) );
  const IntegrationPointVariables&  if_ip_vars( pref.IntegrationPointVariablesAt( INTER_FACE_INTEGRATION_POINT ) );

  // 1. getting mesh manager to build interfaces according to specifications
  // -----------------------------------------------------------------------
  this->elmt_vec_.reserve( ifset.size() );

  // for all the interfaces of the new split boundary
  for ( auto it = ifset.begin(); it != ifset.end(); ++it )
    // getting the element type that the interface shall represent
    // from the first higher dimensional neighbor element
    // InterFaceSet member:   pair<pair<Element<dim>*,size_t>, pair<Element<dim>*,size_t> >
    //                        first high-dim. nbor interface at interface
    // construction with connectivity and number continueing from already existing interfaces
    this->elmt_vec_.push_back( mesh.AddInterFace( ifset.InnerElement(it), ifset.InnerFaceID(it),
                                                  ifset.OuterElement(it), ifset.OuterFaceID(it),
                                                  ifvars, if_ip_vars ) );

  // 2. establising interface neighbor connectivity and interior vs. perimeter includig sorting
  // ---------------------------------------------------------------------------------------------------
  establishNeighborConnectivity( this->elmt_vec_, INSIDE );
  establishNeighborConnectivity( this->elmt_vec_, OUTSIDE );

  // 3. building the interface node vector
  // ---------------------------------------------------------------------------------------------------
  this->node_vec_.reserve( this->elmt_vec_.size() );
  for ( auto& it : this->elmt_vec_ )
    for ( size_t i = 0U; i<it->Nodes(); ++i ) this->node_vec_.push_back( it->N( i ) );
  // sorting node vector and making it unique
  sort( this->node_vec_.begin(), this->node_vec_.end() );
  this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );

  // 4. sorting interfaces and nodes and building the boundary interface vector
  // ---------------------------------------------------------------------------------------------------
  this->IdentifyPerimeter();

  // 5. allocating the storage for subdomain properties
  // --------------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( SPLIT_BOUNDARY ) );

} // end constructor





/// Does not delet interfaces, Delete() has to be called for this
template<size_t dim>
SplitBoundary<dim>::~SplitBoundary()
{
}



template<size_t dim>
IntegrationPointVariables  SplitBoundary<dim>::InterFaceIntegrationPointVariables() const
{ return this->pref_.IntegrationPointVariablesAt( INTER_FACE ); }



template<size_t dim>
LocalVariables  SplitBoundary<dim>::InterFaceVariables() const
{ return this->pref_.LocalVariablesAt( INTER_FACE ); }


// LOCAL VARIABLE STORAGE INTERFACE
template<size_t dim>
bool SplitBoundary<dim>::ValidVariable( const char* variableName ) const
{
  const PLACEMENT p( this->pref_.Placement( variableName ) );
  if ( p == NODE || p == INTER_FACE || p == SPLIT_BOUNDARY )
    return true;
  return false;
}


// VISITORS INTERFACE

/// visitation of a split boundary
template<size_t dim>
void SplitBoundary<dim>::Accept( Visitor<dim>& v )
{
  if ( v.ApplicationLevel() == MODEL || v.ApplicationLevel() == SPLIT_BOUNDARY )
    v.Visit( this );

  switch ( v.ApplicationTarget() ) {
    case MODEL:
      throw csmp::Exception( ERROR, "Region<dim>::Accept",
                             "ApplicationTarget MODEL; Visitor should have never arrived at this SplitBoundary" );
      break;
    case SPLIT_BOUNDARY:
      return;
      // element, interface and interface are treated the same
    case INTER_FACE:
      for ( typename vector<InterFace<dim>*>::iterator
            it = this->ElementsBegin(); it != this->ElementsEnd(); it++ )
        (*it)->Accept( v );
      return;
    case NODE:
      for ( typename vector<csmp::Node<dim>*>::iterator
            nd_it = this->NodesBegin(); nd_it != this->NodesEnd(); nd_it++ )
        (*nd_it)->Accept( v );
      return;
    default:
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::Accept",
                             "ApplicationTarget was not resolved; nothing was done" );
  }
} // end Accept


/// Returns the position of and adjacent region relative to the boundary. Relies on element Idx
template<size_t dim>
INTERFACE_SIDE SplitBoundary<dim>::RegionLocation( const Region<dim>& region )
{
  assert( !this->elmt_vec_.empty() );
  const auto regionElementsEnd( region.ElementsEnd() );
  const auto sbElementsEnd( this->ElementsEnd() );
  for ( auto ifit( this->ElementsBegin() ); ifit != sbElementsEnd; ++ifit )
    {
      const Element<dim>* const innerParent( (*ifit)->Parent( INSIDE ) );
      const Element<dim>* const outerParent( (*ifit)->Parent( OUTSIDE ) );
      for ( auto eit( region.ElementsBegin() ); eit != regionElementsEnd; ++eit )
        {
          if ( (*eit) == innerParent ) return INSIDE;
          if ( (*eit) == outerParent ) return OUTSIDE;
        } // region elements
    } // split boundary interfaces
  throw csmp::Exception( ERROR, "SplitBoundary<dim>::RegionLocation", "Region seems not to be adjacent to split boundary!" );
  // shouldn't get here
  return OUTSIDE;
}


// -----------------------------------------------
//  input/output
// -----------------------------------------------

/**
   for the assignment of properties that are unique to the instance of this subclass
   
      @author SKM
      @date 7/6/2020
*/
template<size_t dim>
template<typename Var>
void SplitBoundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index prop_key(this->pref_.StorageKey(input_prop));
      if ( prop_key.place == SPLIT_BOUNDARY ) {
           if ( sd != COMPLETE )
             csmp_error.notice( WARNING, "Region<dim>::InputPropertyValue",
                                          input_prop, "is a Region property and no distinction between INTERIOR and PERIMETER can be made" );
           this->Store( prop_key, new_value );
           return;
        }
      
      // incorrect applications of method  
      if ( prop_key.place == BOUNDARY || prop_key.place == REGION || prop_key.place == MODEL ||
           prop_key.place == ELEMENT || prop_key.place == FACE )
        csmp_error.notice( ERROR, "Region<dim>::InputPropertyValue",
                           input_prop, "must be a SPLIT_BOUNDARY, INTER_FACE/IP or NODE property for this method call to work" );        
    
     // for any different property placement, the method of the base-class is called
     ModelSubDomain<dim,InterFace>::InputPropertyValue( input_prop, new_value, sd );
      
  } // end InputPropertyValue

template void SplitBoundary<1U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const VectorVariable<1U>&, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const VectorVariable<2U>&, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const VectorVariable<3U>&, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const TensorVariable<1U>&, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const TensorVariable<2U>&, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const TensorVariable<3U>&, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );



template<size_t dim>
template<typename Var>
void SplitBoundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index key(this->pref_.StorageKey(input_prop));
      
      if ( key.place == SPLIT_BOUNDARY ) {
           if ( sd != COMPLETE )
             csmp_error.notice( WARNING, "SplitBoundary<dim>::InputPropertyValue",
                                input_prop, "is a SPLIT_BOUNDARY property and no distinction between INTERIOR and PERIMETER can be made" );
                                
           // only overwriting those variable components / rows that are not flagged 'do_not_overwrite' 
           writeVariableIf( this, key, new_value, do_not_overwrite );
           return;
        }
      
      // incorrect applications of method  
      if ( key.place == REGION  || key.place == BOUNDARY || key.place == MODEL || 
           key.place == ELEMENT || key.place == FACE )
        csmp_error.notice( ERROR, "SplitBoundary<dim>::InputPropertyValue",
                           input_prop, "must be a SPLIT_BOUNDARY, INTER_FACE/IP or NODE property for this method call to work" );        
    
     // for any different property placement, the method of the base-class is called
     ModelSubDomain<dim,InterFace>::InputPropertyValue( input_prop, new_value, do_not_overwrite, sd );
      
  } // end InputPropertyValue

// explicit instantiations
template void SplitBoundary<1U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const VectorVariable<1U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const VectorVariable<2U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const VectorVariable<3U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const TensorVariable<1U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const TensorVariable<2U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const TensorVariable<3U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<1U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<2U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void SplitBoundary<3U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );







// ------------------------------------------------------------------
// Building blocks
// -------------------------------------------------------------------


/**
       Indiscriminately accumulates inside, outside and intervening nodes, if any.
*/
template<size_t dim>
void SplitBoundary<dim>::CreateNodePointerVector()
{
  if ( this->elmt_vec_.empty() )
    throw csmp::Exception( ERROR, "SplitBoundary<dim>::CreateNodePointerVector:",
                           this->Name(), "interface vector is empty; nothing could be done." );

  if ( !this->node_vec_.empty() )
    this->node_vec_.clear();

  // creating the node index vector
  set<csmp::Node<dim>*>  nodes_set;
  for ( typename vector<InterFace<dim>*>::const_iterator it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
    for ( typename vector<Node<dim>*>::size_type i = 0U; i<(*it)->FE()->Nodes(); i++ )
    {
      nodes_set.insert( (*it)->N( i, INSIDE ) );
      nodes_set.insert( (*it)->N( i, OUTSIDE ) );
      if ( (*it)->HasInterveningElement() )
        nodes_set.insert( (*it)->N( i, MIDDLE ) );
    }

  this->node_vec_.assign( nodes_set.begin(), nodes_set.end() );
}


/**
Detects of how many spatial dimensions interface types are contained in model.
It returns a pair: first value gives number of different spatial dimensions contained,
second value returns the highest spatial dimension contained.

@author SKM 1/11/2013
*/
template<size_t dim>
pair<int32_t, int32_t>  SplitBoundary<dim>::InterFaceSpatialDimensions() const
{
  return this->SpatialDimensions();

} // end ElementSpatialDimensions




/**
    Creates a split boundary from a boundary. Requires unique indices.
    @author SKM 1/11/2013
    @author SKM 21/9/2021
*/
template<size_t dim>
bool  SplitBoundary<dim>::CreateFrom( MeshManager<dim>& mesh,
                                      Boundary<dim>& boundary )
{
  //LVS
  const LocalVariables lvsInterFace( InterFaceVariables() );
  const IntegrationPointVariables lvsIntegrationPoint( InterFaceIntegrationPointVariables() );

  // allocating SubDomain element container
  this->elmt_vec_.clear();
  this->elmt_vec_.reserve( boundary.Elements() );

  const typename vector<Face<dim>*>::const_iterator facesEnd( boundary.ElementsEnd() );
  for ( typename vector<Face<dim>*>::const_iterator fit( boundary.ElementsBegin() ); fit != facesEnd; ++fit )
    // the InterFace that is being build from the current interface
    this->elmt_vec_.push_back( mesh.ReplaceFaceByInterFace( (*fit), lvsInterFace, lvsIntegrationPoint ) );

    // free excessive allocated capacity
  vector<InterFace<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  // initialize splitboundary essentials 
  establishNeighborConnectivity( this->elmt_vec_ );
  
  this->IdentifyPerimeter();

  return true;
  
} // CreateFrom









/// CALCULATIONS
/**
Computes length (m) of the SplitBoundary object's perimeter curve.
Operation makes sense only in 33 because the perimeter of a line are just its end points.
*/
template<size_t dim>
double  SplitBoundary<dim>::Perimeter( INTERFACE_SIDE side ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( dim != 3U ) {
    csmp_error.notice( WARNING, "SplitBoundary<>::Perimeter:",
                       "returning 1.0 since perimeter is a point." );
    return 1.;
  }

  double        perimeter_length( 0. );
  vector<size_t>  fnids;
  size_t          n( 0U );

  for ( auto it = this->PerimeterElementsBegin(); it != this->ElementsEnd(); it++, n++ )
    for ( size_t i = 0U; i<this->PerimeterFaces( n ); i++ ) {
      (*it)->FE()->NodesOfFace( this->PerimeterFace( n, i ), fnids );
      perimeter_length += ((*it)->N( fnids[1] )->Coordinate() -
                            (*it)->N( fnids[0] )->Coordinate()).Length();
    }

  return perimeter_length;
}



/**
Is calculated on the basis of the Splitboundary bisector if the split nodes were
moved apart in the simulation process; else a particular side is used.
*/
template<size_t dim>
double  SplitBoundary<dim>::Area( INTERFACE_SIDE side ) const
{
  double  integrated_area( 0. );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if constexpr ( dim == 3U ) {
      for ( auto it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
        if ( (*it)->FE()->IsSurfaceElement() )
          integrated_area += (*it)->Area();
    }
  
  if constexpr ( dim == 2U ) {
      for ( typename vector<InterFace<dim>*>::const_iterator
            it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
        if ( (*it)->FE()->IsLineElement() )
          integrated_area += (*it)->Area();
    }
  
  if constexpr ( dim == 1U ) {
       // TODO: should be a static assert
       csmp_error.notice( ERROR, "SplitBoundary<dim>::Area", "not defined in 1D" );
    }

  return integrated_area;
}





/**
Integrates a scalar variable of interest over the target side of the SplitBoundary using
the current finite element basis functions. Dependent on the placement of the variable it is either
read from the InterFace or the target sidfe of the boundary.
Potential incompatibilities are checked.

@attention for simplex element node and integration point properties are integrated using their value
projected to the element barycentre.

If the variable is vector quantity, it is projected onto the outward pointing normals of the boundary InterFace
objects, multiplying the boundary normal componet with their area.

@author SKM 21/8/2018
*/
template<size_t dim>
double SplitBoundary<dim>::SurfaceIntegral( const PropertyDatabase<dim>& p, const char* property,
                                              INTERFACE_SIDE side ) const
{
  csmp::Index prop_key = p.StorageKey( property );

  if ( prop_key.type == TENSOR ) {
    throw csmp::Exception( ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                           property, "is a tensor property; no implementation for tensor normal projections yet." );
    return std::numeric_limits<double>::quiet_NaN();
  }

  if ( prop_key.place == FACE or prop_key.place == BOUNDARY or prop_key.place == REGION ) {
    throw csmp::Exception( ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                           property, "placed on FACE, BOUNDARY or REGION cannot be assigned on split boundary." );
    return std::numeric_limits<double>::quiet_NaN();
  }

  double  property_integral( 0. );

  // 1. if the property is a scalar various options exist
  // - scalar placed on SplitBoundary
  // - scalar placed on Element
  // - scalar placed on Node
  // - scalar placed on InterFace

  if ( prop_key.type == SCALAR ) {
    if ( prop_key.place == SPLIT_BOUNDARY )
      return this->Area( side ) * this->Read( prop_key );

    if ( prop_key.place == ELEMENT ) {
      // the property is read from ther higher dimensional neighbor element on the target side
      // and integrated over the area of its interface
      double property_integral( 0. ), prop_value;
      for ( auto& ife : this->elmt_vec_ ) {
        double face_area = ife->Parent( side )->FaceArea( ife->ParentFaceID( side ) );
        if ( side != MIDDLE ) prop_value = ife->Parent( side )->Read( prop_key );
        else {
          if ( ife->HasInterveningElement() ) prop_value = ife->InterveningElement()->Read( prop_key );
          else
            throw csmp::Exception( ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                                   property, "placed on Element cannot be retrieved as there is no base element at InterFace MIDDLE." );
        }
        property_integral += prop_value * face_area;
      }
      return property_integral;
    }

    if ( prop_key.place == INTER_FACE ) {
      // the property is read from ther higher dimensional neighbor element on the target side
      // and integrated over the area of its interface
      double property_integral( 0. );
      for ( auto& ife : this->elmt_vec_ )
        property_integral += ife->Area() * ife->Read( prop_key );

      return property_integral;
    }

    if ( prop_key.place == NODE ) {
      // the property value is interpolated to the interface integration points and then integrated
      // using their integration weights
      double        property_integral( 0. );
      ScalarVariable  sc;
      for ( auto& ife : this->elmt_vec_ ) {
        const double interface_area = ife->Area( side );
        for ( size_t i = 0U; i<ife->IntegrationPoints(); ++i ) {
          ife->PropertyValueAtIntegrationPoint( prop_key, i, sc );
          property_integral += interface_area * ife->WeightAtIntegrationPoint( i ) * sc();
        }
      }
      return property_integral;
    }

    if ( prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
      // the property value is integrated using corresponding integration weights
      double property_integral( 0. );
      for ( auto& ife : this->elmt_vec_ ) {
        const double interface_area = ife->Area( side );
        for ( size_t i = 0U; i<ife->IntegrationPoints(); ++i )
          property_integral += interface_area * ife->WeightAtIntegrationPoint( i ) * ife->Read( prop_key );
      }
      return property_integral;
    }

  } // end scalar
  
    // TODO: still needs extra cases" element props. etc.
    // 2. if the property is a vector
  if ( prop_key.type == VECTOR ) {
    // the values projected onto the normal are integrated over the interface.
    if ( prop_key.place == FACE or prop_key.place == INTER_FACE ) {
      VectorVariable<dim>  unrml, vc;
      for ( auto& it : this->elmt_vec_ ) {
        it->UnitNormal( unrml );
        it->Read( prop_key, vc );
        property_integral += dotProduct( unrml, vc ) * it->Area( side );
      }
    }
    // nodal properties are interpolated to the barycentre because this is where the normal is placed
    else if ( prop_key.place == NODE ) { // for nodes on first side of interface
      VectorVariable<dim>  unrml, vc;
      for ( auto& it : this->elmt_vec_ ) {
        it->UnitNormal( unrml );
        it->PropertyValueAtBaryCenter( prop_key, vc );
        property_integral += dotProduct( unrml, vc ) * it->Area( side );
      }
    }
    else {
      throw csmp::Exception( FATAL_ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                             parsePlacement( prop_key.place ), "Property placement not handled yet." );
    }
  }

  return property_integral;

} // end SurfaceIntegral


/**
Property assignment to nodes on either side of the interface or elements colocated with the InterFace objects (MIDDLE).

@todo ugly implementation where the nodes get written too many times as their side of the interface is only known to the InterFace.
*/
template<size_t dim>
template<class Var>
void SplitBoundary<dim>::InputNodePropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART part, INTERFACE_SIDE innerOuter )
{
  Index ipKey( this->pref_.StorageKey( input_prop ) );

  if ( ipKey.place != NODE )
    throw csmp::Exception( ERROR, "SplitBoundary<dim>::InputNodePropertyValue:", "This method applies to node properties only!" );

  if ( part == COMPLETE ) {
    const typename vector<InterFace<dim>*>::const_iterator ifEnd( this->ElementsEnd() );
    for ( typename vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != ifEnd; ++ifit )
      for ( size_t n( 0 ); n < (*ifit)->FE()->Nodes(); ++n )
        (*ifit)->N( n, innerOuter )->Store( ipKey, new_value );
  }
  else if ( part == INTERIOR ) {
    const typename vector<InterFace<dim>*>::const_iterator ifEnd( this->PerimeterElementsBegin() );
    for ( typename vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != ifEnd; ++ifit )
      for ( size_t n( 0 ); n < (*ifit)->FE()->Nodes(); ++n )
        (*ifit)->N( n, innerOuter )->Store( ipKey, new_value );
  }
  else if ( part == PERIMETER ) {
    const typename vector<InterFace<dim>*>::const_iterator ifEnd( this->ElementsEnd() );
    for ( typename vector<InterFace<dim>*>::const_iterator ifit( this->PerimeterElementsBegin() ); ifit != ifEnd; ++ifit )
      for ( size_t n( 0 ); n < (*ifit)->FE()->Nodes(); ++n )
        (*ifit)->N( n, innerOuter )->Store( ipKey, new_value );
  }
  else
    cerr << "\nSplitBoundary<dim>::InputNodePropertyValue: subdomain part not recognized; nothing was done.\n";

} // InputPropertyValue

template void SplitBoundary<1>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<1>::InputNodePropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputNodePropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputNodePropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<1>::InputNodePropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputNodePropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputNodePropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART, INTERFACE_SIDE );


// SCREEN OUTPUT

template<size_t dim>
void SplitBoundary<dim>::Out() const
{
  cout << "\nSplitBoundary<dim>::Out(): ";
  cout << " member interfaces: interior=" << this->InteriorElements();
  cout << ", perimeter=" << this->elmt_vec_.size() - this->InteriorElements() << ": " << endl;

  for ( typename vector<InterFace<dim>*>::const_iterator
        it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ ) {
    if ( (*it) == nullptr )
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::Out",
                             "member interface pointer not initialised" );
  }

  cout << "\n\n edge interfaces and their edges (current local numbering): " << endl;
  vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit( this->bd_face_vec_.begin() );
  for ( size_t i = this->InteriorElements(); i<this->elmt_vec_.size(); i++, bit++ ) {
    cout << "\ninterface " << i << ": edge numbers: ";
    for ( vector<ONE_BYTE_NUMBER>::const_iterator
          ft = (*bit).begin(); ft != (*bit).end(); ft++ ) cout << (*ft) << " ";
  }

  cout << "\n\n edge nodes: " << this->node_vec_.size() - this->first_bd_node_ << " (current local numbering):" << endl;
  for ( size_t i = this->first_bd_node_; i<this->node_vec_.size(); i++ ) {
    if ( this->node_vec_[i] == nullptr )
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::Out", "member node pointer not initialised." );
    else cout << this->node_vec_[i]->Idx() << " ";
  }

  cout << endl;
}


template class SplitBoundary<1>;
template class SplitBoundary<2>;
template class SplitBoundary<3>;

} // end csmp
