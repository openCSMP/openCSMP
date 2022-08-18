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
#include "PL_Utilities.h"
#include "FEM_Data.h"

#include "Visitor.h"
#include "variableOperations.h"

//#define SPLITBOUNDARY_DEBUG

using namespace std;


namespace csmp {


template<uint32_t dim>
SplitBoundary<dim>::SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>& pref )
  : ModelSubDomain<dim, InterFace>( splitboundaryname, pref )
{
  this->ResizePropertyStorage( this->pref_.LocalVariablesAt( Placement() ) );
}


template<uint32_t dim>
SplitBoundary<dim>::SplitBoundary( const SplitBoundary& ed )
  : ModelSubDomain<dim, InterFace>( ed ),
    LocalVariableStorage<dim,SplitBoundary>( ed )
{
}


template<uint32_t dim>
SplitBoundary<dim>::SplitBoundary( SplitBoundary&& ed )
  : ModelSubDomain<dim, InterFace>( move(ed) ),
    LocalVariableStorage<dim,SplitBoundary>( move(ed) )
{
}


template<uint32_t dim>
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
template<uint32_t dim>
SplitBoundary<dim>::SplitBoundary( const PropertyDatabase<dim>& pref,
                                   MeshManager<dim>& mesh,
                                   const SubDomainInfo& info )
  : ModelSubDomain<dim, InterFace>( info.name, pref )
{
  assert( info.interior_nodes.empty() );
  assert( info.perimeter_nodes.empty() );
  
  // building the interface vector (for this particular region)
  // ----------------------------------------------------------
  this->cell_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );
  const size_t n_elements_plus_faces{ mesh.Elements() + mesh.Faces() };
  for ( auto i : info.interior_elmts ) this->cell_vec_.push_back( &(*next(mesh.InterFacesBegin(),i-n_elements_plus_faces)) );
  for ( auto i : info.perimeter_elmts ) this->cell_vec_.push_back( &(*next(mesh.InterFacesBegin(),i-n_elements_plus_faces)) );

  // sorting of the pointers is necessary because the memory addresses of the new pointers will be different than in the last model
  if ( info.interior_elmts.size() == 0 ) {
       sort( this->cell_vec_.begin(), this->cell_vec_.end() );
    }
  else {
       sort( this->cell_vec_.begin(), next(this->cell_vec_.begin(),info.interior_elmts.size()) );
       sort( next(this->cell_vec_.begin(),info.interior_elmts.size()), this->cell_vec_.end() );
    }

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
template<uint32_t dim>
SplitBoundary<dim>::SplitBoundary( std::string splitboundaryname,
                                   const PropertyDatabase<dim>& pref,
                                   const FiniteElementManager& femgr,
                                   MeshManager<dim>& mesh,
                                   const InterFaceParentElements<dim>& ifset )
  : ModelSubDomain<dim, InterFace>( splitboundaryname, pref )
{
  const LocalVariables&             ifvars( pref.LocalVariablesAt( INTER_FACE ) );
  const IntegrationPointVariables&  if_ip_vars( pref.IntegrationPointVariablesAt( INTER_FACE_INTEGRATION_POINT ) );

  // 1. getting mesh manager to build interfaces according to specifications
  // -----------------------------------------------------------------------
  this->cell_vec_.reserve( ifset.size() );

  // for all the interfaces of the new split boundary
  for ( auto it = ifset.begin(); it != ifset.end(); ++it )
    // getting the element type that the interface shall represent
    // from the first higher dimensional neighbor element
    // InterFaceSet member:   pair<pair<Element<dim>*,size_t>, pair<Element<dim>*,size_t> >
    //                        first high-dim. nbor interface at interface
    // construction with connectivity and number continueing from already existing interfaces
    this->cell_vec_.push_back( mesh.AddInterFace( ifset.InnerElement(it), ifset.InnerFaceID(it),
                                                  ifset.OuterElement(it), ifset.OuterFaceID(it),
                                                  ifvars, if_ip_vars ) );

  // 2. establising interface neighbor connectivity and interior vs. perimeter includig sorting
  // ---------------------------------------------------------------------------------------------------
  mesh.BuildInterFaceConnectivity( this->cell_vec_.begin(), this->cell_vec_.end() );

  // 3. building the interface node vector
  // ---------------------------------------------------------------------------------------------------
  this->node_vec_.reserve( this->cell_vec_.size() );
  for ( auto& it : this->cell_vec_ )
    for ( auto i{0U}; i<it->Nodes(); ++i ) this->node_vec_.push_back( it->N( i ) );
  // sorting node vector and making it unique
  sort( this->node_vec_.begin(), this->node_vec_.end() );
  this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );
  
  // assigning the INTERNAL box boundary flag to nodes where they do not already have another flag (like EDGE etc)
  for ( auto& nit : this->node_vec_ )
    if ( nit->AtBoundary() == NOT )
      nit->AtBoundary( INTERNAL );
  

  // 4. sorting interfaces and nodes and building the boundary interface vector
  // ---------------------------------------------------------------------------------------------------
  this->IdentifyPerimeter();

  // 5. allocating the storage for subdomain properties
  // --------------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( SPLIT_BOUNDARY ) );

} // end constructor





/// Does not delet interfaces, Delete() has to be called for this
template<uint32_t dim>
SplitBoundary<dim>::~SplitBoundary()
{
}



template<uint32_t dim>
IntegrationPointVariables  SplitBoundary<dim>::InterFaceIntegrationPointVariables() const
{ return this->pref_.IntegrationPointVariablesAt( INTER_FACE ); }



template<uint32_t dim>
LocalVariables  SplitBoundary<dim>::InterFaceVariables() const
{ return this->pref_.LocalVariablesAt( INTER_FACE ); }


// LOCAL VARIABLE STORAGE INTERFACE
template<uint32_t dim>
bool SplitBoundary<dim>::ValidVariable( const char* variableName ) const
{
  const PLACEMENT p( this->pref_.Placement( variableName ) );
  if ( p == NODE || p == INTER_FACE || p == SPLIT_BOUNDARY )
    return true;
  return false;
}


// VISITORS INTERFACE

/// visitation of a split boundary
template<uint32_t dim>
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
      for ( auto& it : this->CellVector() ) it->Accept( v );
      return;
    case NODE:
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::Accept",
                             "Visitor application target NODE not valid for SplitBoundary; nothing was done" );
      return;
    default:
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::Accept",
                             "ApplicationTarget was not resolved; nothing was done" );
  }
} // end Accept




/// Returns the position of and adjacent region relative to the boundary. Relies on element Idx
template<uint32_t dim>
INTERFACE_SIDE SplitBoundary<dim>::RegionLocation( const Region<dim>& region )
{
  assert( !this->cell_vec_.empty() );
  const auto regionElementsEnd( region.CellsEnd() );
  const auto sbElementsEnd( this->CellsEnd() );
  for ( auto ifit( this->CellsBegin() ); ifit != sbElementsEnd; ++ifit )
    {
      const Element<dim>* const innerParent( (*ifit)->Parent( INSIDE ) );
      const Element<dim>* const outerParent( (*ifit)->Parent( OUTSIDE ) );
      for ( auto eit( region.CellsBegin() ); eit != regionElementsEnd; ++eit )
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
template<uint32_t dim>
template<typename Var>
void SplitBoundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index prop_key(this->pref_.StorageKey(input_prop));
      if ( prop_key.place == SPLIT_BOUNDARY ) {
           if ( sd != COMPLETE )
             csmp_error.Note( WARNING, "SplitBoundary<dim>::InputPropertyValue",
                              input_prop, "is a SplitBoundary property and no distinction between INTERIOR and PERIMETER can be made" );
           this->Store( prop_key, new_value );
           return;
        }
      
      // incorrect applications of method  
      if ( prop_key.place == NODE )
        csmp_error.Note( ERROR, "SplitBoundary<dim>::InputPropertyValue",
                         input_prop, "for node properties, use method InputNodePropertyValue() that allows to distinguish inside from outside" );

      if ( prop_key.place == BOUNDARY || prop_key.place == REGION || prop_key.place == MODEL ||
           prop_key.place == ELEMENT || prop_key.place == FACE )
        csmp_error.Note( ERROR, "SplitBoundary<dim>::InputPropertyValue",
                           input_prop, "must be a SPLIT_BOUNDARY, INTER_FACE/IP property for this method to work" );
    
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



template<uint32_t dim>
template<typename Var>
void SplitBoundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index key(this->pref_.StorageKey(input_prop));
      
      if ( key.place == SPLIT_BOUNDARY ) {
           if ( sd != COMPLETE )
             csmp_error.Note( WARNING, "SplitBoundary<dim>::InputPropertyValue",
                              input_prop, "is a SPLIT_BOUNDARY property and no distinction between INTERIOR and PERIMETER can be made" );
                                
           // only overwriting those variable components / rows that are not flagged 'do_not_overwrite' 
           writeVariableIf( this, key, new_value, do_not_overwrite );
           return;
        }
      
      if ( key.place == NODE )
        csmp_error.Note( ERROR, "SplitBoundary<dim>::InputPropertyValue",
                         input_prop, "for node properties, use method InputNodePropertyValue() that allows to distinguish inside from outside" );

      // incorrect applications of method
      if ( key.place == REGION  || key.place == BOUNDARY || key.place == MODEL || 
           key.place == ELEMENT || key.place == FACE )
        csmp_error.Note( ERROR, "SplitBoundary<dim>::InputPropertyValue",
                         input_prop, "must be a SPLIT_BOUNDARY or INTER_FACE/IP property for this method to work" );
    
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
    Only accumulates the inside nodes as these will later be sorted by ModelSubDomain::PartitionCellVector, i.e.  IdentifyPerimeter
*/
/*
template<uint32_t dim>
void SplitBoundary<dim>::CreateNodePointerVector()
{
  if ( this->cell_vec_.empty() )
    throw csmp::Exception( ERROR, "SplitBoundary<dim>::CreateNodePointerVector:",
                           this->Name(), "interface vector is empty; nothing could be done." );

  if ( !this->node_vec_.empty() ) this->node_vec_.clear();
    this->node_vec_.reserve( this->cell_vec_.size() );

  // creating the node index vector
  for ( auto& it : this->cell_vec_ )
    for ( auto i{0U}; i<it->FE()->Nodes(); i++ )
      this->node_vec_.push_back( it->N( i, INSIDE ) );

   //  making the vector unique and trimming of excess memory
   sort( this->node_vec_.begin(), this->node_vec_.end() );
   this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );
   this->node_vec_.shrink_to_fit();
}
*/

 /// returns a pointer to the node manifold associated with node of the split boundary; on the perimeter, a null pointer might be returned if the node is not a manifold
/*
template<uint32_t dim>
const NodeManifold<dim>* const SplitBoundary<dim>::ManifoldNode( size_t inside_node_idx ) const
 {
    return this->N( inside_node_idx )->Manifold();
 }
    
template<uint32_t dim>
NodeManifold<dim>* const SplitBoundary<dim>::ManifoldNode( size_t inside_node_idx )
 {
    return this->N( inside_node_idx )->Manifold();
 }
*/



/**
       Only those nodes that are manifolds.
*/
template<uint32_t dim>
vector<NodeManifold<dim>*> SplitBoundary<dim>::NodeManifolds() const
 {
    if ( this->cell_vec_.empty() )
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::NodeManifolds:",
                             this->Name(), "interface vector is empty; nothing could be done." );

    vector<NodeManifold<dim>*> manifold_vec;
    manifold_vec.reserve( this->Cells() );
    
    // creating the node index vector
    for ( auto& it : this->cell_vec_ )
      for ( auto i{0U}; i<it->FE()->Nodes(); i++ )
        if ( it->N(i)->IsManifold() )
          manifold_vec.push_back( it->N(i)->Manifold() );
    
    // making the vector unique
    sort( manifold_vec.begin(), manifold_vec.end() );
    manifold_vec.erase( unique( manifold_vec.begin(), manifold_vec.end() ), manifold_vec.end() );
        
    return manifold_vec;

 } // end NodeManifolds



/**
    Sorted Node Vector with perimeter nodes identified by BOX_BOUNDARY flags and non-manifold status.
*/
template<uint32_t dim>
pair<vector<Node<dim>*>,size_t>  SplitBoundary<dim>::InsideNodes() const
 {
    if ( this->cell_vec_.empty() )
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::InsideNodes:",
                             this->Name(), "interface vector is empty; nothing could be done." );
    // return vector
    vector<Node<dim>*>  inside_nodes;
    inside_nodes.reserve( this->Cells() );

    // creating the inside node vector
    for ( auto& it : this->cell_vec_ ) {
         const auto n_nodes{ it->FE()->Nodes() };
         for ( auto i{0U}; i<n_nodes; i++ )
           inside_nodes.push_back( it->N( i, INSIDE ) );
      }

     //  making the vector unique and trimming of excess memory
     sort( inside_nodes.begin(), inside_nodes.end() );
     inside_nodes.erase( unique( inside_nodes.begin(), inside_nodes.end() ), inside_nodes.end() );
     inside_nodes.shrink_to_fit();
     
     // partitioning vector into interior and perimeter nodes
     vector<Node<dim>*> interior_nodes, perimeter_nodes;
     interior_nodes.reserve( inside_nodes.size() );
     //perimeter_nodes.reserve();
     // node that sorting is retained
     for ( const auto& nit : inside_nodes )
       if ( !nit->IsManifold() || nit->AtBoundary() != NOT )
         perimeter_nodes.push_back( nit );
       else
         interior_nodes.push_back( nit );
         
     // overwriting the original vector with the identified ranges
     inside_nodes = interior_nodes;
     inside_nodes.insert( inside_nodes.end(), perimeter_nodes.begin(), perimeter_nodes.end() );
         
     return make_pair( inside_nodes, interior_nodes.size() );
     
 } // end InsideNodes


template<uint32_t dim>
pair<vector<Node<dim>*>,size_t>  SplitBoundary<dim>::OutsideNodes() const
 {
    if ( this->cell_vec_.empty() )
      throw csmp::Exception( ERROR, "SplitBoundary<dim>::OutsideNodes:",
                             this->Name(), "interface vector is empty; nothing could be done." );
    // return vector
    vector<Node<dim>*>  outside_nodes;
    outside_nodes.reserve( this->Cells() );

    // creating the inside node vector
    for ( auto& it : this->cell_vec_ ) {
         const auto n_nodes{ it->FE()->Nodes() };
         for ( auto i{0U}; i<n_nodes; i++ )
           outside_nodes.push_back( it->N( i, INSIDE ) );
      }

     //  making the vector unique and trimming of excess memory
     sort( outside_nodes.begin(), outside_nodes.end() );
     outside_nodes.erase( unique( outside_nodes.begin(), outside_nodes.end() ), outside_nodes.end() );
     outside_nodes.shrink_to_fit();
     
     // partitioning vector into interior and perimeter nodes
     vector<Node<dim>*> interior_nodes, perimeter_nodes;
     interior_nodes.reserve( outside_nodes.size() );
     //perimeter_nodes.reserve();
     // node that sorting is retained
     for ( const auto& nit : outside_nodes )
       if ( !nit->IsManifold() || nit->AtBoundary() != NOT )
         perimeter_nodes.push_back( nit );
       else
         interior_nodes.push_back( nit );
         
     // overwriting the original vector with the identified ranges
     outside_nodes = interior_nodes;
     outside_nodes.insert( outside_nodes.end(), perimeter_nodes.begin(), perimeter_nodes.end() );
         
     return make_pair( outside_nodes, interior_nodes.size() );
     
 } // end OutsideNodes



/**
Detects of how many spatial dimensions interface types are contained in model.
It returns a pair: first value gives number of different spatial dimensions contained,
second value returns the highest spatial dimension contained.

@author SKM 1/11/2013
*/
template<uint32_t dim>
pair<int32_t, int32_t>  SplitBoundary<dim>::InterFaceSpatialDimensions() const
{
  return this->SpatialDimensions();

} // end ElementSpatialDimensions




template<uint32_t dim>
size_t SplitBoundary<dim>::AccumulateByNumber( MeshManager<dim>& mesh,
                                               vector<size_t>& cell_ids )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( cell_ids.empty() )
    csmp_error.Note( ERROR, "SplitBoundary<dim>::AccumulateByNumber",
                      "user-supplied interface-number vector is empty. Nothing is done." );

  if ( !this->cell_vec_.empty() ) {
      csmp_error.Note( WARNING, "SplitBoundary<dim>::AccumulateByNumber",
                         "SplitBoundary is not empty", "erasing all members..." );
      this->cell_vec_.clear();
    }

  // eliminating potential duplicates from element index vector
#ifdef DEBUG
  const size_t n_cells{cell_ids.size()};
  sort( cell_ids.begin(), cell_ids.end() );
  cell_ids.erase( unique( cell_ids.begin(), cell_ids.end() ), cell_ids.end() );
  if ( cell_ids.size() < n_cells )
    csmp_error.Note( WARNING, "SplitBoundary<dim>::AccumulateByNumber",
                      "user-supplied interface ID set contained duplicates which were removed." );
#endif
  if ( cell_ids.size() > mesh.InterFaces() )
    csmp_error.Note( ERROR, "SplitBoundary<dim>::AccumulateByNumber",
                       "user-supplied interface-number vector is larger than range of index-to-element-pointer mapping." );

  // creating the element vector for the region
  // NB: assumes that the Faces are numbered consecutively from 0..n-1, while the supplied IDs start at the number of elements
  const auto offset = mesh.Elements() + mesh.Faces();
  this->cell_vec_.reserve( cell_ids.size() );
  for ( auto& idx : cell_ids ) {
       InterFace<dim>* ifptr = &(*next(mesh.InterFacesBegin(),idx-offset));
       assert( ifptr != nullptr );
       assert( ifptr->Idx() == idx );
       this->cell_vec_.push_back( ifptr );
    }
    
  // NB: a SplitBoundary only has a cell vector, but not a node-pointer vector
  this->IdentifyPerimeter();
  
  return this->cell_vec_.size();

} // end AccumulateByNumber


// helper function for method below
template<uint32_t dim>
static size_t countDisconnectedCells( typename vector<InterFace<dim>*>::const_iterator first,
                                      typename vector<InterFace<dim>*>::const_iterator last )
 {
     size_t disconnected_cells{0U};
     while ( first != last ) {
          if ( (*first)->ConnectedNeighbors() == 0U ) disconnected_cells++;
          first++;
       }
       
     return disconnected_cells;
  }


/**
    For post-processing the results of Divide. Here the assumption is made that the faces are already interconnected.
*/
template<uint32_t dim>
bool SplitBoundary<dim>::CreateFrom( const typename vector<InterFace<dim>*>::const_iterator ifacesBegin,
                                     const typename vector<InterFace<dim>*>::const_iterator ifacesEnd )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
  if ( distance(ifacesBegin,ifacesEnd) == 0U ) {
       csmp_error.Note( ERROR, "SplitBoundary<dim>::CreateFrom", "supplied InterFace range is empty; nothing was done");
       return false;
    }
  size_t disconnected_cells = countDisconnectedCells<dim>( ifacesBegin, ifacesEnd );
  if ( disconnected_cells > 0U ) {
       cout <<"\nSplitBoundary<dim>::CreateFrom: creating '"<< this->Name() <<"'";
       csmp_error.Note( WARNING, "SplitBoundary<dim>::CreateFrom:", to_string(disconnected_cells),
                       "InterFace object(s) do(es) not have neighbors, implying a split-boundary patch consisting of a single InterFace");
    }
    
  this->cell_vec_.assign( ifacesBegin, ifacesEnd );

  // NB: a SplitBoundary only has a cell vector, but not a node-pointer vector
  this->IdentifyPerimeter();

  return true;
}


/**
    Creates a split boundary from a boundary. 
    
    @attention this will prompt the MeshManager to delete the faces that the boundary consists of, i.e., destroy the boundary.
    
    @author SKM 1/11/2013
    @author SKM 21/9/2021
*/
/* DEPRECATED together with the possibility to create SplitBoundaries from Boundaries

template<uint32_t dim>
bool  SplitBoundary<dim>::CreateFrom( const PropertyDatabase<dim>& dbase,
                                      MeshManager<dim>& mesh,
                                      Boundary<dim>& boundary )
{
  //LVS
  const LocalVariables lvsInterFace( InterFaceVariables() );
  const IntegrationPointVariables lvsIntegrationPoint( InterFaceIntegrationPointVariables() );

  // this method already updates the connectivity of all elements, nodes etc.
  this->cell_vec_ =  mesh.ReplaceFacesByInterFaces( dbase, boundary.CellVector().begin(),
                                                    next(boundary.CellVector().begin(),boundary.InteriorCells()),
                                                    boundary.CellVector().end(),
                                                    boundary.PerimeterNodesBegin(),
                                                    boundary.NodesEnd());

  // NB: a SplitBoundary only has a cell vector, but not a node-pointer vector
  this->IdentifyPerimeter();

  return true;
  
} // CreateFrom
*/








// CALCULATIONS

/**
Computes length (m) of the SplitBoundary object's perimeter curve.
Operation makes sense only in 33 because the perimeter of a line are just its end points.
*/
template<uint32_t dim>
double  SplitBoundary<dim>::Perimeter( INTERFACE_SIDE side ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if constexpr ( dim != 3U ) {
    csmp_error.Note( WARNING, "SplitBoundary<>::Perimeter:",
                       "returning 1.0 since perimeter is a point." );
    return 1.;
  }

  double perimeter_length{0.};
  size_t n{ this->InteriorCells() };

  for ( auto it = this->PerimeterCellsBegin(); it != this->CellsEnd(); it++, n++ )
    for ( auto i{0U}; i<this->PerimeterFaces( n ); i++ ) {
      auto fnids = (*it)->FE()->NodesOfFace( this->PerimeterFace( n, i ) );
      perimeter_length += ((*it)->N( fnids[1] )->Coordinate() -
                           (*it)->N( fnids[0] )->Coordinate()).Length();
    }

  return perimeter_length;
}



/**
Is calculated on the basis of the Splitboundary bisector if the split nodes were
moved apart in the simulation process; else a particular side is used.
*/
template<uint32_t dim>
double  SplitBoundary<dim>::Area( INTERFACE_SIDE side ) const
{
  double  integrated_area( 0. );
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if constexpr ( dim == 3U ) {
      for ( auto it = this->cell_vec_.begin(); it != this->cell_vec_.end(); it++ )
        if ( (*it)->FE()->IsSurface() )
          integrated_area += (*it)->Area(side);
    }
  
  if constexpr ( dim == 2U ) {
      for ( typename vector<InterFace<dim>*>::const_iterator
            it = this->cell_vec_.begin(); it != this->cell_vec_.end(); it++ )
        if ( (*it)->FE()->IsLine() )
          integrated_area += (*it)->Area(side);
    }
  
  if constexpr ( dim == 1U ) {
       // TODO: should be a static assert
       csmp_error.Note( ERROR, "SplitBoundary<dim>::Area", "not defined in 1D" );
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
template<uint32_t dim>
double SplitBoundary<dim>::SurfaceIntegral( const PropertyDatabase<dim>& p,
                                            const char* property,
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
      double prop_value;
      for ( auto& ife : this->cell_vec_ ) {
        double face_area = ife->Volume();
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
      for ( auto& ife : this->cell_vec_ )
        property_integral += ife->Area() * ife->Read( prop_key );

      return property_integral;
    }

    if ( prop_key.place == NODE ) {
      // the property value is interpolated to the interface integration points and then integrated
      // using their integration weights
      ScalarVariable  sc;
      for ( auto& ife : this->cell_vec_ ) {
        ife->CurrentSide(side);             //needed so that property at int point calls correct fe nodes
        for ( auto i{0U}; i<ife->IntegrationPoints(); ++i ) {
          double det_j = ife->det_JINV_AtIntegrationPoint(i);
          ife->PropertyValueAtIntegrationPoint( prop_key, i, sc );
          property_integral += det_j * ife->WeightAtIntegrationPoint( i ) * sc();
        }
      }
      return property_integral;
    }

    if ( prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
      // the property value is integrated using corresponding integration weights
      for ( auto& ife : this->cell_vec_ ) {
        ife->CurrentSide(side);             //needed so that property at int point calls correct fe nodes
        for ( auto i{0U}; i<ife->IntegrationPoints(); ++i ){
          double det_j = ife->det_JINV_AtIntegrationPoint(i);
          property_integral += det_j * ife->WeightAtIntegrationPoint( i ) * ife->Read( prop_key );
        }
      }
      return property_integral;
    }

  } // end scalar
  
    // TODO: still needs extra cases" element props. etc.
    // 2. if the property is a vector
  if ( prop_key.type == VECTOR ) {
    // the values projected onto the normal are integrated over the interface.
    if ( prop_key.place == FACE or prop_key.place == INTER_FACE ) {
      VectorVariable<dim>  vc;
      for ( auto& it : this->cell_vec_ ) {
        auto unrml = it->UnitNormal();
        it->Read( prop_key, vc );
        property_integral += vc.DotProduct(unrml) * it->Area( side );
      }
    }
    // nodal properties are interpolated to the barycentre because this is where the normal is placed
    else if ( prop_key.place == NODE ) { // for nodes on first side of interface
      VectorVariable<dim>  vc;
      for ( auto& it : this->cell_vec_ ) {
        auto unrml = it->UnitNormal();
        it->PropertyValueAtBaryCenter( prop_key, vc );
        property_integral += vc.DotProduct(unrml) * it->Area( side );
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
    @attention method overrides each node several times but circumvents the need to store vectors to inside and outside nodes.
*/
template<uint32_t dim>
template<class Var>
void SplitBoundary<dim>::InputNodePropertyValue( const char* input_node_prop, const Var& var, SUBDOMAIN_PART part, INTERFACE_SIDE side )
 {
    ErrorHandler&     csmp_error( ErrorHandler::Instance() );
    const csmp::Index prop_key( this->pref_.StorageKey( input_node_prop ) );

    if ( prop_key.place != NODE ) {
         csmp_error.Note( ERROR, "SplitBoundary<dim>::InputNodePropertyValue:", "This method applies to node properties only!" );
         return;
      }
  
   // dealing with the simplest cases (=whole inside or outside first)
   if ( part == COMPLETE ) {
       for ( auto& it : this->cell_vec_ ) {
             const auto n_nodes{ it->FE()->Nodes() };
             for ( uint32_t i{0U}; i<n_nodes; i++ )
               it->N(i,side)->Store( prop_key, var );
         }
        return;
   }

   if ( part == INTERIOR ) {
        for ( auto it=this->CellsBegin(); it!=this->PerimeterCellsBegin(); ++it ) {
                const auto n_nodes{ (*it)->FE()->Nodes() };
                for ( uint32_t i{0U}; i<n_nodes; i++ )
                  (*it)->N(i,side)->Store( prop_key, var );
            }
        return;
   }

   if ( part == PERIMETER ) {
         const size_t n_sb_cells{ this->Cells() };
         for ( size_t i{ this->InteriorCells() }; i<n_sb_cells; i++ )
           for ( auto j{0U}; j < this->PerimeterFaces(i); ++j ) {
                // ascertaining that we are indeed at the model boundary
                assert( this->E(i)->Neighbor( this->PerimeterFace(i,j) ) == nullptr );
                // getting the nodes
                for ( const auto& nit : this->E(i)->FE()->NodesOfFace( this->PerimeterFace(i,j) ) )
                  this->E(i)->N(nit,side)->Store( prop_key, var );
             }

         return;
         }
      
 } // end InputNodePropertyValue

template void SplitBoundary<3U>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2U>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<3U>::InputNodePropertyValue( const char*, const VectorVariable<3U>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2U>::InputNodePropertyValue( const char*, const VectorVariable<2U>&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<3U>::InputNodePropertyValue( const char*, const TensorVariable<3U>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2U>::InputNodePropertyValue( const char*, const TensorVariable<2U>&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<3U>::InputNodePropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2U>::InputNodePropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<3U>::InputNodePropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2U>::InputNodePropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );






template<uint32_t dim>
void SplitBoundary<dim>::ChangeNodePropertyStatus( const char* property,
                                                   VARIABLE_FLAG new_status_of_scalar,
                                                   SUBDOMAIN_PART part,
                                                   INTERFACE_SIDE side )
 {
    ErrorHandler&     csmp_error( ErrorHandler::Instance() );
    const csmp::Index prop_key( this->pref_.StorageKey( property ) );

    if ( prop_key.place != NODE ) {
         csmp_error.Note( ERROR, "SplitBoundary<dim>::ChangeNodePropertyStatus:", "This method applies to node properties only!" );
         return;
      }
    if ( side == MIDDLE ) {
         csmp_error.Note( ERROR, "SplitBoundary<dim>::ChangeNodePropertyStatus:",
                         "To apply node property values to the intervening mesh, access the corresponding region" );
         return;
      }
  
   // dealing with the simplest cases (=whole inside or outside first)
   if ( part == COMPLETE ) {
        if ( side == INSIDE )
          for ( auto& it : this->cell_vec_ ) {
                const auto n_nodes{ it->FE()->Nodes() };
                for ( uint32_t i{0U}; i<n_nodes; i++ )
                  it->N(i,INSIDE)->Status( prop_key, new_status_of_scalar );
            }
        else if ( side == OUTSIDE )
          for ( auto& it : this->cell_vec_ ) {
                const auto n_nodes{ it->FE()->Nodes() };
                for ( uint32_t i{0U}; i<n_nodes; i++ )
                  it->N(i,OUTSIDE)->Status( prop_key, new_status_of_scalar );
            }
        return;
     }
     
   if ( part == INTERIOR ) {
        if ( side == INSIDE )
          for ( auto it=this->CellsBegin(); it!=this->PerimeterCellsBegin(); ++it ) {
                const auto n_nodes{ (*it)->FE()->Nodes() };
                for ( uint32_t i{0U}; i<n_nodes; i++ )
                  (*it)->N(i,INSIDE)->Status( prop_key, new_status_of_scalar );
            }
        else if ( side == OUTSIDE )
          for ( auto it=this->CellsBegin(); it!=this->PerimeterCellsBegin(); ++it ) {
                const auto n_nodes{ (*it)->FE()->Nodes() };
                for ( uint32_t i{0U}; i<n_nodes; i++ )
                  (*it)->N(i,OUTSIDE)->Status( prop_key, new_status_of_scalar );
            }
     }
   else if ( part == PERIMETER ) {
        if ( side == INSIDE ) {
             const size_t n_sb_cells{ this->Cells() };
             for ( size_t i{ this->InteriorCells() }; i<n_sb_cells; i++ )
               for ( auto j{0U}; j < this->PerimeterFaces(i); ++j ) {
                    // ascertaining that we are indeed at the model boundary
                    assert( this->E(i)->Neighbor( this->PerimeterFace(i,j) ) == nullptr );
                    // getting the nodes
                    for ( auto& nit : this->E(i)->FE()->NodesOfFace( this->PerimeterFace(i,j) ) )
                      this->E(i)->N(nit,INSIDE)->Status( prop_key, new_status_of_scalar );
                 }
          }
        else if ( side == OUTSIDE ) {
             const size_t n_sb_cells{ this->Cells() };
             for ( size_t i{ this->InteriorCells() }; i<n_sb_cells; i++ )
               for ( auto j{0U}; j < this->PerimeterFaces(i); ++j ) {
                    // ascertaining that we are indeed at the model boundary
                    assert( this->E(i)->Neighbor( this->PerimeterFace(i,j) ) == nullptr );
                    // getting the nodes
                    for ( auto& nit : this->E(i)->FE()->NodesOfFace( this->PerimeterFace(i,j) ) )
                      this->E(i)->N(nit,OUTSIDE)->Status( prop_key, new_status_of_scalar );
                 }
          }
    }
 } // end ChangeNodePropertyStatus




/**
      Where the value of the node property (length in case of a vector) is within the specified range its flag is getting changed to the new value
      but only for the INSIDE or OUTSIDE  nodes of the SplitBoundary.
*/
template<uint32_t dim>
void SplitBoundary<dim>::ChangeNodePropertyStatusWhere( const char* property,
                                                        VARIABLE_FLAG status, INTERFACE_SIDE side,
                                                        double min_val, double max_val )
  {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const csmp::Index  prop_key = this->pref_.StorageKey(property);
    string src("SplitBoundary<");
    src += to_string(dim);
    src += ",";
    src += "Interface";
    src +=">::ChangeNodePropertyStatusWhere:";

    if ( prop_key.place != NODE ) {
         csmp_error.Note( ERROR, "SplitBoundary<dim>::ChangeNodePropertyStatusWhere:", "This method applies to node properties only!" );
         return;
      }
    if ( side == MIDDLE ) {
         csmp_error.Note( ERROR, "SplitBoundary<dim>::ChangeNodePropertyStatusWhere:",
                         "To apply node property values to the intervening mesh, access the corresponding region" );
         return;
      }

    if ( prop_key.place == REGION or prop_key.place == BOUNDARY )
      throw csmp::Exception( ERROR, src.c_str(), "Use Status() to change flags of REGION or BOUNDARY variables.");

    if ( prop_key.type == TENSOR or prop_key.type == ARRAY or  prop_key.type == FLAGGEDARRAY )
      throw csmp::Exception( ERROR, src.c_str(), "Method not implemented for tensor or array properties yet");

    if ( prop_key.type == SCALAR )
      {
         for ( auto& it : this->cell_vec_ )
           for ( auto i{0U}; i<it->FE()->Nodes(); i++ ) {
                const double prop_val = it->N(i,side)->Read( prop_key );
                if ( prop_val >= min_val && prop_val <= max_val )
                  it->N(i,side)->Status( prop_key, status );
             }
         return;
      }
    else if ( prop_key.type == VECTOR ) {
         VectorVariable<dim> vc;
         for ( auto& it : this->cell_vec_ )
           for ( auto i{0U}; i<it->FE()->Nodes(); i++ ) {
                it->N(i,side)->Read( prop_key, vc );
                const double prop_val = vc.Length();
                if ( prop_val >= min_val && prop_val <= max_val )
                  it->N(i,side)->Status( prop_key, status );
             }
         return;
      }
      
    // TODO: add support for changing the status of ARRAY variables

    cout <<"\n'"<< property <<"' ";
    throw csmp::Exception( ERROR, src.c_str(), "Property placement not recognized");

 } // end ChangeNodePropertyStatusWhere







/**
    Property assignment to nodes on either side of the interface.
    Relies on the design that the node vector contains the inside nodes.
*/
template<uint32_t dim>
template<class Var>
void SplitBoundary<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART part, INTERFACE_SIDE side )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  Index          ipKey( this->pref_.StorageKey( input_prop ) );

  if ( ipKey.place != NODE ) {
       csmp_error.Note( ERROR, "SplitBoundary<dim>::InputPropertyValue:", "This method applies to node properties only!" );
       return;
    }
  if ( side == MIDDLE ) {
       csmp_error.Note( ERROR, "SplitBoundary<dim>::InputPropertyValue:",
                       "To apply node property values to the intervening mesh, access the corresponding region directly" );
       return;
    }
  
  // PERIMETER NODE INDICATORS
  // - non manifold nodes, model-boundary or internal boundary nodes (!= flagged NOT)
  if ( part == PERIMETER ) {
      for ( auto& it : this->cell_vec_ )
        for ( auto i{0U}; i<it->FE()->Nodes(); i++ )
          if ( it->N(i,side)->AtBoundary() != NOT || !it->N(i,side)->IsManifold() )
            {
               it->N(i,side)->Store( ipKey, new_value );
            }
       return;
    }
    
  // INTERIOR nodes of the SplitBoundary, always manifolds
  if ( part == INTERIOR ) {
      for ( auto& it : this->cell_vec_ )
        for ( auto i{0U}; i<it->FE()->Nodes(); i++ )
          if ( it->N(i,side)->AtBoundary() == NOT && it->N(i,side)->IsManifold() )
            {
               it->N(i,side)->Store( ipKey, new_value );
            }
       return;
    }

  // COMPLETE - all nodes but not all of them manifolds
  if ( part == COMPLETE ) {
      for ( auto& it : this->cell_vec_ )
        for ( auto i{0U}; i<it->FE()->Nodes(); i++ )
          it->N(i,side)->Store( ipKey, new_value );
     }
     
} // InputPropertyValue

template void SplitBoundary<1>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<1>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<1>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART, INTERFACE_SIDE );







/**
 *   Iterates of InterFace objects and displaces nodes by param dist/2 on both sides, opening up the interface. Leaves middle nodes untouched.
 *   Uses InterFace->UnitNormal() from bisector plane or middle element
 *
 *   @param dist is the final aperture the nodes will have. Each side is displaced by dist/2.
 *
 */
template<uint32_t dim>
void SplitBoundary<dim>::PullApartSplitBoundary(double dist){

  std::set<Node<dim>*> operated_in_nodes;
  std::set<Node<dim>*> operated_out_nodes;

  for ( auto& ifp : this->CellVector() ){
    uint32_t n_nodes = ifp->FE()->Nodes();
    for (uint32_t n{0U}; n<n_nodes;++n){
      Node<dim>* in_node  = ifp->MatchingN(n,INSIDE);
      Node<dim>* out_node = ifp->MatchingN(n,OUTSIDE);
      if (in_node != out_node){
        if ( operated_in_nodes.find(in_node) == operated_in_nodes.end() ){
          //displace inside node
          in_node->Coordinate(  in_node->Coordinate()  - dist/2.0 * ifp->UnitNormal()  );
          //adding nodes to operated nodes
          operated_in_nodes.insert(  in_node );
        }
        if (operated_out_nodes.find(out_node) == operated_out_nodes.end() ){
          //displace outside node
          out_node->Coordinate( out_node->Coordinate() + dist/2.0 * ifp->UnitNormal()  );
          //insert into operated outside nodes
          operated_out_nodes.insert( out_node );
        }//end of if
      }
    }//end of node loop
  }//end of interface loop


}//end of PullApartSplitBoundary










// SCREEN OUTPUT

template<uint32_t dim>
void SplitBoundary<dim>::Out() const
{
  ErrorHandler& csmp_err( ErrorHandler::Instance() );
  cout <<"\nSplitBoundary<dim>::Out: '"<< this->Name()<<"'";
  cout <<"\n\t"<<"InterFace objects: interior: "<< this->InteriorCells() <<", perimeter: "<< this->PerimeterCells();
  // cout <<"\n\t"<<"INSIDE Node objects: interior: " << this->InteriorNodes() <<", perimeter: "<< this->PerimeterNodes() << endl;

  size_t iface_idx{0};
  double geom_measure{0};
  for ( const auto& it : this->cell_vec_ ) {
        if ( it == nullptr ) {
             cerr <<" interface pointer "<< iface_idx <<" not valid.";
             csmp_err.Note( ERROR, "SplitBoundary<dim>::Out", "'nullptr' detected" );
          }
        else {
             geom_measure += it->Area();
          }
       iface_idx++;
    }
  if constexpr( dim == 2U ) cout <<"\n\t"<<"split boundary length: " << geom_measure << endl;
  if constexpr( dim == 3U ) cout <<"\n\t"<<"split boundary area: " << geom_measure << endl;

  cout <<"\n\t"<<"perimeter InterFace objects and their boundary face indices (current numbering): " << endl;
  auto  bit( this->bd_face_vec_.begin() );
  for ( auto i = this->InteriorCells(); i<this->cell_vec_.size(); i++, bit++ ) {
      cout << "\n\t\t"<<"interface "<< i <<": edge numbers: ";
      for ( auto ft = (*bit).begin(); ft != (*bit).end(); ft++ ) cout << (*ft) << " ";
    }

  auto manifolds = NodeManifolds();
  cout <<"\n\n\t"<<"NodeManifold objects: "<< manifolds.size() <<":";
  cout <<"\n\n\t";
  for ( const auto& i : manifolds ) {
       i->Out();
    }
  cout << endl;
  
} // end Out



/**
   Special version without nodes.
*/
template<uint32_t dim>
void SplitBoundary<dim>::WriteIndexesToBinaryFile( fstream& fp ) const
 {
    // 1. writing name of the region
    binaryFileWrite( fp, this->Name().c_str() );
   
    // 2. writing the interior cell records of the split boundary
    std::vector<uint32_t> IDs( distance(this->CellsBegin(), this->PerimeterCellsBegin() ) );
    transform( this->CellsBegin(), this->PerimeterCellsBegin(),
               IDs.begin(), []( const InterFace<dim>* const ptr ){ return ptr->Idx(); } ); // tested: OK
    binaryFileWrite( fp, IDs );

    // 3. writing the perimeter cell records of the split boundary
    IDs.resize( distance(this->PerimeterCellsBegin(), this->CellsEnd()) );
    transform( this->PerimeterCellsBegin(), this->CellsEnd(),
               IDs.begin(), []( const InterFace<dim>* const ptr ){ return ptr->Idx(); } );
    binaryFileWrite( fp, IDs );
   
    // NB: the connectivity between the cells is not stored because it is handled by MeshManager
   
 } // end WriteDomainIndexesToBinaryFile



void readIndexesFromBinaryFile( uint32_t dim, fstream& fp, SubDomainInfo& info )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    // 1. reading name of the subdomain
    char name[INFO_STRING];
    binaryFileRead( fp, name );
    info.name = name;
    assert( !info.name.empty() );
   
    // 2. reading the interior cell records of the region
    binaryFileRead( fp, info.interior_elmts );
    if (dim > 2 && info.interior_elmts.empty() ) {
        csmp_error.Note( WARNING, "readDomainIndexesFromBinaryFile:",
                          "Model appears to have a region with no interior cells: ", name );
    }

    // 3. reading the perimeter cell records of the region
    binaryFileRead( fp, info.perimeter_elmts );
    assert( !info.perimeter_elmts.empty() );
   
 } // end readRegionIndexesFromBinaryFile



template class SplitBoundary<1>;
template class SplitBoundary<2>;
template class SplitBoundary<3>;

} // end csmp
