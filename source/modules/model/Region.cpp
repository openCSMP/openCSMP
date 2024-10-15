#include "Region.h"

#include "VSet.h"
#include "PropertyData.h"
#include "FEM_Data.h"
#include "LocalVariableStorage.h"

#include "writeVariableIf.h"
#include "Node.h"
#include "MeshManager.h"
#include "Element.h"
#include "Boundary.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"

#include "Visitor.h"
#include "Interrelation.h"
#include "PropertyDatabase.h"
#include "PropertyConstraints.h"
#include "CopyReplaceVisitor.h"

#include "UnionFind.h"

#include "ErrorHandler.h"
#include "meshManagementUtilities.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
Region<dim>::Region( std::string regionname, const PropertyDatabase<dim>& p )
  : ModelSubDomain<dim, Element>( regionname, p )
{
  this->ResizePropertyStorage( this->pref_.LocalVariablesAt( Placement() ) );
}


template<uint32_t dim>
Region<dim>::Region( const Region& g )
  : ModelSubDomain<dim, Element>( g ),
    LocalVariableStorage<dim,Region>(g)
{
}


template<uint32_t dim>
Region<dim>::Region( Region&& g )
  : ModelSubDomain<dim, Element>( std::move(g) ),
    LocalVariableStorage<dim,Region>( std::move(g) )
{
}


/**
Assignment does not affect the const reference to 'pref'.
It remains the same as that before the assignment
because it refers to the one Model in which all the Regions live.
*/
template<uint32_t dim>
Region<dim>&  Region<dim>::operator=( const Region& g )
{
  if ( &g != this ) {
     ModelSubDomain<dim,Element>::operator=( g ); 
     this->LVS( g.LVS() );
  }
  return *this;
} // end assignment




template<uint32_t dim>
Region<dim>::~Region()
{
}




/**
    RECONSTRUCTOR
    
for regions that were stored in the CSMP nativefile  format.

Using the connectivity indices for the Region as retrieved from the CSMP++ native binary file and stored in SubDomainInfo,
the regions are recreated.

@attention the numbering that is provided through the domain info
must match the current state of the MeshManager.

@attention the indices are just used to connect the element and node
pointers of the region to the mesh storage, but, due to the sorting that
this method will perform, they will not correspond to the indices
used by the region accessors N() and E().

@attention this constructor initialises the data of the ModelSubDomain
from which the region is derived. This cannot be done by the subdomain
itself because it does not know whether it will consist of elements,
faces or interfaces.
*/
template<uint32_t dim>
Region<dim>::Region( const PropertyDatabase<dim>& pref,
                     MeshManager<dim>& mesh, ///< not constant because region shall later be able to modify elements and nodes
                     const SubDomainInfo& info )   ///< information on how to connect pointers to mesh stored in MeshManager 
  : ModelSubDomain<dim, Element>( info.name, pref )
{
  // building the element vector
  // ---------------------------
  this->cell_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );
  if ( this->cell_vec_.capacity() > mesh.Elements() )
    throw csmp::Exception( ERROR, "Region(custom constructor)",
                          "attempt to contruct region with more elements than are contained in MeshManager");

  // pushing back the pointers to the interior elements
  for ( size_t i : info.interior_elmts )
    this->cell_vec_.push_back( &(*next(mesh.ElementsBegin(),i)) );

  // assigning pointers to the perimeter elements
  for ( size_t i : info.perimeter_elmts )
    this->cell_vec_.push_back( &(*next(mesh.ElementsBegin(),i)) );

  // building the node vector
  // ------------------------
  this->first_bd_node_ = info.interior_nodes.size();
  this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );

  // assigning pointers to the interior nodes
  for ( size_t i : info.interior_nodes )
    this->node_vec_.push_back( &(*next(mesh.NodesBegin(),i)) );

  // assigning pointers to the perimeter nodes
  for ( size_t i : info.perimeter_nodes )
    this->node_vec_.push_back( &(*next(mesh.NodesBegin(),i)) );

  this->SortVectors( info.interior_elmts.size(), info.interior_nodes.size() );

  // building the vector of vectors of those faces of the elements that lie on the subdomain perimeter
  // -------------------------------------------------------------------------------------------------
  this->BuildPerimeterFaceVector( info.interior_elmts.size() );

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( REGION ) );

} // end region re-constructor (using MeshManager)







/**
      re-constructor for regions via the nodes and elements which are explored by the MeshManager
 */
template<uint32_t dim>
Region<dim>::Region( const PropertyDatabase<dim>& pref,
                     const deque<Node<dim>*>& nodes,
                     const deque<Element<dim>*>& elmts,
                     const SubDomainInfo& info )  ///< contains correctly partitioned vectors and boundary faces
  : ModelSubDomain<dim, Element>( info.name, pref )
{
  // building the element vector
  // ---------------------------
  this->cell_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );

  // assigning pointers to the interior elements
  for ( size_t i : info.interior_elmts )
    this->cell_vec_.push_back( elmts[i] );

  // assigning pointers to the perimeter elements
  for ( size_t i : info.perimeter_elmts )
    this->cell_vec_.push_back( elmts[i] );
    
  // building the node vector
  // ------------------------
  this->first_bd_node_ = info.interior_nodes.size();
  this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );

  // assigning pointers to the interior nodes
  for ( size_t i : info.interior_nodes )
    this->node_vec_.push_back( nodes[i] );

  // assigning pointers to the perimeter nodes
  for ( size_t i : info.perimeter_nodes )
    this->node_vec_.push_back( nodes[i] );

  this->SortVectors( info.interior_elmts.size(), info.interior_nodes.size() );

  // building the vector of vectors of those faces of the elements that lie on the subdomain perimeter
  // -------------------------------------------------------------------------------------------------
  this->BuildPerimeterFaceVector( info.interior_elmts.size() );

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( REGION ) );

} // end region re-constructor (using the nodes and the elements from the MeshManager)


// LOCAL VARIABLE STORAGE INTERFACE
template<uint32_t dim>
bool Region<dim>::ValidVariable( const char* variableName ) const
{
  const PLACEMENT p( this->pref_.Placement( variableName ) );
  if ( p == NODE || p == ELEMENT || p == REGION )
    return true;
  return false;
}



template<uint32_t dim>
IntegrationPointVariables Region<dim>::ElementIntegrationPointVariables() const
{ return this->pref_.IntegrationPointVariablesAt( ELEMENT ); }

template<uint32_t dim>
LocalVariables Region<dim>::ElementVariables() const
{ return this->pref_.LocalVariablesAt( ELEMENT ); }



/**
    Sets the Region ID to ModelSubDomain::domain_idx_ for all elements of the region.
    For unique non-spatially overlapping regions (including lower-dim node sharing regions) the value is used as is.
    
    @attention This mechanism does not apply to Face or InterFace objects since they have no region ID.
    
    @attention For non-unique regions,  no assignment is made!
*/
template<uint32_t dim>
void Region<dim>::SetRegion_ID( int32_t region_id )
{
   // only if a negative number is inserted, the reference-counted domain index is not used
   region_id = ( region_id < 0 ) ? this->DomainIndex() : region_id;
   for ( auto& it : this->cell_vec_ ) it->Region_ID( region_id );
}



// VISITORS INTERFACE

/**
Recursive visitation of regions (or any other object) starts at application level
and finishes at application target.

@todo SKM provide methods for the new variable placements.

@todo SKM refactor the loops using constants and pre-incrementation to speed them up

@author SKM
*/
template<uint32_t dim>
void Region<dim>::Accept( csmp::Visitor<dim>& v )
{
  if ( v.ApplicationLevel() == MODEL or
       v.ApplicationLevel() == REGION )
    v.Visit( this ); // SKM_FIX else visitor will never be applied to region

  switch ( v.ApplicationTarget() ) {
    case MODEL:
      throw csmp::Exception( ERROR, "Region<dim>::Accept",
                             "ApplicationTarget MODEL; Visitor should have never arrived at this region" );
      break;
    case REGION:
      // nothing else needs to be done
      return;
      // element, face and interface are treated the same
    case ELEMENT:
      for ( auto it = Region<dim>::CellsBegin(); it != Region<dim>::CellsEnd(); it++ )
        (*it)->Accept( v );
      return;
    case NODE:
      for ( auto nd_it = Region<dim>::NodesBegin(); nd_it != Region<dim>::NodesEnd(); nd_it++ )
        (*nd_it)->Accept( v );
      return;
    default:
      throw csmp::Exception( ERROR, "Region<dim>::Accept",
                             "ApplicationTarget was not resolved; nothing was done" );
  }

} // end Accept






// -----------------------------------------------
// Binary and VSet input/output
// -----------------------------------------------


/**
   for the assignment of properties that are unique to the instance of this subclass
   
      @author SKM
      @date 7/6/2020
*/
template<uint32_t dim>
template<typename Var>
void Region<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index prop_key(this->pref_.StorageKey(input_prop));
      if ( prop_key.place == REGION ) {
           if ( sd != COMPLETE )
             csmp_error.Note( WARNING, "Region<dim>::InputPropertyValue",
                                          input_prop, "is a Region property and no distinction between INTERIOR and PERIMETER can be made" );
           this->Store( prop_key, new_value );
           return;
        }
      
      // incorrect applications of method  
      if ( prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY || prop_key.place == MODEL ||
           prop_key.place == FACE || prop_key.place == INTER_FACE )
        csmp_error.Note( ERROR, "Region<dim>::InputPropertyValue",
                           input_prop, "must be a REGION, ELEMENT/IP or NODE property for this method call to work" );        
    
     // for any different property placement, the method of the base-class is called
     ModelSubDomain<dim,Element>::InputPropertyValue( input_prop, new_value, sd );
      
  } // end InputPropertyValue

template void Region<1U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const VectorVariable<1U>&, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const VectorVariable<2U>&, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const VectorVariable<3U>&, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const TensorVariable<1U>&, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const TensorVariable<2U>&, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const TensorVariable<3U>&, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );




template<uint32_t dim>
template<typename Var>
void Region<dim>::InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd )
  {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
      const csmp::Index key(this->pref_.StorageKey(input_prop));
      
      if ( key.place == REGION ) {
           if ( sd != COMPLETE )
             csmp_error.Note( WARNING, "Region<dim>::InputPropertyValue",
                                input_prop, "is a REGION property and no distinction between INTERIOR and PERIMETER can be made" );
                                
           // only overwriting those variable components / rows that are not flagged 'do_not_overwrite' 
           writeVariableIf( this, key, new_value, do_not_overwrite );
           return;
        }
      
      // incorrect applications of method  
      if ( key.place == BOUNDARY  || key.place == SPLIT_BOUNDARY || key.place == MODEL || 
           key.place == FACE || key.place == INTER_FACE )
        csmp_error.Note( ERROR, "Region<dim>::InputPropertyValue",
                           input_prop, "must be a REGION, ELEMENT/IP or NODE property for this method call to work" );        
    
     // for any different property placement, the method of the base-class is called
     ModelSubDomain<dim,Element>::InputPropertyValue( input_prop, new_value, do_not_overwrite, sd );
      
  } // end InputPropertyValue

// explicit instantiations
template void Region<1U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const VectorVariable<1U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const VectorVariable<2U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const VectorVariable<3U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const TensorVariable<1U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const TensorVariable<2U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const TensorVariable<3U>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<1U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<2U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void Region<3U>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );





/**
Outputs geometry and property data from the region to a vset.
The data of the region is copied as is. Boundary conditions will be copied if present.
For unresolved boundaries, the flag IRREGULAR is set.

@attention Relies on consecutively numbered indices

@section implementation Implementation

All data is copied into the vset. IDs are used in the current state.

@section application Application

The idea is to create a vset that contains all the information of the region.
This vset can, for example, later be used to create a new Model.
*/
template<uint32_t dim>
void Region<dim>::OutputTo( VSet<dim>& vset, bool with_properties ) const
{
  if ( this->cell_vec_.empty() ) {
    throw csmp::Exception( ERROR, "Region<dim>::OutputTo",
                           "Attempt to output empty group to VSet. Nothing was done." );
    return;
  }

  if ( vset.Vertices() > 0U ) vset.Erase();

  // 0. resizing the VSet
  size_t         counter( 0U );
  deque<uint32_t>  nodes_per_element( this->cell_vec_.size() );
  deque<uint32_t>  elements_per_element( this->cell_vec_.size() );
  deque<int8_t>  etypes( this->cell_vec_.size() );
  set<int8_t>    n_etypes;

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
  {
    n_etypes.insert( (*eit)->FE_Type() );
    etypes[counter] = (*eit)->FE_Type();
    nodes_per_element[counter] = (*eit)->Nodes();
    elements_per_element[counter] = (*eit)->Neighbors();
    counter++;
  }

  if ( n_etypes.size() == 1U ) etypes.resize( 1U );

  vset.Resize( etypes, nodes_per_element, elements_per_element, this->node_vec_.size(), 0, 0 );

  nodes_per_element.clear();
  elements_per_element.clear();

  cout << "\nOutput of Geometry...";
  this->UpdateMemberIndexes();

  // 1. assigning coordinate values
  for ( size_t i{0U}; i<this->node_vec_.size(); i++ ) {
    vset.Px( i, this->node_vec_[i]->x() );
    if ( dim != 1U ) vset.Py( i, this->node_vec_[i]->y() );
    if ( dim == 3U ) vset.Pz( i, this->node_vec_[i]->z() );
  }

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
    // mapping global to local node numbers
    for ( auto i{0U}; i<(*eit)->Nodes(); i++ ) vset.Plist( (*eit)->Idx(), i, (*eit)->N( i )->Idx() );

  // 3. writing the neighbor per element map (pfverts)
  counter = 0U;
  for ( typename vector<csmp::Element<dim>*>::const_iterator
        eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ ) {
    for ( auto i{0U}; i<(*eit)->Neighbors(); i++ )
      if ( (*eit)->Neighbor( i ) != NULL )
        vset.Pfvert( counter, i, static_cast<int32_t>((*eit)->Neighbor( i )->Idx()) );
      else {
        vset.Pfvert( counter, i, -1 );
      }
      counter++;
  }

  // 4. writing nodal boundary flags
  vector<std::int8_t>  bflags( this->PerimeterNodes(), IRREGULAR_OUTSIDE );
  for ( size_t i=this->first_bd_node_; i<this->node_vec_.size(); i++ )
    bflags[i] = this->node_vec_[i]->AtBoundary();

  vset.AddBFlags( bflags.begin(), bflags.end() );
  bflags.erase( bflags.begin(), bflags.end() );

  // 5. output of material properties
  if ( !with_properties ) {
      cout << "\n\nRegion<" << dim << ">::OutputTo: ";
      cout << "region successfully output to VSet." << endl;
      return;
    }
  else
    {
      OutputDataTo( vset );
      OutputFvDataTo( vset );
    }
  cout << "\n\nRegion<" << dim << ">::OutputTo: ";
  cout << "region successfully output to VSet." << endl;

} // end OutputTo





template<uint32_t dim>
void Region<dim>::OutputDataTo( VSet<dim>& vset ) const
{
  map<string, Index>  properties;
  this->pref_.ListVariables( NODE, properties );
  map<string, Index>  propertiesElement;
  this->pref_.ListVariables( ELEMENT, propertiesElement );
  map<string, Index>  propertiesElementIP;
  this->pref_.ListVariables( ELEMENT_INTEGRATION_POINT, propertiesElementIP );
  properties.insert( propertiesElement.begin(), propertiesElement.end() );
  properties.insert( propertiesElementIP.begin(), propertiesElementIP.end() );

  OutputTo( vset, properties );
}





template<uint32_t dim>
void Region<dim>::OutputFvDataTo( VSet<dim>& vset ) const
{
  map<string, Index>  properties;
  map<string, Index>  propertiesElementSEIP;
  map<string, Index>  propertiesElementFAIP;
  this->pref_.ListVariables( SECTOR_INTEGRATION_POINT, propertiesElementSEIP );
  properties.insert( propertiesElementSEIP.begin(), propertiesElementSEIP.end() );
  this->pref_.ListVariables( FACET_INTEGRATION_POINT, propertiesElementFAIP );
  properties.insert( propertiesElementFAIP.begin(), propertiesElementFAIP.end() );

  OutputTo( vset, properties );
}





/**
Writes the properties whose names are supplied via a map
into the input VSet container.

@author SK 2016
*/
template<uint32_t dim>
void Region<dim>::OutputTo( VSet<dim>& vset, const map<string, Index>& properties ) const
{
  for ( map<string, Index>::const_iterator
        it = properties.begin(); it != properties.end(); it++ )
  {
    PropertyData data = OutputVariableTo( (*it).first.c_str() );
    vset.AddData( (*it).first.c_str(), data );
  }
}






/**
Any variable is moved to inline storage in the provided container.

@author SKM 2016
*/
template<uint32_t dim>
PropertyData  Region<dim>::OutputVariableTo( const char* property ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->pref_.IsDefined( property ) )
    throw csmp::Exception( ERROR, "Region<dim>::OutputVariableTo:", property, "is undefined." );

  const csmp::Index key = this->pref_.StorageKey( property );
  assert( key.place != FACE );
  assert( key.place != INTER_FACE );

  // creating the property storage
  PropertyData  data( key.place, key.type, dim, key.dataDepth );

  // reserving sufficient data memory
  if ( isPlacedOnIntegrationPoint( key.place ) ) {
  }
  else { //
    if ( key.place == REGION ) data.Reserve( 1U, 1U );
    else if ( key.place == NODE ) data.Reserve( this->Nodes(), this->Nodes() );
    else if ( key.place == ELEMENT ) data.Reserve( this->Cells(), this->Cells() );
    else
      csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:", property, "output is not handled yet." );
  }

  // output
  if ( !isPlacedOnIntegrationPoint( key.place ) ) {
    // Region variable
    if ( key.place == REGION ) {
      switch ( key.type ) {
        case SCALAR: {
          ScalarVariable sc;
          this->Read( key, sc );
          pushBack( data, sc );
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> vc;
          this->Read( key, vc );
          pushBack( data, vc );
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> ts;
          this->Read( key, ts );
          pushBack( data, ts );
        }
                     break;
        case ARRAY: {
          ArrayVariable av;
          this->Read( key, av );
          pushBack( data, av );
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable fa;
          this->Read( key, fa );
          pushBack( data, fa );
        }
                           break;
        default:
          csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
                             property, "type of variable not recognized." );
      }
    }
    else if ( key.place == ELEMENT ) {
      data.Reserve( this->Cells() );
      switch ( key.type ) {
        case SCALAR: {
          ScalarVariable sc;
          for ( const auto& it : this->cell_vec_ ) {
            (*it).Read( key, sc );
            pushBack( data, sc );
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> vc;
          for ( const auto& it : this->cell_vec_ ) {
            (*it).Read( key, vc );
            pushBack( data, vc );
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> ts;
          for ( const auto& it : this->cell_vec_ ) {
            (*it).Read( key, ts );
            pushBack( data, ts );
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable av;
          for ( const auto& it : this->cell_vec_ ) {
            (*it).Read( key, av );
            pushBack( data, av );
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable fa;
          for ( const auto& it : this->cell_vec_ ) {
            (*it).Read( key, fa );
            pushBack( data, fa );
          }
        }
                           break;
        default:
          csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
                             property, "type of variable not recognized." );
      }
    }
    else if ( key.place == NODE ) {
      data.Reserve( this->Nodes() );
      switch ( key.type ) {
        case SCALAR: {
          ScalarVariable sc;
          for ( const auto& it : this->node_vec_ ) {
            (*it).Read( key, sc );
            pushBack( data, sc );
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> vc;
          for ( const auto& it : this->node_vec_ ) {
            (*it).Read( key, vc );
            pushBack( data, vc );
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> ts;
          for ( const auto& it : this->node_vec_ ) {
            (*it).Read( key, ts );
            pushBack( data, ts );
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable av;
          for ( const auto& it : this->node_vec_ ) {
            (*it).Read( key, av );
            pushBack( data, av );
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable fa;
          for ( const auto& it : this->node_vec_ ) {
            (*it).Read( key, fa );
            pushBack( data, fa );
          }
        }
                           break;
        default:
          csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
                             property, "type of variable not recognized." );
      }
    }
    // TODO: Element integration points etc
    // finite volume-related data
    // etc.
    else
      csmp_error.Note( ERROR, property, "Region<dim>::OutputVariableTo:", "variable type not handled yet." );

  } // end not integration point

  return data;

} // end OutputVariableTo(PropertyData)












/**
Outputs specific property data to a FEM_Data container.
This function is a nested template:
the outer template provides double = data type and dim = dimension,
and Var the data type of the property
(ScalarVariable, VectorVariable, or TensorVariable).

@note The result is returned into FEM_Data container with the data that corresponds to the given property.

@section implementation Implementation

Only data pertinent to this local group is output.
A new data container containing only group relevant data is returned.

Thus, the data extracted only corresponds to -this- group and is organized
by local id numbers.

@section application Application

The idea is to extract the property data that corresponds to this group only.
*/
template<uint32_t dim>
template<class Var>
void Region<dim>::OutputVariableTo( const char* property, FEM_Data<Var>& data ) const
{
  csmp::Index  idx = this->pref_.StorageKey( property );
  Var          var;
  femDataOutputDispatch::initVariable( idx, var );

  this->UpdateMemberIndexes();

  switch ( idx.place ) {
    case ELEMENT: {
      data.Reset( idx, this->cell_vec_.size(), var );
      for ( const auto& eit : this->cell_vec_ ) {
        eit->Read( idx, var );
        data[ eit->Idx() ] = var;
      }
    }
                  break;
    case NODE: {
      data.Reset( idx, this->node_vec_.size(), var );
      for ( const auto& nit : this->node_vec_ ) {
        nit->Read( idx, var );
        data[nit->Idx()] = var;
      }
    }
               break;
    case ELEMENT_INTEGRATION_POINT: {
      data.Reset( idx, this->IntegrationPoints(), var );
      size_t counter( 0U );
      for ( const auto& eit : this->cell_vec_ )
        for ( auto i{0U}; i<eit->IntegrationPoints(); i++ ) {
          eit->Read( i, idx, var );
          data[counter] = var;
          ++counter;
        }
    }
      break;
    case SECTOR_INTEGRATION_POINT: {
      data.Reset( idx, this->SectorIntegrationPoints(), var );
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
        for ( auto j{0U}; j<(*eit)->Sectors(); j++ )
        {
          for ( auto i{0U}; i<(*eit)->IntegrationPointsPerSector(); i++ ) {
            (*eit)->Read( j, i, idx, var );
            data[counter] = var;
            ++counter;
          }
        }
    }
                                   break;
    case FACET_INTEGRATION_POINT: {
      data.Reset( idx, this->FacetIntegrationPoints(), var );
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
        for ( auto j{0U}; j<(*eit)->Facets(); j++ )
        {
          for ( auto i{0U}; i<(*eit)->IntegrationPointsPerFacet(); i++ ) {
            (*eit)->Read( j, i, idx, var );
            data[counter] = var;
            ++counter;
          }
        }
    }
                                  break;
    case REGION: {
      data.Reset( idx, 1U, var );
      this->Read( idx, var );
      data[0U] = var;
    }
                 break;
    default: {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );
      csmp_error.Note( WARNING, "Region::OutputVariableTo",
                         property, "placement not identified" );
    }
  }

} // end OutputVariableTo


template void Region<1>::OutputVariableTo<ScalarVariable >( const char*, FEM_Data<ScalarVariable >& ) const;
template void Region<1>::OutputVariableTo<VectorVariable<1U> >( const char*, FEM_Data<VectorVariable<1U> >& ) const;
template void Region<1>::OutputVariableTo<TensorVariable<1U> >( const char*, FEM_Data<TensorVariable<1U> >& ) const;
template void Region<1>::OutputVariableTo<ArrayVariable>( const char*, FEM_Data<ArrayVariable>& ) const;
template void Region<1>::OutputVariableTo<FlaggedArrayVariable>( const char*, FEM_Data<FlaggedArrayVariable>& ) const;

template void Region<2>::OutputVariableTo<ScalarVariable >( const char*, FEM_Data<ScalarVariable >& ) const;
template void Region<2>::OutputVariableTo<VectorVariable<2U> >( const char*, FEM_Data<VectorVariable<2U> >& ) const;
template void Region<2>::OutputVariableTo<TensorVariable<2U> >( const char*, FEM_Data<TensorVariable<2U> >& ) const;
template void Region<2>::OutputVariableTo<ArrayVariable>( const char*, FEM_Data<ArrayVariable>& ) const;
template void Region<2>::OutputVariableTo<FlaggedArrayVariable>( const char*, FEM_Data<FlaggedArrayVariable>& ) const;

template void Region<3>::OutputVariableTo<ScalarVariable >( const char*, FEM_Data<ScalarVariable >& ) const;
template void Region<3>::OutputVariableTo<VectorVariable<3U> >( const char*, FEM_Data<VectorVariable<3U> >& ) const;
template void Region<3>::OutputVariableTo<TensorVariable<3U> >( const char*, FEM_Data<TensorVariable<3U> >& ) const;
template void Region<3>::OutputVariableTo<ArrayVariable>( const char*, FEM_Data<ArrayVariable>& ) const;
template void Region<3>::OutputVariableTo<FlaggedArrayVariable>( const char*, FEM_Data<FlaggedArrayVariable>& ) const;

/**
Inputs geometry and property data from an FEM_Data object into this region.

The expectation is that nodes and element numbers in the VData set
correspond to these in the current Region and that elements / nodes
in the model have the same ordering as in the input dataset!

@section arguments Input Arguments

The name of the property which shall be input to VSet.

@section implementation Implementation

The current element and node indices are used to retrieve the data
from from the data container.

@section application Application

Enable to store data to disk to faciliate persistance.
*/
template<uint32_t dim>
template<class Var>
void Region<dim>::InputVariableFrom( const char* property,
                                     const FEM_Data<Var>& vdata )
{
  const csmp::Index  idx = this->pref_.StorageKey( property );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( vdata.Size() == 0 )
  {
    string warningMessage( "Empty FEM_Data provided for " );
    warningMessage.append( property );
    warningMessage.append( ". No data imported." );
    csmp_error.Note( WARNING, "Region<dim>::InputVariableFrom:", warningMessage.c_str() );
    return;
  }

  switch ( idx.place ) {
    case ELEMENT:
      assert( this->Cells() == vdata.Size() );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
        (*eit)->Store( idx, vdata[(*eit)->Idx()] );
      break;
    case NODE:
      assert( this->Nodes() == vdata.Size() );
      for ( typename vector<csmp::Node<dim>*>::const_iterator
            nit = this->node_vec_.begin(); nit != this->node_vec_.end(); nit++ )
        (*nit)->Store( idx, vdata[(*nit)->Idx()] );
      break;
    case ELEMENT_INTEGRATION_POINT: {
      // this will work only if the element ordering in the Region
      // from which the integration point data were written out
      // using OutputDataTo() is exactly the same as in this Region.
      size_t  counter( 0U );
      assert( this->IntegrationPoints() == vdata.Size() );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
      {
        for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
          (*eit)->Store( i, idx, vdata[counter + i] );
        counter += (*eit)->IntegrationPoints();
      }
    }
                                    break;
    case SECTOR_INTEGRATION_POINT: {
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
      {
        for ( auto j{0U}; j<(*eit)->Sectors(); j++ )
        {
          for ( auto i{0U}; i<(*eit)->IntegrationPointsPerSector(); i++ ) {
            (*eit)->Store( j, i, idx, vdata[counter] );
            ++counter;
          }
        }
      }
    }
                                   break;
    case FACET_INTEGRATION_POINT: {
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->cell_vec_.begin(); eit != this->cell_vec_.end(); eit++ )
      {
        for ( auto j{0U}; j<(*eit)->Facets(); j++ )
        {
          for ( auto i{0U}; i<(*eit)->IntegrationPointsPerFacet(); i++ ) {
            (*eit)->Store( j, i, idx, vdata[counter] );
            ++counter;
          }
        }
      }
    }
                                  break;
    case REGION:
      if ( vdata.Size() != 1U )
        csmp_error.Note( WARNING, "Region<dim>::InputVariableFrom:",
                           "don't know which location in FEM_Data I should write the region variable to, using [0]." );

      this->Store( idx, vdata[0U] );
      break;
    default:
      throw csmp::Exception( ERROR, "Region<dim>::InputVariableFrom:",
                             property, "Property placement could not be identified" );
  }

} // end InputVariableFrom


template void Region<1>::InputVariableFrom<ScalarVariable>( const char*, const FEM_Data<ScalarVariable>& );
template void Region<2>::InputVariableFrom<ScalarVariable>( const char*, const FEM_Data<ScalarVariable>& );
template void Region<3>::InputVariableFrom<ScalarVariable>( const char*, const FEM_Data<ScalarVariable>& );


template void Region<1>::InputVariableFrom<ArrayVariable>( const char*, const FEM_Data<ArrayVariable>& );
template void Region<2>::InputVariableFrom<ArrayVariable>( const char*, const FEM_Data<ArrayVariable>& );
template void Region<3>::InputVariableFrom<ArrayVariable>( const char*, const FEM_Data<ArrayVariable>& );

template void Region<1>::InputVariableFrom<FlaggedArrayVariable>( const char*, const FEM_Data<FlaggedArrayVariable>& );
template void Region<2>::InputVariableFrom<FlaggedArrayVariable>( const char*, const FEM_Data<FlaggedArrayVariable>& );
template void Region<3>::InputVariableFrom<FlaggedArrayVariable>( const char*, const FEM_Data<FlaggedArrayVariable>& );

template void Region<1>::InputVariableFrom<VectorVariable<1U> >( const char*, const FEM_Data<VectorVariable<1U> >& );
template void Region<2>::InputVariableFrom<VectorVariable<2U> >( const char*, const FEM_Data<VectorVariable<2U> >& );
template void Region<3>::InputVariableFrom<VectorVariable<3U> >( const char*, const FEM_Data<VectorVariable<3U> >& );

template void Region<1>::InputVariableFrom<TensorVariable<1U> >( const char*, const FEM_Data<TensorVariable<1U> >& );
template void Region<2>::InputVariableFrom<TensorVariable<2U> >( const char*, const FEM_Data<TensorVariable<2U> >& );
template void Region<3>::InputVariableFrom<TensorVariable<3U> >( const char*, const FEM_Data<TensorVariable<3U> >& );










// -------------------------------------------------------------------
// Building blocks
// -------------------------------------------------------------------







/**
Detects of how many spatial dimensions element types are contained in model.
It returns a pair: first value gives number of different spatial dimensions contained,
second value returns the highest spatial dimension contained.

@author SKM 1/11/2013
*/
template<uint32_t dim>
pair<int32_t, int32_t>  Region<dim>::ElementSpatialDimensions() const
{
  return this->SpatialDimensions();

} // end ElementSpatialDimensions






/**
Creates set of lower dimensional elements that have a face at the boundary of the
model subdomain. Pointers to these elements are stored in supplied set,
and their number is returned.

@attention SKM method is not implemented yet.

template<uint32_t dim>
size_t  Region<dim>::IdentifyLowerDimensionalBoundaryElements( const std::pair<int32_t,int32_t>&,
set<Element<dim>*>& ldim_bdry_elmts ) const
{
throw logic_error("Region<dim>::IdentifyLowerDimensionalBoundaryElements: method is not implemented yet.");

ldim_bdry_elmts.clear();

size_t  counter(0U);

for ( typename vector<csmp::Node<dim>*>::const_iterator
it=node_vec_.begin(); it!=node_vec_.end(); it++ )
(*it)->Idx( counter++ );

return ldim_bdry_elmts.size();
}
*/








/**

ANDREW BROMAGE

SKM trying to make sense of Andrew Bromage's undocumented code:
12/08/18

*/
/*
template<uint32_t dim>
size_t Region<dim>::FromLargestComponent( MeshManager<dim>& mesh )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->cell_vec_.empty() )
    csmp_error.Note( WARNING, "Region<dim>::FromLargestComponent:",
                       "Region is not empty; deleting all content." );
  this->cell_vec_.clear();

  // traversal of the existing mesh nodes to find all its elements	
  deque<csmp::Node<dim>*>	   nodes;
  deque<csmp::Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( mesh, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  UnionFind<Node<dim>*> union_find;

  Node<dim>* component_node( nullptr );

  // 1. Loop over all elements, unioning node sets
  for ( auto eit : elmts ) {
    auto fe = eit->FE();
    const size_t iNrNodes = fe->Nodes();
    auto n1 = eit->N( 0u );

    for ( size_t iNode = 1; iNode < iNrNodes; ++iNode ) {
      auto n2 = eit->N( iNode );
      union_find.SameComponent( n1, n2 );
    }
  }

  // 2. Loop over all components, and find the largest
  std::deque<std::pair<size_t, Node<dim>*>> components;
  union_find.Components( components );
  size_t component_size = 0;
  for ( auto c : components ) {
    if ( c.first > component_size ) {
      component_node = c.second;
      component_size = c.first;
    }
  }

  if ( !component_node ) {
    csmp_error.Note( ERROR, "Region<dim>::FromLargestComponent:",
                       "Cannot find representative of largest component" );
  }

  // 3. Find the nodes in the component
  {
    std::vector<Node<dim>*> discovered_nodes_in_component;
    discovered_nodes_in_component.reserve( component_size );

    for ( auto nit : nodes ) {
      auto component = union_find.resolve( nit );
      if ( component == component_node ) {
        discovered_nodes_in_component.push_back( nit );
      }
    }
    this->node_vec_.swap( discovered_nodes_in_component );
  }

  // 4. Find the elements in the component
  {
    std::deque<Element<dim>*> discovered_elements_in_component;
    for ( auto eit : elmts ) {
      auto component = union_find.resolve( eit->N( 0u ) );
      if ( component == component_node ) {
        discovered_elements_in_component.push_back( eit );
      }
    }
    std::vector<Element<dim>*> elmts( discovered_elements_in_component.begin(), discovered_elements_in_component.end() );
    this->cell_vec_.swap( elmts );
  }

  // 5. (re)connecting elements up to their neighbors

  // 6. identifying the boundaries
  this->IdentifyPerimeter();

  return this->cell_vec_.size();
} // end FromLargestComponent
*/















/**
     AcculumateAll - but for all potentially disconnected Element patches that make up the model domain.
     
     @attention assumes that 'indexToPointerMapping' is unique and non empty.
*/
template<uint32_t dim>
size_t Region<dim>::AccumulateAll( MeshManager<dim>& mesh )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   if ( mesh.Elements() == 0 ) {
        csmp_error.Note( ERROR, "Region<dim>::AccumulateAll:", "supplied index-to-pointer mapping is empty.");
        return 0U;
     }
     
   if ( !this->cell_vec_.empty() ) this->cell_vec_.clear();
   if ( !this->node_vec_.empty() ) this->node_vec_.clear();
   
   // obtain pointers to all elements
   this->cell_vec_.reserve( mesh.Elements() );
   for ( auto it=mesh.ElementsBegin(); it!=mesh.ElementsEnd(); ++it )
     this->cell_vec_.push_back( &(*it) );
   
   // obtain pointers to all nodes
   this->node_vec_.reserve( mesh.Nodes() );
   for ( auto it=mesh.NodesBegin(); it!=mesh.NodesEnd(); ++it )
     this->node_vec_.push_back( &(*it) );

   // sorting in interior and perimeter ranges
   this->IdentifyPerimeter();
    
   return this->cell_vec_.size();

 } // end AccumulateAll (disconnected domain version)









/**
    Accumulates range of elements supplied as a vector. Assumes that there are no duplicate Element pointers in vector.
    
        @attention needs the input elements to be interconnected correctly already.
*/
template<uint32_t dim>
size_t Region<dim>::Accumulate( typename vector<csmp::Element<dim>*>::const_iterator start,
                                typename vector<csmp::Element<dim>*>::const_iterator end )
{
  if ( start == end )
    throw csmp::Exception( ERROR, "Region<dim>::Accumulate (vector)",
                           "supplied element range is empty. Nothing is done." );
                           
  this->cell_vec_.assign( start, end );
  
  this->CreateNodePointerVector();
  
  this->IdentifyPerimeter();

  return this->cell_vec_.size();

} // end Accumulate (vector)




template<uint32_t dim>
size_t Region<dim>::Accumulate( typename set<csmp::Element<dim>*>::const_iterator start,
                                typename set<csmp::Element<dim>*>::const_iterator end )
{
  if ( start == end )
    throw csmp::Exception( ERROR, "Region<dim>::Accumulate (set)",
                           "supplied element range is empty. Nothing is done." );

  this->cell_vec_.assign( start, end );

  this->CreateNodePointerVector();

  this->IdentifyPerimeter();

  return this->cell_vec_.size();

} // end Accumulate (set)







/**
All those elements for which all property constraints are met are accumulated
into this region. Returns the number of elements found.

@attention assumes that the iterator range is unique.

*/
template<uint32_t dim>
size_t Region<dim>::AccumulateWithinRange( MeshManager<dim>& mesh, const PropertyConstraints& constraints )
{
  if ( constraints.Constraints() == 0U )
    throw csmp::Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                           "No property constraints in supplied PropertyConstraints object" );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->cell_vec_.empty() )
    csmp_error.Note( WARNING, "Region<dim>::AccumulateWithinRange",
                       "Region<dim> already contains elements, they will be deleted" );

  if ( !this->bd_face_vec_.empty() ) this->bd_face_vec_.clear();

  if ( !this->cell_vec_.empty() ) this->cell_vec_.clear();
  this->cell_vec_.reserve( mesh.Elements()/2 );
  
  set<Node<dim>*> node_set;
  
  for ( auto it=mesh.ElementsBegin(); it!=mesh.ElementsEnd(); ++it )
    if ( constraints.CheckConstraints( &(*it) ) )
      {
        this->cell_vec_.push_back( &(*it) );
        const size_t n_nodes{ (*it).Nodes() };
        for ( auto i{0U}; i<n_nodes; ++i )
          node_set.insert( (*it).N(i) );
      }

  if ( !this->cell_vec_.empty() )
    {
       this->node_vec_.assign( node_set.begin(), node_set.end() );
       this->IdentifyPerimeter();
    }
  else
    csmp_error.Note( ERROR, "Region<dim>::AccumulateWithinRange", "no elements in the desired property range were found.");
  
  this->cell_vec_.shrink_to_fit();
  this->node_vec_.shrink_to_fit();
  
  return this->cell_vec_.size();

} // end AccumulateWithinRange(PropertyConstraints)




/**
Accumulates those property-bearing elements into a region, whose property values are greater or
equal to the lower bound of the specified range and smaller or equal to the
upper bound of it.

@attention if the target variable is a node property, by default at least one of the element's
nodes must have a value inside of the target range.
*/
template<uint32_t dim>
size_t Region<dim>::AccumulateWithinRange( MeshManager<dim>& mesh, const char* feature, double min, double max )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const Index prop_key = this->pref_.StorageKey( feature );
  if ( prop_key.place == MODEL ||prop_key.place == REGION ||
       prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY ||
       prop_key.place == FACE || prop_key.place == INTER_FACE ||
       prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
       prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ||
       prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
       prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT )
    {
       csmp_error.Note( ERROR, "Region<dim>::AccumulateWithinRange",
                          "property placement not handled by this method; use other",
                          parsePlacement(prop_key.place) );
       return 0U;
    }

  const auto end = mesh.ElementsEnd();
  auto     start = mesh.ElementsBegin();
  
  this->cell_vec_.reserve( mesh.Elements() );
  this->node_vec_.reserve( mesh.Nodes() );

  switch ( prop_key.place )
    {
      // if one node of an element is in the target range, it is added
      case NODE:
        while ( start != end ) {
            bool applies = false;
            const auto n_nodes{ (*start).Nodes() };
            for ( auto j{0U}; j<n_nodes; j++ )
              if ( (*start).N( j )->IsWithinRange( prop_key, min, max ) ) {
                applies = true;
                break;
              }
            if ( applies == true ) {
                this->cell_vec_.push_back( &(*start) );
                for ( auto i{0U}; i<n_nodes; i++ )
                  this->node_vec_.push_back( (*start).N(i) );
              }
            start++;
          }
        break;
      case ELEMENT_INTEGRATION_POINT:
        while ( start != end ) {
            bool applies = false;
            for ( auto j{0U}; j<(*start).IntegrationPoints(); j++ )
              if ( (*start).IsWithinRange( j, prop_key, min, max ) ) {
                applies = true;
                break;
              }
            if ( applies == true ) {
                this->cell_vec_.push_back( &(*start) );
                const size_t n_nodes{ (*start).Nodes() };
                for ( auto i{0U}; i<n_nodes; i++ )
                  this->node_vec_.push_back( (*start).N(i) );
              }
            start++;
          }
        break;
      case SECTOR_INTEGRATION_POINT:
        while ( start != end ) {
            bool applies = false;
            assert( (*start).IntegrationPointsPerSector() == 1 );
            auto n_sectors{ (*start).FV()->Sectors() };
            for ( auto j{0U}; j<n_sectors; j++ )
              if ( (*start).IsWithinRange( j, 0, prop_key, min, max ) ) {
                applies = true;
                break;
              }
            if ( applies == true ) {
                this->cell_vec_.push_back( &(*start) );
                const size_t n_nodes{ (*start).Nodes() };
                for ( auto i{0U}; i<n_nodes; i++ )
                  this->node_vec_.push_back( (*start).N(i) );
              }
            start++;
          }
        break;
      case FACET_INTEGRATION_POINT:
        while ( start != end ) {
            bool applies = false;
            assert( (*start).IntegrationPointsPerFacet() == 1 );
            const auto n_facets{ (*start).FV()->Facets() };
            for ( auto j{0U}; j<n_facets; j++ )
              if ( (*start).IsWithinRange( j, 0, prop_key, min, max ) ) {
                applies = true;
                break;
              }
            if ( applies == true ) {
                this->cell_vec_.push_back( &(*start) );
                const size_t n_nodes{ (*start).Nodes() };
                for ( auto i{0U}; i<n_nodes; i++ )
                  this->node_vec_.push_back( (*start).N(i) );
              }
            start++;
          }
        break;
      case ELEMENT:
        while ( start != end ) {
            if ( (*start).IsWithinRange( prop_key, min, max ) ) {
                this->cell_vec_.push_back( &(*start) );
                const auto n_nodes{ (*start).Nodes() };
                for ( auto i{0U}; i<n_nodes; i++ )
                  this->node_vec_.push_back( (*start).N(i) );
              }
            start++;
          }
        break;
      default:
        throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                         feature, "placement could not be identified; REGION is not an option" );
    } // end switch

  // removing duplicate nodes (unique needs vector to be sorted)
  sort( this->node_vec_.begin(), this->node_vec_.end() );
  this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );
  
  if( !this->cell_vec_.empty() ) {
      this->IdentifyPerimeter();
  }

  this->cell_vec_.shrink_to_fit();
  this->node_vec_.shrink_to_fit();
    
  return this->cell_vec_.size();

} // end AccumulateWithinRange






template<uint32_t dim>
size_t Region<dim>::AccumulateRectangularRegion( MeshManager<dim>& mesh,
                                                 const Point<dim>& xyz_min,
                                                 const Point<dim>& xyz_max )
{
  //ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const auto end = mesh.ElementsEnd();
  auto     start = mesh.ElementsBegin();
  
  this->cell_vec_.reserve( mesh.Elements() );
  this->node_vec_.reserve( mesh.Nodes() );

  if ( !this->cell_vec_.empty() ) {
      this->cell_vec_.clear();
      this->node_vec_.clear();
    }

  while ( start != end ) {
      size_t check{0};
      const auto n_nodes{ (*start).Nodes() };
      // all nodes have to be inside for the element selection criterion to be fulfilled
      for ( auto j{0U}; j<n_nodes; j++ ) {
          Point<dim>  p = (*start).N( j )->Coordinate();
          if ( p.IsBetween( xyz_min, xyz_max ) ) check++;
        }
      // if all nodes are inside the rectangular region
      // the element becomes part of the new group
      if ( check == n_nodes ) {
          this->cell_vec_.push_back( &(*start) );
          for ( auto i{0U}; i < n_nodes; i++ )
            this->node_vec_.push_back( (*start).N(i) );
        }
      start++;
    }

  // removing duplicate nodes (unique needs vector to be sorted)
  sort( this->node_vec_.begin(), this->node_vec_.end() );
  this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );

  this->IdentifyPerimeter();
  
  this->cell_vec_.shrink_to_fit();
  this->node_vec_.shrink_to_fit();

  return this->cell_vec_.size();

} // end AccumulateRectangularRegion





/**
       @author SKM revised 22/5/2021
*/
template<uint32_t dim>
size_t  Region<dim>::AccumulateByNumber( MeshManager<dim>& mesh,
                                         vector<size_t>& element_ids )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( element_ids.empty() )
    csmp_error.Note( ERROR, "Region<dim>::AccumulateByNumber",
                      "user-supplied element-number vector is empty. Nothing is done." );

  if ( !this->cell_vec_.empty() ) {
      csmp_error.Note( WARNING, "Region<dim>::AccumulateByNumber",
                         "Region is not empty", "erasing all members..." );
      this->cell_vec_.clear();
    }

  // eliminating potential duplicates from element index vector
#ifdef DEBUG
  const size_t n_elements{element_ids.size()};
  sort( element_ids.begin(), element_ids.end() );
  element_ids.erase( unique( element_ids.begin(), element_ids.end() ), element_ids.end() );
  if ( element_ids.size() < n_elements )
    csmp_error.Note( WARNING, "Region<dim>::AccumulateByNumber",
                      "user-supplied element ID set contained duplicates which were removed." );
#endif
  if ( element_ids.size() > mesh.Elements() )
    csmp_error.Note( ERROR, "Region<dim>::AccumulateByNumber",
                       "user-supplied element-number vector is larger than range of index-to-element-pointer mapping." );

  // creating the element vector for the region
  this->cell_vec_.reserve( element_ids.size() );
  for ( const auto& idx : element_ids ) {
       assert( idx < mesh.Elements() );
       Element<dim>* eptr = &(*next(mesh.ElementsBegin(),static_cast<int64_t>(idx)));
       assert( eptr != nullptr );
       assert( eptr->Idx() == idx );
       this->cell_vec_.push_back( eptr );
    }
    
  // creating node vector
  if ( this->node_vec_.empty() ) this->node_vec_.clear();
  this->node_vec_.reserve( element_ids.size() ); // just a loose measure, asuming that there will always be more elements than nodes
  // filling the vector
  for ( const auto& it : this->cell_vec_ ) {
       const size_t n_nodes{it->Nodes()};
       for ( auto i{0U}; i<n_nodes; ++i ) {
            assert( it->N(i) != nullptr );
            this->node_vec_.push_back( it->N(i) );
         }
     }

  // removing duplicates and trimming excess memory from node vector
  sort( this->node_vec_.begin(), this->node_vec_.end() );
  this->node_vec_.erase( unique( this->node_vec_.begin(), this->node_vec_.end() ), this->node_vec_.end() );

  this->IdentifyPerimeter();
  
  this->cell_vec_.shrink_to_fit();
  this->node_vec_.shrink_to_fit();

  return this->cell_vec_.size();

} // end AccumulateByNumber








/**
Form a group from an array of elements supplied as an STL vector.
The assumption made is that the element numbers are unique and method
will sort the supplied number vector.

@section arguments Input Arguments

The supplied element pointer vector must be unique.
The element numbers must be supplied using their global numbering
somewhere in the range 0..n-1. This numbering is expected to pre-exist.

The input vector is not passed by constant reference because it will
be sorted and duplicates are removed when in debug mode.

@section implementation Implementation

The supplied element id vector must be unique and sorted. In the debug
version, this is tested.

@section messages Messages

Method will detect if the supplied vector<double> is empty or if a group by
that name already exists.
*/
template<uint32_t dim>
size_t  Region<dim>::AccumulateByNumber( typename vector<Element<dim>*>::const_iterator start,
                                         typename vector<Element<dim>*>::const_iterator end,
                                         vector<size_t>& element_ids )
{
  if ( start == end )
    throw Exception( ERROR, "Region<dim>::AccumulateByNumber",
                     "user-supplied iterator range is empty. Nothing is done." );

  if ( element_ids.empty() )
    throw Exception( ERROR, "Region<dim>::AccumulateByNumber",
                     "user-supplied element-number vector is empty. Nothing is done." );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->cell_vec_.empty() ) {
    csmp_error.Note( WARNING, "Region<dim>::AccumulateByNumber",
                       "Region is not empty", "erasing all members..." );
    this->cell_vec_.clear();
  }

  // test whether there are consecutive duplicated elements
  sort( element_ids.begin(), element_ids.end() );
  vector<size_t>::iterator  new_end( unique( element_ids.begin(), element_ids.end() ) );
  if ( new_end != element_ids.end() )
    element_ids.erase( new_end, element_ids.end() );

  if ( element_ids.size() > static_cast<uint32_t>(distance( start, end )) )
    csmp_error.Note( ERROR, "Region<dim>::AccumulateByNumber",
                       "user-supplied element-number vector is larger than iterator range." );

  this->cell_vec_.reserve( element_ids.size() );
  set<Node<dim>*>  node_set;

  // selecting elements and nodes from the selected ID range
  while ( start != end )
  {
    assert( *start != nullptr );
    if ( binary_search( element_ids.begin(), element_ids.end(), (*start)->Idx() ) ) {
      // add element with the correct id to the region
      this->cell_vec_.push_back( (*start) );
      // add its nodes as well
      for ( auto i{0U}; i<(*start)->Nodes(); i++ )
        node_set.insert( (*start)->N( i ) );
    }
    start++;
  }

  this->node_vec_.assign( node_set.begin(), node_set.end() );
  assert( !this->node_vec_.empty() );

  this->IdentifyPerimeter();
  
  this->cell_vec_.shrink_to_fit();
  this->node_vec_.shrink_to_fit();

  return this->cell_vec_.size();

} // end AccumulateByNumber











/**
    Removes target elements from Region, rebuilding it afterwards.
    Reporting the number of removed elements.
    
    @attention the elements are not deleted, but pointers to them are returned into the second argument.
*/
template<uint32_t dim>
size_t Region<dim>::RemoveByNumber( vector<size_t>& element_ids, vector<Element<dim>*>& ptrs_to_removed_elements )
  {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( !ptrs_to_removed_elements.empty() )
      ptrs_to_removed_elements.clear();

    if ( element_ids.empty() ) {
         csmp_error.Note( WARNING, "Region<dim>::RemoveByNumber",  "user-supplied element-number vector is empty. Nothing was done." );
         return 0U;
      }

    if ( this->cell_vec_.empty() ) {
         csmp_error.Note( WARNING, "Region<dim>::RemoveByNumber",  "Region is empty. Nothing was done." );
         return 0U;
      }

    // making the element ID vector unique and searchable
    sort( element_ids.begin(), element_ids.end() );
    element_ids.erase( unique( element_ids.begin(), element_ids.end() ), element_ids.end() );
    
    // getting the elements for removal
    for ( typename vector<Element<dim>*>::const_iterator it=this->CellsBegin(); it!=this->CellsEnd(); ++it )
      if ( binary_search( element_ids.begin(), element_ids.end(), (*it)->Idx() ) )
        ptrs_to_removed_elements.push_back( (*it) );
     
    // finding the difference between the removal and the current element vector
    sort( this->cell_vec_.begin(), this->cell_vec_.end() );
    vector<Element<dim>*> elmts_to_retain;
    set_difference( this->cell_vec_.begin(), this->cell_vec_.end(),
                    ptrs_to_removed_elements.begin(), ptrs_to_removed_elements.end(),
                    inserter(elmts_to_retain, elmts_to_retain.begin()));
     
    // rebuilding the region
    this->cell_vec_.assign( elmts_to_retain.begin(), elmts_to_retain.end() );
    elmts_to_retain.clear();

    this->CreateNodePointerVector();

    this->IdentifyPerimeter();
    
    this->cell_vec_.shrink_to_fit();
    this->node_vec_.shrink_to_fit();

    return ptrs_to_removed_elements.size();

 } // end RemoveByNumber






/**
    Removes those element pointers from the region which match the ones in the supplied range.
*/
template<uint32_t dim>
size_t Region<dim>::RemoveRange( typename vector<csmp::Element<dim>*>::iterator begin,
                                 typename vector<csmp::Element<dim>*>::iterator end )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( begin == end  ) {
         csmp_error.Note( WARNING, "Region<dim>::RemoveRange",  "user-supplied element-pointer vector is empty. Nothing was done." );
         return 0U;
      }

    if ( this->cell_vec_.empty() ) {
         csmp_error.Note( WARNING, "Region<dim>::RemoveRange",  "Region is empty. Nothing was done." );
         return 0U;
      }
      
    const size_t n_elements{ this->Cells() };

    this->cell_vec_.erase( remove_if( this->cell_vec_.begin(), this->cell_vec_.end(),
                                      [&](auto x) { return binary_search( begin, end, x ); }),
                           this->cell_vec_.end() );
    
    this->CreateNodePointerVector();

    this->IdentifyPerimeter();
      
    this->cell_vec_.shrink_to_fit();
    this->node_vec_.shrink_to_fit();

    return n_elements - this->cell_vec_.size();
 
  } // end RemoveRange






/**
Adds supplied region to the current region.
*/
template<uint32_t dim>
void  Region<dim>::Add( const Region<dim>& grp )
{
  // if the added region is empty nothing needs to be done
  if ( grp.Empty() ) {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.Note( WARNING, "Region<dim>::Add", "Region to add was empty; nothing was done" );
    return;
  }

  const size_t n_cells_old{ this->Cells() };
  
  // appending the elements of the second region at the end
  this->cell_vec_.reserve( this->Cells() + grp.Cells() );
  for ( auto it=grp.CellsBegin(); it!=grp.CellsEnd(); ++it )
    this->cell_vec_.push_back( (*it) );
    
  // removing potential duplicates
  sort( this->cell_vec_.begin(), this->cell_vec_.end() );
  this->cell_vec_.erase( unique( this->cell_vec_.begin(), this->cell_vec_.end() ), this->cell_vec_.end() );
  vector<csmp::Element<dim>*>( this->cell_vec_ ).swap( this->cell_vec_ );

  const size_t n_cells_new{ this->Cells() };

  this->CreateNodePointerVector();
  
  cout <<"\n"<<"Region<dim>::Add: successfully added region '"<< grp.Name() <<"' to region '"<< this->Name();
  cout <<"', growing its size from "<< n_cells_old <<" to "<< n_cells_new <<" cells."<< endl;
  this->IdentifyPerimeter();

} // end Add










/**

Forms the union of group a and group b and returns it into the result
group 'res'.

@section arguments Input Arguments
The two groups which will be turned into one.

@attention this method cannot be called union because 'union' is a keyword
in the C++ language.

*/
template<uint32_t dim>
size_t  groupUnion( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.Note( WARNING, "union", "first region is empty" );
  if ( b.Empty() )
    csmp_error.Note( WARNING, "union", "second region is empty" );

  res.CellVector().reserve( a.Cells() + b.Cells() );
  for ( auto it=a.CellsBegin(); it!=a.CellsEnd(); ++it )
    res.CellVector().push_back( const_cast<Element<dim>* const>(*it) );
  for ( auto it=b.CellsBegin(); it!=b.CellsEnd(); ++it )
    res.CellVector().push_back( const_cast<Element<dim>* const>(*it) );
  
  sort( res.CellVector().begin(), res.CellVector().end() );
  res.CellVector().erase( unique(res.CellVector().begin(), res.CellVector().end()), res.CellVector().end() );
  res.CellVector().shrink_to_fit();

  res.CreateNodePointerVector();
  res.IdentifyPerimeter();

  return res.Cells();

} // end groupUnion







/**
Forms the intersection group 'res' of the groups 'a' and 'b'.
*/
template<uint32_t dim>
size_t  intersection( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.Note( WARNING, "intersection", "first region is empty" );
  if ( b.Empty() )
    csmp_error.Note( WARNING, "intersection", "second region is empty" );
    
  // since both regions are sorted already, a will be searched
  if ( !res.CellVector().empty() ) res.CellVector().clear();
  res.CellVector().reserve( min(a.Cells(),b.Cells()) );
  for ( auto it=b.CellsBegin(); it!=b.CellsEnd(); ++it )
    if ( binary_search( a.CellsBegin(), a.PerimeterCellsBegin(), (*it) ) ||
         binary_search( a.PerimeterCellsBegin(), a.CellsEnd(), (*it) ) )
      res.CellVector().push_back( const_cast<Element<dim>* const>(*it) );

  res.CellVector().shrink_to_fit();

  if ( !res.Empty() ) {
      res.CreateNodePointerVector();
      res.IdentifyPerimeter();
    }

  return res.Cells();

} // end intersection





/**
Returns into Region 'res' those elements of Region 'a' which are not contained
in Regions  'b'.
*/
template<uint32_t dim>
size_t  difference( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.Note( WARNING, "difference", "first region is empty" );
  if ( b.Empty() )
    csmp_error.Note( WARNING, "difference", "second region is empty" );

  // since both regions are sorted already, a will be searched first
  if ( !res.CellVector().empty() ) res.CellVector().clear();
  res.CellVector().reserve( min(a.Cells(),b.Cells()) );
  for ( auto it=a.CellsBegin(); it!=a.CellsEnd(); ++it )
    if ( !binary_search( b.CellsBegin(), b.PerimeterCellsBegin(), (*it) ) &&
         !binary_search( b.PerimeterCellsBegin(), b.CellsEnd(), (*it) ) )
      res.CellVector().push_back( const_cast<Element<dim>* const>(*it) );

  res.CellVector().shrink_to_fit();

  if ( !res.Empty() ) {
      res.CreateNodePointerVector();
      res.IdentifyPerimeter();
    }

  return res.Cells();

} // end difference






/**
Returns into Region 'res' those elements of Regions 'a' which are not in 'b' and those in 'b' which are
not in 'a'.
*/
template<uint32_t dim>
size_t  symmetricDifference( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.Note( WARNING, "symmetricDifference", "first region is empty" );
  if ( b.Empty() )
    csmp_error.Note( WARNING, "symmetricDifference", "second region is empty" );

  if ( !res.CellVector().empty() ) res.CellVector().clear();
  res.CellVector().reserve( min(a.Cells(),b.Cells()) );
  for ( auto it=b.CellsBegin(); it!=b.CellsEnd(); ++it )
    if ( !binary_search( a.CellsBegin(), a.PerimeterCellsBegin(), (*it) ) &&
         !binary_search( a.PerimeterCellsBegin(), a.CellsEnd(), (*it) ) )
      res.CellVector().push_back( const_cast<Element<dim>* const>(*it) );

  for ( auto it=a.CellsBegin(); it!=a.CellsEnd(); ++it )
    if ( !binary_search( b.CellsBegin(), b.PerimeterCellsBegin(), (*it) ) &&
         !binary_search( b.PerimeterCellsBegin(), b.CellsEnd(), (*it) ) )
      res.CellVector().push_back( const_cast<Element<dim>* const>(*it) );

  // removing potential duplicates
  sort( res.CellVector().begin(), res.CellVector().end() );
  res.CellVector().erase( unique(res.CellVector().begin(), res.CellVector().end()), res.CellVector().end() );
  res.CellVector().shrink_to_fit();

  res.CellVector().shrink_to_fit();

  if ( !res.Empty() ) {
    res.CreateNodePointerVector();
    res.IdentifyPerimeter();
  }
  return res.Cells();

} // end symmetricDifference







/**
    Counts and returns the number of Elements shared between the two regions.
*/
template<uint32_t dim>
size_t  sharedElements( const Region<dim>& g1, const Region<dim>& g2 )
{
  // ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
  auto  first1( g1.CellsBegin() );
  auto  first2( g2.CellsBegin() );
  size_t shared_elements( 0U );

  // comparing the interior nodes
  while ( first1 != g1.PerimeterCellsBegin() and first2 != g2.PerimeterCellsBegin() )
    {
      if ( *first1 < *first2 ) ++first1;
      else if ( *first2 < *first1 ) ++first2;
      else {
        shared_elements++;
        first1++;
        first2++;
      }
    }

  // comparing the boundary nodes
  first1 = g1.PerimeterCellsBegin();
  first2 = g2.PerimeterCellsBegin();

  while ( first1 != g1.CellsEnd() and first2 != g2.CellsEnd() )
    {
      if ( *first1 < *first2 ) ++first1;
      else if ( *first2 < *first1 ) ++first2;
      else {
        shared_elements++;
        first1++;
        first2++;
      }
    }

  return shared_elements;

} // end sharedElements


template size_t  sharedElements( const Region<1>&, const Region<1>& );
template size_t  sharedElements( const Region<2>&, const Region<2>& );
template size_t  sharedElements( const Region<3>&, const Region<3>& );








// TOPOLOGICAL INFORMATION

/**
If Region contains all the elements of Region 'g', Includes() returns true.
*/
template<uint32_t dim>
bool Region<dim>::Includes( const Region<dim>& g ) const
{
  // because region vectors are sorted into 2 seperate ranges,
  // new sorted vectors spanning the whole ranges need to be established first
  // before includes() can be called
  vector<Element<dim>*> tmp_elmt_vec1( this->cell_vec_ );
  vector<Element<dim>*> tmp_elmt_vec2( g.cell_vec_ );
  // sorting the vectors
  sort( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end() );
  sort( tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

  return includes( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end(),
                   tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

} // end Includes

static bool  hasLowerDimensionalRepresentation( const Element<3U>& element )
{
  return !element.FE()->IsVolume();
}

static bool  hasLowerDimensionalRepresentation( const Element<2U>& element )
{
  return !element.FE()->IsSurface();
}

template<>
bool  hasLowerDimensionalRepresentation<1>( const Region<1>& )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  csmp_error.Note( csmp::ERROR, "hasLowerDimensionalRepresentation", "na for 1D" );
  return false;
}

template<uint32_t dim>
bool  hasLowerDimensionalRepresentation( const Region<dim>& region )
{
  const auto elementsEnd( region.CellsEnd() );
  for ( auto it = region.CellsBegin(); it != elementsEnd; ++it )
  {
    if ( !hasLowerDimensionalRepresentation( *(*it) ) )
      return false;
  }
  return true;
}



/// if all BOX_BOUDARY flags != NOT region is on an external boundary; method highlights potential inconsistencies with topology
template<uint32_t dim>
bool  formsPartOfExternalBoundary( const Region<dim>& region, bool check_topo_as_well )
  {
    bool all_nodes_are_at_boundary{ true };
    
    if ( check_topo_as_well ) {
         for ( const auto& nit : region.NodeVector() ) {
             if ( nit->AtBoundary() == NOT )
               all_nodes_are_at_boundary = false;
             if ( nit->Attribute() == MESH_VERTEX ||
                  nit->Attribute() == INTERSECTION_POINT ||
                  nit->Attribute() == INTERIOR_LINE ||
                  nit->Attribute() == INTERIOR_SURFACE ) {
                  all_nodes_are_at_boundary = false;
                  cout <<"\n"<<"formsPartOfExternalBoundary(check_topo_as_well): detected vertex flagged '";
                  cout << parseTopology( nit->Attribute() ) <<"'"<< endl;
               }
              
          }
      }
    else {
        for ( const auto& it : region.NodeVector() ) {
             if ( it->AtBoundary() == NOT )
               all_nodes_are_at_boundary = false;
          }
      }
    return all_nodes_are_at_boundary;
 }
 
template bool formsPartOfExternalBoundary( const Region<2U>&, bool );
template bool formsPartOfExternalBoundary( const Region<3U>&, bool );

template<>
bool formsPartOfExternalBoundary( const Region<1U>&, bool ) {
     throw logic_error("ERROR, formsPartOfExternalBoundary: boundary objects do not exist in one-dimensional models\n");
  }

 


template<uint32_t dim>
bool  containsVolumeElements( const Region<dim>& region )
{
  const auto elementsEnd( region.CellsEnd() );
  for ( auto it = region.CellsBegin(); it != elementsEnd; ++it )
  {
    if ( (*it)->FE()->IsVolume() )
      return true;
  }
  return false;
}

template<uint32_t dim>
bool  containsSurfaceElements( const Region<dim>& region )
{
  const auto elementsEnd( region.CellsEnd() );
  for ( auto it = region.CellsBegin(); it != elementsEnd; ++it )
  {
    if ( (*it)->FE()->IsSurface() )
      return true;
  }
  return false;
}

template<uint32_t dim>
bool  containsLineElements( const Region<dim>& region )
{
  const auto elementsEnd( region.CellsEnd() );
  for ( auto it = region.CellsBegin(); it != elementsEnd; ++it )
  {
    if ( (*it)->FE()->IsLine() )
      return true;
  }
  return false;
}


// auxiliary function instantiations
template size_t  groupUnion<1U>( const Region<1>&, const Region<1>&, Region<1>& );
template size_t  groupUnion<2U>( const Region<2>&, const Region<2>&, Region<2>& );
template size_t  groupUnion<3U>( const Region<3>&, const Region<3>&, Region<3>& );

template size_t  intersection<1U>( const Region<1>&, const Region<1>&, Region<1>& );
template size_t  intersection<2U>( const Region<2>&, const Region<2>&, Region<2>& );
template size_t  intersection<3U>( const Region<3>&, const Region<3>&, Region<3>& );

template size_t  difference<1U>( const Region<1>&, const Region<1>&, Region<1>& );
template size_t  difference<2U>( const Region<2>&, const Region<2>&, Region<2>& );
template size_t  difference<3U>( const Region<3>&, const Region<3>&, Region<3>& );

template size_t  symmetricDifference<1U>( const Region<1>&, const Region<1>&, Region<1>& );
template size_t  symmetricDifference<2U>( const Region<2>&, const Region<2>&, Region<2>& );
template size_t  symmetricDifference<3U>( const Region<3>&, const Region<3>&, Region<3>& );

template bool hasLowerDimensionalRepresentation( const Region<3>& );
template bool hasLowerDimensionalRepresentation( const Region<2>& );

template bool containsVolumeElements( const Region<3>& );
template bool containsVolumeElements( const Region<2>& );
template bool containsVolumeElements( const Region<1>& );
template bool containsSurfaceElements( const Region<3>& );
template bool containsSurfaceElements( const Region<2>& );
template bool containsSurfaceElements( const Region<1>& );
template bool containsLineElements( const Region<3>& );
template bool containsLineElements( const Region<2>& );
template bool containsLineElements( const Region<1>& );



// GEOMETRY MANIPILATIONS


/**
Displaces the node coordinates by displacements supplied as through the
VECTOR variable (dim = dimensions of the model). For each element, the
nodes of which were displaced, the private boolean variable shape_to_date
is set to false such that its volume is newly calculated
once it is requested.

@param vector_variable The name of the nodal vector<double> variable which holds the node
coordinate displacement.

@section implementation Implementation

After the assignment of the nodal displacements, MoveNodeCoordinatesBy()
will set the Element state variable 'shape_to_data' to false. This will
prompt re-calculation of elemental volumes if these are queried in
subsequent computations.

@section application Application

If a mesh shall be deformed using the diplacements of a deformation
calculation stored in a vector<double> variable, ChangeNodeCoordinatesTo() can
be used displace the node coordinates by these displacements.

@section messages Messages

Due to the total garbage results that may arise,
MoveNodeCoordinatesBy() will halt the simulation reporting a fatal
error, if the target property is node a node or vector<double> type property.
*/
template<uint32_t dim>
void Region<dim>::MoveNodeCoordinatesBy( const char* vector_variable )
{
  csmp::Index          prop_key = this->pref_.StorageKey( vector_variable );
  VectorVariable<dim>  vc;

  if ( prop_key.type != VECTOR )
    throw csmp::Exception( ERROR, "Model<dim>::MoveNodeCoordinatesBy",
                           "Only vector<double> variables can be added to coordinates" );

  if ( prop_key.place != NODE )
    throw csmp::Exception( ERROR, "Model<dim>::MoveNodeCoordinatesBy",
                           "Only NODE variables can be added to NODE coordinates" );

  for ( auto nit = this->NodesBegin(); nit != this->NodesEnd(); nit++ ) {
      (*nit)->Read( prop_key, vc );
      (*nit)->Coordinate( (*nit)->Coordinate() + vc );
    }

} // end MoveNodeCoordinatesBy



// Correct Low Dimensional Region Orientation if necessary
// Method chooses refence elemnt and correct the orientation of neighbor elements correspondingly to that element
// NOTE: extend for non-contiguos Regions
// @author Roman M.
// @date   27/07/2014
/*
template<uint32_t dim>
void Region<dim>::CorrectLowDimRegionOrientation(  ) const
{

// pair contains: 1) number of spatial element dimensions in region, 2)  highest contained dimension
// TODO: improved version should include more information about perimeter elements ( stickin out elements )
pair<int32_t,int32_t>  elmt_dim = this->ElementSpatialDimensions();

// Correct Orientation of low dim region

std::set<Element<dim>* > elements_considered;
std::set<Element<dim>* > element_neighbors_to_be_considered;
std::set<Element<dim>* > remaining_elements( this->CellsBegin(), this->CellsEnd() );

Element<dim>*   lowDimElement         ( NULL );
Element<dim>*   lowDimNeighborElement ( NULL );

VectorVariable<dim> low_dim_unit_normal               ( ANY, 0.0 );
VectorVariable<dim> low_dim_neighbor_unit_normal      ( ANY, 0.0 );
VectorVariable<dim> low_dim_edge_unit_normal          ( ANY, 0.0 );
VectorVariable<dim> low_dim_neighbor_edge_unit_normal ( ANY, 0.0 );

// Find first refence element
typename std::vector<csmp::Element<dim>* >::const_iterator eit = this->PerimeterElementsBegin();
lowDimElement = NULL;
while ( eit!=this->CellsEnd() )
{
if( (*eit)->IsLine() && elmt_dim.second == 2 ) // do not consider line element in surface region in 3D
{
elements_considered.insert( *eit );
remaining_elements.erase( *eit );
eit++;
}
else
{
lowDimElement = (*eit);
break;
}
}
if( lowDimElement == NULL )
throw csmp::Exception( ERROR, "Region<dim>::CorrectLowDimRegionOrientation()", "All perimeter elements of surface regions are lines! Cannot apply this method!" );

while( elements_considered.size() != this->Cells() )
{
element_neighbors_to_be_considered.insert( lowDimElement );
elements_considered.insert( lowDimElement );
remaining_elements.erase( lowDimElement );

while( !element_neighbors_to_be_considered.empty() )
{
//Mesh orientation diagnostic
//cerr<<"\nN1 = "<<lowDimElement->N(0)->x() << "; N2 = "<<lowDimElement->N(1)->x() << endl;

// unit normal to low dimensional element
lowDimElement = *element_neighbors_to_be_considered.begin();
low_dim_unit_normal = lowDimElement->UnitNormal( );

for ( size_t face=0U; face<lowDimElement->Neighbors(); face++ )
{
lowDimNeighborElement = lowDimElement->Neighbor( face );

if( lowDimNeighborElement!= NULL )
{
if( ( elements_considered.find( lowDimNeighborElement ) == elements_considered.end() ) )
{
// unit normal to low dimensional neighbor element
low_dim_neighbor_unit_normal = lowDimNeighborElement->UnitNormal( );

bool face_found( false );
size_t neighbor_face = 0;
// unit normal to face that elements share
for( ; neighbor_face< lowDimNeighborElement->Neighbors(); neighbor_face++ )
if( lowDimNeighborElement->Neighbor( neighbor_face ) == lowDimElement )
{
face_found = true;
break;
}
if( !face_found )
throw csmp::Exception( ERROR, "Region<dim>::CorrectLowDimRegionOrientation()", "The neighbor connectivity in low dimensional region is wrong!" );

low_dim_edge_unit_normal          = lowDimElement->UnitNormalToFace( face );
low_dim_neighbor_edge_unit_normal = lowDimNeighborElement->UnitNormalToFace( neighbor_face );

if( ( low_dim_unit_normal.CrossProduct( low_dim_edge_unit_normal ) & low_dim_neighbor_unit_normal.CrossProduct( low_dim_neighbor_edge_unit_normal ) ) > 0.0 )
{
lowDimNeighborElement->LowDimRevertNodeNumbering();
lowDimNeighborElement->LowDimRevertNeighborNumbering();
}

element_neighbors_to_be_considered.insert( lowDimNeighborElement );
elements_considered.insert( lowDimNeighborElement );
remaining_elements.erase( lowDimElement );
}
}
}
element_neighbors_to_be_considered.erase( lowDimElement );
}

if( elements_considered.size() != this->Cells() )
{
typename std::set<csmp::Element<dim>* >::const_iterator eit = remaining_elements.begin();
lowDimElement = NULL;
while ( eit!=remaining_elements.end() )
{
if( (*eit)->IsLine() && elmt_dim.second == 2 ) // do not consider line element in surface region in 3D
{
elements_considered.insert( *eit );
remaining_elements.erase( *eit );
eit++;
}
else
{
lowDimElement = (*eit);
break;
}
}
if( lowDimElement == NULL )
break;
}
}
elements_considered.clear();
remaining_elements.clear();

}
*/




// CALCULATIONS



/**
Returns the area (2D) or the volume (3D) which is occupied by the sum of
elements which constitute the Region.

If regions consist of elements from multiple space dimensions, the highest
ones are used only: Thus, the volume of a 3D group "Model" will exclude
potential surface elements.

@section implementation Implementation

Sums elemental areas/volumes as calculated via the element interpolation
functions of the current finite elements.

@section application Application

For instance, if you want to ask your model which rock volume has been
heated above 400oC at a certain timestep.
*/
template<uint32_t dim>
double  Region<dim>::Volume( bool multiply_with_porosity ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( this->Empty() ) {
    csmp_error.Note( WARNING, "Region<dim>::Volume:", "region is empty; returning NaN." );
    return std::numeric_limits<double>::signaling_NaN();
  }
  double  volume( 0. ), area( 0. ), length( 0. );

  // recording contributions of different types of elements making up the group
  if ( multiply_with_porosity ) {
      csmp::Index  phi_key = this->pref_.StorageKey( "porosity" );
      for ( const auto& it : this->cell_vec_ ) {
          if      ( it->IsVolume() )  volume += it->Volume() * it->Read( phi_key );
          else if ( it->IsSurface() ) area   += it->Volume() * it->Read( phi_key );
          else if ( it->IsLine() )    length += it->Volume() * it->Read( phi_key );
        }
    }
  else
    for ( const auto& it : this->cell_vec_ ) {
        if      ( it->IsVolume() )  volume += it->Volume();
        else if ( it->IsSurface() ) area   += it->Volume();
        else if ( it->IsLine() )    length += it->Volume();
      }

  // counting only the contributions of the elements with the highest spatial dimension
  if constexpr ( dim == 3U ) {
      if ( volume > 0. ) return volume;
      if ( area >   0. ) return area;
      if ( length > 0. ) return length;
    }
  
  if constexpr ( dim == 2U )
    return (area > 0.) ? area : length;

  return length;

} // end Volume



// STUB FOR 1D CALCULATION
inline double triangleArea( const Point<1U>&, const Point<1U>&, const Point<1U>& ) {
     return 1.;
  }

/**
Computes the perimeter (2D) or the surface area (3D) of a model subdomain object.
To compute the surface area, only elements of dim-1 are considered.

@attention method assumes that the first nodes in higher-order elements are the corner nodes.

@todo SKM: method will produce nonsense if line elements intersect the boundary.

*/
template<uint32_t dim>
double  Region<dim>::SurfaceArea() const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  // in 1D there is no meaningful surface area
  if ( dim == 1U ) {
    csmp_error.Note( WARNING, "Region<1U>::SurfaceArea:",
                       "is not defined in one-dimensional model; returning NaN" );
    return std::numeric_limits<double>::signaling_NaN();
  }

  auto    bit( this->bd_face_vec_.begin() );
  double  area{0.};

  if constexpr ( dim == 3U ) {
      // for all elements located on the region boundary
      for ( auto i = this->InteriorCells(); i<this->Cells(); ++i, ++bit )
        // since each element can have multiple boundary faces
        for ( uint32_t j{0U}; j<(*bit).size(); j++ ) {
          const uint32_t face( (*bit)[j] );
          auto fnids = this->cell_vec_[i]->FE()->NodesOfFace( face );
          const CSMP_FEM_TYPE etype( this->cell_vec_[i]->FE()->ElementTypeOfFace( face ) );
          // triangular face
          if ( isTriangular(etype) )
            area += triangleArea( this->cell_vec_[i]->N( fnids[0U] )->Coordinate(),
                                  this->cell_vec_[i]->N( fnids[1U] )->Coordinate(),
                                  this->cell_vec_[i]->N( fnids[2U] )->Coordinate() );
          // quadrilateral face
          else if ( isQuadrilateral(etype) )
            area += facetArea4( this->cell_vec_[i]->N( fnids[0U] )->Coordinate(),
                                this->cell_vec_[i]->N( fnids[1U] )->Coordinate(),
                                this->cell_vec_[i]->N( fnids[2U] )->Coordinate(),
                                this->cell_vec_[i]->N( fnids[3U] )->Coordinate() );
          // linear face
          else if ( isLineElement(etype) ) {
            area += distance( this->cell_vec_[i]->N( fnids[0U] )->Coordinate(), this->cell_vec_[i]->N( fnids[1U] )->Coordinate() );
            //csmp_error.Note( WARNING, "Region<dim>::SurfaceArea", "line-element thickness on perimeter is assumed to be one." );
          }
          // additional case of point face where a line-element is perpendicular to a boundary node
        }
    }
  // in 2D the face is a segment the length of which has to be used
  else if constexpr ( dim == 2U ) {
      // for all surface elements on the region boundary (excluding line elements)
      for ( size_t i = this->InteriorCells(); i<this->Cells(); i++, bit++ )
        if ( this->cell_vec_[i]->FE()->IsSurface() )
          for ( size_t j{0U}; j<(*bit).size(); j++ ) {
            auto fnids = this->cell_vec_[i]->FE()->NodesOfFace( (*bit)[j] );
            area += distance( this->cell_vec_[i]->N( fnids[0U] )->Coordinate(), this->cell_vec_[i]->N( fnids[1U] )->Coordinate() );
          }
    }

  return area;

} // end SurfaceArea


/**
Integrates variable over all elements contained in the Region and returns
the result as volume integral of the variable.

Optionally, the integrated variable can be multiplied with the 'porosity',
assuming that this property is piecewise contant from element to element,
i.e. an element property.

@attention method should only be applied to volumetric regions.
*/
template<uint32_t dim>
double  Region<dim>::VolumeIntegral( const char* prop, bool multiply_with_porosity, bool verbose ) const
{
  csmp::Index prop_key = this->pref_.StorageKey( prop );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( prop_key.type == TENSOR ) {
    csmp_error.Note( ERROR, "Region<dim>::VolumeIntegral:",
                       "No rule to integrate tensor properties. Nothing was done" );
    return std::numeric_limits<double>::signaling_NaN();
  }
  if ( prop_key.type == VECTOR and prop_key.place != ELEMENT ) {
    csmp_error.Note( ERROR, "Region<dim>::VolumeIntegral:",
                       "Vector properties can only be integrated if they are placed on the element. Nothing was done" );
    return std::numeric_limits<double>::signaling_NaN();
  }
  if ( this->cell_vec_.empty() ) {
    csmp_error.Note( ERROR, "Region<dim>::VolumeIntegral", "Region is empty; returning NaN." );
    return std::numeric_limits<double>::signaling_NaN();
  }

  if ( verbose ) {
#ifndef NDEBUG
    if ( dim == 2U && !(*this->cell_vec_.begin())->FE()->IsSurface() ) {
      string info( this->Name() ); info += " ('"; info += prop; info += "')";
      csmp_error.Note( WARNING, "Region<2>::VolumeIntegral:", info, "region is not a surface; integral may not be correct." );
    }
    if ( dim == 3U && !(*this->cell_vec_.begin())->FE()->IsVolume() ) {
      string info( this->Name() ); info += " ('"; info += prop; info += "')";
      csmp_error.Note( WARNING, "Region<3>::VolumeIntegral:", info, "region is not a volume, integral may not be correct." );
    }
#endif
  }
  // region properties
  if ( prop_key.place == REGION ) {
    if ( prop_key.type == SCALAR )
      return this->Read( prop_key ) * Volume( multiply_with_porosity );
    else if ( prop_key.type == VECTOR ) {
      VectorVariable<dim>  vc;
      this->Read( prop_key, vc );
      return vc.Length() * Volume( multiply_with_porosity );
    }
    return std::numeric_limits<double>::signaling_NaN();
  }

  double  integral( 0. );

  // node or integration point properties
  if ( prop_key.place != ELEMENT )
  {
    csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
    // integration of nodal properties using the FEM framework
    if ( prop_key.place == NODE or prop_key.place == ELEMENT_INTEGRATION_POINT ) {
      if ( multiply_with_porosity ) {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key ) * (*iter)->Read( phi_key );
      }
      else {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key );
      }
      return integral;
    }
    else {
      throw csmp::Exception( ERROR, "Region<dim>::VolumeIntegral",
                             "only node, integration point or element properties can be integrated over the region" );
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  else // element property
  {
    if ( prop_key.type == SCALAR ) {
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume() * (*iter)->Read( phi_key );
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume();
    }
    else { // VECTOR property
      VectorVariable<dim>  vc;
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume() * (*iter)->Read( phi_key );
        }
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume();
        }
    }
    return integral;
  }

  return std::numeric_limits<double>::signaling_NaN();

} // end VolumeIntegral





/**
Integrates property over the elements contained in the Region and returns
the volume integral of the variable multiplied with the local element thickness.
Element thickness must be equal to one for volumetric elements.
*/
template<uint32_t dim>
double  Region<dim>::VolumeIntegral_x_Thickness( const char* prop, bool multiply_with_porosity ) const
{
  csmp::Index prop_key = this->pref_.StorageKey( prop );
  csmp::Index thic_key = this->pref_.StorageKey( "thickness" );
  assert( thic_key.place == ELEMENT );
  assert( thic_key.type == SCALAR );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( prop_key.type == TENSOR ) {
    csmp_error.Note( ERROR, "Region<dim>::VolumeIntegral_x_Thickness",
                       "No rule to integrate tensor properties. Nothing was done" );
    return std::numeric_limits<double>::signaling_NaN();
  }
  if ( prop_key.type == VECTOR and prop_key.place != ELEMENT ) {
    csmp_error.Note( ERROR, "Region<dim>::VolumeIntegral_x_Thickness",
                       "Vector properties can only be integrated if they are placed on the element. Nothing was done" );
    return std::numeric_limits<double>::signaling_NaN();
  }
  if ( this->cell_vec_.empty() ) {
    csmp_error.Note( ERROR, "Region<dim>::VolumeIntegral_x_Thickness", "Region is empty" );
    return std::numeric_limits<double>::signaling_NaN();
  }

  // region properties (no thickness multiplier is accounted for)
  if ( prop_key.place == REGION ) {
    csmp_error.Note( WARNING, "Region<dim>::VolumeIntegral_x_Thickness", "there is no thickness attribute for region variable, assuming t=1" );
    if ( prop_key.type == SCALAR )
      return this->Read( prop_key ) * Volume( multiply_with_porosity );
    else if ( prop_key.type == VECTOR ) {
      VectorVariable<dim>  vc;
      this->Read( prop_key, vc );
      return vc.Length() * Volume( multiply_with_porosity );
    }
    return std::numeric_limits<double>::signaling_NaN();
  }

  double  integral( 0. );

  // node or integration point properties
  if ( prop_key.place != ELEMENT )
  {
    csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
    // integration of nodal properties using the FEM framework
    if ( prop_key.place == NODE or prop_key.place == ELEMENT_INTEGRATION_POINT ) {
      if ( multiply_with_porosity ) {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key ) * (*iter)->Read( phi_key ) * (*iter)->Read( thic_key );
      }
      else {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key ) * (*iter)->Read( thic_key );
      }
      return integral;
    }
    else {
      throw csmp::Exception( ERROR, "Region<dim>::VolumeIntegral_x_Thickness",
                             "only node, integration point or element properties can be integrated over the region" );
      return std::numeric_limits<double>::signaling_NaN();
    }
  }

  else // element property
  {
    if ( prop_key.type == SCALAR ) {
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume() * (*iter)->Read( phi_key ) * (*iter)->Read( thic_key );
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume() * (*iter)->Read( thic_key );
    }
    else { // VECTOR property
      VectorVariable<dim>  vc;
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume() * (*iter)->Read( phi_key ) * (*iter)->Read( thic_key );
        }
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->cell_vec_.begin(); iter != this->cell_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume() * (*iter)->Read( thic_key );
        }
    }
    return integral;
  }

  return std::numeric_limits<double>::signaling_NaN();

} // end VolumeIntegral_x_Thickness







// SCREEN OUTPUT



/*
template<uint32_t dim>
void Region<dim>::Out() const
{
cout <<"\n\n\nRegion<"<< dim <<">::Out: ";
cout <<" member elements: interior="<< this->InteriorCells();
cout <<", boundary="<< this->cell_vec_.size()-this->InteriorCells() <<": "<< endl;

for ( typename vector<csmp::Element<dim>*>::const_iterator
it=this->cell_vec_.begin(); it!=this->cell_vec_.end(); it++ ) {
if ( (*it) == NULL )
throw csmp::Exception( ERROR, "Region<dim>::Out",
"member element pointer not initialised");
//         else (*it)->Out();
}

cout <<"\n\n boundary elements and their boundary faces (current local numbering): "<< endl;
vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit(this->bd_face_vec_.begin());
for ( size_t i=this->InteriorCells(); i<this->cell_vec_.size(); i++, bit++ ) {
cout <<"\nelement "<< i <<": boundary face numbers: ";
for ( vector<ONE_BYTE_NUMBER>::const_iterator
ft=(*bit).begin(); ft!=(*bit).end(); ft++ ) cout << (*ft) <<" ";
}

cout <<"\n\n boundary nodes: "<< this->node_vec_.size() - this->first_bd_node_ <<" (current local numbering):"<< endl;
for ( size_t i=this->first_bd_node_; i<this->node_vec_.size(); i++ ) {
if ( this->node_vec_[i] == NULL )
throw csmp::Exception( ERROR, "Region<dim>::Out", "member node pointer not initialised.");
else cout << this->node_vec_[i]->Idx() <<" ";
}

cout << endl;

} // end Out
*/


// explicit instantiations

template class Region<1>;
template class Region<2>;
template class Region<3>;


} // end namespace csmp
