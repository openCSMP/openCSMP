#include "Region.h"

#include "VSet.h"
#include "PropertyData.h"
#include "FEM_Data.h"
#include "LocalVariableStorage.h"
#include "Node.h"
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
#include "CSMP_highLevelUtilities.h"

// #define REGION_DEBUG

using namespace std;

namespace csmp {

template<size_t dim>
Region<dim>::Region( std::string regionname, const PropertyDatabase<dim>& p )
  : ModelSubDomain<dim, Element>( regionname, p )
{
  this->ResizePropertyStorage( this->pref_.LocalVariablesAt( Placement() ) );
}


template<size_t dim>
Region<dim>::Region( const Region& g )
  : ModelSubDomain<dim, Element>( g )
{
}


template<size_t dim>
Region<dim>::Region( Region&& g )
  : ModelSubDomain<dim, Element>( g )
{
}


/**
Assignment does not affect the const reference to 'pref'.
It remains the same as that before the assignment
because it refers to the one Model in which all the Regions live.
*/
template<size_t dim>
Region<dim>&  Region<dim>::operator=( const Region& g )
{
  if ( &g != this ) {
    *this = g;
  }
  return *this;
} // end assignment

/**
The destructor does not reduce the counter of groups n_groups_.
This is avoided because the group which gets deleted is not necessarily
the group with the highest number. Thus, if the counter would be decremented
after each deletion, groups created thereafter might share the same
group ID.

As a consequence of the implemented approach, the group numbering is
no longer contiguous.
*/
template<size_t dim>
Region<dim>::~Region()
{
}




/**
Re-constructor for regions that were stored in the CSMP native
file format.

Using the indices retrieved from binary file and stored in SubDomainInfo,
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
template<size_t dim>
Region<dim>::Region( const PropertyDatabase<dim>& pref,
                     MeshManager<dim>& mesh, ///< not constant because region shall later be able to modify elements and nodes
                     const SubDomainInfo& info )   ///< information on how to connect pointers to mesh stored in MeshManager 
  : ModelSubDomain<dim, Element>( info.name, pref )
{
  // traversal of the existing mesh nodes to find all its elements
  deque<csmp::Node<dim>*>		nodes;
  deque<csmp::Element<dim>*>	elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // building the element vector
  // ---------------------------
  this->elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );

  // assigning pointers to the interior elements
  for ( size_t i : info.interior_elmts )
    this->elmt_vec_.push_back( elmts[i] );

  // assigning pointers to the perimeter elements
  for ( size_t i : info.perimeter_elmts )
    this->elmt_vec_.push_back( elmts[i] );

  // building the vector of vectors of those faces of the elements that lie on the subdomain perimeter
  // -------------------------------------------------------------------------------------------------
  const auto perimeterElementsBegin( next( this->elmt_vec_.begin(), info.interior_elmts.size() ) );
  this->bd_face_vec_.reserve( info.perimeter_faces.size() );
  const auto elementsEnd( this->elmt_vec_.end() );
  for ( auto it = perimeterElementsBegin; it != elementsEnd; ++it ) {
    assert( (*it)->Faces() == (*it)->Neighbors() );
    // for all the faces of the element that are located on the model boundary
    vector<ONE_BYTE_NUMBER>  boundary_faces;
    const size_t faces( (*it)->Faces() );
    boundary_faces.reserve( faces );
    for ( size_t face = 0U; face<faces; ++face )
      if ( (*it)->Neighbor( face ) == nullptr )
        boundary_faces.push_back( static_cast<ONE_BYTE_NUMBER>(face) );
    // storing the boundary face vector for the current element
    this->bd_face_vec_.push_back( boundary_faces );
  }

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

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( REGION ) );

} // end region re-constructor (using MeshManager)


/// re-constructor for regions via the nodes and elements which are explored by the MeshManager
template<size_t dim>
Region<dim>::Region( const PropertyDatabase<dim>& pref,
                     const deque<Node<dim>*>& nodes,
                     const deque<Element<dim>*>& elmts,
                     const SubDomainInfo& info )  ///< contains correctly partitioned vectors and boundary faces
  : ModelSubDomain<dim, Element>( info.name, pref )
{
  // building the element vector
  // ---------------------------
  this->elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );

  // assigning pointers to the interior elements
  for ( size_t i : info.interior_elmts )
    if ( i < elmts.size() ) this->elmt_vec_.push_back( elmts[i] );

  // assigning pointers to the perimeter elements
  for ( size_t i : info.perimeter_elmts )
    if ( i < elmts.size() ) this->elmt_vec_.push_back( elmts[i] );

  // building the vector of vectors of those faces of the elements that lie on the subdomain perimeter
  // -------------------------------------------------------------------------------------------------
  const auto perimeterElementsBegin( next( this->elmt_vec_.begin(), info.interior_elmts.size() ) );
  this->bd_face_vec_.reserve( info.perimeter_faces.size() );
  const auto elementsEnd( this->elmt_vec_.end() );
  for ( auto it = perimeterElementsBegin; it != elementsEnd; ++it ) {
    assert( (*it)->Faces() == (*it)->Neighbors() );
    // for all the faces of the element that are located on the model boundary
    vector<ONE_BYTE_NUMBER>  boundary_faces;
    const size_t faces( (*it)->Faces() );
    boundary_faces.reserve( faces );
    for ( size_t face = 0U; face<faces; ++face )
      if ( (*it)->Neighbor( face ) == nullptr )
        boundary_faces.push_back( static_cast<ONE_BYTE_NUMBER>(face) );
    // storing the boundary face vector for the current element
    this->bd_face_vec_.push_back( boundary_faces );
  }

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

  // allocating the storage for boundary properties
  // ----------------------------------------------
  this->ResizePropertyStorage( pref.LocalVariablesAt( REGION ) );

} // end region re-constructor (using the nodes and the elements from the MeshManager)


// LOCAL VARIABLE STORAGE INTERFACE
template<size_t dim>
bool Region<dim>::ValidVariable( const char* variableName ) const
{
  const PLACEMENT p( this->pref_.Placement( variableName ) );
  if ( p == NODE || p == ELEMENT || p == REGION )
    return true;
  return false;
}



template<size_t dim>
IntegrationPointVariables Region<dim>::ElementIntegrationPointVariables() const
{ return this->pref_.IntegrationPointVariablesAt( ELEMENT ); }

template<size_t dim>
LocalVariables Region<dim>::ElementVariables() const
{ return this->pref_.LocalVariablesAt( ELEMENT ); }



// VISITORS INTERFACE

/**
Recursive visitation of regions (or any other object) starts at application level
and finishes at application target.

@todo SKM provide methods for the new variable placements.

@todo SKM refactor the loops using constants and pre-incrementation to speed them up

@author SKM
*/
template<size_t dim>
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
      for ( typename vector<Element<dim>*>::iterator
            it = Region<dim>::ElementsBegin(); it != Region<dim>::ElementsEnd(); it++ )
        (*it)->Accept( v );
      return;
    case NODE:
      for ( typename vector<csmp::Node<dim>*>::iterator
            nd_it = Region<dim>::NodesBegin(); nd_it != Region<dim>::NodesEnd(); nd_it++ )
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
Outputs geometry and property data from the region to a vset.
The data of the region is copied as is. Boundary conditions will be copied if present.
For unresolved boundaries, the flag IRREGULAR is set.

@attention Relies on consecutively numbered indices

@return The vset which has the geometry that corresponds to the group.

@section implementation Implementation

All data is copied into the vset. IDs are used in the current state.

@section application Application

The idea is to create a vset that contains all the information of the region.
This vset can, for example, later be used to create a new Model.
*/
template<size_t dim>
void Region<dim>::OutputTo( VSet<dim>& vset, bool with_properties ) const
{
  if ( this->elmt_vec_.empty() ) {
    throw csmp::Exception( ERROR, "Region<dim>::OutputTo",
                           "Attempt to output empty group to VSet. Nothing was done." );
    return;
  }

  if ( vset.Vertices() > 0U ) vset.Erase();

  // 0. resizing the VSet
  size_t         counter( 0U );
  deque<size_t>  nodes_per_element( this->elmt_vec_.size() );
  deque<size_t>  elements_per_element( this->elmt_vec_.size() );
  deque<int32>   etypes( this->elmt_vec_.size() );
  set<int32>     n_etypes;

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
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

  // 1. assigning coordinate values
  for ( size_t i = 0U; i<this->node_vec_.size(); i++ ) {
    vset.Px( i, this->node_vec_[i]->x() );
    if ( dim != 1U ) vset.Py( i, this->node_vec_[i]->y() );
    if ( dim == 3U ) vset.Pz( i, this->node_vec_[i]->z() );
  }

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
    // mapping global to local node numbers
    for ( size_t i = 0U; i<(*eit)->Nodes(); i++ ) vset.Plist( (*eit)->Idx(), i, (*eit)->N( i )->Idx() );

  // 3. writing the neighbor per element map (pfverts)
  counter = 0U;
  for ( typename vector<csmp::Element<dim>*>::const_iterator
        eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ ) {
    for ( size_t i = 0U; i<(*eit)->Neighbors(); i++ )
      if ( (*eit)->Neighbor( i ) != NULL )
        vset.Pfvert( counter, i, static_cast<int32>((*eit)->Neighbor( i )->Idx()) );
      else {
        vset.Pfvert( counter, i, -1 );
      }
      counter++;
  }

  // 4. writing nodal boundary flags
  unordered_map<size_t, long64>  bflags;

  for ( size_t i = this->first_bd_node_; i<this->node_vec_.size(); i++ )
    bflags.insert( make_pair( i, this->node_vec_[i]->AtBoundary() ) );

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





template<size_t dim>
void Region<dim>::OutputDataTo( VSet<dim>& vset ) const
{
  map<string, Index>  properties;
  this->pref_.ListProperties( NODE, properties );
  map<string, Index>  propertiesElement;
  this->pref_.ListProperties( ELEMENT, propertiesElement );
  map<string, Index>  propertiesElementIP;
  this->pref_.ListProperties( ELEMENT_INTEGRATION_POINT, propertiesElementIP );
  properties.insert( propertiesElement.begin(), propertiesElement.end() );
  properties.insert( propertiesElementIP.begin(), propertiesElementIP.end() );

  OutputTo( vset, properties );
}





template<size_t dim>
void Region<dim>::OutputFvDataTo( VSet<dim>& vset ) const
{
  map<string, Index>  properties;
  map<string, Index>  propertiesElementSEIP;
  map<string, Index>  propertiesElementFAIP;
  this->pref_.ListProperties( SECTOR_INTEGRATION_POINT, propertiesElementSEIP );
  properties.insert( propertiesElementSEIP.begin(), propertiesElementSEIP.end() );
  this->pref_.ListProperties( FACET_INTEGRATION_POINT, propertiesElementFAIP );
  properties.insert( propertiesElementFAIP.begin(), propertiesElementFAIP.end() );

  OutputTo( vset, properties );
}





/**
Writes the properties whose names are supplied via a map
into the input VSet container.

@author SK 2016
*/
template<size_t dim>
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
template<size_t dim>
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
    else if ( key.place == ELEMENT ) data.Reserve( this->Elements(), this->Elements() );
    else
      csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:", property, "output is not handled yet." );
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
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             property, "type of variable not recognized." );
      }
    }
    else if ( key.place == ELEMENT ) {
      data.Reserve( this->Elements() );
      switch ( key.type ) {
        case SCALAR: {
          ScalarVariable sc;
          for ( const auto it : this->elmt_vec_ ) {
            (*it).Read( key, sc );
            pushBack( data, sc );
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> vc;
          for ( const auto it : this->elmt_vec_ ) {
            (*it).Read( key, vc );
            pushBack( data, vc );
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> ts;
          for ( const auto it : this->elmt_vec_ ) {
            (*it).Read( key, ts );
            pushBack( data, ts );
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable av;
          for ( const auto it : this->elmt_vec_ ) {
            (*it).Read( key, av );
            pushBack( data, av );
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable fa;
          for ( const auto it : this->elmt_vec_ ) {
            (*it).Read( key, fa );
            pushBack( data, fa );
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             property, "type of variable not recognized." );
      }
    }
    else if ( key.place == NODE ) {
      data.Reserve( this->Nodes() );
      switch ( key.type ) {
        case SCALAR: {
          ScalarVariable sc;
          for ( const auto it : this->node_vec_ ) {
            (*it).Read( key, sc );
            pushBack( data, sc );
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> vc;
          for ( const auto it : this->node_vec_ ) {
            (*it).Read( key, vc );
            pushBack( data, vc );
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> ts;
          for ( const auto it : this->node_vec_ ) {
            (*it).Read( key, ts );
            pushBack( data, ts );
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable av;
          for ( const auto it : this->node_vec_ ) {
            (*it).Read( key, av );
            pushBack( data, av );
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable fa;
          for ( const auto it : this->node_vec_ ) {
            (*it).Read( key, fa );
            pushBack( data, fa );
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             property, "type of variable not recognized." );
      }
    }
    // TODO: Element integration points etc
    // finite volume-related data
    // etc.
    else
      csmp_error.notice( ERROR, property, "Region<dim>::OutputVariableTo:", "variable type not handled yet." );

  } // end not integration point

  return data;

} // end OutputVariableTo(PropertyData)












/**
Outputs specific property data to a FEM_Data container.
This function is a nested template:
the outer template provides double64 = data type and dim = dimension,
and Var the data type of the property
(ScalarVariable, VectorVariable, or TensorVariable).

@return The FEM_Data container with the data that corresponds to the given property.

@section implementation Implementation

Only data pertinent to this local group is output.
A new data container containing only group relevant data is returned.

Thus, the data extracted only corresponds to -this- group and is organized
by local id numbers.

@section application Application

The idea is to extract the property data that corresponds to this group only.
*/
template<size_t dim>
template<class Var>
void Region<dim>::OutputVariableTo( const char* property, FEM_Data<Var>& data ) const
{
  csmp::Index  idx = this->pref_.StorageKey( property );
  Var          var;
  femDataOutputDispatch::initVariable( idx, var );

  switch ( idx.place ) {
    case ELEMENT: {
      data.Reset( idx, this->elmt_vec_.size(), var );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ ) {
        (*eit)->Read( idx, var );
        data[(*eit)->Idx()] = var;
      }
    }
                  break;
    case NODE: {
      data.Reset( idx, this->node_vec_.size(), var );
      for ( typename vector<csmp::Node<dim>*>::const_iterator
            nit = this->node_vec_.begin(); nit != this->node_vec_.end(); nit++ ) {
        (*nit)->Read( idx, var );
        data[(*nit)->Idx()] = var;
      }
    }
               break;
    case ELEMENT_INTEGRATION_POINT: {
      data.Reset( idx, this->IntegrationPoints(), var );
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
        for ( size_t i = 0U; i<(*eit)->IntegrationPoints(); i++ ) {
          (*eit)->Read( i, idx, var );
          data[counter] = var;
          ++counter;
        }
    }
                                    break;
    case SECTOR_INTEGRATION_POINT: {
      data.Reset( idx, this->SectorIntegrationPoints(), var );
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
        for ( size_t j = 0U; j<(*eit)->Sectors(); j++ )
        {
          for ( size_t i = 0U; i<(*eit)->IntegrationPointsPerSector(); i++ ) {
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
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
        for ( size_t j = 0U; j<(*eit)->Facets(); j++ )
        {
          for ( size_t i = 0U; i<(*eit)->IntegrationPointsPerFacet(); i++ ) {
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
      csmp_error.notice( WARNING, "Region::OutputVariableTo",
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
template<size_t dim>
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
    csmp_error.notice( WARNING, "Region<dim>::InputVariableFrom:", warningMessage.c_str() );
    return;
  }

  switch ( idx.place ) {
    case ELEMENT:
      assert( this->Elements() == vdata.Size() );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
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
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
      {
        for ( size_t i = 0U; i<(*eit)->IntegrationPoints(); i++ )
          (*eit)->Store( i, idx, vdata[counter + i] );
        counter += (*eit)->IntegrationPoints();
      }
    }
                                    break;
    case SECTOR_INTEGRATION_POINT: {
      size_t counter( 0U );
      for ( typename vector<csmp::Element<dim>*>::const_iterator
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
      {
        for ( size_t j = 0U; j<(*eit)->Sectors(); j++ )
        {
          for ( size_t i = 0U; i<(*eit)->IntegrationPointsPerSector(); i++ ) {
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
            eit = this->elmt_vec_.begin(); eit != this->elmt_vec_.end(); eit++ )
      {
        for ( size_t j = 0U; j<(*eit)->Facets(); j++ )
        {
          for ( size_t i = 0U; i<(*eit)->IntegrationPointsPerFacet(); i++ ) {
            (*eit)->Store( j, i, idx, vdata[counter] );
            ++counter;
          }
        }
      }
    }
                                  break;
    case REGION:
      if ( vdata.Size() != 1U )
        csmp_error.notice( WARNING, "Region<dim>::InputVariableFrom:",
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


template<size_t dim>
void Region<dim>::CreateNodePointerVector()
{
  assert( !this->elmt_vec_.empty() );

  if ( !this->node_vec_.empty() )
    this->node_vec_.clear();

  // creating the node index vector
  set<csmp::Node<dim>*>  nodes_set;
  for ( typename vector<Element<dim>*>::const_iterator it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
    for ( typename vector<Node<dim>*>::size_type i = 0U; i<(*it)->Nodes(); i++ )
      nodes_set.insert( (*it)->N( i ) );

  this->node_vec_.assign( nodes_set.begin(), nodes_set.end() );
}





/*
In order to use the region node/element flags 'INTERIOR' or 'PERIMETER',
boundary nodes and elements must be identified first.
This is done every time when a new region is formed.

@attention convention: elements are considered 'PERIMETER' elements
of a Region only if they have at least one edge or face at the
region boundary. This is the case irrespective of the dimensionality
of such elements. In the extreme case of line elements inside a volume,
they are flagged boundary, if they have one node on the perimeter of
the volume region.
Another way of looking at this is to consider all elements boundary
elements, that share a face with the boundary.

@attention convention: If a Region consists of elements of different
dimensionality, the highest dimensional ones define the position of the boundary, i.e.
all lower-dimensional elements that stick out of the region (and including all their
nodes) are flagged 'PERIMETER'.
Lower-dimensional elements inside a higher dimensional region are considered
'INTERIOR'.

@section application Application

IdentifyPerimeter() is called as part of the region-forming process.
The perimeter flagging is used to access boundaries selectively.

@todo (1) Does not work for quadratic regions
*/

/* ORIGINAL VERSION

template<size_t dim>
void Region<dim>::IdentifyPerimeter()
{
// 1. putting the boundary elements at the end of the elmt_vec and sorting
//    interior and boundary element ranges subsequently
// ----------------------------------------------------
const size_t first_bd_elmt = PartitionElementVector( "Region");
assert( first_bd_elmt <= this->elmt_vec_.size() );

// removing any potentially pre-existing boundary face information (O.K.)
if ( !this->bd_face_vec_.empty() )
this->bd_face_vec_.clear();
this->bd_face_vec_.reserve( this->elmt_vec_.size() - first_bd_elmt );

// 2. finding boundary faces of boundary elements
//    and storing them in vector with same ordering
// ------------------------------------------------
vector<ONE_BYTE_NUMBER>  bfaces;
enum { MAX_FACES=6 }; // hexahedron
// for all elements on the region boundary
for ( size_t i=first_bd_elmt; i<this->elmt_vec_.size(); i++ ) {
bfaces.reserve(MAX_FACES);
// identify their boundary faces as those who have no neighbor (because on model boundary)
// and those whose neighbor elements do not belong to the region.
for ( size_t j=0U; j<this->elmt_vec_[i]->Neighbors(); j++ ) {
csmp::Element<dim>* const eptr(this->elmt_vec_[i]->Neighbor(j));
// searching the pointers of the group
if ( eptr == NULL or // both range searches (interior and boundary) do not give hits
(!binary_search( this->elmt_vec_.begin(), this->elmt_vec_.begin() + static_cast<long>(first_bd_elmt), eptr ) and
!binary_search( this->elmt_vec_.begin() + static_cast<long>(first_bd_elmt), this->elmt_vec_.end(), eptr )) )
bfaces.push_back( static_cast<ONE_BYTE_NUMBER>(j) );
}
//         assert( !bfaces.empty() );
if ( !bfaces.empty() ) this->bd_face_vec_.push_back( bfaces );
bfaces.clear();
}
vector<vector<ONE_BYTE_NUMBER> >( this->bd_face_vec_ ).swap( this->bd_face_vec_ );

// 3. identifying the boundary nodes
// ---------------------------------
vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit(this->bd_face_vec_.begin());
set<csmp::Node<dim>*>                             boundary_nodes;
vector<size_t>                                    fnids;

// extracting a set of boundary nodes in order to partition the node vector into
// interior- and exterior nodes
for ( typename vector<Element<dim>*>::const_iterator
it=PerimeterElementsBegin(); it!=ElementsEnd(); it++, bit++ )
for ( size_t j=0U; j<(*bit).size(); j++ ) {
(*it)->FE()->NodesOfFace( (*bit)[j], fnids );
for ( size_t k=0U; k<fnids.size(); k++ )
boundary_nodes.insert( (*it)->N(fnids[k]) );
}
this->first_bd_node_ = this->node_vec_.size() - boundary_nodes.size();

// rebuilding and sorting the node vector (set nodes are already sorted)
vector<csmp::Node<dim>*>  temp;
temp.reserve(this->node_vec_.size());
// first, the interior elements are inserted
for ( typename vector<csmp::Node<dim>*>::const_iterator
nit=this->node_vec_.begin(); nit!=this->node_vec_.end(); nit++ )
if ( boundary_nodes.find(*nit) == boundary_nodes.end() )
temp.push_back( *nit );
// now the vector is sorted
sort( temp.begin(), temp.end() );
// second, the already sorted boundary nodes are appended
for ( typename set<csmp::Node<dim>*>::const_iterator
nit=boundary_nodes.begin(); nit!=boundary_nodes.end(); nit++ )
temp.push_back( *nit );
// now the temporary vector is assigned to the permanent one
this->node_vec_ = temp;

} // end IdentifyPerimeter
*/




/**
Detects of how many spatial dimensions element types are contained in model.
It returns a pair: first value gives number of different spatial dimensions contained,
second value returns the highest spatial dimension contained.

@author SKM 1/11/2013
*/
template<size_t dim>
pair<int32, int32>  Region<dim>::ElementSpatialDimensions() const
{
  return this->SpatialDimensions();

} // end ElementSpatialDimensions






/**
Creates set of lower dimensional elements that have a face at the boundary of the
model subdomain. Pointers to these elements are stored in supplied set,
and their number is returned.

@attention SKM method is not implemented yet.

template<size_t dim>
size_t  Region<dim>::IdentifyLowerDimensionalBoundaryElements( const std::pair<int32,int32>&,
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
SKM trying to make sense of Andrew Bromage's undocumented code:
12/08/18
*/
template<size_t dim>
size_t Region<dim>::FromLargestComponent( MeshManager<dim>& mesh,
                                          bool reestablishNeighborConnectivity )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->elmt_vec_.empty() )
    csmp_error.notice( WARNING, "Region<dim>::FromLargestComponent:",
                       "Region is not empty; deleting all content." );
  this->elmt_vec_.clear();

  // traversal of the existing mesh nodes to find all its elements	
  deque<csmp::Node<dim>*>	nodes;
  deque<csmp::Element<dim>*>	elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
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
    csmp_error.notice( ERROR, "Region<dim>::FromLargestComponent:",
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
    this->elmt_vec_.swap( elmts );
  }

  // 5. (re)connecting elements up to their neighbors
  //    TODO: this is a very time-consuming step; is there a speed-up?
  if ( reestablishNeighborConnectivity )
    this->EstablishNeighborConnectivity();

  // 6. identifying the boundaries
  this->IdentifyPerimeter();

  return this->elmt_vec_.size();
} // end FromLargestComponent













/**
Visit all elements without relying on their storage in a container.
The Idx numbering of elements and nodes is not altered by this method.

@return number of elements that were discovered.

@attention, this method always gets called when a region is first formed.

@attention nodes must have been assigned their parent elements for this method to work.

@author SKM 9/20/2008, CSMP Castasegna workshop, Switzerland.
*/
template<size_t dim>
size_t Region<dim>::AccumulateAll( const csmp::Node<dim>* root_node,
                                   bool reestablishNeighborConnectivity )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( root_node == NULL )
    csmp_error.notice( FATAL_ERROR, "Region<dim>::AccumulateAll:", "Root node pointer is dangling!" );

  if ( root_node->Parent( 0 ) == NULL )
    csmp_error.notice( FATAL_ERROR, "Region<dim>::AccumulateAll:",
                       "Root node must have been assigned parent elements; else this method cannot operate." );

  if ( !this->elmt_vec_.empty() )
    csmp_error.notice( WARNING, "Region<dim>::AccumulateAll:",
                       "Region is not empty; deleting all content." );
  this->elmt_vec_.clear();

  // 1. traversal of the existing mesh nodes to find all its elements
  set<csmp::Element<dim>*>       explored_elements;
  set<const csmp::Node<dim>*>    discovered_nodes;
  deque<const csmp::Node<dim>*>  current_nodes;
  // starting at the root element
  discovered_nodes.insert( root_node );
  current_nodes.push_back( root_node );

  // MESH TRAVERSAL
  while ( !current_nodes.empty() ) {
    const csmp::Node<dim>*  n_ptr( *current_nodes.begin() );
    // for all parent elements of the current node
    for ( size_t i = 0U; i<n_ptr->Parents(); i++ ) {
      // for all the nodes of each parent element
      for ( size_t j = 0U; j<n_ptr->Parent( i )->Nodes(); j++ )
        // if this node is not the one from which we started
        if ( j != n_ptr->ParentNodeNumber( i ) ) {
          pair<typename set<const csmp::Node<dim>*>::iterator, bool>
            new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
          if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
        }
      // storing the explored element
      explored_elements.insert( n_ptr->Parent( i ) );
    }
    // removing the node from the discovered (but not yet explored) deque
    current_nodes.pop_front();
  }

  // 2. assigning and trimming excess storage from the element pointer vector
  this->elmt_vec_.assign( explored_elements.begin(), explored_elements.end() );
  vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  // 3. creating pointers to the nodes of the identified elements
  set<csmp::Node<dim>*>  node_set;
  for ( typename vector<csmp::Element<dim>*>::const_iterator
        it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
    for ( typename vector<csmp::Node<dim>*>::const_iterator
          nit = (*it)->NodesBegin(); nit != (*it)->NodesEnd(); nit++ )
      node_set.insert( (*nit) );

  // TODO: use emplace here?
  this->node_vec_.assign( node_set.begin(), node_set.end() );

  // 4. (re)connecting elements up to their neighbors
  //    TODO: this is a very time-consuming step; is there a speed-up?
  if ( reestablishNeighborConnectivity )
    this->EstablishNeighborConnectivity();

  // 5. identifying the boundaries
  this->IdentifyPerimeter();

  return this->elmt_vec_.size();

} // end AccumulateAll






/**
Accumulates range of elements supplied.
*/

template<size_t dim>
void  Region<dim>::Accumulate( typename deque<csmp::Element<dim> >::iterator start,
                               typename deque<csmp::Element<dim> >::iterator end )
{
  if ( start == end )
    throw csmp::Exception( ERROR, "Region<dim>::Accumulate (deque)",
                           "supplied element range is empty. Nothing is done." );

  this->elmt_vec_.clear();
  this->elmt_vec_.reserve( static_cast<size_t>(distance( start, end )) );

  set<csmp::Node<dim>*>  node_set;
  for ( typename deque<csmp::Element<dim> >::iterator
        it = start; it != end; it++ )
  {
    this->elmt_vec_.push_back( &(*it) );
    for ( typename vector<csmp::Node<dim>*>::const_iterator
          nit = (*it).NodesBegin(); nit != (*it).NodesEnd(); nit++ )
      node_set.insert( (*nit) );
  }

  vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  this->node_vec_.assign( node_set.begin(), node_set.end() );

  this->IdentifyPerimeter();

} // end Accumulate (deque)


template<size_t dim>
void  Region<dim>::Accumulate( typename vector<csmp::Element<dim>*>::const_iterator start,
                               typename vector<csmp::Element<dim>*>::const_iterator end )
{
  if ( start == end )
    throw csmp::Exception( ERROR, "Region<dim>::Accumulate (vector)",
                           "supplied element range is empty. Nothing is done." );

  this->elmt_vec_.clear();
  this->elmt_vec_.reserve( static_cast<size_t>(distance( start, end )) );

  // making sure that there a no duplicate element pointers
  unique_copy( start, end, back_inserter( this->elmt_vec_ ) );

  vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  set<csmp::Node<dim>*>  node_set;

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
    for ( typename vector<csmp::Node<dim>*>::const_iterator
          nit = (*it)->NodesBegin(); nit != (*it)->NodesEnd(); nit++ )
      node_set.insert( (*nit) );

  this->node_vec_.assign( node_set.begin(), node_set.end() );
  this->IdentifyPerimeter();

} // end Accumulate (vector)




template<size_t dim>
void  Region<dim>::Accumulate( typename set<csmp::Element<dim>*>::const_iterator start,
                               typename set<csmp::Element<dim>*>::const_iterator end )
{
  if ( start == end )
    throw csmp::Exception( ERROR, "Region<dim>::Accumulate (set)",
                           "supplied element range is empty. Nothing is done." );
  this->elmt_vec_.clear();
  this->elmt_vec_.reserve( static_cast<size_t>(distance( start, end )) );

  while ( start != end ) {
    this->elmt_vec_.push_back( *start );
    start++;
  }

  vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  set<csmp::Node<dim>*>  node_set;

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ )
    for ( typename vector<csmp::Node<dim>*>::const_iterator
          nit = (*it)->NodesBegin(); nit != (*it)->NodesEnd(); nit++ )
      node_set.insert( (*nit) );

  this->node_vec_.assign( node_set.begin(), node_set.end() );
  this->IdentifyPerimeter();

} // end Accumulate (set)



/**
All those elements for which all property constraints are met are accumulated
into this region.

Use this also to modify an existing group.

@code
mesh<Element<dim>*>::iterator
@endcode
*/
template<size_t dim>
void Region<dim>::AccumulateWithinRange( typename vector<csmp::Element<dim>*>::const_iterator start,
                                         typename vector<csmp::Element<dim>*>::const_iterator end,
                                         const PropertyConstraints& constraints )
{
  if ( start == end )
    throw csmp::Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                           "supplied element range is empty. Nothing is done." );

  if ( constraints.Constraints() == 0U )
    throw csmp::Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                           "No property constraints are supplied" );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->elmt_vec_.empty() )
    csmp_error.notice( WARNING, "Region<dim>::AccumulateWithinRange",
                       "Region<dim> already contains elements, they will be deleted" );

  if ( !this->elmt_vec_.empty() ) this->elmt_vec_.clear();
  if ( !this->node_vec_.empty() ) this->node_vec_.clear();
  if ( !this->bd_face_vec_.empty() ) this->bd_face_vec_.clear();

  set<csmp::Element<dim>*>  element_set;
  set<csmp::Node<dim>*>     node_set;

  while ( start != end ) {
    if ( constraints.CheckConstraints( *(*start) ) )
    {
      element_set.insert( (*start) );
      for ( size_t i = 0U; i<(*start)->Nodes(); i++ )
        node_set.insert( (*start)->N( i ) );
    }
    start++;
  }

  if ( !element_set.empty() )
  {
    this->node_vec_.assign( node_set.begin(), node_set.end() );
    this->elmt_vec_.assign( element_set.begin(), element_set.end() );
    this->IdentifyPerimeter();
  }

} // end AccumulateWithinRange(PropertyConstraints)


/**
Accumulates those property-bearing elements into a region, whose property values are greater or
equal to the lower bound of the specified range and smaller or equal to the
upper bound of it.

@attention if the target variable is a node property, by default at least one of the element's
nodes must have a value inside of the target range.
*/
template<size_t dim>
void Region<dim>::AccumulateWithinRange( typename vector<Element<dim>*>::const_iterator start,
                                         typename vector<Element<dim>*>::const_iterator end,
                                         const char* feature, double64 min, double64 max )
{
  if ( start == end )
    throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                     "supplied element range is empty. Nothing is done." );

  set<Element<dim>*>  element_set;
  set<Node<dim>*>     node_set;
  Index                     prop_key = this->pref_.StorageKey( feature );
  bool                      applies;

  switch ( prop_key.place )
  {
    case NODE:
      while ( start != end ) {
        applies = false;
        for ( size_t j = 0U; j<(*start)->Nodes(); j++ )
          if ( (*start)->N( j )->IsWithinRange( prop_key, min, max ) ) {
            applies = true;
            break;
          }

        if ( applies == true ) {
          element_set.insert( (*start) );
          for ( size_t i = 0U; i<(*start)->Nodes(); i++ ) {
            assert( (*start)->N( i ) != nullptr );
            node_set.insert( (*start)->N( i ) );
          }
        }
        start++;
      }
      break;
    case ELEMENT_INTEGRATION_POINT:
      while ( start != end ) {
        applies = false;
        for ( size_t j = 0U; j<(*start)->IntegrationPoints(); j++ )
          if ( (*start)->IsWithinRange( j, prop_key, min, max ) ) {
            applies = true;
            break;
          }

        if ( applies == true ) {
          element_set.insert( (*start) );
          for ( size_t i = 0U; i<(*start)->Nodes(); i++ ) {
            assert( (*start)->N( i ) != nullptr );
            node_set.insert( (*start)->N( i ) );
          }
        }
        start++;
      }
      break;
    case ELEMENT:
      while ( start != end ) {
        if ( (*start)->IsWithinRange( prop_key, min, max ) ) {
          element_set.insert( (*start) );
          for ( size_t i = 0U; i<(*start)->Nodes(); i++ ) {
            assert( (*start)->N( i ) != nullptr );
            node_set.insert( (*start)->N( i ) );
          }
        }
        start++;
      }
      break;
    default:
      throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                       feature, "placement could not be identified; REGION is not an option" );
  } // end switch

  if ( !element_set.empty() ) {
    this->node_vec_.assign( node_set.begin(), node_set.end() );
    this->elmt_vec_.assign( element_set.begin(), element_set.end() );
    this->IdentifyPerimeter();
  }

} // end AccumulateWithinRange


template<size_t dim>
void Region<dim>::AccumulateWithinRange( const MeshManager<dim>& mesh, const PropertyConstraints& constraints )
{
  if ( mesh.Elements() < 1U )
    throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                     "supplied element range is empty. Nothing is done." );

  if ( constraints.Constraints() == 0U )
    throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                     "No property constraints are supplied" );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->elmt_vec_.empty() )
    csmp_error.notice( WARNING, "Region<dim>::AccumulateWithinRange",
                       "Region<dim> already contains elements, they will be deleted" );

  if ( !this->elmt_vec_.empty() ) this->elmt_vec_.clear();
  if ( !this->node_vec_.empty() ) this->node_vec_.clear();
  if ( !this->bd_face_vec_.empty() ) this->bd_face_vec_.clear();

  // traversal of the existing mesh nodes to find all its elements	
  deque<const csmp::Node<dim>*>	nodes;
  deque<csmp::Element<dim>*>		elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  set<Element<dim>*>  element_set;
  set<Node<dim>*>     node_set;

  for ( auto start : elmts ) {
    if ( constraints.CheckConstraints( *start ) )
    {
      element_set.insert( start );
      for ( size_t i = 0U; i<start->Nodes(); i++ )
        node_set.insert( start->N( i ) );
    }
  }

  if ( !element_set.empty() )
  {
    this->node_vec_.assign( node_set.begin(), node_set.end() );
    this->elmt_vec_.assign( element_set.begin(), element_set.end() );
    this->IdentifyPerimeter();
  }

} // end AccumulateWithinRange(PropertyConstraints)



template<size_t dim>
void Region<dim>::AccumulateRectangularRegion( typename vector<Element<dim>*>::const_iterator start,
                                               typename vector<Element<dim>*>::const_iterator end,
                                               const Point<dim>& xyz_min,
                                               const Point<dim>& xyz_max )
{
  if ( start == end )
    throw Exception( ERROR, "Region<dim>::AccumulateRectangularRegion",
                     "supplied element range is empty. Nothing is done." );

  size_t                    check;
  set<Element<dim>*>  element_set;
  set<Node<dim>*>     node_set;

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->elmt_vec_.empty() ) {
    csmp_error.notice( WARNING, "Region<dim>::AccumulateRectangularRegion",
                       "Region already contains elements" );
    this->elmt_vec_.clear();
    this->node_vec_.clear();
  }

  while ( start != end ) {
    check = 0U;
    // all nodes have to be inside for the element selection criterion to be fulfilled
    for ( size_t j = 0U; j<(*start)->Nodes(); j++ ) {
      Point<dim>  p = (*start)->N( j )->Coordinate();
      if ( p.IsBetween( xyz_min, xyz_max ) ) check++;
    }
    // if all nodes are inside the rectangular region
    // the element becomes part of the new group
    if ( check == (*start)->Nodes() ) {
      element_set.insert( (*start) );
      for ( size_t i = 0U; i < (*start)->Nodes(); i++ ) {
        node_set.insert( (*start)->N( i ) );
      }
    }
    start++;
  }

  if ( !element_set.empty() ) {
    this->node_vec_.assign( node_set.begin(), node_set.end() );
    this->elmt_vec_.assign( element_set.begin(), element_set.end() );
    this->IdentifyPerimeter();
  }

} // end AccumulateRectangularRegion

template<size_t dim>
void Region<dim>::AccumulateRectangularRegion( const MeshManager<dim>& mesh, const Point<dim>& xyz_min, const Point<dim>& xyz_max )
{
  if ( mesh.Elements() < 1U )
    throw Exception( ERROR, "Region<dim>::AccumulateRectangularRegion",
                     "supplied element range is empty. Nothing is done." );

  size_t                    check;
  set<Element<dim>*>  element_set;
  set<Node<dim>*>     node_set;

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->elmt_vec_.empty() ) {
    csmp_error.notice( WARNING, "Region<dim>::AccumulateRectangularRegion",
                       "Region already contains elements" );
    this->elmt_vec_.clear();
    this->node_vec_.clear();
  }

  // traversal of the existing mesh nodes to find all its elements	
  deque<const csmp::Node<dim>*>	nodes;
  deque<csmp::Element<dim>*>		elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  for ( auto e : elmts ) {
    check = 0U;
    // all nodes have to be inside for the element selection criterion to be fulfilled
    for ( size_t j = 0U; j<e->Nodes(); j++ ) {
      Point<dim>  p = e->N( j )->Coordinate();
      if ( p.IsBetween( xyz_min, xyz_max ) ) check++;
    }
    // if all nodes are inside the rectangular region
    // the element becomes part of the new group
    if ( check == e->Nodes() ) {
      element_set.insert( e );
      for ( size_t i = 0U; i < e->Nodes(); i++ ) {
        node_set.insert( e->N( i ) );
      }
    }
  }

  if ( !element_set.empty() ) {
    this->node_vec_.assign( node_set.begin(), node_set.end() );
    this->elmt_vec_.assign( element_set.begin(), element_set.end() );
    this->IdentifyPerimeter();
  }

} // end AccumulateRectangularRegion

template<size_t dim>
void Region<dim>::AccumulateWithinRange( const MeshManager<dim>& mesh, const char* feature, double64 min, double64 max )
{
  if ( mesh.Elements() < 1U )
    throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                     "supplied element range is empty. Nothing is done." );

  // traversal of the existing mesh nodes to find all its elements	
  deque<const csmp::Node<dim>*>	nodes;
  deque<csmp::Element<dim>*>		elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  set<Element<dim>*>  element_set;
  set<Node<dim>*>     node_set;
  Index                     prop_key = this->pref_.StorageKey( feature );
  bool                      applies;

  switch ( prop_key.place )
  {
    case NODE:
      for ( auto e : elmts ) {
        applies = false;
        for ( size_t j = 0U; j<e->Nodes(); j++ )
          if ( e->N( j )->IsWithinRange( prop_key, min, max ) ) {
            applies = true;
            break;
          }

        if ( applies == true ) {
          element_set.insert( e );
          for ( size_t i = 0U; i<e->Nodes(); i++ ) {
            assert( e->N( i ) != nullptr );
            node_set.insert( e->N( i ) );
          }
        }
      }
      break;
    case ELEMENT_INTEGRATION_POINT:
      for ( auto e : elmts ) {
        applies = false;
        for ( size_t j = 0U; j<e->IntegrationPoints(); j++ )
          if ( e->IsWithinRange( j, prop_key, min, max ) ) {
            applies = true;
            break;
          }

        if ( applies == true ) {
          element_set.insert( e );
          for ( size_t i = 0U; i<e->Nodes(); i++ ) {
            assert( e->N( i ) != nullptr );
            node_set.insert( e->N( i ) );
          }
        }
      }
      break;
    case ELEMENT:
      for ( auto e : elmts ) {
        if ( e->IsWithinRange( prop_key, min, max ) ) {
          element_set.insert( e );
          for ( size_t i = 0U; i<e->Nodes(); i++ ) {
            assert( e->N( i ) != nullptr );
            node_set.insert( e->N( i ) );
          }
        }
      }
      break;
    default:
      throw Exception( ERROR, "Region<dim>::AccumulateWithinRange",
                       feature, "placement could not be identified; REGION is not an option" );
  } // end switch

  if ( !element_set.empty() ) {
    this->node_vec_.assign( node_set.begin(), node_set.end() );
    this->elmt_vec_.assign( element_set.begin(), element_set.end() );
    this->IdentifyPerimeter();
  }

} // end AccumulateWithinRange


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

Method will detect if the supplied vector<double64> is empty or if a group by
that name already exists.
*/
template<size_t dim>
void  Region<dim>::AccumulateByNumber( typename vector<Element<dim>*>::const_iterator start,
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

  if ( !this->elmt_vec_.empty() ) {
    csmp_error.notice( WARNING, "Region<dim>::AccumulateByNumber",
                       "Region is not empty", "erasing all members..." );
    this->elmt_vec_.clear();
  }

  sort( element_ids.begin(), element_ids.end() );

  // if in debug mode, tests whether there are consecutive duplicated elements
  vector<size_t>::iterator  new_end( unique( element_ids.begin(), element_ids.end() ) );
  if ( new_end != element_ids.end() )
    element_ids.erase( new_end, element_ids.end() );

  if ( element_ids.size() > static_cast<size_t>(distance( start, end )) )
    csmp_error.notice( ERROR, "Region<dim>::AccumulateByNumber",
                       "user-supplied element-number vector is larger than iterator range." );

  this->elmt_vec_.reserve( element_ids.size() );
  set<Node<dim>*>  node_set;

  // selecting elements and nodes from the selected ID range
  while ( start != end )
  {
    assert( *start != NULL );
    if ( binary_search( element_ids.begin(), element_ids.end(), (*start)->Idx() ) ) {
      // add element with the correct id to the region
      this->elmt_vec_.push_back( (*start) );
      // add its nodes as well
      for ( size_t i = 0U; i<(*start)->Nodes(); i++ )
        node_set.insert( (*start)->N( i ) );
    }
    start++;
  }

  // TODO: use emplace here?
  this->node_vec_.assign( node_set.begin(), node_set.end() );
  assert( !this->node_vec_.empty() );

  this->IdentifyPerimeter();

} // end AccumulateByNumber

template<size_t dim>
void  Region<dim>::AccumulateByNumber( const MeshManager<dim>& mesh, std::vector<size_t>& element_ids )
{
  if ( mesh.Elements() < 1U )
    throw Exception( ERROR, "Region<dim>::AccumulateByNumber",
                     "user-supplied iterator range is empty. Nothing is done." );

  if ( element_ids.empty() )
    throw Exception( ERROR, "Region<dim>::AccumulateByNumber",
                     "user-supplied element-number vector is empty. Nothing is done." );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !this->elmt_vec_.empty() ) {
    csmp_error.notice( WARNING, "Region<dim>::AccumulateByNumber",
                       "Region is not empty", "erasing all members..." );
    this->elmt_vec_.clear();
  }

  sort( element_ids.begin(), element_ids.end() );

  // if in debug mode, tests whether there are consecutive duplicated elements
  vector<size_t>::iterator  new_end( unique( element_ids.begin(), element_ids.end() ) );
  if ( new_end != element_ids.end() )
    element_ids.erase( new_end, element_ids.end() );

  if ( element_ids.size() > mesh.Elements() )
    csmp_error.notice( ERROR, "Region<dim>::AccumulateByNumber",
                       "user-supplied element-number vector is larger than iterator range." );

  this->elmt_vec_.reserve( element_ids.size() );
  set<Node<dim>*>  node_set;

  // selecting elements and nodes from the selected ID range
  // traversal of the existing mesh nodes to find all its elements	
  deque<const csmp::Node<dim>*>	nodes;
  deque<csmp::Element<dim>*>		elmts;
  exploreNodesAndElementsFromMesh( &mesh, nodes, elmts );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  for ( auto e : elmts ) {
    assert( e != NULL );
    if ( binary_search( element_ids.begin(), element_ids.end(), e->Idx() ) ) {
      // add element with the correct id to the region
      this->elmt_vec_.push_back( e );
      // add its nodes as well
      for ( size_t i = 0U; i<e->Nodes(); i++ )
        node_set.insert( e->N( i ) );
    }
  }

  this->node_vec_.assign( node_set.begin(), node_set.end() );
  assert( !this->node_vec_.empty() );

  this->IdentifyPerimeter();

} // end AccumulateByNumber




/**
Adds supplied region to the current region.
*/
template<size_t dim>
void  Region<dim>::Add( const Region<dim>&  grp )
{
  // if the added region is empty nothing needs to be done
  if ( grp.Empty() ) {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( WARNING, "Region<dim>::Add", "Region to add was empty; nothing was done" );
    return;
  }

  set<Element<dim>*>  elmt_set( this->elmt_vec_.begin(), this->elmt_vec_.end() );

  for ( typename vector<csmp::Element<dim>*>::const_iterator
        it = grp.ElementsBegin(); it != grp.ElementsEnd(); it++ )
    elmt_set.insert( *it );

  this->elmt_vec_.assign( elmt_set.begin(), elmt_set.end() );

  vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  this->CreateNodePointerVector();
  this->IdentifyPerimeter();

} // end Add




/**
Method builds a new region of elements between the argument regions.

@author P. Lang
@author R. Manasipov
@date 22/8/2014

This method assumes that elements are uniquely and throughgoingly numbered.
The Elements are build so that their normals point from region 1 to region 2.
Variable storage is assigned for both, the elements and region itself.
*/
template<size_t dim>
bool Region<dim>::CreateBetween( MeshManager<dim>& meshManager,
                                 const FiniteElementManager& finiteElementManager,
                                 const Region<dim>& region1,
                                 const Region<dim>& region2 )
{

  // LVS
  const LocalVariables lvsElements( ElementVariables() );
  const IntegrationPointVariables lvsIntegrationPoints( ElementIntegrationPointVariables() );

  // pointer to running element, fem type of new faces
  Element<dim>*       ePtr( NULL );
  Element<dim>*       ePtrNeighbor( NULL );
  FiniteElement*      femPtr( NULL );
  std::vector<size_t> face_node_ids;

  // reserving storage for boundary elements (logic: the number of faces created cannot be
  // larger than the minimum number boundary elements of the two neighboring groups)
  this->elmt_vec_.reserve( std::min( region1.PerimeterElements(),
                           region2.PerimeterElements() ) );

  // searching for elements of region1 that are neighbors of ones in group1.
  // If so, there is a shared boundary and faces or interfaces are constructed.
  const size_t n_elements( region1.Elements() );
  for ( size_t i = region1.InteriorElements(); i < n_elements; ++i )
  {
    ePtr = region1.E( i );
    const size_t perimeter_faces( region1.PerimeterFaces( i ) );
    for ( size_t j = 0U; j < perimeter_faces; ++j )
    {
      const size_t face = region1.PerimeterFace( i, j );
      ePtrNeighbor = ePtr->Neighbor( face );

      // checking whether neighbor element forms part of the boundary of group2
      if ( ePtrNeighbor != NULL )
        if ( region2.IsPerimeterElement( ePtrNeighbor ) )
        {
          // if the neighbor is in the boundary, the new Face is build
          femPtr = finiteElementManager.E( ePtr->FE()->ElementTypeOfFace( face ) );

          // getting the mesh manager to construct a new element
          Element<dim> new_elmt( femPtr, NULL, lvsElements, lvsIntegrationPoints );
          Element<dim>* elmtObj = meshManager.AddIfUnique( new_elmt );

          // the new face is connected to the elements it is sandwiched between
          // this assignment also includes connecting the element to its nodes
          //                 inner        outer  element w.r.t. to normal of face
          face_node_ids.clear();
          face_node_ids.resize( ePtr->FE()->NodesPerFace( face ) );
          ePtr->FE()->NodesOfFace( face, face_node_ids );
          for ( size_t nid = 0U; nid < face_node_ids.size(); ++nid )
            elmtObj->Assign( nid, ePtr->N( face_node_ids[i] ) );

          // added to boundary
          this->elmt_vec_.emplace_back( elmtObj );

        } // neighboring elements

    } // perimeter faces

  } // perimeter elements

    // free
  vector<csmp::Element<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

  this->CreateNodePointerVector();
  establishNeighborConnectivity( this->elmt_vec_ );
  this->IdentifyPerimeter();
  this->UpdateMemberIndexes();

  //done
  return true;

} // CreateBetween







/**
Forms the union of group a and group b and returns it into the result
group 'res'.

@section arguments Input Arguments
The two groups which will be turned into one.

@attention this method cannot be called union because 'union' is a keyword
in the C++ language.
*/
template<size_t dim>
size_t  groupUnion( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.notice( WARNING, "groupUnion", "first region is empty" );
  if ( b.Empty() )
    csmp_error.notice( WARNING, "groupUnion", "second region is empty" );

  merge( a.ElementsBegin(), a.PerimeterElementsBegin(),
         b.ElementsBegin(), b.PerimeterElementsBegin(),
         back_inserter( res.ElementVector() ) );

  merge( a.PerimeterElementsBegin(), a.ElementsEnd(),
         b.PerimeterElementsBegin(), b.ElementsEnd(),
         back_inserter( res.ElementVector() ) );

  vector<csmp::Element<dim>*>( res.ElementVector() ).swap( res.ElementVector() );

  res.CreateNodePointerVector();
  res.IdentifyPerimeter();

  return res.Elements();

} // end groupUnion




/**
Forms the intersection group 'res' of the groups 'a' and 'b'.
*/
template<size_t dim>
size_t  intersection( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.notice( WARNING, "intersection", "first region is empty" );
  if ( b.Empty() )
    csmp_error.notice( WARNING, "intersection", "second region is empty" );

  // because region vectors are sorted into 2 seperate ranges,
  // new sorted vectors spanning the whole ranges need to be established first
  vector<Element<dim>*> tmp_elmt_vec1( a.ElementsBegin(), a.ElementsEnd() );
  vector<Element<dim>*> tmp_elmt_vec2( b.ElementsBegin(), b.ElementsEnd() );
  // sorting the vectors
  sort( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end() );
  sort( tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

  // see chapter 16, Nelson STL Programmer's guide, p. 559
  set_intersection( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end(),
                    tmp_elmt_vec2.begin(), tmp_elmt_vec2.end(),
                    back_inserter( res.ElementVector() ) );

  if ( !res.Empty() ) {
    vector<csmp::Element<dim>*>( res.ElementVector() ).swap( res.ElementVector() );
    res.CreateNodePointerVector();
    res.IdentifyPerimeter();
  }

  return res.Elements();

} // end intersection





/**
Returns into Region 'res' those elements of Region 'a' which are not shared
by Regions 'a' and 'b'.
*/
template<size_t dim>
size_t  difference( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.notice( WARNING, "difference", "first region is empty" );
  if ( b.Empty() )
    csmp_error.notice( WARNING, "difference", "second region is empty" );

  // because region vectors are sorted into 2 seperate ranges,
  // new sorted vectors spanning the whole ranges need to be established first
  vector<Element<dim>*> tmp_elmt_vec1( a.ElementsBegin(), a.ElementsEnd() );
  vector<Element<dim>*> tmp_elmt_vec2( b.ElementsBegin(), b.ElementsEnd() );
  // sorting the vectors
  sort( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end() );
  sort( tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

  // see chapter 16, Nelson STL Programmer's guide, p. 559
  set_difference( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end(),
                  tmp_elmt_vec2.begin(), tmp_elmt_vec2.end(),
                  back_inserter( res.ElementVector() ) );

  if ( !res.Empty() ) {
    vector<csmp::Element<dim>*>( res.ElementVector() ).swap( res.ElementVector() );
    res.CreateNodePointerVector();
    res.IdentifyPerimeter();
  }

  return res.Elements();

} // end difference






/**
Returns into Region 'res' those elements of Regions 'a' and 'b' which are
not shared by Regions 'a' and 'b'.
*/
template<size_t dim>
size_t  symmetricDifference( const Region<dim>& a, const Region<dim>& b, Region<dim>& res )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( a.Empty() )
    csmp_error.notice( WARNING, "symmetricDifference", "first region is empty" );
  if ( b.Empty() )
    csmp_error.notice( WARNING, "symmetricDifference", "second region is empty" );

  // because region vectors are sorted into 2 seperate ranges,
  // new sorted vectors spanning the whole ranges need to be established first
  vector<Element<dim>*> tmp_elmt_vec1( a.ElementsBegin(), a.ElementsEnd() );
  vector<Element<dim>*> tmp_elmt_vec2( b.ElementsBegin(), b.ElementsEnd() );
  // sorting the vectors
  sort( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end() );
  sort( tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

  set_symmetric_difference( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end(),
                            tmp_elmt_vec2.begin(), tmp_elmt_vec2.end(),
                            back_inserter( res.ElementVector() ) );
  if ( !res.Empty() ) {
    vector<csmp::Element<dim>*>( res.ElementVector() ).swap( res.ElementVector() );
    res.CreateNodePointerVector();
    res.IdentifyPerimeter();
  }
  return res.Elements();

} // end symmetricDifference







/**
Counts and returns the number of Elements shared between the two regions.
*/
template<size_t dim>
size_t  sharedElements( const Region<dim>& g1, const Region<dim>& g2 )
{
  typename vector<Element<dim>*>::const_iterator  first1( g1.ElementsBegin() );
  typename vector<Element<dim>*>::const_iterator  first2( g2.ElementsBegin() );
  size_t shared_elements( 0U );

  // comparing the interior nodes
  while ( first1 != g1.PerimeterElementsBegin() and first2 != g2.PerimeterElementsBegin() )
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
  first1 = g1.PerimeterElementsBegin();
  first2 = g2.PerimeterElementsBegin();

  while ( first1 != g1.ElementsEnd() and first2 != g2.ElementsEnd() )
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
template<size_t dim>
bool Region<dim>::Includes( const Region<dim>& g ) const
{
  // because region vectors are sorted into 2 seperate ranges,
  // new sorted vectors spanning the whole ranges need to be established first
  // before includes() can be called
  vector<Element<dim>*> tmp_elmt_vec1( this->elmt_vec_ );
  vector<Element<dim>*> tmp_elmt_vec2( g.elmt_vec_ );
  // sorting the vectors
  sort( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end() );
  sort( tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

  return includes( tmp_elmt_vec1.begin(), tmp_elmt_vec1.end(),
                   tmp_elmt_vec2.begin(), tmp_elmt_vec2.end() );

} // end Includes

bool  isOfLowerDimensionalRepresentation( const Element<3U>& element )
{
  return !element.FE()->IsVolumeElement();
}

bool  isOfLowerDimensionalRepresentation( const Element<2U>& element )
{
  return !element.FE()->IsSurfaceElement();
}

template<>
bool  isOfLowerDimensionalRepresentation<1>( const Region<1>& )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  csmp_error.notice( csmp::ERROR, "isOfLowerDimensionalRepresentation", "na for 1D" );
  return false;
}

template<size_t dim>
bool  isOfLowerDimensionalRepresentation( const Region<dim>& region )
{
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator it = region.ElementsBegin(); it != elementsEnd; ++it )
  {
    if ( !isOfLowerDimensionalRepresentation( *(*it) ) )
      return false;
  }
  return true;
}

template<size_t dim>
bool  containsVolumeElements( const Region<dim>& region )
{
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator it = region.ElementsBegin(); it != elementsEnd; ++it )
  {
    if ( (*it)->FE()->IsVolumeElement() )
      return true;
  }
  return false;
}

template<size_t dim>
bool  containsSurfaceElements( const Region<dim>& region )
{
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator it = region.ElementsBegin(); it != elementsEnd; ++it )
  {
    if ( (*it)->FE()->IsSurfaceElement() )
      return true;
  }
  return false;
}

template<size_t dim>
bool  containsLineElements( const Region<dim>& region )
{
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator it = region.ElementsBegin(); it != elementsEnd; ++it )
  {
    if ( (*it)->FE()->IsLineElement() )
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

template bool isOfLowerDimensionalRepresentation( const Region<3>& );
template bool isOfLowerDimensionalRepresentation( const Region<2>& );

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

@param vector_variable The name of the nodal vector<double64> variable which holds the node
coordinate displacement.

@section implementation Implementation

After the assignment of the nodal displacements, MoveNodeCoordinatesBy()
will set the Element state variable 'shape_to_data' to false. This will
prompt re-calculation of elemental volumes if these are queried in
subsequent computations.

@section application Application

If a mesh shall be deformed using the diplacements of a deformation
calculation stored in a vector<double64> variable, ChangeNodeCoordinatesTo() can
be used displace the node coordinates by these displacements.

@section messages Messages

Due to the total garbage results that may arise,
MoveNodeCoordinatesBy() will halt the simulation reporting a fatal
error, if the target property is node a node or vector<double64> type property.
*/
template<size_t dim>
void Region<dim>::MoveNodeCoordinatesBy( const char* vector_variable )
{
  csmp::Index          prop_key = this->pref_.StorageKey( vector_variable );
  VectorVariable<dim>  vc;

  if ( prop_key.type != VECTOR )
    throw csmp::Exception( ERROR, "Model<dim>::MoveNodeCoordinatesBy",
                           "Only vector<double64> variables can be added to coordinates" );

  if ( prop_key.place != NODE )
    throw csmp::Exception( ERROR, "Model<dim>::MoveNodeCoordinatesBy",
                           "Only NODE variables can be added to NODE coordinates" );

  for ( typename vector<csmp::Node<dim>*>::iterator
        nit = this->NodesBegin(); nit != this->NodesEnd(); nit++ ) {
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
template<size_t dim>
void Region<dim>::CorrectLowDimRegionOrientation(  ) const
{

// pair contains: 1) number of spatial element dimensions in region, 2)  highest contained dimension
// TODO: improved version should include more information about perimeter elements ( stickin out elements )
pair<int32,int32>  elmt_dim = this->ElementSpatialDimensions();

// Correct Orientation of low dim region

std::set<Element<dim>* > elements_considered;
std::set<Element<dim>* > element_neighbors_to_be_considered;
std::set<Element<dim>* > remaining_elements( this->ElementsBegin(), this->ElementsEnd() );

Element<dim>*   lowDimElement         ( NULL );
Element<dim>*   lowDimNeighborElement ( NULL );

VectorVariable<dim> low_dim_unit_normal               ( ANY, 0.0 );
VectorVariable<dim> low_dim_neighbor_unit_normal      ( ANY, 0.0 );
VectorVariable<dim> low_dim_edge_unit_normal          ( ANY, 0.0 );
VectorVariable<dim> low_dim_neighbor_edge_unit_normal ( ANY, 0.0 );

// Find first refence element
typename std::vector<csmp::Element<dim>* >::const_iterator eit = this->PerimeterElementsBegin();
lowDimElement = NULL;
while ( eit!=this->ElementsEnd() )
{
if( (*eit)->IsLineElement() && elmt_dim.second == 2 ) // do not consider line element in surface region in 3D
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

while( elements_considered.size() != this->Elements() )
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

if( elements_considered.size() != this->Elements() )
{
typename std::set<csmp::Element<dim>* >::const_iterator eit = remaining_elements.begin();
lowDimElement = NULL;
while ( eit!=remaining_elements.end() )
{
if( (*eit)->IsLineElement() && elmt_dim.second == 2 ) // do not consider line element in surface region in 3D
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

#ifdef REGION_DEBUG

VectorVariable<dim> un( ANY, 0.0 );
typename std::vector<csmp::Element<dim>* >::const_iterator reit = this->ElementsBegin();
while ( reit != this->ElementsEnd() )
{
(*reit)->UnitNormal( un );
std::cerr<<" UN = (";
for( size_t i=0; i<dim; i++)
std::cerr<< " "<< un[i];
std::cerr<<" )\n";
reit++;
}
#endif

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
template<size_t dim>
double64  Region<dim>::Volume( bool multiply_with_porosity ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( this->Empty() ) {
    csmp_error.notice( WARNING, "Region<dim>::Volume:", "region is empty; returning NaN." );
    return std::numeric_limits<double64>::quiet_NaN();
  }
  double64  volume( 0. ), area( 0. ), length( 0. );

  // recording contributions of different types of elements making up the group
  if ( multiply_with_porosity ) {
    csmp::Index  phi_key = this->pref_.StorageKey( "porosity" );
    for ( typename vector<csmp::Element<dim>*>::const_iterator
          it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ ) {
      if ( (*it)->FE()->IsVolumeElement() )  volume += (*it)->Volume() * (*it)->Read( phi_key );
      else if ( (*it)->FE()->IsSurfaceElement() ) area += (*it)->Volume() * (*it)->Read( phi_key );
      else if ( (*it)->FE()->IsLineElement() )    length += (*it)->Volume() * (*it)->Read( phi_key );
    }
  }
  else
    for ( typename vector<csmp::Element<dim>*>::const_iterator
          it = this->elmt_vec_.begin(); it != this->elmt_vec_.end(); it++ ) {
      if ( (*it)->FE()->IsVolumeElement() )  volume += (*it)->Volume();
      else if ( (*it)->FE()->IsSurfaceElement() ) area += (*it)->Volume();
      else if ( (*it)->FE()->IsLineElement() )    length += (*it)->Volume();
    }

  // counting only the contributions of the elements with the highest spatial dimension
  if ( dim == 3U ) {
    if ( volume > 0. ) return volume;
    if ( area >   0. ) return area;
    if ( length > 0. ) return length;
  }
  else if ( dim == 2U )
    return (area > 0.) ? area : length;

  return length;

} // end Volume


/**
Computes the perimeter (2D) or the surface area (3D) of a model subdomain object.
To compute the surface area, only elements of dim-1 are considered.

@attention method assumes that the first nodes in higher-order elements are the corner nodes.

@todo SKM: method will produce nonsense if line elements intersect the boundary.

TODO: SKM: implement and use virtual void FiniteElement::AreaOfFace() rather than iffy statement in area calculation
*/
template<size_t dim>
double64  Region<dim>::SurfaceArea() const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  // in 1D there is no meaningful surface area
  if ( dim == 1U ) {
    csmp_error.notice( WARNING, "Region<1U>::SurfaceArea:",
                       "is not defined in one-dimensional model; returning NaN" );
    return std::numeric_limits<double64>::quiet_NaN();
  }

  vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit( this->bd_face_vec_.begin() );
  vector<size_t>  fnids;
  double64        area( 0. );

  if ( dim == 3U ) {
    // for all elements located on the region boundary
    for ( size_t i = this->InteriorElements(); i<this->Elements(); ++i, ++bit )
      // since each element can have multiple boundary faces
      for ( size_t j = 0U; j<(*bit).size(); j++ ) {
        const size_t face( (*bit)[j] );
        this->elmt_vec_[i]->FE()->NodesOfFace( face, fnids );
        const CSMP_FEM_TYPE etype( this->elmt_vec_[i]->FE()->ElementTypeOfFace( face ) );
        // triangular face
        if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE or
             etype == LINEAR_TRIANGLE3D or
             etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE )
          area += triangleArea( this->elmt_vec_[i]->N( fnids[0U] )->Coordinate(),
                                this->elmt_vec_[i]->N( fnids[1U] )->Coordinate(),
                                this->elmt_vec_[i]->N( fnids[2U] )->Coordinate() );
        // quadrilateral face
        if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL or etype == LINEAR_RECTANGLE or
             etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL )
          area += facetArea4( this->elmt_vec_[i]->N( fnids[0U] )->Coordinate(),
                              this->elmt_vec_[i]->N( fnids[1U] )->Coordinate(),
                              this->elmt_vec_[i]->N( fnids[2U] )->Coordinate(),
                              this->elmt_vec_[i]->N( fnids[3U] )->Coordinate() );
        // linear face
        if ( etype == ISOPARAMETRIC_LINEAR_BAR or
             etype == ISOPARAMETRIC_QUADRATIC_BAR ) {
          area += this->elmt_vec_[i]->N( fnids[0U] )->Coordinate().DistanceTo( this->elmt_vec_[i]->N( fnids[1U] )->Coordinate() );
          //csmp_error.notice( WARNING, "Region<dim>::SurfaceArea", "line-element thickness on perimeter is assumed to be one." );
        }
        // additional case of point face where a line-element is perpendicular to a boundary node
      }
  }
  // in 2D the face is a segment the length of which has to be used
  else if ( dim == 2U ) {
    // for all surface elements on the region boundary (excluding line elements)
    for ( size_t i = this->InteriorElements(); i<this->Elements(); i++, bit++ )
      if ( this->elmt_vec_[i]->FE()->IsSurfaceElement() )
        for ( size_t j = 0U; j<(*bit).size(); j++ ) {
          this->elmt_vec_[i]->FE()->NodesOfFace( (*bit)[j], fnids );
          area += this->elmt_vec_[i]->N( fnids[0U] )->Coordinate().DistanceTo( this->elmt_vec_[i]->N( fnids[1U] )->Coordinate() );
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
template<size_t dim>
double64  Region<dim>::VolumeIntegral( const char* prop, bool multiply_with_porosity, bool verbose ) const
{
  csmp::Index prop_key = this->pref_.StorageKey( prop );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( prop_key.type == TENSOR ) {
    csmp_error.notice( ERROR, "Region<dim>::VolumeIntegral:",
                       "No rule to integrate tensor properties. Nothing was done" );
    return std::numeric_limits<double64>::quiet_NaN();
  }
  if ( prop_key.type == VECTOR and prop_key.place != ELEMENT ) {
    csmp_error.notice( ERROR, "Region<dim>::VolumeIntegral:",
                       "Vector properties can only be integrated if they are placed on the element. Nothing was done" );
    return std::numeric_limits<double64>::quiet_NaN();
  }
  if ( this->elmt_vec_.empty() ) {
    csmp_error.notice( ERROR, "Region<dim>::VolumeIntegral", "Region is empty; returning NaN." );
    return std::numeric_limits<double64>::quiet_NaN();
  }

  if ( verbose ) {
#ifndef NDEBUG
    if ( dim == 2U && !(*this->elmt_vec_.begin())->FE()->IsSurfaceElement() ) {
      string info( this->Name() ); info += " ('"; info += prop; info += "')";
      csmp_error.notice( WARNING, "Region<2>::VolumeIntegral:", info, "region is not a surface; integral may not be correct." );
    }
    if ( dim == 3U && !(*this->elmt_vec_.begin())->FE()->IsVolumeElement() ) {
      string info( this->Name() ); info += " ('"; info += prop; info += "')";
      csmp_error.notice( WARNING, "Region<3>::VolumeIntegral:", info, "region is not a volume, integral may not be correct." );
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
    return std::numeric_limits<double64>::quiet_NaN();
  }

  double64  integral( 0. );

  // node or integration point properties
  if ( prop_key.place != ELEMENT )
  {
    csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
    // integration of nodal properties using the FEM framework
    if ( prop_key.place == NODE or prop_key.place == ELEMENT_INTEGRATION_POINT ) {
      if ( multiply_with_porosity ) {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key ) * (*iter)->Read( phi_key );
      }
      else {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key );
      }
      return integral;
    }
    else {
      throw csmp::Exception( ERROR, "Region<dim>::VolumeIntegral",
                             "only node, integration point or element properties can be integrated over the region" );
      return std::numeric_limits<double64>::quiet_NaN();
    }
  }

  else // element property
  {
    if ( prop_key.type == SCALAR ) {
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume() * (*iter)->Read( phi_key );
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume();
    }
    else { // VECTOR property
      VectorVariable<dim>  vc;
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume() * (*iter)->Read( phi_key );
        }
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume();
        }
    }
    return integral;
  }

  return std::numeric_limits<double64>::quiet_NaN();

} // end VolumeIntegral





/**
Integrates property over the elements contained in the Region and returns
the volume integral of the variable multiplied with the local element thickness.
Element thickness must be equal to one for volumetric elements.
*/
template<size_t dim>
double64  Region<dim>::VolumeIntegral_x_Thickness( const char* prop, bool multiply_with_porosity ) const
{
  csmp::Index prop_key = this->pref_.StorageKey( prop );
  csmp::Index thic_key = this->pref_.StorageKey( "thickness" );
  assert( thic_key.place == ELEMENT );
  assert( thic_key.type == SCALAR );

  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( prop_key.type == TENSOR ) {
    csmp_error.notice( ERROR, "Region<dim>::VolumeIntegral_x_Thickness",
                       "No rule to integrate tensor properties. Nothing was done" );
    return std::numeric_limits<double64>::quiet_NaN();
  }
  if ( prop_key.type == VECTOR and prop_key.place != ELEMENT ) {
    csmp_error.notice( ERROR, "Region<dim>::VolumeIntegral_x_Thickness",
                       "Vector properties can only be integrated if they are placed on the element. Nothing was done" );
    return std::numeric_limits<double64>::quiet_NaN();
  }
  if ( this->elmt_vec_.empty() ) {
    csmp_error.notice( ERROR, "Region<dim>::VolumeIntegral_x_Thickness", "Region is empty" );
    return std::numeric_limits<double64>::quiet_NaN();
  }

  // region properties (no thickness multiplier is accounted for)
  if ( prop_key.place == REGION ) {
    csmp_error.notice( WARNING, "Region<dim>::VolumeIntegral_x_Thickness", "there is no thickness attribute for region variable, assuming t=1" );
    if ( prop_key.type == SCALAR )
      return this->Read( prop_key ) * Volume( multiply_with_porosity );
    else if ( prop_key.type == VECTOR ) {
      VectorVariable<dim>  vc;
      this->Read( prop_key, vc );
      return vc.Length() * Volume( multiply_with_porosity );
    }
    return std::numeric_limits<double64>::quiet_NaN();
  }

  double64  integral( 0. );

  // node or integration point properties
  if ( prop_key.place != ELEMENT )
  {
    csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
    // integration of nodal properties using the FEM framework
    if ( prop_key.place == NODE or prop_key.place == ELEMENT_INTEGRATION_POINT ) {
      if ( multiply_with_porosity ) {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key ) * (*iter)->Read( phi_key ) * (*iter)->Read( thic_key );
      }
      else {
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->PropertyIntegral( prop_key ) * (*iter)->Read( thic_key );
      }
      return integral;
    }
    else {
      throw csmp::Exception( ERROR, "Region<dim>::VolumeIntegral_x_Thickness",
                             "only node, integration point or element properties can be integrated over the region" );
      return std::numeric_limits<double64>::quiet_NaN();
    }
  }

  else // element property
  {
    if ( prop_key.type == SCALAR ) {
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume() * (*iter)->Read( phi_key ) * (*iter)->Read( thic_key );
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ )
          integral += (*iter)->Read( prop_key ) * (*iter)->Volume() * (*iter)->Read( thic_key );
    }
    else { // VECTOR property
      VectorVariable<dim>  vc;
      if ( multiply_with_porosity ) {
        csmp::Index phi_key = this->pref_.StorageKey( "porosity" );
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume() * (*iter)->Read( phi_key ) * (*iter)->Read( thic_key );
        }
      }
      else
        for ( typename vector<csmp::Element<dim>*>::const_iterator
              iter = this->elmt_vec_.begin(); iter != this->elmt_vec_.end(); iter++ ) {
          (*iter)->Read( prop_key, vc );
          integral += vc.Length() * (*iter)->Volume() * (*iter)->Read( thic_key );
        }
    }
    return integral;
  }

  return std::numeric_limits<double64>::quiet_NaN();

} // end VolumeIntegral_x_Thickness




// SCREEN OUTPUT



/*
template<size_t dim>
void Region<dim>::Out() const
{
cout <<"\n\n\nRegion<"<< dim <<">::Out: ";
cout <<" member elements: interior="<< this->InteriorElements();
cout <<", boundary="<< this->elmt_vec_.size()-this->InteriorElements() <<": "<< endl;

for ( typename vector<csmp::Element<dim>*>::const_iterator
it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ ) {
if ( (*it) == NULL )
throw csmp::Exception( ERROR, "Region<dim>::Out",
"member element pointer not initialised");
//         else (*it)->Out();
}

cout <<"\n\n boundary elements and their boundary faces (current local numbering): "<< endl;
vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit(this->bd_face_vec_.begin());
for ( size_t i=this->InteriorElements(); i<this->elmt_vec_.size(); i++, bit++ ) {
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
