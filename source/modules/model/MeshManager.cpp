#include "MeshManager.h"
#include "PropertyDatabase.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "VSet.h"
#include "Exception.h"
#include "AP_BoolVector.h"
#include "CSMP_highLevelUtilities.h"
#include "ErrorHandler.h"
#include "PropertyData.h"
#include "Box.h"
#include "ModelTopology.h"

using namespace std;

namespace csmp {

// refactored
template<size_t dim>
MeshManager<dim>::MeshManager()
  : hybrid_element_mesh_( false ), n_nodes_( 0 ), n_elmts_( 0 ), n_faces_( 0 ), n_interfaces_( 0 )
{
}


// refactored
template<size_t dim>
MeshManager<dim>::MeshManager( const PropertyDatabase<dim>& pref, const FiniteElementManager& fem_manager, VSet<dim>& vset )
  : hybrid_element_mesh_( vset.HybridElementTypeMesh() )
{
  Initialize( pref, fem_manager, vset );
}


// Deallocate all the dynamic storage
template<size_t dim>
MeshManager<dim>::~MeshManager()
{
  // traversal of the existing mesh root nodes to find all nodes and elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );

  // traversal of the existing mesh root faces to find all faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );

  // traversal of the existing mesh root interfaces to find all interfaces	
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );

  // deallocate the memories and the root pointers
  for ( size_t i = 0U; i < nodes.size(); i++ ) {
    delete nodes[i];
    nodes[i] = NULL;
  }

  for ( size_t i = 0U; i < elmts.size(); i++ ) {
    delete elmts[i];
    elmts[i] = NULL;
  }

  for ( size_t i = 0U; i < faces.size(); i++ ) {
    delete faces[i];
    faces[i] = NULL;
  }

  for ( size_t i = 0U; i < interfaces.size(); i++ ) {
    delete interfaces[i];
    interfaces[i] = NULL;
  }

  for ( auto f : root_node_group_ )		 f = NULL;
  for ( auto f : root_elmt_group_ )		 f = NULL;
  for ( auto f : root_face_group_ )		 f = NULL;
  for ( auto f : root_interface_group_ ) f = NULL;

  root_node_group_.clear();
  root_elmt_group_.clear();
  root_face_group_.clear();
  root_interface_group_.clear();

  n_nodes_ = 0;
  n_elmts_ = 0;
  n_faces_ = 0;
  n_interfaces_ = 0;
}


/**

Assigns a copy of the supplied MeshManager to this MeshManager manager and
re-creates the connectivity for the new storage locations.

@param mmgr A MeshManager object.


@section implementation Implementation

Apart from copying the storage vectors for Nodes, IntegrationPoints and
Elements, the Elements must be initialized by updating their connections
to the new storage locations of their member Nodes, IntegrationPoints and
neighbor elements. This extra step is required since the
connections are pointers that are not re-addressed when the vectors are
copied.
*/
template<size_t dim>
MeshManager<dim>&  MeshManager<dim>::operator=( const MeshManager<dim>& mmgr )
{
  if ( this == &mmgr ) return *this;

  throw Exception( WARNING, "MeshManager<dim>::operator=",
                   "operator has not been tested yet" );

  hybrid_element_mesh_ = mmgr.hybrid_element_mesh_;
  n_nodes_ = mmgr.n_nodes_;
  n_elmts_ = mmgr.n_elmts_;
  n_faces_ = mmgr.n_faces_;
  n_interfaces_ = mmgr.n_interfaces_;
  root_node_group_ = mmgr.root_node_group_;
  root_elmt_group_ = mmgr.root_elmt_group_;
  root_face_group_ = mmgr.root_face_group_;
  root_interface_group_ = mmgr.root_interface_group_;

  return *this;

} // end operator=


/**
Configures the distributed variable storage, inialises the finite element container
with the element types that are contained in the current mesh and builds the
mesh including the connectivity among its elements.

@note the neighbors of each element include only the elements of the same type, i.e.
a line element only has line neighbors, a surface element surface element neighbors
and so forth.

@attention lower dimensional elements may have multiple neighbors for each of their
faces. Yet only one of them will be assigned.

@todo SKM: create manifolds to deal with lower-dimensional elements that have
multiple neighbors per face.
*/
template<size_t dim>
bool  MeshManager<dim>::Initialize( const PropertyDatabase<dim>& phys_vars,
                                    const FiniteElementManager& fem_manager,
                                    const VSet<dim>& vset )
{
  assert( n_nodes_ == 0 );
  assert( n_elmts_ == 0 );
  assert( n_faces_ == 0 );
  assert( n_interfaces_ == 0 );

  // temporary deques to handle primitives
  deque<Node<dim>*>			 node_connector;
  deque<Element<dim>*>	 elmt_connector;
  deque<Face<dim>*>			 face_connector;
  deque<InterFace<dim>*> interface_connector;

  hybrid_element_mesh_ = vset.HybridElementTypeMesh();

  // 1. checking the availability of the necessary finite element types, and the valid model topology
  // ------------------------------------------------------------------
  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: checking the availability of the necessary finite element types..." << endl;

  set<CSMP_FEM_TYPE> input_etypes;
  input_etypes.insert( parseFiniteElementTypeEnum( vset.ElementType( 0U ) ) );

  if ( vset.HybridElementTypeMesh() )
    for ( size_t i = 0U; i<vset.TotalNumberOfCells(); ++i )
      input_etypes.insert( parseFiniteElementTypeEnum( vset.ElementType( i ) ) );

  cout << "\nMeshManager<" << dim << ">::Initialize: ";
  cout << "input VSet contains the following finite element types:\n\t";
  for ( typename set<CSMP_FEM_TYPE>::const_iterator
        iit = input_etypes.begin(); iit != input_etypes.end(); iit++ ) {
    cerr << parseFiniteElementType( (*iit) ) << "  ";
    if ( !fem_manager.ContainsElementType( *iit ) ) {
      cerr << "\n\n\tFinite element type not available: " << parseFiniteElementType( *iit ) << endl;
      fem_manager.Out();
      throw Exception( FATAL_ERROR,
                       "MeshManager<dim>::Initialize(VSet):",
                       "'FiniteElementManager' lacks finite-element type required by VSet." );
    }
  }
  cout << endl;

  // 2. constructing the nodes and the elements using the VSet element type information
  // --------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: building storage and assigning nodes to elements..." << endl;

  size_t idx( 0U );

  /// Variable storage for nodes
   {
      vector<double64> coord( dim );
      const LocalVariables nvars( phys_vars.LocalVariablesAt( NODE ) );
      for ( size_t i = 0U; i < vset.Vertices(); ++i )
        {
          for ( size_t j = 0U; j<dim; ++j ) coord[j] = vset.P( j, idx );
          Node<dim>* node = new Node<dim>( idx, Point<dim>( coord ), nvars, NOT );
          node_connector.push_back( node );
          idx++;
        }
      n_nodes_ = node_connector.size();
    }
     
  /// Variable storage for elements and integration points
   {
      const LocalVariables evars( phys_vars.LocalVariablesAt( ELEMENT ) );
      const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( ELEMENT ) );
      typename deque<vector<size_t>>::const_iterator first( vset.PlistElmtsBegin() ), last( vset.PlistElmtsEnd() );

      idx = 0U;
      // 2.1 If the MeshManager contains only one element type
      if ( !vset.HybridElementTypeMesh() ) {
          const int32 csmpElementType = vset.ElementType( 0U );

          while ( first != last )
            {
              // TODO: this Element and all other primitives should be constructed directly in the MeshManager
              Element<dim>* elmt = new Element<dim>( idx, fem_manager.E( csmpElementType ), evars, cvars, NOT );
              const size_t nodes( fem_manager.E( csmpElementType )->Nodes() );
              for ( size_t j = 0U; j < nodes; j++ )
              {
                elmt->Assign( j, node_connector[vset.Plist( elmt->Idx(), j )] );
              }
              elmt_connector.push_back( elmt );
              idx++;
              first++;
            }
        }
      // 2.2 If there are multiple element types
      else {
          while ( first != last )
            {
              const int32 csmpElementType = vset.ElementType( idx );
              Element<dim>* elmt = new Element<dim>( idx, fem_manager.E( csmpElementType ), evars, cvars, NOT );
              const size_t nodes( fem_manager.E( csmpElementType )->Nodes() );
              for ( size_t j = 0U; j < nodes; j++ )
              {
                elmt->Assign( j, node_connector[vset.Plist( elmt->Idx(), j )] );
              }
              elmt_connector.push_back( elmt );
              idx++;
              first++;
            }
        }
      n_elmts_ = elmt_connector.size();
   }
  
  // ----------------------------------------------------------------
  // 2. Assigning neighbor elements to elements, faces and interfaces
  // ----------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: assigning neighbors to elements..." << endl;
  size_t eidx( 0U );
  if ( vset.WithNeighbourConnectivity() ) {
    if ( !vset.HybridElementTypeMesh() ) {
      const int32 csmpElementType = vset.ElementType( 0U );
      for ( auto& e : elmt_connector ) {
        const size_t neighbors( fem_manager.E( csmpElementType )->Neighbors() );
        for ( size_t j = 0U, nidx = 0U; j < neighbors; j++ ) {
          // if there is a neighbor (as is the case if the stored index is greater than zero)			
          if ( j < vset.PfvertsSize( e->Idx() ) ) {
            const int32 index( static_cast<int32>(vset.Pfvert( e->Idx(), j )) );
            if ( index >= 0 && index < n_elmts_ )
              e->Assign( nidx++, elmt_connector[static_cast<size_t>(index)] );
            else
              e->Assign( nidx++, static_cast<Element<dim>*>(NULL) );
          }
        }
      }
    }
    else {
      for ( auto& e : elmt_connector ) {
        const int32 csmpElementType = vset.ElementType( eidx++ );
        const size_t neighbors( fem_manager.E( csmpElementType )->Neighbors() );
        for ( size_t j = 0U, nidx = 0U; j < neighbors; j++ ) {
          // if there is a neighbor (as is the case if the stored index is greater than zero)			
          if ( j < vset.PfvertsSize( e->Idx() ) ) {
            const int32 index( static_cast<int32>(vset.Pfvert( e->Idx(), j )) );
            if ( index >= 0 && index < n_elmts_ )
              e->Assign( nidx++, elmt_connector[static_cast<size_t>(index)] );
            else
              e->Assign( nidx++, static_cast<Element<dim>*>(NULL) );
          }
        }
      }
    }
  }
  else
    csmp_error.notice( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity; nothing was done." );

  // 3. constructing the faces using the VSet nodes information
  // ----------------------------------------------------------
  // (continuous running index 'idx' will be used so that face-IDs start with n-elements)
  // different face types are intrinsic to the VSet if so initialized
  if ( vset.Faces() > 0 )
  {
    if ( csmp_error.Verbose() )
      cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to faces..." << endl;
    idx = n_elmts_;
    // since faces are lower-dimensional, the mesh must contain different element types
    assert( vset.HybridElementTypeMesh() );
    const LocalVariables evars( phys_vars.LocalVariablesAt( FACE ) );
    const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( FACE ) );

    typename deque<vector<size_t> >::const_iterator  first( vset.PlistFacesBegin() ), last( vset.PlistFacesEnd() );
    while ( first != last ) {
      const int32 csmpElementType = vset.ElementType( idx );
      if ( csmpElementType == UNKNOWN ) { ++first; continue;}

      Face<dim>* face = new Face<dim>( idx, fem_manager.E( csmpElementType ), evars, cvars );
      const size_t nodes( face->Nodes() );
      for ( size_t j = 0U; j<nodes; ++j ) {
        // assigning node indices
        if( vset.Plist( face->Idx(), j ) < node_connector.size() )
          face->Assign( j, node_connector[vset.Plist( face->Idx(), j )] );
      }
      face_connector.push_back( face );
      ++first;
      ++idx;
    }
    n_faces_ = face_connector.size();

    if ( csmp_error.Verbose() )
      cout << "\nMeshManager<" << dim << ">::Initialize: connecting faces to their higher-dimensional neighbors..." << endl;
    const size_t elements( elmt_connector.size() );
    const size_t faces( face_connector.size() );
    for ( auto& e : face_connector ) {
      // equidimensional neighbors first
      // -------------------------------		
      const size_t neighbors( e->Neighbors() );
      for ( size_t j = 0U; j<neighbors; ++j ) {
        // if there is a neighbor (as is the case if the stored index is greater than zero)
        // (e->Idx() starts with elements=first face)
        if ( j < vset.PfvertsSize( e->Idx() ) ) {
          const long64 index( vset.Pfvert( e->Idx(), j ) );
          if ( index >= 0 ) {
            assert( index >= elements );
            assert( index < elements + faces ); // (-) elements because face container is numbered from 0..n-1                   
            if ( index >= elements && index < (elements + faces) )
              e->Assign( j, face_connector[static_cast<size_t>(index) - elements] );
            else
              e->Assign( j, static_cast<Face<dim>*>(NULL) );
          }
        }
      }
      // higher-dimensional neighbors
      // ----------------------------
      // (are stored in VSet 'pfverts' record after the equidimensional neighbors)
      // index of inner neighbor element i which is always there

      const long64 index1( vset.Pfvert( e->Idx(), neighbors ) );
      // index of outer neighbor element which may be there		
      const long64 index2( vset.Pfvert( e->Idx(), neighbors + 1U ) );

      Element<dim>* const innerElement = (index1 < 0) ? NULL : elmt_connector[index1];
      Element<dim>* const outerElement = (index2 < 0) ? NULL : elmt_connector[index2];
      // assigning inner and outer higher-dimensional neighbors
      //              inner element             outer element
      e->Assign( innerElement, outerElement );
    }

    // traversal of the existing mesh to find all its faces
    deque<set<Face<dim>*>>	explored_face_groups;
    set<Face<dim>*>			discovered_faces;
    deque<Face<dim>*>		current_faces;
    map<size_t, Face<dim>*> faces_map;
    for ( auto f : face_connector )
      faces_map[f->Idx()] = f;

    // starting at the first face		
    current_faces.push_back( faces_map.begin()->second );

    size_t group_idx = 0U;
    while ( !faces_map.empty() )
    {
      while ( !current_faces.empty() ) {
        Face<dim>*  n_ptr( *current_faces.begin() );
        // for all neighbor faces of the current face
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_face = discovered_faces.insert( n_ptr->Neighbor( i ) );
          if ( new_face.second ) {
            current_faces.push_back( n_ptr->Neighbor( i ) );
            faces_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the face from the discovered (but not yet explored) deque
        current_faces.pop_front();
      }
      if ( discovered_faces.size() == 0 ) {
        auto first_face = faces_map.begin()->second;
        discovered_faces.insert( first_face );
        explored_face_groups.push_back( discovered_faces );
        faces_map.erase( first_face->Idx() );
        discovered_faces.clear();
        group_idx++;
        continue;
      }
      explored_face_groups.push_back( discovered_faces );

      // if there is only one single face left, it creates a new group of faces.
      if ( faces_map.size() == 1 ) {
        discovered_faces.clear();
        auto first_face = faces_map.begin()->second;
        discovered_faces.insert( first_face );
        explored_face_groups.push_back( discovered_faces );
        break;
      }
      if ( faces_map.size() > 0 ) {
        auto first_face = faces_map.begin()->second;
        current_faces.push_back( first_face );
      }

      discovered_faces.clear();
      group_idx++;
    }

    for ( auto group : explored_face_groups )
      root_face_group_.push_back( (*group.begin()) );
  } // end faces

    // 4. constructing the interfaces using the VSet nodes information
    // ---------------------------------------------------------------
    // (continuous running index 'idx' will be used)
  if ( vset.InterFaces() > 0 )
  {
    if ( csmp_error.Verbose() )
      cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to interfaces..." << endl;
    assert( vset.HybridElementTypeMesh() );
    const LocalVariables evars( phys_vars.LocalVariablesAt( INTER_FACE ) );
    const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( INTER_FACE ) );

    typename deque<vector<size_t> >::const_iterator  first( vset.PlistInterFacesBegin() ), last( vset.PlistInterFacesEnd() );
    while ( first != last ) {
      const int32 csmpElementType = vset.ElementType( idx );
      InterFace<dim>* inter_face = new InterFace<dim>( idx, fem_manager.E( csmpElementType ), evars, cvars );
      const size_t nodes( inter_face->Nodes() );
      for ( size_t j = 0U; j<nodes; ++j ) {
        // assigning node indices
        inter_face->Assign( j, node_connector[ vset.Plist( inter_face->Idx(), j )], INSIDE );
      }
      interface_connector.push_back( inter_face );
      ++idx;
      ++first;
    }
    n_interfaces_ = interface_connector.size();

    if ( csmp_error.Verbose() )
      cout << "\nMeshManager<" << dim << ">::Initialize: connecting interfaces to their higher-dimensional neighbors..." << endl;
    const size_t elements( elmt_connector.size() );
    const size_t faces( face_connector.size() );
    // connecting interfaces to their higher-dimensional neighbors
    for ( auto& e : interface_connector ) {
      // equidimensional neighbors first
      const size_t neighbors( e->Neighbors() );
      for ( size_t j = 0U; j<neighbors; ++j ) {
        // if there is a neighbor (as is the case if the stored index is greater than zero)
        const long64 index( vset.Pfvert( e->Idx(), j ) );
        if ( index >= 0 ) {
          assert( index >= elements + faces );
          assert( index < elements + faces + interfaces ); // (-) because interface container is numbered from 0..n-1
          e->Assign( j, interface_connector[static_cast<size_t>(index) - elements - faces] );
        }
      }
      // higher-dimensional neighbors
      // ----------------------------
      // (the 2 sides will always be present because interfaces exist only on internal boundaries)
      if ( vset.Pfvert( e->Idx(), neighbors ) < 0 || vset.Pfvert( e->Idx(), neighbors + 1U) < 0 || vset.Pfvert( e->Idx(), neighbors ) >= elmt_connector.size() || vset.Pfvert( e->Idx(), neighbors + 1U ) >= elmt_connector.size() )
        continue;
      e->Assign( elmt_connector[vset.Pfvert( e->Idx(), neighbors )], elmt_connector[vset.Pfvert( e->Idx(), neighbors + 1U )] );
    }

    // traversal of the existing mesh nodes to find all its interfaces
    deque<set<InterFace<dim>*>>	explored_interface_groups;
    set<InterFace<dim>*>		discovered_interfaces;
    deque<InterFace<dim>*>		current_interfaces;
    map<size_t, InterFace<dim>*> interfaces_map;
    for ( auto f : interface_connector )
      interfaces_map[f->Idx()] = f;

    // starting at the first interface		
    current_interfaces.push_back( interfaces_map.begin()->second );

    size_t group_idx = 0U;
    while ( !interfaces_map.empty() )
    {
      while ( !current_interfaces.empty() ) {
        InterFace<dim>*  n_ptr( *current_interfaces.begin() );
        // for all neighbor interfaces of the current interface
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == nullptr ) continue;

          // if this neighbor is a new one					
          auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
          if ( new_interface.second ) {
            current_interfaces.push_back( n_ptr->Neighbor( i ) );
            interfaces_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the interface from the discovered (but not yet explored) deque
        current_interfaces.pop_front();
      }
      if ( discovered_interfaces.size() == 0 ) {
        auto first_interface = interfaces_map.begin()->second;
        discovered_interfaces.insert( first_interface );
        explored_interface_groups.push_back( discovered_interfaces );
        interfaces_map.erase( first_interface->Idx() );
        discovered_interfaces.clear();
        group_idx++;
        continue;
      }
      explored_interface_groups.push_back( discovered_interfaces );

      // if there is only one single interface left, it creates a new group of interfaces.
      if ( interfaces_map.size() == 1 ) {
        discovered_interfaces.clear();
        auto first_interface = interfaces_map.begin()->second;
        discovered_interfaces.insert( first_interface );
        explored_interface_groups.push_back( discovered_interfaces );
        break;
      }
      if ( interfaces_map.size() > 0 ) {
        auto first_interface = interfaces_map.begin()->second;
        current_interfaces.push_back( first_interface );
      }

      discovered_interfaces.clear();
      group_idx++;
    }

    for ( auto group : explored_interface_groups )
      root_interface_group_.push_back( (*group.begin()) );
  } // end interfaces

    // ---------------------------------------------------------------------
    // 3. Flagging nodes located at the model boundary
    // ---------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: flagging boundary objects..." << endl;
  if ( vset.BFlags() > 0 )
  {
    // nodes were initially constructed as not located at the model boundary
    for ( auto bit = vset.BFlagsBegin(); bit != vset.BFlagsEnd(); bit++ )
      node_connector[(*bit).first]->AtBoundary( intToBOX_BOUNDARY( (*bit).second ) );
  }

  // ---------------------------------------------------------------------
  // 4. Flagging the elements using boundary flags from the nodes
  // ---------------------------------------------------------------------
  flagElementUsingNodal_BOX_BOUNDARY_Flags<dim>( elmt_connector.begin(), elmt_connector.end() );

  // ------------------------------------------------------------------------------
  // 5. Assigning parent elements (these are the elements that share the node) and
  // their respective internal node-id numbers to the nodes
  // -------------------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: assigning parent element information to nodes..." << endl;
  vector<size_t>  parent_elmts_per_node( node_connector.size(), 0U );

  // counting how many parent elements each node has
  for ( auto& e : elmt_connector ) {
    for ( auto nit = e->NodesBegin(); nit != e->NodesEnd(); nit++ )
      parent_elmts_per_node[(*nit)->Idx()]++;
  }

  // reserving the memory for the parent storage and zeroing parent vector for next step
  size_t i( 0 );
  for ( auto& n : node_connector ) {
    n->ResizeParentStorage( parent_elmts_per_node[i++] );
  }

  // assigning the parent element information to the nodes
  for ( auto& e : elmt_connector ) {
    for ( size_t j = 0U; j<e->Nodes(); j++ )
      e->N( j )->Assign( j, e );
  }

  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: forming regions for contiguous subdomains..." << endl;

  // forming root pointers of the nodes for contiguous regions
  if ( node_connector.size()>0 ) {
    // assgining the root pointers
    // traversal of the existing mesh nodes to find all its nodes
    deque<set<Node<dim>*>>	explored_node_groups;
    set<Node<dim>*>			discovered_nodes;
    deque<Node<dim>*>		current_nodes;
    map<size_t, Node<dim>*> nodes_map;
    for ( auto n : node_connector )
      nodes_map[n->Idx()] = n;

    // starting at the first node
    current_nodes.push_back( nodes_map.begin()->second );

    size_t group_idx = 0U;
    while ( !nodes_map.empty() )
    {
      while ( !current_nodes.empty() ) {
        auto n_ptr( *current_nodes.begin() );
        // for all neighbor nodes of the current node
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_node = discovered_nodes.insert( n_ptr->Neighbor( i ) );
          if ( new_node.second ) {
            current_nodes.push_back( n_ptr->Neighbor( i ) );
            nodes_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the node from the discovered (but not yet explored) deque
        current_nodes.pop_front();
      }
      if ( discovered_nodes.size() == 0 ) {
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        nodes_map.erase( first_node->Idx() );
        discovered_nodes.clear();
        group_idx++;
        continue;
      }
      // avoid the infinite iterative looping search
      bool infinite_looping( false );
      for ( size_t i = 0U; i < explored_node_groups.size(); i++)
      {
        if ( explored_node_groups[i] == discovered_nodes ) {
          infinite_looping = true;
          if ( csmp_error.Verbose() )
            cout << "\nMeshManager<" << dim << ">::Initialize: node connections might be wrong..." << endl;
        }

      }

      if ( infinite_looping ) {
        discovered_nodes.clear();
        if ( nodes_map.size() > 0 ) {
          auto first_node = nodes_map.begin()->second;
          discovered_nodes.insert( first_node );
          nodes_map.erase( first_node->Idx() );
        }        
      }      
      
      explored_node_groups.push_back( discovered_nodes );      

      // if there is only one single node left, it creates a new group of nodes.
      if ( nodes_map.size() == 1 ) {
        discovered_nodes.clear();
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        if ( csmp_error.Verbose() )
          cout << "\nMeshManager<" << dim << ">::Initialize: the node (" << first_node->Idx() << ") doesn't have any parents and neighbours..." << endl;
        break;
      }
      if ( nodes_map.size() > 0 ) {
        auto first_node = nodes_map.begin()->second;
        current_nodes.push_back( first_node );
      }

      discovered_nodes.clear();
      group_idx++;
    }

    for ( auto group : explored_node_groups )
      root_node_group_.push_back( (*group.begin()) );

  } // end nodes

  set<Element<dim>*>	explored_elmt_groups;
  if ( elmt_connector.size()>0 ) {
    // assgining the root pointers
    for ( auto root_node : root_node_group_ )
    {
      for ( size_t i = 0; i < root_node->Parents(); i++ ) {
        auto parent = root_node->Parent( i );
        if ( parent ) {
          explored_elmt_groups.insert( parent );
          break;
        }
      }
    }
  } // end elements

  for ( auto root_elmt : explored_elmt_groups )
    root_elmt_group_.push_back( root_elmt );

  if ( root_node_group_.size() != root_elmt_group_.size() ) {
    root_node_group_.clear();
    for ( auto root_elmt : root_elmt_group_ )
      root_node_group_.push_back( root_elmt->N( 0 ) );
  }

  // 5. verify the constructed elements according to the VSet and assgin the root node and element
  if ( node_connector.size() != vset.Vertices() || elmt_connector.size() != vset.Elements() ) {
    cout << "\nMeshManager<" << dim << ">::Initialize: ";
    cout << " The constructed elements have a different number of nodes than are stored in the VSet" << vset.Elements() << endl;
    return false;
  }
  return true;
  
} // end Initialise






/**
    Update  root pointers of the mesh after modifying the mesh.
*/
template<size_t dim>
void MeshManager<dim>::Update( std::deque<Node<dim>*> nodes, std::deque<Element<dim>*> elmts )
{
  // forming root pointers of the nodes for contiguous subdomains
  if ( nodes.size()>0 ) {
    // finding the root pointers		
    deque<set<Node<dim>*>>	explored_node_groups;
    set<Node<dim>*>			discovered_nodes;
    deque<Node<dim>*>		current_nodes;
    map<size_t, Node<dim>*> nodes_map;
    for ( auto n : nodes )
      nodes_map[n->Idx()] = n;

    // starting at the first node
    current_nodes.push_back( nodes_map.begin()->second );

    size_t group_idx = 0U;
    while ( !nodes_map.empty() )
    {
      while ( !current_nodes.empty() ) {
        Node<dim>*  n_ptr( *current_nodes.begin() );
        // for all neighbor nodes of the current node
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_node = discovered_nodes.insert( n_ptr->Neighbor( i ) );
          if ( new_node.second ) {
            current_nodes.push_back( n_ptr->Neighbor( i ) );
            nodes_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the node from the discovered (but not yet explored) deque
        current_nodes.pop_front();
      }
      if ( discovered_nodes.size() == 0 ) {
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        nodes_map.erase( first_node->Idx() );
        discovered_nodes.clear();
        group_idx++;
        continue;
      }
      explored_node_groups.push_back( discovered_nodes );

      // if there is only one single node left, it creates a new group of nodes.
      if ( nodes_map.size() == 1 ) {
        discovered_nodes.clear();
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        break;
      }
      if ( nodes_map.size() > 0 ) {
        auto first_node = nodes_map.begin()->second;
        current_nodes.push_back( first_node );
      }

      discovered_nodes.clear();
      group_idx++;
    }

    //update the node root pointers
    std::deque<Node<dim>*> updated_root_node_group;
    for ( auto group : explored_node_groups )
    {
      deque<Node<dim>*> nodes_vec;
      nodes_vec.assign( group.begin(), group.end() );
      sort( nodes_vec.begin(), nodes_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_node_group.push_back( nodes_vec.front() );
      nodes_vec.clear();
    }

    root_node_group_.swap( updated_root_node_group );
  } // end nodes	

  if ( elmts.size()>0 ) {
    // assgining the root pointers
    set<Element<dim>*>	explored_elmt_groups;
    for ( auto root_node : root_node_group_ )
    {
      for ( size_t i = 0; i < root_node->Parents(); i++ ) {
        auto parent = root_node->Parent( i );
        if ( parent ) {
          explored_elmt_groups.insert( parent );
          break;
        }
      }
    }

    //update the element root pointers
    std::deque<Element<dim>*> updated_root_elmt_group;
    for ( auto root_elmt : explored_elmt_groups )
      updated_root_elmt_group.push_back( root_elmt );

    root_elmt_group_.swap( updated_root_elmt_group );

    if ( root_node_group_.size() != root_elmt_group_.size() ) {
      root_node_group_.clear();
      for ( auto root_elmt : root_elmt_group_ )
        root_node_group_.push_back( root_elmt->N( 0 ) );
    }
  } // end elements
}





template<size_t dim>
void MeshManager<dim>::Update()
{
  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );
  sort( interfaces.begin(), interfaces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // forming root pointers of the nodes for contiguous subdomains
  if ( nodes.size()>0 ) {
    // finding the root pointers		
    deque<set<Node<dim>*>>	explored_node_groups;
    set<Node<dim>*>			discovered_nodes;
    deque<Node<dim>*>		current_nodes;
    map<size_t, Node<dim>*> nodes_map;
    for ( auto n : nodes )
      nodes_map[n->Idx()] = n;

    // starting at the first node
    current_nodes.push_back( nodes_map.begin()->second );

    size_t group_idx = 0U;
    while ( !nodes_map.empty() )
    {
      while ( !current_nodes.empty() ) {
        Node<dim>*  n_ptr( *current_nodes.begin() );
        // for all neighbor nodes of the current node
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_node = discovered_nodes.insert( n_ptr->Neighbor( i ) );
          if ( new_node.second ) {
            current_nodes.push_back( n_ptr->Neighbor( i ) );
            nodes_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the node from the discovered (but not yet explored) deque
        current_nodes.pop_front();
      }
      if ( discovered_nodes.size() == 0 ) {
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        nodes_map.erase( first_node->Idx() );
        discovered_nodes.clear();
        group_idx++;
        continue;
      }
      explored_node_groups.push_back( discovered_nodes );

      // if there is only one single node left, it creates a new group of nodes.
      if ( nodes_map.size() == 1 ) {
        discovered_nodes.clear();
        auto first_node = nodes_map.begin()->second;
        discovered_nodes.insert( first_node );
        explored_node_groups.push_back( discovered_nodes );
        break;
      }
      if ( nodes_map.size() > 0 ) {
        auto first_node = nodes_map.begin()->second;
        current_nodes.push_back( first_node );
      }

      discovered_nodes.clear();
      group_idx++;
    }

    //update the node root pointers
    std::deque<Node<dim>*> updated_root_node_group;
    for ( auto group : explored_node_groups )
    {
      deque<Node<dim>*> nodes_vec;
      nodes_vec.assign( group.begin(), group.end() );
      sort( nodes_vec.begin(), nodes_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_node_group.push_back( nodes_vec.front() );
      nodes_vec.clear();
    }

    root_node_group_.swap( updated_root_node_group );
  } // end nodes	

  if ( elmts.size()>0 ) {
    // assgining the root pointers
    set<Element<dim>*>	explored_elmt_groups;
    for ( auto root_node : root_node_group_ )
    {
      for ( size_t i = 0; i < root_node->Parents(); i++ ) {
        auto parent = root_node->Parent( i );
        if ( parent ) {
          explored_elmt_groups.insert( parent );
          break;
        }
      }
    }

    //update the element root pointers
    std::deque<Element<dim>*> updated_root_elmt_group;
    for ( auto root_elmt : explored_elmt_groups )
      updated_root_elmt_group.push_back( root_elmt );

    root_elmt_group_.swap( updated_root_elmt_group );

    if ( root_node_group_.size() != root_elmt_group_.size() ) {
      root_node_group_.clear();
      for ( auto root_elmt : root_elmt_group_ )
        root_node_group_.push_back( root_elmt->N( 0 ) );
    }
  } // end elements

    // forming root pointers of the faces for contiguous regions
  if ( faces.size() > 0 ) {
    // finding the root pointers
    // traversal of the existing mesh nodes to find all its faces
    deque<set<Face<dim>*>>		explored_face_groups;
    set<Face<dim>*>				discovered_faces;
    deque<Face<dim>*>			current_faces;
    map<size_t, Face<dim>*>		faces_map;
    for ( auto f : faces )
      faces_map[f->Idx()] = f;

    // starting at the first face		
    current_faces.push_back( faces.front() );

    size_t group_idx = 0U;
    while ( !faces_map.empty() )
    {
      while ( !current_faces.empty() ) {
        Face<dim>*  n_ptr( *current_faces.begin() );
        // for all neighbor faces of the current face
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_face = discovered_faces.insert( n_ptr->Neighbor( i ) );
          if ( new_face.second ) {
            current_faces.push_back( n_ptr->Neighbor( i ) );
            faces_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the face from the discovered (but not yet explored) deque
        current_faces.pop_front();
      }
      if ( discovered_faces.size() == 0 ) {
        auto first_face = faces_map.begin()->second;
        discovered_faces.insert( first_face );
        explored_face_groups.push_back( discovered_faces );
        faces_map.erase( first_face->Idx() );
        discovered_faces.clear();
        group_idx++;
        continue;
      }
      explored_face_groups.push_back( discovered_faces );

      // if there is only one single face left, it creates a new group of faces.
      if ( faces_map.size() == 1 ) {
        discovered_faces.clear();
        auto first_face = faces_map.begin()->second;
        discovered_faces.insert( first_face );
        explored_face_groups.push_back( discovered_faces );
        break;
      }
      if ( faces_map.size() > 0 ) {
        auto first_face = faces_map.begin()->second;
        current_faces.push_back( first_face );
      }

      discovered_faces.clear();
      group_idx++;
    }

    // update the face root pointers
    std::deque<Face<dim>*> updated_root_face_group;
    for ( auto group : explored_face_groups )
    {
      deque<Face<dim>*> faces_vec;
      faces_vec.assign( group.begin(), group.end() );
      sort( faces_vec.begin(), faces_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_face_group.push_back( faces_vec.front() );
      faces_vec.clear();
    }
    root_face_group_.swap( updated_root_face_group );
  } // end faces

    // forming root pointers of the interfaces for contiguous regions
  if ( interfaces.size() > 0 ) {
    // finding the root pointers
    // traversal of the existing mesh nodes to find all its interfaces
    deque<set<InterFace<dim>*>>		explored_interface_groups;
    set<InterFace<dim>*>			discovered_interfaces;
    deque<InterFace<dim>*>			current_interfaces;
    map<size_t, InterFace<dim>*>	interfaces_map;
    for ( auto f : interfaces )
      interfaces_map[f->Idx()] = f;

    // starting at the first interface		
    current_interfaces.push_back( interfaces.front() );

    size_t group_idx = 0U;
    while ( !interfaces_map.empty() )
    {
      while ( !current_interfaces.empty() ) {
        InterFace<dim>*  n_ptr( *current_interfaces.begin() );
        // for all neighbor interfaces of the current interface
        for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
          if ( n_ptr->Neighbor( i ) == NULL ) continue;

          // if this neighbor is new one					
          auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
          if ( new_interface.second ) {
            current_interfaces.push_back( n_ptr->Neighbor( i ) );
            interfaces_map.erase( n_ptr->Neighbor( i )->Idx() );
          }
        }
        // removing the interface from the discovered (but not yet explored) deque
        current_interfaces.pop_front();
      }
      if ( discovered_interfaces.size() == 0 ) {
        auto first_interface = interfaces_map.begin()->second;
        discovered_interfaces.insert( first_interface );
        explored_interface_groups.push_back( discovered_interfaces );
        interfaces_map.erase( first_interface->Idx() );
        discovered_interfaces.clear();
        group_idx++;
        continue;
      }
      explored_interface_groups.push_back( discovered_interfaces );

      // if there is only one single interface left, it creates a new group of interfaces.
      if ( interfaces_map.size() == 1 ) {
        discovered_interfaces.clear();
        auto first_interface = interfaces_map.begin()->second;
        discovered_interfaces.insert( first_interface );
        explored_interface_groups.push_back( discovered_interfaces );
        break;
      }
      if ( interfaces_map.size() > 0 ) {
        auto first_interface = interfaces_map.begin()->second;
        current_interfaces.push_back( first_interface );
      }

      discovered_interfaces.clear();
      group_idx++;
    }

    // update the face root pointers
    std::deque<InterFace<dim>*> updated_root_face_group;
    for ( auto group : explored_interface_groups )
    {
      deque<InterFace<dim>*> interfaces_vec;
      interfaces_vec.assign( group.begin(), group.end() );
      sort( interfaces_vec.begin(), interfaces_vec.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
      updated_root_face_group.push_back( interfaces_vec.front() );
      interfaces_vec.clear();
    }
    root_interface_group_.swap( updated_root_face_group );
  } // end interfaces
}





/**
Rebuilds the storage of parents of each element from the model and assigns their relationships.
*/
template<size_t dim>
void MeshManager<dim>::RebuildParentRelationships( typename vector<Node<dim>*>::iterator begin, typename vector<Node<dim>*>::iterator end )
{
  deque< pair<size_t, Element<dim>*>> parents;
  for ( auto it = begin; it != end; ++it ) {
    auto n = *it;
    for ( size_t i = 0; i < n->Parents(); ++i ) {
      auto e = n->Parent( i );
      parents.push_back( make_pair( n->ParentNodeNumber( i ), e ) );
    }
    n->EraseParents();
    n->ResizeParentStorage( parents.size() );
    for ( auto& p : parents ) {
      n->Assign( p.first, p.second );
    }
    parents.clear();
  }
}





/**
Counts elements the nodes of which are all located on the model boundary.

@attention such elements typically give rise to problems with the assignment of
boundary conditions and should be eliminated.

@attention method will work only if the elements are stored in elmt_connector deque.

@author SKM
*/
template<size_t dim>
size_t MeshManager<dim>::DetectElementsWithAllNodesOnBoundary( set<size_t>& belmts ) const
{
  belmts.clear();

  assert( Elements() > 0 );
  if ( Elements() == 0 ) return 0U;

  // traversal of the existing mesh nodes to find all its elements	
  deque<const Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );

  size_t boundary_only_elements( 0U );

  // counting how many parent elements each node has
  for ( auto& e : elmts ) {
    size_t counter( 0U );
    for ( auto nit = e->NodesBegin(); nit != e->NodesEnd(); nit++ )
      if ( (*nit)->AtBoundary() != NOT )
        counter++;
    if ( counter == e->Nodes() ) {
      belmts.insert( e->Idx() );
      boundary_only_elements++;
    }
  }

  return boundary_only_elements;
}


/**
Assigns fem_manager to fvs_manager. Assigns the finite volume stencil pointers of the elements to the
corresponding finite volume stencils.

@attention If a stencil is assigned already, noting is done. Remove stencil first (NULL ptr in elements)

@todo (2-C) This should not also assign the femManager to the stencilManager - method does too much
*/
template<size_t dim>
void MeshManager<dim>::InitializeFiniteVolumeStencils( const PropertyDatabase<dim>& pref,
                                                       const FiniteElementManager& fem_manager,
                                                       FiniteVolumeStencilManager<dim>& fvs_manager )
{
  // 1. find all nodes and elements from the mesh
  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );

  if ( Elements() <= 1U )
    throw Exception( FATAL_ERROR, "MeshManager<dim>::InitializeFiniteVolumeStencils",
                     "Currently no model exists to which stencils could be assigned." );

  // 2. initialize the stencil manager (the stencils are build and assigned the correct properties
  fvs_manager.Initialize( fem_manager );

  // 3. Now the stencil pointers in each finite element are connected to the correct corresponding stencils and update variable
  // storage for fv integration (sector/facet) point properties
  const LocalVariables lvs( pref.LocalVariablesAt( ELEMENT ) );
  const IntegrationPointVariables ipvs( pref.IntegrationPointVariablesAt( ELEMENT ) );

  for ( auto e : elmts ) {
    if ( e->FV() == nullptr )
    {
      e->AssignFiniteVolume( fvs_manager.Stencil( e->FE_Type() ) );
      e->ResizePropertyStorage( lvs, ipvs );
    }
  }

} // end InitializeFiniteVolumeStencils


  /**
  Inserts the corresponding node
  */
template<size_t dim>
Node<dim>* MeshManager<dim>::AddIfUnique( Node<dim>& node )
{
  set<Element<dim>*>		discovered_elements;
  set<const Node<dim>*>	discovered_nodes;
  deque<const Node<dim>*>	current_nodes;
  for ( auto root_node : root_node_group_ ) {
    discovered_nodes.insert( root_node );
    current_nodes.push_back( root_node );
    while ( !current_nodes.empty() ) {
      const Node<dim>*  n_ptr( *current_nodes.begin() );
      for ( size_t i = 0U; i < n_ptr->Parents(); i++ ) {
        for ( size_t j = 0U; j < n_ptr->Parent( i )->Nodes(); j++ ) {
          if ( j != n_ptr->ParentNodeNumber( i ) ) {
            if ( *(n_ptr->Parent( i )->N( j )) == node )
              return n_ptr->Parent( i )->N( j );
            pair<typename set<const Node<dim>*>::const_iterator, bool>
              new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
            if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
          }
        }
        discovered_elements.insert( n_ptr->Parent( i ) );
      }
      current_nodes.pop_front();
    }
  }

  Node<dim>* new_node = new Node<dim>( node );
  *new_node = node;
  n_nodes_++;

  return new_node;
}


/**
Inserts the corresponding element
*/
template<size_t dim>
Element<dim>* MeshManager<dim>::AddIfUnique( Element<dim>& elmt )
{
  set<Element<dim>*>		discovered_elements;
  set<const Node<dim>*>	discovered_nodes;
  deque<const Node<dim>*>	current_nodes;
  for ( auto root_node : root_node_group_ ) {
    // starting at the root element
    discovered_nodes.insert( root_node );
    current_nodes.push_back( root_node );
    while ( !current_nodes.empty() ) {
      const Node<dim>*  n_ptr( *current_nodes.begin() );
      for ( size_t i = 0U; i < n_ptr->Parents(); i++ ) {
        for ( size_t j = 0U; j < n_ptr->Parent( i )->Nodes(); j++ ) {
          if ( j != n_ptr->ParentNodeNumber( i ) ) {
            pair<typename set<const Node<dim>*>::const_iterator, bool>
              new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
            if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
          }
        }
        if ( *(n_ptr->Parent( i )) == elmt ) return n_ptr->Parent( i );
        discovered_elements.insert( n_ptr->Parent( i ) );
      }
      current_nodes.pop_front();
    }
  }

  Element<dim>* new_elmt = new Element<dim>( elmt );
  *new_elmt = elmt;
  n_elmts_++;

  return new_elmt;
}


/**
Inserts the corresponding face
*/
template<size_t dim>
Face<dim>* MeshManager<dim>::AddIfUnique( Face<dim>& face )
{
  set<Face<dim>*>		discovered_faces;
  deque<Face<dim>*>	current_faces;

  Face<dim>* found_face = NULL;

  for ( auto root_face : root_face_group_ ) {
    discovered_faces.insert( root_face );
    current_faces.push_back( root_face );
    while ( !current_faces.empty() || found_face != NULL ) {
      Face<dim>*  n_ptr( *current_faces.begin() );
      for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
        if ( n_ptr->Neighbor( i ) == NULL ) continue;
        if ( *(n_ptr->Neighbor( i )) == face )
          return n_ptr->Neighbor( i );
        auto new_face = discovered_faces.insert( n_ptr->Neighbor( i ) );
        if ( new_face.second ) current_faces.push_back( n_ptr->Neighbor( i ) );
      }
      current_faces.pop_front();
    }
  }

  Face<dim>* new_face = new Face<dim>( face );
  *new_face = face;
  n_faces_++;
  return new_face;
}


/**
Inserts the corresponding interface
*/
template<size_t dim>
InterFace<dim>* MeshManager<dim>::AddIfUnique( InterFace<dim>& interface )
{
  set<InterFace<dim>*>	discovered_interfaces;
  deque<InterFace<dim>*>	current_interfaces;

  for ( auto root_interface : root_interface_group_ ) {
    discovered_interfaces.insert( root_interface );
    current_interfaces.push_back( root_interface );
    while ( !current_interfaces.empty() ) {
      InterFace<dim>*  n_ptr( *current_interfaces.begin() );
      for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
        if ( n_ptr->Neighbor( i ) == NULL ) continue;
        if ( *(n_ptr->Neighbor( i )) == interface ) return n_ptr->Neighbor( i );
        auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
        if ( new_interface.second ) current_interfaces.push_back( n_ptr->Neighbor( i ) );
      }
      current_interfaces.pop_front();
    }
  }

  InterFace<dim>* new_interface = new InterFace<dim>( interface );
  *new_interface = interface;
  n_interfaces_++;
  return new_interface;
}


/**
   Duplicates existing node inside of the MeshManager
*/
template<size_t dim>
Node<dim>* MeshManager<dim>::Duplicate( const Node<dim>& node )
  {
    Node<dim>* new_node = new Node<dim>( node );
    n_nodes_++;

    if ( root_node_group_.size() == 0 )
      root_node_group_.push_back( new_node );

    return new_node;
  }


/**
   Duplicates the corresponding element
*/
template<size_t dim>
Element<dim>* MeshManager<dim>::Duplicate( const Element<dim>& elmt )
{
  Element<dim>* new_elmt = new Element<dim>( elmt );
  n_elmts_++;

  // delete the element's nodes which are not connected to any other elements
  deque<Node<dim>*> neighbor_nodes;
  for ( auto n : new_elmt->NodeVector() ) {
    neighbor_nodes.push_back( n );
  }

  // update node-to-element pointers in remaining node objects
  for ( auto n : neighbor_nodes )
    if ( n != NULL ) {
      n->ResizeParentStorage( n->Parents() + 1 );
      n->Assign( n->Parents() - 1, new_elmt );
    }
  neighbor_nodes.clear();

  if ( root_elmt_group_.size() == 0 )
    root_elmt_group_.push_back( new_elmt );

  return new_elmt;
}


/**
Inserts the corresponding node
*/
template<size_t dim>
Node<dim>* MeshManager<dim>::Add( Node<dim>&& node )
{
  Node<dim>* new_node = new Node<dim>( node );
  n_nodes_++;

  if ( root_node_group_.size() == 0 )
    root_node_group_.push_back( new_node );

  return new_node;
}


/**
Inserts the corresponding element
*/
template<size_t dim>
Element<dim>* MeshManager<dim>::Add( Element<dim>&& elmt )
{
  Element<dim>* new_elmt = new Element<dim>( elmt );
  n_elmts_++;

  // delete the element's nodes which are not connected to any other elements
  deque<Node<dim>*> neighbor_nodes;
  for ( auto n : new_elmt->NodeVector() ) {
    neighbor_nodes.push_back( n );
  }

  // update node-to-element pointers in remaining node objects
  for ( auto n : neighbor_nodes )
    if ( n != NULL ) {
      n->ResizeParentStorage( n->Parents() + 1 );
      n->Assign( n->Parents() - 1, new_elmt );
    }
  neighbor_nodes.clear();

  if ( root_elmt_group_.size() == 0 )
    root_elmt_group_.push_back( new_elmt );

  return new_elmt;
}


/**
Inserts the corresponding face
*/
template<size_t dim>
Face<dim>* MeshManager<dim>::Add( Face<dim>&& face )
{
  Face<dim>* new_face = new Face<dim>( face );
  n_faces_++;

  if ( root_face_group_.size() == 0 )
    root_face_group_.push_back( new_face );

  return new_face;
}


/**
Inserts the corresponding interface
*/
template<size_t dim>
InterFace<dim>* MeshManager<dim>::Add( InterFace<dim>&& interface )
{
  InterFace<dim>* new_interface = new InterFace<dim>( interface );
  n_interfaces_++;

  if ( root_interface_group_.size() == 0 )
    root_interface_group_.push_back( new_interface );

  return new_interface;
}


/**
Deletes the corresponding node
*/
template<size_t dim>
void MeshManager<dim>::Erase( Node<dim>& node )
{
  set<Node<dim>*>		discovered_nodes;
  deque<Node<dim>*>	current_nodes;

  Node<dim>* found_node = NULL;
  for ( auto root_node : root_node_group_ ) {
    discovered_nodes.insert( root_node );
    current_nodes.push_back( root_node );
    while ( !current_nodes.empty() ) {
      const Node<dim>*  n_ptr( *current_nodes.begin() );
      for ( size_t i = 0U; i < n_ptr->Parents(); i++ ) {
        for ( size_t j = 0U; j < n_ptr->Parent( i )->Nodes(); j++ ) {
          if ( j != n_ptr->ParentNodeNumber( i ) ) {
            if ( *(n_ptr->Parent( i )->N( j )) == node ) {
              found_node = n_ptr->Parent( i )->N( j ); break;
            }
            pair<typename set<Node<dim>*>::const_iterator, bool>
              new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
            if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
          }
        }
      }
      if ( found_node != NULL ) break;
      current_nodes.pop_front();
    }
    if ( found_node != NULL ) break;
  }
  discovered_nodes.clear();
  current_nodes.clear();

  // if the found node is the root node
  // change the root into one of its neighbors and delete the found node,		
  Node<dim>* new_root_node( NULL );
  Node<dim>* root_node( NULL );
  if ( root_elmt_group_.size() > 0 ) {
    for ( auto root : root_node_group_ )
      if ( found_node == root ) {
        root_node = root; break;
      }

    if ( root_node ) {
      for ( size_t i = 0U; i < root_node->Neighbors(); i++ ) {
        auto n = root_node->Neighbor( i );
        if ( new_root_node != n ) {
          new_root_node = n; break;
        }
      }
      if ( new_root_node == NULL )
        cerr << "MeshManager::Erase(Node): this is the root node which cannot be deleted. \n";
      else
        root_node = new_root_node;
    }
  }

  // delete the node
  if ( found_node ) {
    delete found_node;
    found_node = NULL;
    n_nodes_--;
  }
}


/**
Deletes the corresponding node
*/
template<size_t dim>
void MeshManager<dim>::Erase( Node<dim>* node )
{
  // if the found node is the root node
  // change the root into one of its neighbors and delete the found node,		
  Node<dim>* new_root_node( NULL );
  Node<dim>* root_node( NULL );
  if ( root_elmt_group_.size() > 0 ) {
    for ( auto root : root_node_group_ )
      if ( node == root ) {
        root_node = root; break;
      }

    if ( root_node ) {
      for ( size_t i = 0U; i < root_node->Neighbors(); i++ ) {
        auto n = root_node->Neighbor( i );
        if ( new_root_node != n ) {
          new_root_node = n; break;
        }
      }
      if ( new_root_node == NULL )
        cerr << "MeshManager::Erase(Node): this is the root node which cannot be deleted. \n";
      else
        root_node = new_root_node;
    }
  }

  // delete the node
  if ( node ) {
    delete node;
    node = NULL;
    n_nodes_--;
  }
}


/**
Deletes the corresponding element
*/
template<size_t dim>
void MeshManager<dim>::Erase( Element<dim>& elmt )
{
  set<Node<dim>*>		discovered_nodes;
  deque<Node<dim>*>	current_nodes;

  Element<dim>* found_elmt = NULL;
  for ( auto root_node : root_node_group_ ) {
    discovered_nodes.insert( root_node );
    current_nodes.push_back( root_node );
    while ( !current_nodes.empty() ) {
      const Node<dim>*  n_ptr( *current_nodes.begin() );
      for ( size_t i = 0U; i < n_ptr->Parents(); i++ ) {
        for ( size_t j = 0U; j < n_ptr->Parent( i )->Nodes(); j++ ) {
          if ( j != n_ptr->ParentNodeNumber( i ) ) {
            pair<typename set<Node<dim>*>::const_iterator, bool>
              new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
            if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
          }
        }
        if ( *n_ptr->Parent( i ) == elmt ) {
          found_elmt = n_ptr->Parent( i ); break;
        }
      }
      if ( found_elmt != NULL ) break;
      current_nodes.pop_front();
    }
    if ( found_elmt != NULL ) break;
  }
  discovered_nodes.clear();
  current_nodes.clear();

  // if the found element is the root element
  // change the root element into one of its neighbors and delete the found element,		
  Element<dim>* new_root_elmt( NULL );
  if ( root_elmt_group_.size() > 0 ) {
    size_t group_idx( 0U );
    for ( auto root : root_elmt_group_ ) {
      if ( found_elmt == root ) { // if the found element is the root element				
        for ( size_t i = 0U; i < root_node_group_[group_idx]->Neighbors(); i++ ) {
          Node<dim>* n = root_node_group_[group_idx]->Neighbor( i );
          for ( size_t j = 0U; j < n->Parents(); j++ ) {
            auto e = n->Parent( j );
            if ( found_elmt != e ) {
              new_root_elmt = e;
              break;
            }
          }
        }
        if ( new_root_elmt == NULL )
          cerr << "MeshManager::Erase(Element): this is the root element which cannot be deleted. \n";
        else
          root_elmt_group_[group_idx] = new_root_elmt;
        break;
      }
      group_idx++;
    }
  }

  // delete the element	
  if ( found_elmt ) {
    // delete the element's nodes which are not connected to any other elements
    deque<Node<dim>*> neighbor_nodes;
    for ( auto n : found_elmt->NodeVector() ) {
      neighbor_nodes.push_back( n );
    }

    // update node-to-element pointers in remaining node objects
    for ( auto n : neighbor_nodes )
      if ( n != NULL )
        n->Unassign( found_elmt );
    neighbor_nodes.clear();

    delete found_elmt;
    found_elmt = NULL;
    n_elmts_--;
  }
}


/**
Deletes the corresponding element
*/
template<size_t dim>
void MeshManager<dim>::Erase( Element<dim>* elmt )
{
  // if the found element is the root element
  // change the root element into one of its neighbors and delete the found element,		
  Element<dim>* new_root_elmt( NULL );
  if ( root_elmt_group_.size() > 0 ) {
    size_t group_idx( 0U );
    for ( auto root : root_elmt_group_ ) {
      if ( elmt == root ) { // if the found element is the root element				
        for ( size_t i = 0U; i < root_node_group_[group_idx]->Neighbors(); i++ ) {
          Node<dim>* n = root_node_group_[group_idx]->Neighbor( i );
          for ( size_t j = 0U; j < n->Parents(); j++ ) {
            auto e = n->Parent( j );
            if ( elmt != e ) {
              new_root_elmt = e;
              break;
            }
          }
        }
        if ( new_root_elmt == NULL )
          cerr << "MeshManager::Erase(*Element): this is the root element which cannot be deleted. \n";
        else
          root_elmt_group_[group_idx] = new_root_elmt;
        break;
      }
      group_idx++;
    }
  }

  // delete the element	
  if ( elmt ) {
    // delete the element's nodes which are not connected to any other elements
    deque<Node<dim>*> neighbor_nodes;
    for ( auto n : elmt->NodeVector() ) {
      neighbor_nodes.push_back( n );
    }

    // update node-to-element pointers in remaining node objects
    for ( auto n : neighbor_nodes )
      if ( n != NULL )
        n->Unassign( elmt );
    neighbor_nodes.clear();

    delete elmt;
    elmt = NULL;
    n_elmts_--;
  }
}


/**
Deletes the corresponding face
*/
template<size_t dim>
void MeshManager<dim>::Erase( Face<dim>& face )
{
  // traversal of the existing mesh nodes to find all its faces
  set<Face<dim>*>		discovered_faces;
  deque<Face<dim>*>	current_faces;

  Face<dim>* found_face = NULL;
  Face<dim>* found_root = NULL;

  for ( auto root_face : root_face_group_ ) {
    // starting at the first face
    discovered_faces.insert( root_face );
    current_faces.push_back( root_face );
    found_root = root_face;
    while ( !current_faces.empty() ) {
      Face<dim>*  n_ptr( *current_faces.begin() );
      // for all neighbor faces of the current face
      for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
        if ( n_ptr->Neighbor( i ) == NULL ) continue;
        if ( *n_ptr->Neighbor( i ) == face ) found_face = n_ptr->Neighbor( i );

        // if this neighbor is new one					
        auto new_face = discovered_faces.insert( n_ptr->Neighbor( i ) );
        if ( new_face.second ) current_faces.push_back( n_ptr->Neighbor( i ) );
      }
      if ( found_face != NULL ) break;
      // removing the face from the discovered (but not yet explored) deque
      current_faces.pop_front();
    }
    if ( found_face != NULL ) break;
  }

  // change the root face into one of its neighbors and delete the found face,
  // if the found face is the root face
  Face<dim>* new_root_face( NULL );
  if ( found_face == found_root ) {
    for ( auto f : discovered_faces ) {
      if ( found_root != f ) {
        new_root_face = f;
        break;
      }
    }
    if ( new_root_face == NULL )
      cerr << "MeshManager::Erase(Face): this is the root face which cannot be deleted. \n";
    else
      found_root = new_root_face;
  }

  if ( found_face ) {
    // unassign the connections of its neighbors
    for ( auto n : found_face->NeighborElementVector() )
      if ( n != NULL ) n->DisconnectNeighbor( found_face );

    delete found_face;
    found_face = NULL;
    n_faces_--;
  }

  discovered_faces.clear();
  current_faces.clear();
}


/**
Deletes the corresponding face
*/
template<size_t dim>
void MeshManager<dim>::Erase( Face<dim>* face )
{
  // if the found face is the root face
  // change the root face into one of its neighbors and delete the found face.
  Face<dim>* new_root_face( NULL );
  Face<dim>* root_face( NULL );
  if ( root_face_group_.size() > 0 ) {
    for ( auto root : root_face_group_ )
      if ( face == root ) root_face = root;

    if ( root_face ) {
      for ( auto e : root_face->NeighborElementVector() ) {
        if ( root_face != e ) {
          new_root_face = e;
          break;
        }
      }
      if ( new_root_face == NULL )
        cerr << "MeshManager::Erase(Face): this is the root face which cannot be deleted. \n";
      else
        root_face = new_root_face;
    }
  }

  // delete the face
  if ( face ) {
    // unassign the connections of its neighbors
    for ( auto n : face->NeighborElementVector() )
      if ( n != NULL ) n->DisconnectNeighbor( face );

    delete face;
    face = NULL;
    n_faces_--;
  }
}


/**
Deletes the corresponding interface
*/
template<size_t dim>
void MeshManager<dim>::Erase( InterFace<dim>& interface )
{
  // traversal of the existing mesh nodes to find all its interfaces
  set<InterFace<dim>*>	discovered_interfaces;
  deque<InterFace<dim>*>	current_interfaces;

  InterFace<dim>* found_interface = NULL;
  InterFace<dim>*	found_root = NULL;

  for ( auto root_interface : root_interface_group_ ) {
    // starting at the first interface
    discovered_interfaces.insert( root_interface );
    current_interfaces.push_back( root_interface );
    while ( !current_interfaces.empty() ) {
      InterFace<dim>*  n_ptr( *current_interfaces.begin() );
      // for all neighbor interfaces of the current interface
      for ( size_t i = 0U; i < n_ptr->Neighbors(); i++ ) {
        if ( n_ptr->Neighbor( i ) == NULL ) continue;
        if ( *n_ptr->Neighbor( i ) == interface ) found_interface = n_ptr->Neighbor( i );

        // if this neighbor is new one					
        auto new_interface = discovered_interfaces.insert( n_ptr->Neighbor( i ) );
        if ( new_interface.second ) current_interfaces.push_back( n_ptr->Neighbor( i ) );
      }
      if ( found_interface != NULL ) break;
      // removing the interface from the discovered (but not yet explored) deque			
      current_interfaces.pop_front();
    }
    if ( found_interface != NULL ) break;
  }

  // change the root face into one of its neighbors and delete the found face,
  // if the found face is the root face
  InterFace<dim>* new_root_face( NULL );
  if ( found_interface == found_root ) {
    for ( auto f : discovered_interfaces ) {
      if ( found_root != f ) {
        new_root_face = f;
        break;
      }
    }
    if ( new_root_face == NULL )
      cerr << "MeshManager::Erase(InterFace): this is the root interface which cannot be deleted. \n";
    else
      found_root = new_root_face;
  }

  if ( found_interface ) {
    // unassign the connections of its neighbors
    for ( auto n : found_interface->NeighborElementVector() )
      if ( n != NULL ) n->DisconnectNeighbor( found_interface );

    delete found_interface;
    found_interface = NULL;
    n_interfaces_--;
  }

  discovered_interfaces.clear();
  current_interfaces.clear();
}


/**
Deletes the corresponding interface
*/
template<size_t dim>
void MeshManager<dim>::Erase( InterFace<dim>* interface )
{
  // if the found interface is the root interface
  // change the root interface into one of its neighbors and delete the found interface.
  InterFace<dim>* new_root_interface( NULL );
  InterFace<dim>* root_interface( NULL );
  if ( root_interface_group_.size() > 0 ) {
    for ( auto root : root_interface_group_ )
      if ( interface == root ) root_interface = root;

    if ( root_interface ) {
      for ( auto e : root_interface->NeighborElementVector() ) {
        if ( root_interface != e ) {
          new_root_interface = e;
          break;
        }
      }
      if ( new_root_interface == NULL )
        cerr << "MeshManager::Erase(InterFace): this is the root interface which cannot be deleted. \n";
      else
        root_interface = new_root_interface;
    }
  }

  // delete the interface
  if ( interface ) {
    // unassign the connections of its neighbors
    for ( auto n : interface->NeighborElementVector() )
      if ( n != NULL ) n->DisconnectNeighbor( interface );

    delete interface;
    interface = NULL;
    n_interfaces_--;
  }
}


/**
erases all nodes, discerning those that do not have any parent element connections.
*/
template<size_t dim>
bool MeshManager<dim>::EraseNodes()
{
  bool emptyNodes( false );

  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );

  for ( auto n : nodes ) {
    if ( n ) {
      delete n;
      n = NULL;
      n_nodes_--;
    }
  }

  if ( n_nodes_ == 0 ) {
    root_node_group_.clear();
    emptyNodes = true;
  }

  return emptyNodes;
}


/**
erases all elements, discerning those that have any connections.
*/
template<size_t dim>
bool MeshManager<dim>::EraseElements()
{
  bool emptyElements( false );

  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );

  for ( auto e : elmts ) {
    if ( e ) {
      delete e;
      e = NULL;
      n_elmts_--;
    }
  }

  for ( auto n : nodes ) {
    if ( n ) {
      delete n;
      n = NULL;
      n_nodes_--;
    }
  }

  if ( n_elmts_ == 0 && n_nodes_ == 0 ) {
    root_node_group_.clear();
    root_elmt_group_.clear();
    emptyElements = true;
  }

  return emptyElements;
}


/**
erases all faces, discerning those that have any connections.
*/
template<size_t dim>
bool MeshManager<dim>::EraseFaces()
{
  bool emptyFaces( false );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );

  // 2. assigning and trimming excess storage from the face pointer vector
  for ( auto f : faces ) {
    if ( f ) {
      delete f;
      f = NULL;
      n_faces_--;
    }
  }

  if ( n_faces_ == 0 ) {
    root_face_group_.clear();
    emptyFaces = true;
  }

  return emptyFaces;
}


/**
erases all interfaces, discerning those that have any connections.
*/
template<size_t dim>
bool MeshManager<dim>::EraseInterFaces()
{
  bool emptyInterFaces( false );

  // traversal of the existing mesh root interfaces to find all its interfaces
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );

  // 2. assigning and trimming excess storage from the interface pointer vector
  for ( auto f : interfaces ) {
    if ( f ) {
      delete f;
      f = NULL;
      n_interfaces_--;
    }
  }

  if ( n_interfaces_ == 0 ) {
    root_interface_group_.clear();
    emptyInterFaces = true;
  }

  return emptyInterFaces;
}


/** (Re)number all cells; either continuous for all cells or seperate ranges for all entity types

@note this member function is constant because the idx_ is a mutable variable in the cell classes
*/
template<size_t dim>
void MeshManager<dim>::AssignUniqueNumbers( bool in_a_single_sequence ) const
{
  // traversal of the existing mesh nodes to find all its elements	
  deque<const Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<const Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<const InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );
  sort( interfaces.begin(), interfaces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  size_t n( 0U );

  for_each( nodes.begin(), nodes.end(), [&n]( const Node<dim>* o ) { o->Idx( n++ ); return o; } );

  n = 0U; // resetting the counter
  for_each( elmts.begin(), elmts.end(), [&n]( const Element<dim>* o ) { o->Idx( n++ ); return o; } );

  if ( !in_a_single_sequence ) n = 0U;
  for_each( faces.begin(), faces.end(), [&n]( const Face<dim>* o ) { o->Idx( n++ ); return o; } );

  if ( !in_a_single_sequence ) n = 0U;
  for_each( interfaces.begin(), interfaces.end(), [&n]( const InterFace<dim>* o ) { o->Idx( n++ ); return o; } );

} // end AssignUniqueNumbers


/** Writes MeshManager to VSet

Writes the elements, faces and interfaces
stored in the current MeshManager
to the supplied VSet.
This also includes the connectivity information,
i.e. the connections between these entities.

@param vset A VSet preferably empty, which is
resized first, if necessary and into which the
MeshManager connectivity information is input.

@section conventions Conventions

Connections will exist between elements of the
topology, i.e., between:

- volume elements
- surface elements
- line elements

these connections are made via pointers and exist
across the faces, edges, and end end point of these entities.

Where there is no neighbor, like for a face that sits on a model
boundary, the neighbor pointer will be a 'NULL'.

Special rules apply to Face and Interface objects that
are essentially elements, but they are lower-dimensional
(for instance surfaces in a volumetric model)

Faces and Interfaces also have extra connections
to their higher-dimensional "parent" elements.

All these entities are stored in the VSet in a specific order:

1. highest dimensional elements (volumes in 3D)
2. dim-1 elements
3. dim-2 elements, if any
4. Face objects
5. InterFace objects

Since pointers cannot be written to file as such,
the first step in the output is to give all the entities
a unique numbering: 0..n-1 elements, then faces, last interfaces.

@section application Application

The data of the MeshManager is output to a VSet, for instance, when the
whole model is to be saved to a binary file.

@attention If the VSet is initialized successfully
the method will report this.

SKM 15/9/2014: fixed major bug in the output code blocks for
plist and pfverts.

@attention renumber the elements of the mesh first, using

@code
// renumbering the mesh; face numbers follow those of the elements
// last are the interface objects
// NB: we remember how many faces and interfaces there were to deduce the offsets at a later point
const bool in_a_single_sequence(true);
AssignUniqueNumbers( in_a_single_sequence );
@endcode

*/
template<size_t dim>
void MeshManager<dim>::OutputMeshTo( VSet<dim>& vset ) const
{
  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  // 1. resizing the VSet
  // --------------------
  const size_t higherDimParents( 2U );
  const size_t interfaceMultiplier( 2U );

  // traversal of the existing mesh nodes to find all its elements	
  deque<const Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
    
  // traversal of the existing mesh root faces to find all its faces	
  deque<const Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<const InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );
  sort( interfaces.begin(), interfaces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  if ( HybridElementMesh() || faces.size() > 0 || interfaces.size() > 0 ) {
    deque<size_t>  nodes_per_element;
    deque<size_t>  neighbors_per_element;
    deque<int32>   csp_fem_types;

    // 1.1 identifying how many nodes and neighbors there are per element		
    for ( auto e : elmts ) {
      nodes_per_element.push_back( e->Nodes() );
      neighbors_per_element.push_back( e->Neighbors() );
      csp_fem_types.push_back( e->FE_Type() );
    }

    // 1.2 adding Face information after the elements
    for ( auto f : faces ) {
      nodes_per_element.push_back( f->Nodes() );
      neighbors_per_element.push_back( f->Neighbors() + higherDimParents );
      csp_fem_types.push_back( f->FE_Type() );
    }

    // 1.3 adding InterFace information after the faces
    for ( auto f : interfaces ) {
      // multiplier takes care of the multiplicated interface nodes that the InterFace will be connected to
      nodes_per_element.push_back( f->Nodes() * interfaceMultiplier );
      neighbors_per_element.push_back( f->Neighbors() * interfaceMultiplier + higherDimParents );
      csp_fem_types.push_back( f->FE_Type() );
    }
    // 1.4 resizing the VSet
    vset.Resize( csp_fem_types, nodes_per_element, neighbors_per_element,
                 nodes.size(), faces.size(), interfaces.size() );
  }
  else { // if there is only a single element type
    auto el = RootElement( 0 );
    auto fe = el->FE();
    assert( fe != NULL );
    vset.Resize( fe->Nodes(),
                 fe->Neighbors(),
                 fe->ElementType(),
                 Nodes(), Elements() );

    vset.ElementType( 0, el->FE_Type() );
  }

  // 2. adding node coordinates and boundary flags (BOX_BOUNDARY)
  // ------------------------------------------------------------
  size_t i( 0U );
  if ( dim == 1U ) {
    for ( auto n : nodes ) {
      vset.Px( i, n->x() );
      ++i;
    }
  }
  else if ( dim == 2U ) {
    for ( auto n : nodes ) {
      vset.Px( i, n->x() );
      vset.Py( i, n->y() );
      ++i;
    }
  }
  else {
    for ( auto n : nodes ) {
      vset.Px( i, n->x() );
      vset.Py( i, n->y() );
      vset.Pz( i, n->z() );
      ++i;
    }
  }
  
  // 3. adding 'plist' connectivity list
  // -----------------------------------
  size_t eidx = 0;
  // elements
  for ( auto e : elmts ) {
    for ( size_t j = 0U; j<e->Nodes(); ++j )
      vset.Plist( eidx, j, (e->N( j )->Idx()) );
    ++eidx;
  }

  // faces
  for ( auto f : faces ) {
    for ( size_t j = 0U; j<f->Nodes(); ++j )
      vset.Plist( eidx, j, (f->N( j )->Idx()) );
    ++eidx;
  }

  // interfaces
  for ( auto f : interfaces ) {
    for ( size_t j = 0U; j<f->Nodes(); ++j )
      vset.Plist( eidx, j, (f->N( j )->Idx()) );
    ++eidx;
  }


  // 4. adding 'pfverts' neighbors per element list
  // ----------------------------------------------
  eidx = 0;

  // 'pfverts' elements
  for ( auto e : elmts ) {
    for ( size_t j = 0U; j<e->Neighbors(); ++j ) {
      Element<dim>* const ptr( e->Neighbor( j ) );
      if ( ptr != NULL )
        vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
      else
        vset.Pfvert( eidx, j, e->AtBoundary() );
    }
    ++eidx;
  }

  // 'pfverts' faces
  // ---------------
  // the faces are stored after the elements including connections to their higher-dimensional neighbors
  // add the end of the pfverts entries
  for ( auto f : faces ) {
    // equidimensional neighbors first
    const size_t neighbors( f->Neighbors() );
    for ( size_t j = 0U; j<neighbors; ++j ) {
      Face<dim>* const ptr( f->Neighbor( j ) );
      // if the neighbor exists (which it must on the inside of the Face)
      if ( ptr != NULL )
        vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
      else
        vset.Pfvert( eidx, j, f->InnerParent()->AtBoundary() );
    }
    // higher-dimensional neighbors second
    // inner neighbor
    if ( dim == 3U ) assert( f->InnerParent()->IsVolumeElement() );
    else if ( dim == 2U ) assert( f->InnerParent()->IsSurfaceElement() );
    assert( f->InnerParent()->Idx() < Elements() );
    vset.Pfvert( eidx, neighbors, static_cast<int32>(f->InnerParent()->Idx()) );
    // outer neighbor
    if ( f->OuterParent() != NULL && dim == 3U ) assert( f->OuterParent()->IsVolumeElement() );
    else if ( f->OuterParent() != NULL && dim == 2U ) assert( f->OuterParent()->IsSurfaceElement() );
    if ( f->OuterParent() != NULL ) {
      assert( f->OuterParent()->Idx() < Elements() );
      vset.Pfvert( eidx, neighbors + 1U, static_cast<int32>(f->OuterParent()->Idx()) );
    }
    else {
      // if there is no neighbor, the inner element parent should be at the model boundary
      if ( f->InnerParent()->AtBoundary() == NOT ) {
        f->Out();
        csmp_error.notice( WARNING, "MeshManager<dim>::OutputMeshTo (face neighbors):",
                           "inner dim+1 element should be at model boundary because Face has no outer element." );
        f->InnerParent()->AtBoundary( IRREGULAR );
      }
      vset.Pfvert( eidx, neighbors + 1U, f->InnerParent()->AtBoundary() );
    }
    ++eidx;
  }

  // 'pfverts' interfaces
  // --------------------
  for ( auto f : interfaces ) {
    // equidimensional neighbors (=other interfaces) first
    const size_t neighbors( f->Neighbors() );
    for ( size_t j = 0U; j<neighbors; ++j ) {
      InterFace<dim>* const ptr( f->Neighbor( j ) );
      if ( ptr != NULL ) {
        // building search maps that we will use to find the shared interfaces
        // key=pointset   face iD
        map<set<Point<dim> >, pair<INTERFACE_SIDE, size_t> >   inner_elmt_faces, outer_elmt_faces;
        vector<size_t>  nids;
        // first element
        Element<dim>* e1 = f->InnerParent();
        for ( size_t face = 0U; face<e1->Faces(); ++face ) {
          e1->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j = 0U; j<nids.size(); ++j )
            face_key.insert( e1->N( nids[j] )->Coordinate() );
          outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
        }
        // second element
        Element<dim>* e2 = f->OuterParent();
        for ( size_t face = 0U; face<e2->Faces(); ++face ) {
          e2->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j = 0U; j<nids.size(); ++j )
            face_key.insert( e2->N( nids[j] )->Coordinate() );
          inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
        }

        // 2. finding the shared faces
        bool found( false );
        long64 inner_face_id( -1 ), outer_face_id( -1 );
        for ( auto inner_face : inner_elmt_faces ) {
          for ( auto outer_face : outer_elmt_faces ) {
            if ( inner_face.first == outer_face.first ) {
              inner_face_id = inner_face.second.second;
              outer_face_id = outer_face.second.second;
              found = true;
              break;
            }
          }
          if ( found ) break;
        }

        vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
      }
      else
        vset.Pfvert( eidx, j, REGION_BOUNDARY );
    }
    // inner and outer higher-dimensional neighbors
    // ( they must always exist because SplitBoundaries are internal model boundaries)
    assert( f->InnerParent() != NULL );
    assert( f->OuterParent() != NULL );
    assert( f->InnerParent()->Idx() < elmts.size() );
    assert( f->OuterParent()->Idx() < elmts.size() );

    vset.Pfvert( eidx, neighbors, static_cast<int32>(f->InnerParent()->Idx()) );
    vset.Pfvert( eidx, neighbors + 1U, static_cast<int32>(f->OuterParent()->Idx()) );    
    if ( f->InnerParent()->Idx() < 0 || f->OuterParent()->Idx() < 0 )
      cout << "\n\t Negative Interface Idx in Element #" << eidx;
    ++eidx;
  }

  // 5. adding boundary flags
  // ------------------------
  for ( auto n : nodes ) {
    if ( n->AtBoundary() != NOT ) vset.AddBFlag( n->Idx(), n->AtBoundary() );
  }

  cout << "\nMeshManager<" << dim << ">::OutputMeshTo: MeshManager successfully output to VSet..." << endl;

} // end OutputMeshTo


/**
Storing distributed variables associated with the mesh in the VSet

For all finite volumes and elements, their integration points and nodes,
but not for any of the variables stored on the model, region, boundaries or splitboundaries
this method stores the current values in the VSet.

To store the properties in the VSet, they are first written to PropertyData objects.
These are then added to the VSet property storage.

@attention  a continuous numbering of elements, faces, interfaces and nodes has to be
created with AssignUniqueNumbers() before this method is called.

@note region properties are stored together with the regions in respective binary files

@note this method is anything, but nice. Yet it will be quite a challenge to come up with a better
design; hopefully there will be no more extra variable placements or types in the future.

@attention in the VSet, the variables are identified only by their (unique) names. The property
database is therefore essential to retrieve all other variable related information.

@author SKM 5/5/2016

*/
template<size_t dim>
void MeshManager<dim>::OutputStoredVariablesTo( const PropertyDatabase<dim>& database, VSet<dim>& vset ) const
{
  ErrorHandler&		csmp_error( ErrorHandler::Instance() );
  map<string, Index>  properties;

  // traversal of the existing mesh nodes to find all its elements	
  deque<const Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<const Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<const InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );
  sort( interfaces.begin(), interfaces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  database.ListProperties( NODE, properties );
  const size_t n_nodes( nodes.size() );
  // for all node properties
  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    // setting the specifications for the property storage (no memory allocation yet)
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    // for the given property type
    const size_t flag_capacity( n_nodes * (*pit).second.flagDepth );
    const size_t data_capacity( n_nodes * (*pit).second.dataDepth );
    // allocating memory to store the property flags and values
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        // estimating the storage required
        ScalarVariable value;
        for ( auto it : nodes ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( auto it : nodes ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( auto it : nodes ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        for ( auto it : nodes ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        for ( auto it : nodes ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of node variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }

  // -------------------------------------------------
  // -------------------------------------------------
  // element properties (including integration points)
  // -------------------------------------------------
  // -------------------------------------------------
  database.ListProperties( ELEMENT, properties );
  const size_t elements( elmts.size() );
  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    // creating the property storage
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    // for the given property type
    const size_t flag_capacity( elements * (*pit).second.flagDepth );
    const size_t data_capacity( elements * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( auto it : elmts ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( auto it : elmts ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( auto it : elmts ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        for ( auto it : elmts ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        for ( auto it : elmts ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of element variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }

  // element integration point properties
  // ------------------------------------
  if ( database.ListProperties( ELEMENT_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_ips( RootElement( 0 )->IntegrationPoints() ); // just an estimate

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      // creating the property storage
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( elements * elmt_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * elmt_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : elmts ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : elmts ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : elmts ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : elmts ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : elmts ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of element integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // element sector integration point properties
  // -------------------------------------------
  if ( database.ListProperties( SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_sector_ips( RootElement( 0 )->IntegrationPointsPerSector() );
    const size_t sectors_per_element( RootElement( 0 )->Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : elmts ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : elmts ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : elmts ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : elmts ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : elmts ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of element sector integraton point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // element facet integration point properties
  // ------------------------------------------
  if ( database.ListProperties( FACET_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_facet_ips( RootElement( 0 )->IntegrationPointsPerFacet() );
    const size_t facets_per_element( RootElement( 0 )->Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : elmts ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : elmts ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : elmts ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : elmts ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : elmts ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // ----------------------------------------------
  // ----------------------------------------------
  // face properties (including integration points)
  // ----------------------------------------------
  // ----------------------------------------------
  database.ListProperties( FACE, properties );
  const size_t n_faces( Faces() );

  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    const size_t flag_capacity( n_faces * (*pit).second.flagDepth );
    const size_t data_capacity( n_faces  * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( auto it : faces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( auto it : faces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( auto it : faces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        for ( auto it : faces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        for ( auto it : faces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }


  // face integration point properties
  // ---------------------------------
  if ( database.ListProperties( FACE_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_ips( RootFace( 0 )->IntegrationPoints() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * face_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * face_ips  * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : faces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : faces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : faces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : faces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : faces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // face sector integration point properties
  // ----------------------------------------
  if ( database.ListProperties( FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_sector_ips( RootFace( 0 )->IntegrationPointsPerSector() );
    const size_t sectors_per_face( RootFace( 0 )->Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * sectors_per_face * face_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * sectors_per_face * face_sector_ips  * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : faces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : faces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : faces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : faces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : faces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face-sector integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // face facet integration point properties
  // ---------------------------------------
  if ( database.ListProperties( FACE_FACET_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_facet_ips( RootFace( 0 )->IntegrationPointsPerFacet() );
    const size_t facets_per_face( RootFace( 0 )->Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * facets_per_face * face_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * facets_per_face * face_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : faces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : faces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : faces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : faces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : faces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // ---------------------------------------------------
  // ---------------------------------------------------
  // interface properties (including integration points)
  // ---------------------------------------------------
  // ---------------------------------------------------
  database.ListProperties( INTER_FACE, properties );
  const size_t n_interfaces( InterFaces() );

  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    const size_t flag_capacity( n_interfaces * (*pit).second.flagDepth );
    const size_t data_capacity( n_interfaces  * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( auto it : interfaces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( auto it : interfaces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( auto it : interfaces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        for ( auto it : interfaces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        for ( auto it : interfaces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }

  // interface integration point properties
  // --------------------------------------
  if ( database.ListProperties( INTER_FACE_INTEGRATION_POINT, properties ) > 0 && InterFaces() > 0 ) {
    const size_t interface_ips( RootInterFace( 0 )->IntegrationPoints() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * interface_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * interface_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : interfaces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : interfaces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : interfaces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : interfaces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : interfaces ) {
            const size_t integration_points( (*it).IntegrationPoints() );
            for ( size_t i = 0U; i<integration_points; ++i ) {
              (*it).Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of interface integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }

  // interface sector integration point properties
  // ---------------------------------------------
  if ( database.ListProperties( INTER_FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 && InterFaces() > 0 ) {
    assert( InterFaces() > 0 );
    const size_t interface_sector_ips( RootInterFace( 0 )->IntegrationPointsPerSector() );
    const size_t sectors_per_interface( RootInterFace( 0 )->Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : interfaces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : interfaces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : interfaces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : interfaces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : interfaces ) {
            const size_t sectors( (*it).Sectors() );
            for ( size_t i = 0U; i<sectors; ++i ) {
              const size_t ips_per_sector( (*it).IntegrationPointsPerSector() );
              for ( size_t j = 0U; j<ips_per_sector; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of interface sector integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }

  // interface facet integration point properties
  // --------------------------------------------
  if ( database.ListProperties( INTER_FACE_FACET_INTEGRATION_POINT, properties ) > 0 && InterFaces() > 0 ) {
    assert( InterFaces() > 0 );
    const size_t interface_facet_ips( RootInterFace( 0 )->IntegrationPointsPerFacet() );
    const size_t facets_per_interface( RootInterFace( 0 )->Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * facets_per_interface * interface_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * facets_per_interface * interface_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( auto it : interfaces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( auto it : interfaces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( auto it : interfaces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                     break;
        case ARRAY: {
          ArrayVariable value;
          for ( auto it : interfaces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                    break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value;
          for ( auto it : interfaces ) {
            const size_t facets( (*it).Facets() );
            for ( size_t i = 0U; i<facets; ++i ) {
              const size_t ips_per_facet( (*it).IntegrationPointsPerFacet() );
              for ( size_t j = 0U; j<ips_per_facet; ++j ) {
                (*it).Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
                           break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of interface facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }

} // end OutputStoredVariablesTo


  /**
  Read and assign property values written exactly by OutputStoredVariablesTo()

  @ test SKM 28/6/2016
  */
template<size_t dim>
void MeshManager<dim>::InputStoredVariablesFrom( const PropertyDatabase<dim>& database, const VSet<dim>& vset )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );
  sort( interfaces.begin(), interfaces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // -------------------
  // assigns properties
  // -------------------
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    // apart from the name string key in the map, PropertyData contains the most important variable specifications
    if ( (*pit).second.Placement() != NODE ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == NODE );
    assert( (*pit).second.Size() == nodes.size() * key.dataDepth );

    // for the given property type
    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t i( 0U );
        for ( auto n : nodes ) {
          read( (*pit).second, i++, value );
          n->Store( key, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t i( 0U );
        for ( auto n : nodes ) {
          read( (*pit).second, i++, value );
          n->Store( key, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t i( 0U );
        for ( auto n : nodes ) {
          read( (*pit).second, i++, value );
          n->Store( key, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t i( 0U );
        for ( auto n : nodes ) {
          read( (*pit).second, i++, value );
          n->Store( key, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t i( 0U );
        for ( auto n : nodes ) {
          read( (*pit).second, i++, value );
          n->Store( key, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of node variable not recognized." );
    }
  }

  // -------------------------------------------------
  // -------------------------------------------------
  // element properties (including integration points)
  // -------------------------------------------------
  // -------------------------------------------------
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != ELEMENT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == ELEMENT );
    assert( (*pit).second.Size() / key.dataDepth == Elements() );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t i( 0U );
        for ( auto e : elmts ) {
          read( (*pit).second, e->Idx(), value );
          e->Store( key, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t i( 0U );
        for ( auto e : elmts ) {
          read( (*pit).second, i++, value );
          e->Store( key, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t i( 0U );
        for ( auto e : elmts ) {
          read( (*pit).second, i++, value );
          e->Store( key, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t i( 0U );
        for ( auto e : elmts ) {
          read( (*pit).second, i++, value );
          e->Store( key, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t i( 0U );
        for ( auto e : elmts ) {
          read( (*pit).second, i++, value );
          e->Store( key, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of element variable not recognized." );
    }
  }

  // element integration point properties
  // ------------------------------------
  // ELEMENT_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != ELEMENT_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == ELEMENT_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U ); // running index
        for ( auto e : elmts ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of element integration point variable not recognized." );
    }
  }

  // element sector integration point properties
  // -------------------------------------------
  // SECTOR_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != SECTOR_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == SECTOR_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of element sector integraton point variable not recognized." );
    }
  }

  // element facet integration point properties
  // ------------------------------------------
  // FACET_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACET_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACET_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : elmts ) {
          if ( e->FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of facet integration point variable not recognized." );
    }
  }



  // ----------------------------------------------
  // ----------------------------------------------
  // face properties (including integration points)
  // ----------------------------------------------
  // ----------------------------------------------
  // FACE
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE );
    assert( (*pit).second.Size() / key.dataDepth == Faces() );

    // for the given property type
    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t i( 0U );
        for ( auto f : faces ) {
          read( (*pit).second, i++, value );
          f->Store( key, value );
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t i( 0U );
        for ( auto f : faces ) {
          read( (*pit).second, i++, value );
          f->Store( key, value );
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t i( 0U );
        for ( auto f : faces ) {
          read( (*pit).second, i++, value );
          f->Store( key, value );
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t i( 0U );
        for ( auto f : faces ) {
          read( (*pit).second, i++, value );
          f->Store( key, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t i( 0U );
        for ( auto f : faces ) {
          read( (*pit).second, i++, value );
          f->Store( key, value );
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face variable not recognized." );
    }
  }

  // face integration point properties
  // ---------------------------------
  // FACE_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face integration point variable not recognized." );
    }
  }

  // face sector integration point properties
  // ----------------------------------------
  // FACE_SECTOR_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE_SECTOR_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE_SECTOR_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face-sector integration point variable not recognized." );
    }
  }

  // face facet integration point properties
  // ---------------------------------------
  // FACE_FACET_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE_FACET_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE_FACET_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : faces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face facet integration point variable not recognized." );
    }
  }



  // ---------------------------------------------------
  // ---------------------------------------------------
  // interface properties (including integration points)
  // ---------------------------------------------------
  // ---------------------------------------------------
  // INTER_FACE
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE );
    assert( (*pit).second.Size() / key.dataDepth == InterFaces() );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t i( 0U );
        for ( auto e : interfaces ) {
          read( (*pit).second, i, value );
          e->Store( key, value );
          i++;
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t i( 0U );
        for ( auto e : interfaces ) {
          read( (*pit).second, i, value );
          e->Store( key, value );
          i++;
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t i( 0U );
        for ( auto e : interfaces ) {
          read( (*pit).second, i, value );
          e->Store( key, value );
          i++;
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t i( 0U );
        for ( auto e : interfaces ) {
          read( (*pit).second, i, value );
          e->Store( key, value );
          i++;
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t i( 0U );
        for ( auto e : interfaces ) {
          read( (*pit).second, i, value );
          e->Store( key, value );
          i++;
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface variable not recognized." );
    }
  }

  // interface integration point properties
  // --------------------------------------
  // INTER_FACE_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t integration_points( e->IntegrationPoints() );
          for ( size_t i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e->Store( i, key, value );
            entry++;
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface integration point variable not recognized." );
    }
  }

  // interface sector integration point properties
  // ---------------------------------------------
  // INTER_FACE_SECTOR_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE_SECTOR_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE_SECTOR_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t sectors( e->Sectors() );
          for ( size_t i = 0U; i<sectors; ++i ) {
            const size_t ips_per_sector( e->IntegrationPointsPerSector() );
            for ( size_t j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface sector integration point variable not recognized." );
    }
  }

  // interface facet integration point properties
  // --------------------------------------------
  // INTER_FACE_FACET_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE_FACET_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE_FACET_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value;
        size_t entry( 0U );
        for ( auto e : interfaces ) {
          const size_t facets( e->Facets() );
          for ( size_t i = 0U; i<facets; ++i ) {
            const size_t ips_per_facet( e->IntegrationPointsPerFacet() );
            for ( size_t j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e->Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
                         break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface facet integration point variable not recognized." );
    }
  }

} // end InputStoredVariablesFrom


template<size_t dim>
bool  MeshManager<dim>::HybridElementMesh() const
{
  return hybrid_element_mesh_;
}


template<size_t dim>
size_t  MeshManager<dim>::Nodes() const
{
  return this->n_nodes_;
}

template<size_t dim>
size_t  MeshManager<dim>::NodeGroups() const
{
  return this->root_node_group_.size();
}

template<size_t dim>
size_t  MeshManager<dim>::Elements() const
{
  return this->n_elmts_;
}

template<size_t dim>
size_t  MeshManager<dim>::ElementGroups() const
{
  return this->root_elmt_group_.size();
}

template<size_t dim>
size_t  MeshManager<dim>::Faces() const
{
  return this->n_faces_;
}


template<size_t dim>
size_t  MeshManager<dim>::FaceGroups() const
{
  return this->root_face_group_.size();
}


template<size_t dim>
size_t  MeshManager<dim>::InterFaces() const
{
  return this->n_interfaces_;
}


template<size_t dim>
size_t  MeshManager<dim>::InterFaceGroups() const
{
  return this->root_interface_group_.size();
}


template<size_t dim>
Node<dim>*  MeshManager<dim>::RootNode( size_t group_idx )
{
  assert( group_idx < this->root_node_group_.size() );

  return this->root_node_group_[group_idx];
}


template<size_t dim>
void  MeshManager<dim>::SetRootNode( Node<dim>* root_node )
{
  this->root_node_group_.push_back( root_node );
}

template<size_t dim>
Element<dim>*  MeshManager<dim>::RootElement( size_t group_idx )
{
  assert( group_idx < this->root_elmt_group_.size() );

  return this->root_elmt_group_[group_idx];
}


template<size_t dim>
void  MeshManager<dim>::SetRootElement( Element<dim>* root_elmt )
{
  this->root_elmt_group_.push_back( root_elmt );
}

template<size_t dim>
Face<dim>*  MeshManager<dim>::RootFace( size_t group_idx )
{
  assert( group_idx < this->root_face_group_.size() );

  return this->root_face_group_[group_idx];
}


template<size_t dim>
void  MeshManager<dim>::SetRootFace( Face<dim>* root_face )
{
  this->root_face_group_.push_back( root_face );
}


template<size_t dim>
InterFace<dim>*  MeshManager<dim>::RootInterFace( size_t group_idx )
{
  assert( group_idx < this->root_interface_group_.size() );

  return this->root_interface_group_[group_idx];
}


template<size_t dim>
void  MeshManager<dim>::SetRootInterFace( InterFace<dim>* root_interface )
{
  this->root_interface_group_.push_back( root_interface );
}

template<size_t dim>
const Node<dim>*  MeshManager<dim>::RootNode( size_t group_idx ) const
{
  assert( group_idx < this->root_node_group_.size() );

  return this->root_node_group_[group_idx];
}


template<size_t dim>
const Element<dim>*  MeshManager<dim>::RootElement( size_t group_idx ) const
{
  assert( group_idx < this->root_elmt_group_.size() );

  return this->root_elmt_group_[group_idx];
}


template<size_t dim>
const Face<dim>*  MeshManager<dim>::RootFace( size_t group_idx ) const
{
  assert( group_idx < this->root_face_group_.size() );

  return this->root_face_group_[group_idx];
}


template<size_t dim>
const InterFace<dim>*  MeshManager<dim>::RootInterFace( size_t group_idx ) const
{
  assert( group_idx < this->root_interface_group_.size() );

  return this->root_interface_group_[group_idx];
}


template<size_t dim>
void MeshManager<dim>::Out() const
{
  cout << "\nMeshManager<" << dim << ">::Out: " << endl;

  // traversal of the existing mesh nodes to find all its elements	
  deque<const Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( this, nodes, elmts );

  // traversal of the existing mesh root faces to find all its faces	
  deque<const Face<dim>*> faces;
  exploreFacesFromMesh( this, faces );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<const InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( this, interfaces );

  // nodes
  cout << "\nNODES: " << endl;
  for ( auto n : nodes ) {
    string bound = parseBoundary( n->AtBoundary() );
    cout << "\nNode ID: " << n->Idx() << " ";
    cout << n->Coordinate();
    cout << " Boundary flag: " << bound << endl;
  }

  // elements
  cout << "\nELEMENTS: " << endl;
  for ( auto e : elmts ) {
    string bound = parseBoundary( e->AtBoundary() );
    cout << "\nElement ID: " << e->Idx() << " Boundary flag: " << bound << endl;
    cout << "Member Nodes: " << endl;
    for ( size_t i = 0U; i < e->Nodes(); i++ )
      cout << e->N( i )->Idx() << "\t";
    cout << "\nNeighbor elements: " << endl;
    for ( size_t i = 0U; i < e->Neighbors(); i++ )
      if ( e->Neighbor( i ) != NULL )
        cout << e->Neighbor( i )->Idx() << "\t";
      else
        cout << "NO NEIGHBOR\t";

    cout << endl;
  }

  // faces
  cout << "\nFACES: " << endl;
  for ( auto f : faces ) {
    cout << "\nFace ID: " << f->Idx() << endl;
    cout << "Member Nodes: " << endl;
    for ( size_t i = 0U; i < f->Nodes(); i++ )
      cout << f->N( i )->Idx() << "\t";
    cout << "\nNeighbor faces: " << endl;
    for ( size_t i = 0U; i < f->Neighbors(); i++ )
      if ( f->Neighbor( i ) != NULL )
        cout << f->Neighbor( i )->Idx() << "\t";
      else
        cout << "NO NEIGHBOR\t";

    cout << endl;
  }

  // inter faces
  cout << "\nINTERFACES: " << endl;
  for ( auto f : interfaces ) {
    cout << "\nInterFace ID: " << f->Idx() << endl;
    cout << "Member Nodes: " << endl;
    for ( size_t i = 0U; i < f->Nodes(); i++ )
      cout << f->N( i )->Idx() << "\t";
    cout << "\nNeighbor faces: " << endl;
    for ( size_t i = 0U; i < f->Neighbors(); i++ )
      if ( f->Neighbor( i ) != NULL )
        cout << f->Neighbor( i )->Idx() << "\t";
      else
        cout << "NO NEIGHBOR\t";

    cout << endl;
  }

  // parent elements ID's for each node
  cout << endl << endl;
  cout << "PARENT ELEMENT INFORMATION FOR ALL NODES: " << endl;
  for ( auto n : nodes ) {
    cout << "\nNode: " << n->Idx() << ", parent elements: " << endl;
    for ( size_t i = 0u; i < n->Parents(); i++ )
      cout << n->Parent( i )->Idx() << " ";
    cout << endl;
  }
} // end Out

template class MeshManager<1U>;
template class MeshManager<2U>;
template class MeshManager<3U>;

} // end namespace csmp 
