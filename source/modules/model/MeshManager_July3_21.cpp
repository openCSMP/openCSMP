#include "MeshManager.h"
#include "PropertyDatabase.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "VSet.h"
#include "AP_BoolVector.h"
#include "ErrorHandler.h"
#include "PropertyData.h"
#include "Box.h"
#include "ModelTopology.h"
#include "MeshIterator.h"
#include "IndexToPointerMapping.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp {



// ===================================================================================
//
//                 NON-MEMBER FUNCTIONS
//
// ===================================================================================


/**
     Traverses a mesh of elements, returning the discovered nodes and elements into the supplied deques.
     
        Method uses node-to-parent element connectivity to explore the patch.
        
        TODO: output set rather than deque!
        
        Reasons
        - to build region "Model" OK, Model is then used to build any other region 
        - to delete mesh partches and all mesh -> OK
        - to delete nodes, faces, elements etc. selectively OK
        - to update mesh after insertions, deletions - DO LOCALLY WHEN MESH IS MODIFIED
        - used to number elements, nodes etc. OK
        - used to output mesh to VSet
        - used to input variables to model      
        
        TODO: understand whether there is a difference between between "Model" region and Mesh 
        
*/
template<size_t dim>
void exploreNodesAndElementsFromMesh( const MeshManager<dim>& mesh, std::deque<const Node<dim>*>& nodes, std::deque<const Element<dim>*>& elmts )
{
	// traversal of the existing mesh nodes to find all nodes and elements
	set<const Element<dim>*> explored_elements;
	set<const Node<dim>*>	   discovered_nodes;
	deque<const Node<dim>*>	 current_nodes;
	for (size_t g = 0U; g < mesh.NodeGroups(); g++) {
		auto root_node = mesh.RootNode(g);
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

template void exploreNodesAndElementsFromMesh( const MeshManager<1U>&, std::deque<const Node<1U>*>&, std::deque<const Element<1U>*>& );
template void exploreNodesAndElementsFromMesh( const MeshManager<2U>&, std::deque<const Node<2U>*>&, std::deque<const Element<2U>*>& );
template void exploreNodesAndElementsFromMesh( const MeshManager<3U>&, std::deque<const Node<3U>*>&, std::deque<const Element<3U>*>& );




template<size_t dim>
void exploreFacesFromMesh( const MeshManager<dim>& mesh, std::deque<const Face<dim>*>& faces)
{
	// traversal of the existing mesh root faces to find all faces	
	set<const Face<dim>*>	discovered_faces;
	deque<const Face<dim>*>	current_faces;

	for (size_t g = 0U; g < mesh.FaceGroups(); g++) {
		auto root_face = mesh.RootFace(g);
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

template void exploreFacesFromMesh(const MeshManager<1U>&, std::deque<const Face<1U>*>&);
template void exploreFacesFromMesh(const MeshManager<2U>&, std::deque<const Face<2U>*>&);
template void exploreFacesFromMesh(const MeshManager<3U>&, std::deque<const Face<3U>*>&);



template<size_t dim>
void exploreInterFacesFromMesh( const MeshManager<dim>& mesh, std::deque<const InterFace<dim>*>& interfaces)
{
	// traversal of the existing mesh root interfaces to find all interfaces	
	set<const InterFace<dim>*>		discovered_interfaces;
	deque<const InterFace<dim>*>	current_interfaces;

	for (size_t g = 0U; g < mesh.InterFaceGroups(); g++) {
		auto root_interface = mesh.RootInterFace(g);
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

template void exploreInterFacesFromMesh(const MeshManager<1U>&, std::deque<const InterFace<1U>*>&);
template void exploreInterFacesFromMesh(const MeshManager<2U>&, std::deque<const InterFace<2U>*>&);
template void exploreInterFacesFromMesh(const MeshManager<3U>&, std::deque<const InterFace<3U>*>&);




// ===================================================================================
//
//                 MESH MANAGER METHODS 
//
// ===================================================================================





// refactored
template<size_t dim>
MeshManager<dim>::MeshManager()
  : hybrid_element_mesh_( false ), n_nodes_( 0 ), n_elmts_( 0 ), n_faces_( 0 ), n_interfaces_( 0 )
{
}


/**
       
*/
template<size_t dim>
MeshManager<dim>::MeshManager( const PropertyDatabase<dim>& pref, const FiniteElementManager& fem_manager,
                               const VSet<dim>& vset, IndexToPointerMapping<dim>& mapping )
  : hybrid_element_mesh_( vset.HybridElementTypeMesh() )
{
   Initialize( pref, fem_manager, vset, mapping );
}


/**
     Deallocate all the dynamically allocated nodes, elements, faces, interfaces.
     
     @note Only the nodes of the volumetric elements are deleted because they are shared  with the surface mesh.
     @note Any pointer used to delete an object is set to null after the deletion.
     
     Sequence of operations
     
     1. Delete the nodes
     2. For each of the mesh trees, do a deletion traversal
     
 */
template<size_t dim>
MeshManager<dim>::~MeshManager()
  {
    n_nodes_ = 0;
    n_elmts_ = 0;
    n_faces_ = 0;
    n_interfaces_ = 0;
    
 } // end destructor


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

  hybrid_element_mesh_   = mmgr.hybrid_element_mesh_;
  n_nodes_               = mmgr.n_nodes_;
  n_elmts_               = mmgr.n_elmts_;
  n_faces_               = mmgr.n_faces_;
  n_interfaces_          = mmgr.n_interfaces_;
  elmt_tree_roots_       = mmgr.elmt_tree_roots_;
  face_tree_roots_       = mmgr.face_tree_roots_;
  itfc_tree_roots_       = mmgr.itfc_tree_roots_;
  node_manifold_manager_ = mmgr.node_manifold_manager_;

  return *this;

} // end operator=


/**
Configures the distributed variable storage, inialises the finite element container
with the element types that are contained in the current mesh and builds the
mesh including the connectivity among its elements.

While the method detects insconsistencies / omissions in the element/face/interface neighbor arrays it will continue execution assigning nullptr where incorrect neighbors are encountered.
However it will report when and where this happened.

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
                                    const VSet<dim>& vset,
                                    IndexToPointerMapping<dim>& indexToPtrMapping )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  assert( n_nodes_ == 0 );
  assert( n_elmts_ == 0 );
  assert( n_faces_ == 0 );
  assert( n_interfaces_ == 0 );

  // temporary deques to handle primitives
  indexToPtrMapping.nodeConnector_.reserve( vset.Vertices() );
  indexToPtrMapping.elementConnector_.reserve( vset.Elements() );

  hybrid_element_mesh_ = vset.HybridElementTypeMesh();
  cout <<"\nMeshManager<"<< dim <<">::Initialize: building mesh with "<< vset.Elements() <<" elements, "<< vset.Vertices() <<" nodes, ";
  cout << vset.Faces() <<" faces, and "<< vset.InterFaces() <<" interfaces.\n";
  if ( vset.HybridElementTypeMesh() ) cout <<"mesh consists of multiple element types.\n";
  if ( vset.Faces() > 0 ) cout <<"mesh contains 'Boundary' objects.\n";
  if ( vset.InterFaces() > 0 ) cout <<"mesh contains 'SplitBoundary' objects.\n";
  cout << endl;

  // ------------------------------------------------------------------------------------
  // 1. check availability of necessary finite element types, and valid model topology
  // ------------------------------------------------------------------------------------

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

  // --------------------------------------------------------------------------
  // 2. construct nodes and elements using the VSet element type information
  // --------------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: building storage and assigning nodes to elements..." << endl;

  // Storage for nodes
   {
      vector<double64> coord( dim );
      const LocalVariables nvars( phys_vars.LocalVariablesAt( NODE ) );
      for ( size_t idx = 0U; idx < vset.Vertices(); ++idx ) {
          for ( size_t j = 0U; j<dim; ++j ) coord[j] = vset.P( j, idx );
          indexToPtrMapping.nodeConnector_.push_back( new Node<dim>( idx, Point<dim>( coord ), nvars, NOT ) );
        }
      n_nodes_ = indexToPtrMapping.nodeConnector_.size();
    }
     
  // storage for elements
   {
      const LocalVariables evars( phys_vars.LocalVariablesAt( ELEMENT ) );
      const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( ELEMENT ) );
      typename deque<vector<size_t>>::const_iterator first( vset.PlistElmtsBegin() ), last( vset.PlistElmtsEnd() );

      size_t elmt_idx(0U);
      // 2.1 If the MeshManager contains only one element type
      if ( !vset.HybridElementTypeMesh() ) {
          const int32 csmpElementType = vset.ElementType( 0U );
          while ( first != last )
            {
              indexToPtrMapping.elementConnector_.push_back( new Element<dim>( elmt_idx, fem_manager.E( csmpElementType ), evars, cvars, NOT ) );
              const size_t nodes( fem_manager.E( csmpElementType )->Nodes() );
              for ( size_t j = 0U; j < nodes; ++j )
                indexToPtrMapping.elementConnector_[elmt_idx]->Assign( j, indexToPtrMapping.nodeConnector_[ vset.Plist( elmt_idx, j )] );
              elmt_idx++;
              first++;
            }
        }
      // 2.2 If there are multiple element types
      else {
          while ( first != last )
            {
              const int32 csmpElementType = vset.ElementType( elmt_idx );
              indexToPtrMapping.elementConnector_.push_back( new Element<dim>( elmt_idx, fem_manager.E( csmpElementType ), evars, cvars, NOT ) );
              const size_t nodes( fem_manager.E( csmpElementType )->Nodes() );
              for ( size_t j = 0U; j < nodes; j++ )
                indexToPtrMapping.elementConnector_[elmt_idx]->Assign( j, indexToPtrMapping.nodeConnector_[vset.Plist( elmt_idx, j )] );
              elmt_idx++;
              first++;
            }
        }
      assert( indexToPtrMapping.elementConnector_.size() == vset.Elements() );
      n_elmts_ = indexToPtrMapping.elementConnector_.size();
   }
  
  // 2.3 Assign neighbor elements to elements
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: assigning neighbors to elements..." << endl;

  if ( vset.WithNeighbourConnectivity() ) {
       for ( auto& e : indexToPtrMapping.elementConnector_ ) {
            const int32 csmpElementType = (!hybrid_element_mesh_) ? vset.ElementType( 0U ) : vset.ElementType( e->Idx() );
            const size_t neighbors( fem_manager.E( csmpElementType )->Neighbors() );

            for ( size_t j = 0U, nidx = 0U; j < neighbors; ++j ) {
                  // if there is a neighbor (as is the case if the stored index is greater than zero)
                  if ( j < vset.PfvertsSize( e->Idx() ) ) {
                      const int32 index( static_cast<int32>(vset.Pfvert( e->Idx(), j )) );
                      if ( index >= n_elmts_ ) {
                           cerr <<"\n\t"<< index <<" vs. number of elements = "<< n_elmts_ << endl;
                           csmp_error.notice( ERROR, "MeshManager::Initialise: ", "element ID in 'pfverts' out of range.");
                        }
                      else if ( index >= 0 )
                        e->Assign( nidx++, indexToPtrMapping.elementConnector_[static_cast<size_t>(index)] );
                      else
                        e->Assign( nidx++, static_cast<Element<dim>*>(nullptr) );
                   }
              }
         }
    }
  else
  csmp_error.notice( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity; nothing was done." );

// TODO: find root Elements and insert them into deque/set or whatever mesh patch container


  // ----------------------------------------------------------
  // 3. constructing the Faces using the VSet node information
  // ----------------------------------------------------------
  // (continuous running index 'idx' will be used so that face-IDs start with n-elements)
  // different face types are intrinsic to the VSet if so initialized
  if ( vset.Faces() > 0 )
    {
       assert( vset.HybridElementTypeMesh() );
       indexToPtrMapping.faceConnector_.reserve( vset.Faces() );
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to faces..." << endl;
      
       // since faces are lower-dimensional, the mesh must contain different element types
       const LocalVariables evars( phys_vars.LocalVariablesAt( FACE ) );
       const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( FACE ) );

       // the faces are numbered  elements to (elements + faces - 1), but they are stored in connector at Face 0..n-1
       size_t  face_idx(vset.Elements());
       typename deque<vector<size_t> >::const_iterator  first( vset.PlistFacesBegin() ), last( vset.PlistFacesEnd() );
       while ( first != last ) {
            const int32 csmpElementType = vset.ElementType( face_idx );
            if ( csmpElementType == UNKNOWN ) {
                 cerr <<"\n\t"<< parseFiniteElementType(csmpElementType) <<" encountered for Face "<< face_idx <<"\n";
                 csmp_error.notice( FATAL_ERROR, "MeshManager::Initialise:", "encountered UNKNOWN Face element type." );
              }
            indexToPtrMapping.faceConnector_.push_back( new Face<dim>( face_idx, fem_manager.E( csmpElementType ), evars, cvars ) );
            const size_t nodes( indexToPtrMapping.faceConnector_[face_idx]->Nodes() );
            // assigning nodes to faces
            for ( size_t j = 0U; j<nodes; ++j ) {
                const size_t node = vset.Plist( face_idx, j );
                assert( node < n_nodes_ );
                indexToPtrMapping.faceConnector_[face_idx]->Assign( j, indexToPtrMapping.nodeConnector_[node] );
              }
            ++face_idx;
            ++first;
          }
       assert( indexToPtrMapping.faceConnector_.size() == vset.Faces() );
       n_faces_ = indexToPtrMapping.faceConnector_.size();

       // connecting the faces to their equi- and higher-dimensional neighbors
       // --------------------------------------------------------------------
       // necessary info is stored in 'pfverts' record for Face
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: connecting faces to their higher-dimensional neighbors..." << endl;
       if ( !vset.WithNeighbourConnectivity() )
         csmp_error.notice( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity for Face objects; nothing was done." );
       else
         {
            for ( auto& e : indexToPtrMapping.faceConnector_ )
              {
                 // equidimensional neighbor Faces first
                 // ------------------------------------
                 const size_t neighbors( e->Neighbors() );
                 for ( size_t j = 0U; j<neighbors; ++j )
                   {
                      // if there is a neighbor (as is the case if the stored index is greater than zero)
                      // (e->Idx() starts with elements=first face)
                      const long64 index( vset.Pfvert( e->Idx(), j ) );
                      if ( index >= n_elmts_+n_faces_ ) {
                           cerr <<"\n\t"<< index <<" vs. number of elements+faces = "<< n_elmts_+n_faces_ << endl;
                           csmp_error.notice( ERROR, "MeshManager::Initialise: ", "face ID in 'pfverts' out of range.");
                        }
                      assert( index >= vset.Elements() );
                      assert( index < vset.Elements() + vset.Faces() ); // (-) elements because face container is numbered from 0..n-1
                      if ( index >= n_elmts_ )
                        e->Assign( j, indexToPtrMapping.faceConnector_[ static_cast<size_t>(index) - n_elmts_] );
                      else
                        e->Assign( j, static_cast<Face<dim>*>(nullptr) );
                   }
              
                 // Two higher-dimensional Element neighbors of Faces
                 // -------------------------------------------------
                 // (are stored in VSet 'pfverts' record after the equidimensional neighbors)
                 // index of inner neighbor element i which is always there
                 // Inside neighbor 1
                 const long64 index1( vset.Pfvert( e->Idx(), neighbors ) );
                 if ( index1 >= n_elmts_ ) {
                      cerr <<"\n\tFace "<< e->Idx() <<": "<< index1 <<" vs. "<< n_elmts_ <<" elements.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of first higher-dimensional element of Face out of range.");
                   }
                 else if ( index1 < 0 ) {
                      cerr <<"\n\tFace "<< e->Idx() <<": "<< index1 <<"\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "First higher-dimensional element out of range.");
                   }

                 // index of outer neighbor element which will be there only if the Face is located within the model
                 const long64 index2( vset.Pfvert( e->Idx(), neighbors + 1U ) );
                 if ( index2 >= n_elmts_ ) {
                      cerr <<"\n\tFace "<< e->Idx() <<": "<< index2 <<" vs. "<< n_elmts_ <<" elements.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of second higher-dimensional element of Face out of range.");
                   }
                 Element<dim>* const innerElement = (index1 < 0) ? nullptr : indexToPtrMapping.elementConnector_[index1];
                 Element<dim>* const outerElement = (index2 < 0) ? nullptr : indexToPtrMapping.elementConnector_[index2];
                 // assigning inner and outer higher-dimensional neighbors
                 indexToPtrMapping.faceConnector_[index1 - n_elmts_]->Assign( innerElement, outerElement );
                 
               } // end face loop
              
             } // end else

// TODO: find root Faces and insert them into deque/set or whatever mesh patch container
        
    } // end construction and initialisation of Face objects


  // ------------------------------------------------------------------
  // 4. constructing Interfaces objects using the VSet node information
  // ------------------------------------------------------------------
  // (continuous running index 'idx' will also be used for interfaces)
  if ( vset.InterFaces() > 0 )
    {
       assert( vset.HybridElementTypeMesh() );
       indexToPtrMapping.interFaceConnector_.reserve( vset.InterFaces() );
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to interfaces..." << endl;
       
       const LocalVariables evars( phys_vars.LocalVariablesAt( INTER_FACE ) );
       const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( INTER_FACE ) );

       typename deque<vector<size_t> >::const_iterator  first( vset.PlistInterFacesBegin() ),
                                                        last( vset.PlistInterFacesEnd() );

       size_t interface_idx(vset.Elements() + vset.Faces());
       while ( first != last )
         {
            const int32 csmpElementType = vset.ElementType( interface_idx );
            indexToPtrMapping.interFaceConnector_.push_back( new InterFace<dim>( interface_idx, fem_manager.E( csmpElementType ), evars, cvars ) );
            const size_t nodes( indexToPtrMapping.interFaceConnector_[interface_idx]->Nodes() );
            // assigning nodes
            for ( size_t j = 0U; j<nodes; ++j ) {
                 const size_t node(vset.Plist( interface_idx, j ));
                 if ( node >= n_nodes_ ) {
                      cerr <<"\n\tInterFace "<< interface_idx <<": node j "<< node <<" vs. "<< n_nodes_ <<" nodes.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of InterFace node out of range.");
                   }
                 indexToPtrMapping.interFaceConnector_[interface_idx]->Assign( j, indexToPtrMapping.nodeConnector_[node], INSIDE );
              }
            ++interface_idx;
            ++first;
          }
       assert( indexToPtrMapping.interFaceConnector_.size() == vset.InterFaces() );
       n_interfaces_ = indexToPtrMapping.interFaceConnector_.size();

       // connecting the faces to their equi- and higher-dimensional neighbors
       // --------------------------------------------------------------------
       // necessary info is stored in 'pfverts' record for Face
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: connecting interfaces to their higher-dimensional neighbors..." << endl;
       const size_t cells(n_elmts_+n_faces_+n_interfaces_);
       // connecting interfaces to their higher-dimensional neighbors
       for ( auto& e : indexToPtrMapping.interFaceConnector_ )
         {
            // equidimensional neighbors first
            const size_t neighbors( e->Neighbors() );
            for ( size_t j = 0U; j<neighbors; ++j )
              {
                 // if there is a neighbor (as is the case if the stored index is greater than zero)
                 const long64 index( vset.Pfvert( e->Idx(), j ) );
                 if ( index >= cells ) {
                      cerr <<"\n\t"<< index <<" vs. number of elements+faces+interfaces = "<< cells << endl;
                      csmp_error.notice( ERROR, "MeshManager::Initialise: ", "interface ID in 'pfverts' out of range.");
                   }
                 assert( index >= vset.Elements() + vset.Faces() );
                 assert( index < vset.Elements() + vset.Faces() + vset.InterFaces() ); // (-) elements because face container is numbered from 0..n-1
                
                 // the number of the interface in the container is the number from the VSet - elements and faces
                 // because the interface container is counts from 0..n-1
                 e->Assign( j, indexToPtrMapping.interFaceConnector_[ static_cast<size_t>(index) - n_elmts_ - n_faces_ ] );
              }

           // higher-dimensional neighbors
           // ----------------------------
           // (the higher-dimensional neighbors on both sides must always be present because interfaces can exist only inside of a model)
           if ( vset.Pfvert( e->Idx(), neighbors ) < 0 || vset.Pfvert( e->Idx(), neighbors + 1U) < 0 )
             {
                cerr <<"\n\tInterFace "<< e->Idx() <<": inner neighbor "<< vset.Pfvert( e->Idx(), neighbors ) <<" and outer "<< vset.Pfvert( e->Idx(), neighbors+1U ) <<"\n";
                csmp_error.notice( ERROR, "MeshManager::Initialise: ", "Higher dimensional neighbor of InterFace not defined in 'pfverts'.");
             }
           if ( vset.Pfvert( e->Idx(), neighbors ) >= cells || vset.Pfvert( e->Idx(), neighbors + 1U ) >= cells )
             {
                cerr <<"\n\tInterFace "<< e->Idx() <<": inner neighbor "<< vset.Pfvert( e->Idx(), neighbors ) <<" and outer "<< vset.Pfvert( e->Idx(), neighbors+1U ) <<"\n";
                csmp_error.notice( ERROR, "MeshManager::Initialise: ", "Higher dimensional neighbor of InterFace out of range.");
             }
           const long64 index1 = vset.Pfvert( e->Idx(), neighbors );
           const long64 index2 = vset.Pfvert( e->Idx(), neighbors+1U );
           Element<dim>* const innerElement = (index1 < 0) ? nullptr : indexToPtrMapping.elementConnector_[index1];
           Element<dim>* const outerElement = (index2 < 0) ? nullptr : indexToPtrMapping.elementConnector_[index2];
           // assigning higher dimensional neighbors or null pointers
           e->Assign( innerElement, outerElement );
        }

// TODO: find root InterFaces and insert them into deque/set or whatever mesh patch container
// TODO: create / restore manifolds where collocated nodes were encountered (perhaps these need to be stored in VSet).
        
    } // end construction of Interface objects


  // ---------------------------------------------------------------------
  // 5. Flagging nodes at model boundary with BOX_BOUNDARY flags
  // ---------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: flagging boundary objects..." << endl;
  if ( vset.BFlags() > 0 ) {
       // nodes were initially constructed as not located at the model boundary
       // TODO: why are the bflags long64?!
       for ( unordered_map<size_t,long64>::const_iterator bit = vset.BFlagsBegin(); bit != vset.BFlagsEnd(); bit++ )
         indexToPtrMapping.nodeConnector_[(*bit).first]->AtBoundary( intToBOX_BOUNDARY( (*bit).second ) );
    }
    
    
  // ------------------------------------------------------------------------------
  // 6. Assigning parent elements (these are the elements that share the node) and
  // their respective internal node-id numbers to the nodes
  // -------------------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: assigning parent element information to nodes..." << endl;
  vector<size_t>  parent_elmts_per_node( vset.Vertices(), 0U );

  // counting how many parent elements each node has
  for ( auto& e : indexToPtrMapping.elementConnector_ ) {
      assert( e != nullptr );
      for ( auto nit = e->NodesBegin(); nit != e->NodesEnd(); nit++ ) {
           assert( (*nit) != nullptr );
           assert( (*nit)->Idx() < n_nodes_ );
           parent_elmts_per_node[ (*nit)->Idx() ]++;
        }
    }

  // reserving the memory for the parent storage and zeroing parent vector for next step
  size_t i( 0 );
  for ( auto& n : indexToPtrMapping.nodeConnector_ )
    n->ResizeParentStorage( parent_elmts_per_node[i++] );

   // assigning the parent element information to the nodes
   for ( auto& e : indexToPtrMapping.elementConnector_ )
     for ( size_t j = 0U; j<e->Nodes(); ++j )
       e->N( j )->Assign( j, e );

   if ( csmp_error.Verbose() )
     cout << "\nMeshManager<" << dim << ">::Initialize: forming regions for contiguous subdomains..." << endl;


   // ------------------------------------------------------------------------------
   // 7. forming root-node and root-element pointers for contiguous regions
   // ------------------------------------------------------------------------------
   size_t elmt_patches = initialisePointersToStandAloneMeshPatches( indexToPtrMapping.elementConnector_.begin(),
                                                                    indexToPtrMapping.elementConnector_.end(),
                                                                    elmt_tree_roots_ );
   cout <<"\nMeshManager<"<< dim <<">::Intialise: ";
   cout <<" created "<< elmt_patches <<" contiguous mesh trees.\n";

   if ( !indexToPtrMapping.faceConnector_.empty() ) {
   size_t face_patches = initialisePointersToStandAloneMeshPatches( indexToPtrMapping.faceConnector_.begin(),
                                                                    indexToPtrMapping.faceConnector_.end(),
                                                                    face_tree_roots_ );
        cout <<"\nMeshManager<"<< dim <<">::Intialise: ";
        cout <<" mesh contains Face objects (Boundaries), ";
        cout <<" created "<< face_patches <<" contiguous boundary mesh trees.\n";
     }
     
   if ( !indexToPtrMapping.interFaceConnector_.empty() ) {
   size_t itfc_patches = initialisePointersToStandAloneMeshPatches( indexToPtrMapping.interFaceConnector_.begin(),
                                                                    indexToPtrMapping.interFaceConnector_.end(),
                                                                    itfc_tree_roots_ );
        cout <<"\nMeshManager<"<< dim <<">::Intialise: ";
        cout <<" mesh contains InterFace objects (SplitBoundaries), ";
        cout <<" created "<< itfc_patches <<" contiguous split-boundary mesh trees.\n";
     }

   return true;
  
} // end Initialise





  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================
  /*
      Think about use cases:
      - what are the individual steps
      - who in the hierarchy of nodes, elements is reponsible for what (Element is master?)
      - Do we need reference counting for nodes?


    // To erase elements the following steps are necessary
    // ---------------------------------------------------
    // 1. set the neighbor pointers to the element to zero
    // 2. remove the pointers from the connected nodes to this parent-element
    // 3. disconnect the nodes, deleting them when they have no parent elements anymore
    // 4. disconnect the neighbor elements
    // 5. delete element
 
  */








/**
   Rebuilds the node to element parent relationship. JCK 2018.
*/
template<size_t dim>
void MeshManager<dim>::RebuildParentRelationships( typename vector<Node<dim>*>::iterator begin,
                                                   typename vector<Node<dim>*>::iterator end )
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
  size_t boundary_only_elements( 0U );
  for ( size_t g=0U; g<NodeGroups(); ++g ) 
    {
		   const MeshIterator<dim,Element>  mit( const_cast<Node<dim>*>(RootNode(g)) );
       for ( typename set<Element<dim>*>::const_iterator 
             it=mit.CellsBegin(); it!=mit.CellsEnd(); ++it ) {
            const size_t nodes((*it)->Nodes());
            size_t       counter(0U);
            for ( size_t i=0U; i<nodes; ++i )
              if ( (*it)->N(i)->AtBoundary() != NOT ) counter++;
            if ( counter == nodes ) {
                  belmts.insert( (*it)->Idx() );
                  boundary_only_elements++;
              } 
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
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );

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
   Duplicates existing node inside of the MeshManager
*/
template<size_t dim>
Node<dim>* MeshManager<dim>::Duplicate( const Node<dim>* const node )
  {
    Node<dim>* new_node = new Node<dim>( node );
    n_nodes_++;

    return new_node;
  }


/**
   Duplicates the corresponding element
*/
template<size_t dim>
Element<dim>* MeshManager<dim>::Duplicate( const Element<dim>* const elmt )
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
    if ( n != nullptr ) {
      n->ResizeParentStorage( n->Parents() + 1 );
      n->Assign( n->Parents() - 1, new_elmt );
    }
  neighbor_nodes.clear();

  return new_elmt;
}


/**
Inserts the corresponding node
*/
template<size_t dim>
Node<dim>* MeshManager<dim>::AddNodeAt( const Point<dim>& location )
{
  Node<dim>* new_node = new Node<dim>( location );
  n_nodes_++;

  return new_node;
}








/**
Deletes the corresponding node
*/
template<size_t dim>
void MeshManager<dim>::Erase( Node<dim>* node )
  {
     assert( node != nullptr );
     assert( node->Parents() == 0U );
     
     if ( node == nullptr ) return;
     if ( node->Parents() > 0U ) return;
     
     // delete the node
     delete node;
     node = nullptr;
     n_nodes_--;
  }




/**
    Deletes the  element pointed to by the pointer:

    // To erase elements the following steps are necessary
    // ---------------------------------------------------
    // 1. set the neighbor pointers to the element to zero
    // 2. remove the pointers from the connected nodes to this parent-element
    // 3. disconnect the nodes, deleting them when they have no parent elements anymore
         // 3.2 if a node gets deleted and is part of a Manifold the latter must be modified
    // 4. disconnect the neighbor elements
    // 5. delete element
    
    @author SKM
    @date 22/5/2021

*/
template<size_t dim>
void MeshManager<dim>::Erase( Element<dim>* elmt )
  {
     // has the element already been deleted?
     assert( elmt != nullptr );
     if ( elmt == nullptr ) return;
     
     // 0. if the element is a root element, a new root element must be set
     // TODO: this code is problematic because the root element searches are costly
     const typename map<Element<dim>*,MeshPatchAttributes>::iterator eit=elmt_tree_roots_.find( elmt );
     // if it is a root element
     if ( eit != elmt_tree_roots_.end() ) {
          // the root element is substituted with its first non-null neighbor
          for ( size_t i=0U; i<elmt->Neighbors(); ++i )
            if ( elmt->Neighbor(i) != nullptr ) {
                (*eit).first = elmt->Neighbor(i);
                break;
             }
          // changing entry in map
          elmt_tree_roots_.insert( make_pair( (*eit), (*eit).second ) );
          elmt_tree_roots_.Erase( (*eit) );
       }
     
     // 1. set the neighbor pointers to the element to zero
     for ( size_t i=0U; i<elmt->Neighbors(); ++i )
       if ( elmt->Neighbor(i) != nullptr )
         for ( size_t j=0U; j<elmt->Neighbor(i)->Neighbors(); ++j )
           if ( elmt->Neighbor(i)->Neighbor(j) != nullptr )
            if ( elmt->Neighbor(i)->Neighbor(j) == elmt->Neighbor(i) )
              elmt->Neighbor(i)->Neighbor(j) = nullptr;

     // 2. remove the pointers from the connected nodes to this parent-element
     for ( size_t i=0U; elmt->Nodes(); ++i ) {
          assert( elmt->N(i) != nullptr );
          if ( elmt->N(i) != nullptr )
            elmt->N(i)->Unassign( elmt );
       }
     
     // 3. disconnect the nodes, deleting them when they have no parent elements anymore
     for ( size_t i=0U; elmt->Nodes(); ++i )
       if ( elmt->N(i) != nullptr ) {
            if ( elmt->N(i)->ParentManifold() != nullptr ) {
                 // removing itself from the manifold
                 elmt->N(i)->ParentManifold().Delete( elmt->N(i) );
                 // disconnecting the node from the manifold
                 elmt->N(i)->AssignParentManifold( nullptr );
              }
            // deleting the node
            if ( elmt->N(i)->Parents() == 0U ) {
                 Erase( elmt->N(i) );
                 elmt->N(i) = nullptr;
              }
         }
     
     // 4. disconnect the neighbor elements
     for ( size_t i=0U; i<elmt->Neighbors(); ++i )
       elmt->Unassign( elmt->Neighbor(i) );
     
     // 5. delete element
     // -----------------
     delete elmt;
     elmt = nullptr;

} // end Erase(element-ptr)








/**
Deletes the corresponding face
*/
template<size_t dim>
void MeshManager<dim>::Erase( Face<dim>* face )
{
  // if the found face is the root face
  // change the root face into one of its neighbors and delete the found face.
  Face<dim>* new_root_face( nullptr );
  Face<dim>* root_face( nullptr );
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
      if ( new_root_face == nullptr )
        cerr << "MeshManager::Erase(Face): this is the root face which cannot be deleted. \n";
      else
        root_face = new_root_face;
    }
  }

  // delete the face
  if ( face ) {
    // unassign the connections of its neighbors
    for ( auto n : face->NeighborElementVector() )
      if ( n != nullptr ) n->Unassign( face );

    delete face;
    face = nullptr;
    n_faces_--;
  }
  
} // end Erase( Face pointer )







/**
Deletes the corresponding interface
*/
template<size_t dim>
void MeshManager<dim>::Erase( InterFace<dim>* interface )
{
  // if the found interface is the root interface
  // change the root interface into one of its neighbors and delete the found interface.
  InterFace<dim>* new_root_interface( nullptr );
  InterFace<dim>* root_interface( nullptr );
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
      if ( new_root_interface == nullptr )
        cerr << "MeshManager::Erase(InterFace): this is the root interface which cannot be deleted. \n";
      else
        root_interface = new_root_interface;
    }
  }

  // delete the interface
  if ( interface ) {
    // unassign the connections of its neighbors
    for ( auto n : interface->NeighborElementVector() )
      if ( n != NULL ) n->Unassign( interface );

    delete interface;
    interface = nullptr;
    n_interfaces_--;
  }
}





/**
    erases the supplied sequence of elements returning the number of erasures, the pointers to the erased elements are nulled. The connectivity of the affected mesh neighborhood will get fixed.
*/
template<size_t dim>
size_t MeshManager<dim>::EraseElements( typename vector<Element<dim>*>::iterator first_elmt,
                                        typename vector<Element<dim>*>::iterator last_elmt )
 {
    size_t erased_elements(0U);
    
    // To erase elements the following steps are necessary
    // ---------------------------------------------------
    // 1. set the neighbor pointers to the element to zero
    // 2. remove the pointers from the connected nodes to this parent-element
    // 3. disconnect the nodes, deleting them when they have no parent elements anymore
    // 4. disconnect the neighbor elements
    // 5. delete element
 
    return erased_elements;
    
 } // end EraseElements





/**
erases all elements, discerning those that have any connections.
TODO: traverse and erase at the same time
*/
template<size_t dim>
bool MeshManager<dim>::EraseElements()
{
  bool emptyElements( false );

  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );

  for ( auto e : elmts ) {
    if ( e ) {
      delete e;
      e = nullptr;
      n_elmts_--;
    }
  }

  for ( auto n : nodes ) {
    if ( n ) {
      delete n;
      n = nullptr;
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
  exploreFacesFromMesh( *this, faces );

  // 2. assigning and trimming excess storage from the face pointer vector
  for ( auto f : faces ) {
    if ( f ) {
      delete f;
      f = nullptr;
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
  exploreInterFacesFromMesh( *this, interfaces );

  // 2. assigning and trimming excess storage from the interface pointer vector
  for ( auto f : interfaces ) {
    if ( f ) {
      delete f;
      f = nullptr;
      n_interfaces_--;
    }
  }

  if ( n_interfaces_ == 0 ) {
    root_interface_group_.clear();
    emptyInterFaces = true;
  }

  return emptyInterFaces;
}





/**
      Converts supplied range of faces into interfaces. Prior to InterFace creation, perimeter Faces of the Face region are identified, receiving a special treatment that
      avoids that the nodes along the tiplines are duplicated, see source.
      
      Supplied Face objects are replaced by InterFace objects (Face objects are deleted); returns pointer to first perimeter interface.
  Perimeter nodes are not duplicated, see CSMP User' s guide/

 @author SKM modified original function
 @date 16/05/2021
 
Duplicates nodes as necessary. 
@note Originally this was the Split() function from JC 1/7/2019

TODO: needs complete rewrite

*/
/*
template<size_t dim>
typename vector<InterFace<dim>*>::iterator  MeshManager<dim>::ConvertFacesToInterFaces( typename vector<Face<dim>*>::iterator faces_begin,
                                                                                        typename vector<Face<dim>*>::iterator faces_end,
                                                                                        typename vector<Node<dim>*>::iterator interior_nodes_begin,
                                                                                        typename vector<Node<dim>*>::iterator interior_nodes_end,
                                                                                        typename vector<InterFace<dim>*>::iterator interfaces_begin,
                                                                                        typename vector<InterFace<dim>*>::iterator interfaces_end )
{
  // generating InterFace objects
  
  std::set<Node<dim>*> outsideElementNodes, insideElementNodes;
  std::map<Node<dim>*, Node<dim>*> manyfoldNodes;
  
  // finding higher-dimensional neighbors on inside
/*
  for ( typename vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != this->ElementsEnd(); ++ifit )
    {
        Element<dim>* eit = (*ifit)->OuterParent();
        for ( size_t en( 0 ); en < eit->Nodes(); ++en )
        {
          bool found( false );
          for ( size_t ifn( 0 ); ifn < (*ifit)->Nodes(); ++ifn )
          {
            if ( (*ifit)->N( ifn )->Idx() == eit->N( en )->Idx() )
              found = true;
          }
          if( found )
            outsideElementNodes.insert( eit->N( en ) );
        }
    }

  // finding higher-dimensional neighbors on outside
  for ( typename std::vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != this->ElementsEnd(); ++ifit )
    {
        Element<dim>* eit = (*ifit)->InnerParent();
        for ( size_t en( 0 ); en < eit->Nodes(); ++en )
        {
          bool found( false );
          for ( size_t ifn( 0 ); ifn < (*ifit)->Nodes(); ++ifn )
          {
            if ( (*ifit)->N( ifn )->Idx() == eit->N( en )->Idx() )
              found = true;
          }
          if ( found )
          insideElementNodes.insert( eit->N( en ) );
        }
    }
*/
  // TODO: here the new nodes should be added to the mesh tree
  /*
  for ( auto oen : outsideElementNodes ) {
    Node<dim> new_node( *oen );
    Node<dim>* duplicatedNode = model.Mesh().Add( new_node );
    duplicatedNode->Idx( model.Mesh().Nodes() );
    manyfoldNodes.insert( make_pair( oen, duplicatedNode ) );
  }
  */
  
/* dealing with the Node manifolds

  Region<dim>* outer_region = nullptr;
  for ( typename std::vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != this->ElementsEnd(); ++ifit )
  {
    Element<dim>* oeit = (*ifit)->OuterParent();
    for ( typename std::map<std::string, csmp::Region<dim> >::iterator
          it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); ++it ) {
      if ( (*it).second.Contains( oeit ) ) {
        outer_region = &(*it).second;
        break;
      }
    }
  }
  for ( typename std::vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != this->ElementsEnd(); ++ifit )
  {
    Element<dim>* oeit = (*ifit)->OuterParent();
    for ( size_t en( 0 ); en < oeit->Nodes(); ++en ) {
      if ( manyfoldNodes.find( oeit->N( en ) ) != manyfoldNodes.end() ) {
        oeit->Assign( en, manyfoldNodes[oeit->N( en )] );
      }
    }

    if ( model.UniqueRegions() > 1 ) {
      Region<dim>&  mref( model.Region( "Model" ) );
      for ( typename vector<Element<dim>*>::iterator eit( mref.ElementsBegin() ); eit != mref.ElementsEnd(); ++eit )
        for ( size_t en( 0 ); en < (*eit)->Nodes(); ++en )
          if ( outer_region->Contains( (*eit) ) )
            if ( manyfoldNodes.find( (*eit)->N( en ) ) != manyfoldNodes.end() )
              (*eit)->Assign( en, manyfoldNodes[(*eit)->N( en )] );
    }
  }

} // end ConvertFacesToInterFaces

*/








/** (Re)number all cells; either continuous for all cells or seperate ranges for all entity types

@note this member function is constant because the idx_ is a mutable variable in the cell classes

TODO: replace this with an efficient algorithm
*/
template<size_t dim>
void MeshManager<dim>::AssignUniqueNumbers( bool in_a_single_sequence )
{
  // traversal of the existing mesh nodes to find all its elements	
  deque<Node<dim>*> nodes;
  deque<Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( *this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( *this, interfaces );
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

TODO: recode this method or MOVE METHOD TO MODEL - so that it can use the existing regions, boundaries and split boundaries

*/
template<size_t dim>
void MeshManager<dim>::OutputMeshTo( VSet<dim>& vset,
                                     deque<const Node<dim>*>& nodes, deque<const Element<dim>*>& elmts, 
                                     deque<const Face<dim>*>& faces, deque<const InterFace<dim>*>& interfaces ) const
{
  ErrorHandler& csmp_error( ErrorHandler::Instance() );

  // 1. resizing the VSet
  // --------------------
  const size_t higherDimParents( 2U );
  const size_t interfaceMultiplier( 2U );

  // traversal of the existing mesh nodes to find all its elements	
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
    
  // traversal of the existing mesh root faces to find all its faces	
  exploreFacesFromMesh( *this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  exploreInterFacesFromMesh( *this, interfaces );
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
          csmp_error.notice( WARNING, "MeshManager<dim>::OutputMeshTo (face neighbors):",
                            "inner dim+1 neighbor element of Face should be flagged as model boundary because Face has no outer element; flagging element as irregular" );
          cerr <<"\nDiagnostics:";
          f->InnerParent()->Out();
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
established, for instance with AssignUniqueNumbers() before this method is called.

@note region properties are stored together with the regions in respective binary files

@note this method is anything, but nice. Yet it will be quite a challenge to come up with a better
design; hopefully there will be no more extra variable placements or types in the future.

@attention in the VSet, the variables are identified only by their (unique) names. The property
database is therefore essential to retrieve all other variable related information.

@author SKM 5/5/2016

*/
template<size_t dim>
void MeshManager<dim>::OutputStoredVariablesTo( const PropertyDatabase<dim>& database, 
                                                const deque<const Node<dim>*>& nodes, 
                                                const deque<const Element<dim>*>& elmts, 
                                                const deque<const Face<dim>*>& faces, 
                                                const deque<const InterFace<dim>*>& interfaces, 
                                                VSet<dim>& vset ) const
{
  ErrorHandler&		csmp_error( ErrorHandler::Instance() );
  map<string, Index>  properties;

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
          ArrayVariable value( (*pit).second.dataDepth );
          for ( auto it : nodes ) {
                (*it).Read( (*pit).second, value );
                pushBack( data, value );
              }
            }
          break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
        ArrayVariable value( (*pit).second.dataDepth );
        for ( auto it : elmts ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
        ArrayVariable value( (*pit).second.dataDepth );
        for ( auto it : faces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
        ArrayVariable value( (*pit).second.dataDepth );
        for ( auto it : interfaces ) {
          (*it).Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
          ArrayVariable value( (*pit).second.dataDepth );
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
          FlaggedArrayVariable value( (*pit).second.dataDepth );
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
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );
  sort( nodes.begin(), nodes.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );
  sort( elmts.begin(), elmts.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root faces to find all its faces	
  deque<Face<dim>*> faces;
  exploreFacesFromMesh( *this, faces );
  sort( faces.begin(), faces.end(), []( auto& lhs, auto& rhs ) {return lhs->Idx() < rhs->Idx(); } );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( *this, interfaces );
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
        ArrayVariable value( key.dataDepth );
        size_t i( 0U );
        for ( auto n : nodes ) {
          read( (*pit).second, i++, value );
          n->Store( key, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
        size_t i( 0U );
        for ( auto e : elmts ) {
          read( (*pit).second, i++, value );
          e->Store( key, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
        size_t i( 0U );
        for ( auto f : faces ) {
          read( (*pit).second, i++, value );
          f->Store( key, value );
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
        size_t i( 0U );
        for ( auto e : interfaces ) {
          read( (*pit).second, i, value );
          e->Store( key, value );
          i++;
        }
      }
                  break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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
        ArrayVariable value( key.dataDepth );
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
        FlaggedArrayVariable value( key.dataDepth );
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

/**
    Only use these const object methods when you DO NOT want to modify the state of the Node
*/
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


//Qi added
template<size_t dim>
NodeManifoldManager<dim>&  MeshManager<dim>::Manifold() { return node_manifold_manager_; }

template<size_t dim>
const NodeManifoldManager<dim>&  MeshManager<dim>::Manifold() const { return node_manifold_manager_; }


template<size_t dim>
void MeshManager<dim>::Out() const
{
  cout << "\nMeshManager<" << dim << ">::Out: " << endl;

  // traversal of the existing mesh nodes to find all its elements	
  deque<const Node<dim>*>    nodes;
  deque<const Element<dim>*> elmts;
  exploreNodesAndElementsFromMesh( *this, nodes, elmts );

  // traversal of the existing mesh root faces to find all its faces	
  deque<const Face<dim>*> faces;
  exploreFacesFromMesh( *this, faces );

  // traversal of the existing mesh root interfaces to find all its interfaces	
  deque<const InterFace<dim>*> interfaces;
  exploreInterFacesFromMesh( *this, interfaces );

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


// ======================================================================================
//
//         Non-Member Functions
//
// ======================================================================================


/**  findContiguousMeshPatch

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

@author SKM
@date 21/5/2021

*/
template<size_t dim,template<size_t> class CELL>
void findContiguousMeshPatch( CELL<dim>* const eptr, set<CELL<dim>*>& elements_contiguous_subset )
 {
    assert( eptr != nullptr );
    // identifying the neighbors of the first element to be looked at
    deque<CELL<dim>*>  neighbor_elements;
    const size_t  neighbors(eptr->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( eptr->Neighbor(i) != nullptr )
        neighbor_elements.push_back( eptr->Neighbor(i) );
 
    // performing the floodfill, starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( static_cast<CELL<dim>*>(eptr) );
      
    // element set for subsequent passes
    deque<CELL<dim>*>  new_neighbor_elements;
   
     while( !neighbor_elements.empty() )
       {
          // 1. loop over those neighbors that are not already part of the deque
          for ( typename deque<CELL<dim>*>::const_iterator
                nit=neighbor_elements.begin(); nit!=neighbor_elements.end(); ++nit )
            // if the element has not already been dealt with
            if ( elements_contiguous_subset.find( (*nit) ) == elements_contiguous_subset.end() ) {
                const size_t  neighbors((*nit)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( (*nit)->Neighbor(j) != nullptr )
                    new_neighbor_elements.push_back( (*nit)->Neighbor(j) );
                elements_contiguous_subset.insert( (*nit) );
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_elements.clear();
      }
     
 } // end findContiguousMeshPatch


template void findContiguousMeshPatch( Element<1U>* const, set<Element<1U>*>& );
template void findContiguousMeshPatch( Element<2U>* const, set<Element<2U>*>& );
template void findContiguousMeshPatch( Element<3U>* const, set<Element<3U>*>& );

template void findContiguousMeshPatch( Face<1U>* const, set<Face<1U>*>& );
template void findContiguousMeshPatch( Face<2U>* const, set<Face<2U>*>& );
template void findContiguousMeshPatch( Face<3U>* const, set<Face<3U>*>& );

template void findContiguousMeshPatch( InterFace<1U>* const, set<InterFace<1U>*>& );
template void findContiguousMeshPatch( InterFace<2U>* const, set<InterFace<2U>*>& );
template void findContiguousMeshPatch( InterFace<3U>* const, set<InterFace<3U>*>& );





/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these a pointer is inserted into the argument deque (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
      TODO: implement version that does not rely on std::set and uses binary_search on vector
      TODO: think of the attributes that mesh patches should be endowed with (number of elements, surface, volume, mixed, inside outside etc.
*/
template<size_t dim, template<size_t> class CELL>
size_t  initialisePointersToStandAloneMeshPatches( std::vector<CELL<dim>*>::const_iterator begin,
                                                   std::vector<CELL<dim>*>::const_iterator end,
                                                   std::map<CELL<dim>*,MeshPatchAttributes>& root_pointers )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( begin == end ) {
             csmp_error.notice( WARNING, "initialisePointersToIndependentMeshPatches:",
                                         "input CELL pointer range is empty." );
            return 0U;
        }
     
      // getting a copy of the element pointers of the mesh
      vector<CELL<dim>*>  elements( begin, end );
      sort( elements.begin(), elements.end() );

      // detecting via a flood-fill whether the group can be partitioned, else nothing is done
      set<CELL<dim>*>  elements_contiguous_subset;
      findContiguousMeshPatch( elements.begin(), elements_contiguous_subset );
      
      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( elements.size() <= elements_contiguous_subset.size() ) {
           std::cout <<"\nModel<" << dim << ">::initialisePointersToStandAloneMeshPatches: ";
           std::cout <<"mesh is already contiguous, nothing was done."<< std::endl;
           return 0U;
        }

      // else partitions can be created
      std::string  group_name("meshPatch");
      std::string  subgroup_name;
      char         num[128];
      size_t       n_subgroups(1);
    
      // creating new contiguous group from the element subset
      while ( !elements.empty() )
        {
           // creating name of contiguous subgroup
           sprintf( num, "%lu", n_subgroups );
           subgroup_name = group_name + num;
           if ( n_subgroups == 1U ) {
                 std::cout <<"\ninitialisePointersToStandAloneMeshPatches: ";
                 std::cout <<"mesh is divided into the subregion(s):\n";
             }
           std::cout <<"\t\t\t'"<< subgroup_name <<"'";
           std::cout <<" ("<< elements_contiguous_subset.size() <<" elmts)"<< std::endl;
         
           // subtracting the elements that constitute the new group from the remaining element list
           for ( typename std::set<Element<dim>*>::const_iterator
                 sit=elements_contiguous_subset.begin(); sit!=elements_contiguous_subset.end(); ++sit )
             elements.erase( (*sit) );

           // computing the next subset
           if ( elements.empty() ) break;
           else {
                MeshPatchAttributes attributes( elements_contiguous_subset.size(),
                                                parseFiniteElementDimension((*elements.begin()).FE_Type()) );
                
                pair<map<CELL<dim>*,MeshPatchAttributes>::iterator,bool>
                  insertion = root_pointers.insert( make_pair( (*elements.begin()), attributes ) );
                if ( insertion.second == false ) {
                     csmp_error.notice( WARNING, "initialisePointersToIndependentMeshPatches:",
                                                 "mesh patch could not be inserted into root cell map. Does it already exist?" );
                  }
                findContiguousMeshPatch( (*elements.begin()), elements_contiguous_subset );
             }
           n_subgroups++;
        }

      return n_subgroups;
    
   } // end initialisePointersToStandAloneMeshPatches

// 3D version
template size_t  initialisePointersToStandAloneMeshPatches( vector<Element<3U>*>::const_iterator,
                                                            vector<Element<3U>*>::const_iterator,
                                                            map<Element<3U>*,MeshPatchAttributes>& );

template size_t  initialisePointersToStandAloneMeshPatches( vector<Face<3U>*>::const_iterator,
                                                            vector<Face<3U>*>::const_iterator,
                                                            map<Face<3U>*,MeshPatchAttributes>& );

template size_t  initialisePointersToStandAloneMeshPatches( vector<InterFace<3U>*>::const_iterator,
                                                            vector<InterFace<3U>*>::const_iterator,
                                                            map<InterFace<3U>*,MeshPatchAttributes>& );


} // end namespace csmp 
