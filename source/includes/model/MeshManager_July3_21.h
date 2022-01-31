#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "MeshPatchAttributes.h"
#include "NodeManifoldManager.h" 

namespace csmp {

class ModelTopology;
class FiniteElementManager;
template<size_t> struct IndexToPointerMapping;
template<size_t> class PropertyDatabase;
template<size_t> class VSet;
template<size_t> class FiniteVolumeStencilManager;

/**
@brief Helper class of the Model which takes care of the storage of Element, Face and InterFace objects;
internal application is hidden and may vary between models (tree-storage is default).

@author S.K. Matthai
@date 2007

@remark gain access via Mesh() public interface of Model.

@attention the MeshManager takes care of the creation and destruction of Elements, Faces or Interfaces.
Region or Boundary objects merely contain pointers to these.
*/
template<size_t dim>
class MeshManager {
public:
  MeshManager();
  MeshManager( const PropertyDatabase<dim>&, const FiniteElementManager&, const VSet<dim>&, IndexToPointerMapping<dim>& );
  MeshManager&  operator=( const MeshManager& );

  ~MeshManager();

  /// sets up distributed storage for variables, finite elements, and mesh connectivity, returns vectors of pointers remembering index-pointer mapping
  bool Initialize( const PropertyDatabase<dim>&,
                   const FiniteElementManager&,
                   const VSet<dim>&,
                   IndexToPointerMapping<dim>& );

  /// updates the root pointers of the mesh after modifying the mesh
  void Update();
  void Update( std::deque<Node<dim>*> nodes, std::deque<Element<dim>*> elmts );

  /// assigns finite volume stencils to the FV pointers stored in each element
  void InitializeFiniteVolumeStencils( const PropertyDatabase<dim>&,
                                       const FiniteElementManager&,
                                       FiniteVolumeStencilManager<dim>& );

  /// returns true if the mesh consists of multiple element types
  bool HybridElementMesh() const;
  
  /// returns whether the mesh consists out of disconnected mesh patches; TODO: analyse this properly
  bool IsContiguous() const { return (elmt_tree_roots_.size() == 1); }

  /// counts and returns current indices of elements that may give rise to problems during the assignment of boundary conditions
  size_t DetectElementsWithAllNodesOnBoundary( std::set<size_t>& ) const;

  /// returns number of nodes=vertices in the current mesh
  size_t Nodes() const;

  /// returns number of node groups
  size_t NodeGroups() const;

  /// returns number of elements in the current mesh
  size_t Elements() const;

  /// returns number of element groups
  size_t ElementGroups() const;

  /// returns number of Faces=lower-dimensional elements in current mesh
  size_t Faces() const;

  /// returns number of face groups
  size_t FaceGroups() const;

  /// returns number of interfaces=faces with multiplicated nodes
  size_t InterFaces() const;

  /// returns number of interface groups
  size_t InterFaceGroups() const;

  /// the node's root pointer at the basis of node tree for contiguous regions
  Node<dim>*		 RootNode( size_t mesh_patch_root_node_index );

  /// assigns node's root pointer at the basis of node tree for contiguous regions
  void				   SetRootNode( Node<dim>* );

  /// the element's root pointer at the basis of element tree for contiguous regions
  Element<dim>*	 RootElement( size_t group_idx );

  /// assigns element's root pointer at the basis of element tree for contiguous regions
  void				   SetRootElement( Element<dim>* );

  /// the face's root pointer at the basis of face tree for multiple boundaries
  Face<dim>*		 RootFace( size_t group_idx );

  /// assigns face's root pointer at the basis of face tree for multiple boundaries
  void				   SetRootFace( Face<dim>* );

  /// the interface's root pointer at the basis of interface tree for multiple split boundaries
  InterFace<dim>*  RootInterFace( size_t group_idx );

  /// assigns interface's root pointer at the basis of interface tree for multiple split boundaries
  void				     SetRootInterFace( InterFace<dim>* );

  /// the node's root pointer at the basis of node tree for contiguous regions
  const Node<dim>*	RootNode( size_t group_idx ) const;

  /// the element's root pointer at the basis of element tree for contiguous regions
  const Element<dim>*	  RootElement( size_t group_idx ) const;

  /// the face's root pointer at the basis of face tree for multiple boundaries
  const Face<dim>*	    RootFace( size_t group_idx ) const;

  /// the interface's root pointer at the basis of interface tree for multiple split boundaries
  const InterFace<dim>* RootInterFace( size_t group_idx ) const;
  
  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================

  /// inserts new  object if it does not already exist in the tree, otherwise returns pointer to existing one.
  Node<dim>*      Duplicate( const Node<dim>* const );
  Element<dim>*   Duplicate( const Element<dim>* const );

  Node<dim>*			AddNodeAt( const Point<dim>& ); // give connections as well
  Element<dim>*		AddElment( std::vector<Node<dim>* const>& );
  
  /// replace lower-dimensional element with Face object
  Face<dim>*      ReplaceElementWithFace( const Element<dim>& );

  /// removes primitive without checking its connectivities since its pointer indicates its original primitive.
  void Erase( Node<dim>* );
  void Erase( Element<dim>* );
  void Erase( Face<dim>* );
  void Erase( InterFace<dim>* );

  /// erase all objects of the given type
  bool EraseNodes();
  bool EraseElements();
  bool EraseFaces();
  bool EraseInterFaces();
  
  // SKM new 22/5/21
  /// erases the supplied sequence of elements returning the number of erasures, the pointers to the erased elements are nulled. The connectivity of the affected mesh neighborhood will get fixed.
  size_t EraseElements( typename std::vector<Element<dim>*>::iterator first_elmt,
                        typename std::vector<Element<dim>*>::iterator last_elmt );


  /// supplied Face objects are replaced by InterFace objects (Face objects are deleted); returns pointer to first perimeter interface.
  typename std::vector<InterFace<dim>*>::iterator  ConvertFacesToInterFaces( typename std::vector<Face<dim>*>::iterator faces_begin,
                                                                             typename std::vector<Face<dim>*>::iterator faces_end,
                                                                             typename std::vector<Node<dim>*>::iterator perimeter_nodes_begin,
                                                                             typename std::vector<Node<dim>*>::iterator perimeter_nodes_end,
                                                                             typename std::vector<InterFace<dim>*>::iterator interfaces_begin,
                                                                             typename std::vector<InterFace<dim>*>::iterator interfaces_end );

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence );

  /// computes deques of numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>&, 
                     std::deque<const Node<dim>*>&, std::deque<const Element<dim>*>&, std::deque<const Face<dim>*>&, std::deque<const InterFace<dim>*>&  ) const;

  /// adds distributed variables to the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, 
                                const std::deque<const Node<dim>*>&, 
                                const std::deque<const Element<dim>*>&, const std::deque<const Face<dim>*>&, const std::deque<const InterFace<dim>*>&, 
                                VSet<dim>& ) const;
  
  /// reads distributed variables from VSet
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

  /// Rebuild parent relationships, for example after a region was removed
  void RebuildParentRelationships( typename std::vector<Node<dim>*>::iterator begin, typename std::vector<Node<dim>*>::iterator end );

  // TODO: why is there a need to expose the manager?
  const NodeManifoldManager<dim>&  Manifold() const;
  NodeManifoldManager<dim>&  Manifold();  
  
  /// prints stored objects and their connectivity to screen
  void Out() const;

private:
  /// access is via root node or element only    
  bool								hybrid_element_mesh_;	///< true if the mesh consists of different FE types

  size_t							n_nodes_;				      ///< total number of nodes
  size_t							n_elmts_;				      ///< total number of elements
  size_t							n_faces_;				      ///< total number of faces
  size_t							n_interfaces_;        ///< total number of interfaces

  // root pointers to contiguous mesh patches
  std::map<Element<dim>*,MeshPatchAttributes>   elmt_tree_roots_;       ///< pointers to root elements for contiguous regions
  std::map<Face<dim>*,MeshPatchAttributes>      face_tree_roots_;       ///< pointers to root faces for multiple boundaries
  std::map<InterFace<dim>*,MeshPatchAttributes> itfc_tree_roots_;       ///< pointers to root interfaces for multiple split boundaries
  NodeManifoldManager<dim>                      node_manifold_manager_; ///< stores which nodes are co-located at SplitBoundaries
};


// NON-MEMBER FUNCTIONS

template<size_t dim, template<size_t> class CELL>
size_t  initialisePointersToStandAloneMeshPatches( typename std::vector<CELL<dim>*>::const_iterator begin,
                                                   typename std::vector<CELL<dim>*>::const_iterator end,
                                                   std::map<Element<dim>*,MeshPatchAttributes>& );

template<size_t dim,template<size_t> class CELL>
void findContiguousMeshPatch( CELL<dim>* const eptr, std::set<CELL<dim>*>& elements_contiguous_subset );


} // end namespace csmp

#endif
