#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"

namespace csmp {

class ModelTopology;
class FiniteElementManager;
template<size_t> class PropertyDatabase;
template<size_t> class VSet;
template<size_t> class FiniteVolumeStencilManager;

/**
@brief Helper class of the Model which takes care of the storage of Element, Face and InterFace objects;
internal application is hidden and may vary between models (tree-storage is default).

@author S.K. Matthaei
@date 2007

@remark gain access via Mesh() public interface of Model.

@attention the MeshManager takes care of the creation and destruction of Elements, Faces or Interfaces.
Region or Boundary objects merely contain pointers to these.
*/
template<size_t dim>
class MeshManager {
public:
  MeshManager();
  MeshManager( const PropertyDatabase<dim>&, const FiniteElementManager&, VSet<dim>& );
  MeshManager&  operator=( const MeshManager& );

  ~MeshManager();

  /// sets up distributed storage for variables, finite elements, and mesh connectivity
  bool Initialize( const PropertyDatabase<dim>&,
                   const FiniteElementManager&,
                   const VSet<dim>& );

  /// updates the root pointers of the mesh after modifying the mesh
  void Update();
  void Update( std::deque<Node<dim>*> nodes, std::deque<Element<dim>*> elmts );

  /// assigns finite volume stencils to the FV pointers stored in each element
  void InitializeFiniteVolumeStencils( const PropertyDatabase<dim>&,
                                       const FiniteElementManager&,
                                       FiniteVolumeStencilManager<dim>& );

  /// returns true if the mesh consists of multiple element types
  bool HybridElementMesh() const;

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
  Node<dim>*			RootNode( size_t group_idx );

  /// assigns node's root pointer at the basis of node tree for contiguous regions
  void				SetRootNode( Node<dim>* );

  /// the element's root pointer at the basis of element tree for contiguous regions
  Element<dim>*		RootElement( size_t group_idx );

  /// assigns element's root pointer at the basis of element tree for contiguous regions
  void				SetRootElement( Element<dim>* );

  /// the face's root pointer at the basis of face tree for multiple boundaries
  Face<dim>*			RootFace( size_t group_idx );

  /// assigns face's root pointer at the basis of face tree for multiple boundaries
  void				SetRootFace( Face<dim>* );

  /// the interface's root pointer at the basis of interface tree for multiple split boundaries
  InterFace<dim>*		RootInterFace( size_t group_idx );

  /// assigns interface's root pointer at the basis of interface tree for multiple split boundaries
  void				SetRootInterFace( InterFace<dim>* );

  /// the node's root pointer at the basis of node tree for contiguous regions
  const Node<dim>*	RootNode( size_t group_idx ) const;

  /// the element's root pointer at the basis of element tree for contiguous regions
  const Element<dim>*	RootElement( size_t group_idx ) const;

  /// the face's root pointer at the basis of face tree for multiple boundaries
  const Face<dim>*	RootFace( size_t group_idx ) const;

  /// the interface's root pointer at the basis of interface tree for multiple split boundaries
  const InterFace<dim>* RootInterFace( size_t group_idx ) const;

  /// inserts new primitive if it does not already exist in the tree, otherwise returns the pointer of the existing one.
  Node<dim>*			AddIfUnique( Node<dim>& );
  Element<dim>*		AddIfUnique( Element<dim>& );
  Face<dim>*			AddIfUnique( Face<dim>& );
  InterFace<dim>*		AddIfUnique( InterFace<dim>& );

  /// inserts new primitive without checking whether it does not already exist in the tree and its connectivity is valid
  Node<dim>*			Add( Node<dim>& );
  Element<dim>*		Add( Element<dim>& );
  Face<dim>*			Add( Face<dim>& );
  InterFace<dim>*		Add( InterFace<dim>& );

  /// removes primitive after checking its connectivities
  void Erase( Node<dim>& );
  void Erase( Element<dim>& );
  void Erase( Face<dim>& );
  void Erase( InterFace<dim>& );

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

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence ) const;

  /// outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>& ) const;

  /// storing distributed variables associated with the mesh in the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, VSet<dim>& ) const;
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

  /// Rebuild parent relationships, for example after a region was removed
  void RebuildParentRelationships( typename std::vector<Node<dim>*>::iterator begin, typename std::vector<Node<dim>*>::iterator end );

  /// prints stored objects and their connectivity to screen
  void Out() const;

private:
  /// access is via root node or element only    
  bool								hybrid_element_mesh_;	///< true if the mesh consists of different FE types

  size_t								n_nodes_;				///< total number of nodes
  size_t								n_elmts_;				///< total number of elements
  size_t								n_faces_;				///< total number of faces
  size_t								n_interfaces_;			///< total number of interfaces

  std::deque<Node<dim>*>				root_node_group_;		///< pointers to root nodes for contiguous regions
  std::deque<Element<dim>*>			root_elmt_group_;		///< pointers to root elements for contiguous regions
  std::deque<Face<dim>*>				root_face_group_;		///< pointers to root faces for multiple boundaries
  std::deque<InterFace<dim>*>			root_interface_group_;	///< pointers to root interfaces for multiple split boundaries
};

} // end namespace csmp

#endif
