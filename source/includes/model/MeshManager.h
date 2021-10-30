#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"

namespace csmp {

class MeshManager_Test;
class ModelTopology;
template<size_t> struct IndexToPointerMapping;
template<size_t> class PropertyDatabase;
template<size_t> class VSet;
template<size_t> class NodeManifoldManager;
enum class ManifoldType : int8_t;

/**
@brief Helper class of the Model which takes care of the storage of Element, Face and InterFace objects;
internal application is hidden and may vary between models (tree-storage is default).

@author S.K. Matthai
@date 2021 (complete rewrite)

@remark gain access via Mesh() public interface of Model.

@attention the MeshManager takes care of the creation and destruction of Elements, Faces or Interfaces.
Region or Boundary objects merely contain pointers to these.

TODO: which kind of mesh error diagnostics should the MeshManager implement? - should these be in a separate compilation unit?
*/
template<size_t dim>
class MeshManager {
public:
  MeshManager();
  MeshManager( const PropertyDatabase<dim>&, const VSet<dim>& );
  
  /// MeshManager is not copy constructible
  MeshManager( const MeshManager& ) = delete;
  MeshManager&  operator=( const MeshManager& ) = delete;

  ~MeshManager();

  /// sets up distributed storage for variables, finite elements, and mesh connectivity, returns vectors of pointers remembering index-pointer mapping
  bool Initialize( const PropertyDatabase<dim>&, const VSet<dim>& );

  // ==============================================================
  //
  // MESH DIAGNOSTICS & ACCESS
  //
  // ==============================================================

  /// returns true if the mesh consists of multiple element types
  bool HybridElementMesh() const;
  
  /// uses a floodfill on the highest-dimensional elements in the mesh to identify whether the model consists  of disconnected mesh patches
  bool IsContiguous() const; 

  /// returns number of nodes=vertices in the current mesh
  size_t Nodes() const;

  /// returns number of elements in the current mesh
  size_t Elements() const;

  /// returns number of Faces=lower-dimensional elements in current mesh
  size_t Faces() const;

  /// returns number of InterFaces=lower-dimensional elements in current mesh
  size_t InterFaces() const;
  
  typename std::deque<Node<dim>*>::iterator      NodesBegin();
  typename std::deque<Node<dim>*>::iterator      NodesEnd();

  typename std::deque<Element<dim>*>::iterator   ElementsBegin();
  typename std::deque<Element<dim>*>::iterator   ElementsEnd();

  typename std::deque<Face<dim>*>::iterator      FacesBegin();
  typename std::deque<Face<dim>*>::iterator      FacesEnd();

  typename std::deque<InterFace<dim>*>::iterator InterFacesBegin();
  typename std::deque<InterFace<dim>*>::iterator InterFacesEnd();

  // const versions
  typename std::deque<Node<dim>*>::const_iterator      NodesBegin() const;
  typename std::deque<Node<dim>*>::const_iterator      NodesEnd() const;

  typename std::deque<Element<dim>*>::const_iterator   ElementsBegin() const;
  typename std::deque<Element<dim>*>::const_iterator   ElementsEnd() const;

  typename std::deque<Face<dim>*>::const_iterator      FacesBegin() const;
  typename std::deque<Face<dim>*>::const_iterator      FacesEnd() const;

  typename std::deque<InterFace<dim>*>::const_iterator InterFacesBegin() const;
  typename std::deque<InterFace<dim>*>::const_iterator InterFacesEnd() const;
  
  Node<dim>* const      N( size_t ) const;
  Element<dim>* const   E( size_t ) const;
  Face<dim>* const      F( size_t ) const;
  InterFace<dim>* const I( size_t ) const;
  
  /// direct access for backward compatibility
  const FiniteElementManager& FiniteElements() const { return fem_manager_; }

  /// direct access for backward compatibility
  const FiniteVolumeStencilManager<dim>& FiniteVolumes() const { return fvm_manager_; }

  
  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================

  /// by location only, no parent element  gets connected
  Node<dim>* const		 AddNodeAt( const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT );

  /// only if there is not already a node at this location, else a pointer to that node is returned, no parent element  gets connected
  Node<dim>* const		 AddNodeAtUniqueLocation( const Point<dim>&, size_t nearby_node,
                                                const LocalVariables& node_variables,
                                                BOX_BOUNDARY = NOT );

  /// method tries to find neighbors through the parent connectivity of the nodes
  Element<dim>*	const AddElement( CSMP_FEM_TYPE,
                                  const LocalVariables& element_variables,
                                  const IntegrationPointVariables& element_integration_point_variables,
                                  const std::vector<Node<dim>*>& nodes, int32 material_id );

  /// puts lower-dimensional element inside of an InterFace, connecting it to its base pointer; the neighbors are not connected yet
  Element<dim>*	const AddInterveningElement( csmp::InterFace<dim>* const,
                                             const LocalVariables&,
                                             const IntegrationPointVariables&,
                                             const std::vector<Node<dim>*>& nodes,
                                             int32 material_id );

  /// compatibility checks are performed
  Face<dim>* const ReplaceElementByFace( csmp::Element<dim>* eptr,
                                         csmp::Element<dim>* inner_eptr,
                                         csmp::Element<dim>* outer_eptr,
                                         size_t adjacent_face_of_inner_element,
                                         size_t adjacent_face_of_outer_element,
                                         const LocalVariables& face_variables,
                                         const IntegrationPointVariables& face_integration_point_variables );
     
  /// the neighbor element pointers are not assigned; @note node pointers must be supplied in CCW order from outside looking in; deduces element type
  Face<dim>* const AddFace( Element<dim>* const inner_parent, size_t inner_parent_face_id,
                            Element<dim>* const outer_parent, size_t outer_parent_face_id,
                            const LocalVariables&,
                            const IntegrationPointVariables& );

  /// creates a lower-dimensional face with parents of adjacent higher-dimensional Face objects as parents; deduces element type
  Face<dim>* const AddEdgeFace( Face<dim>* const adjacent_face1, size_t parent_elmt1_segm_id,
                                Face<dim>* const adjacent_face2, size_t parent_elmt2_segm_id,
                                const LocalVariables&,
                                const IntegrationPointVariables&,
                                const std::vector<Node<dim>*>& nodes );

  /// adds Face that caps a higher-dimensional Element at the model boundary
  Face<dim>* const AddBoundaryFace( csmp::Element<dim>* const innerParent,
                                    size_t local_face_id,
                                    const LocalVariables&,
                                    const IntegrationPointVariables& ); ///< optional

   /// assuming that the nodes on either side of the interface are already there, the face gets replaced
  InterFace<dim>* const ReplaceFaceByInterFace( csmp::Face<dim>* eptr,
                                                const LocalVariables&,
                                                const IntegrationPointVariables& );
 
  /// For connecting node-matched mesh patches, creating / updating their node manifolds
  InterFace<dim>*	const	AddInterFace( Element<dim>* const inner_parent, size_t inner_element_face_id,
                                      Element<dim>* const outer_parent, size_t outer_element_face_id,
                                      const LocalVariables& interface_variables,
                                      const IntegrationPointVariables& interface_integration_point_variables );

   /// duplicates Node, automatically creating a node manifold or adding it to an existing one.
  Node<dim>* const      Duplicate( Node<dim>* const nptr_inside,
                                   INTERFACE_SIDE new_node_side,
                                   ManifoldType geometry );

  /// Rebuild node-to-element parent relationships, for example after a region was removed
  void RebuildNodeParentElementRelationships( typename std::vector<Element<dim>*>::iterator begin,
                                              typename std::vector<Element<dim>*>::iterator end );

  /// updates all connectivity (elements, faces, interfaces, nodes to parents); however, node manifolds are not reconstructed
  void UpdateConnectivity();
  // TODO: create version of method that permits selective update of cells

  /// re-establishes the neighbor connectivity between cells of the same dimensionality (Elements & Faces)
  /// @todo disambiguate connectivity between Face and InterFace object at manifolds
  template<template<size_t> class CELL>
  void BuildConnectivity( typename std::deque<CELL<dim>*>::iterator first,
                          typename std::deque<CELL<dim>*>::iterator last );


  // DELETIONS & MAINTANANCE OF MESH CONNECTIVITY
  // --------------------------------------------
  // NB: elements are responsible for the nodes, nodes for their manifolds

  /// after disconnecting the nodes from potential manifolds, the supplied range of nodes is deleted
  size_t Delete( typename std::deque<Node<dim>*>::iterator first,
                 typename std::deque<Node<dim>*>::iterator last );

  /// deletes the supplied range of elements, and singly owned nodes if any; pointers are nulled 
  size_t Delete( typename std::deque<Element<dim>*>::iterator first,
                 typename std::deque<Element<dim>*>::iterator last );

  size_t Delete( typename std::deque<InterFace<dim>*>::iterator first,
                 typename std::deque<InterFace<dim>*>::iterator last );

  size_t Delete( typename std::deque<Face<dim>*>::iterator first,
                 typename std::deque<Face<dim>*>::iterator last );
                
  /// for any type of cells using a vector iterator
  template<template<size_t> class CELL>
  size_t Delete( typename std::vector<CELL<dim>*>::iterator first,
                 typename std::vector<CELL<dim>*>::iterator last );

  /// JCK's method to test the connectivity of a mesh after it had been read from binary file
  int32 CheckElementConnectivity() const;

  /// prints stored objects and their connectivity to screen
  void Out() const;
  
  
private:

  template<template<size_t> class CELL>
  void BuildVolumeElementConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                                       typename std::vector<CELL<dim>*>::iterator last );
                                       
  template<template<size_t> class CELL>
  void BuildSurfaceElementConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                                        typename std::vector<CELL<dim>*>::iterator last );
                                        
  template<template<size_t> class CELL>
  void BuildLineElementConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                                     typename std::vector<CELL<dim>*>::iterator last );


  /// compacts deques, first filling in deleted cells with cells from the back; then erasing cells at the back
  size_t EraseNullPointerCells();

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence=false );
  
  /// puts nodes, elements, faces, and interfaces into the order given by Idx() variables; removes nullptr cells first
  void ReorderObjectsByIndexes();

  /// computes deques of numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>&, bool get_indices_from_stored_variables=false );

  /// adds distributed variables to the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, VSet<dim>& ) const;
  
  /// reads distributed variables from VSet
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );
  
  /// detecting and counting potentially empty cells or nodes in storage for prompting an update
  std::pair<std::array<size_t,4>,bool>  NullPointersInStorage() const;

private:

  FiniteElementManager             fem_manager_;
  FiniteVolumeStencilManager<dim>  fvm_manager_; ///< current finite volume specifications // TODO: make this a trait class because it needs no dynamic data!

  /// access is via root node or element only
  bool hybrid_element_mesh_;	///< true if the mesh consists of different FE types

  // root pointers to contiguous mesh patches; mutable to allow for behind scene updates
  std::deque<Node<dim>*>      nodes_;          ///<  nodes
  std::deque<Element<dim>*>   elements_;       ///<  pointers elements
  std::deque<Face<dim>*>      faces_;          ///<  pointers faces making up the boundaries
  std::deque<InterFace<dim>*> interfaces_;     ///<  pointers to interfaces making up the split boundaries
  // only used in models that contain node SplitBoundaries / IterFace objects
  NodeManifoldManager<dim>*   node_manifold_manager_ = nullptr; ///<  node manifolds of SplitBoundaries

  friend class MeshManager_Test; ///< so that private methods can be tested
  friend class Model<dim>;       ///<  exclusive access to private member functions
};

// POTENTIAL METHODS?


   /// replaces face, constructing new nodes & manifolds where indicated by vector
//  InterFace<dim>* const ReplaceFaceByInterFace( csmp::Face<dim>* eptr,
//                                                const std::vector<bool>&  nodes_to_duplicate,
//                                                const LocalVariables&,
//                                                const IntegrationPointVariables& );



} // end namespace csmp

#endif
