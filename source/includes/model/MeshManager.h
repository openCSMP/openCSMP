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
template<size_t> struct IndexToPointerMapping;
template<size_t> class PropertyDatabase;
template<size_t> class VSet;
template<size_t> class FiniteVolumeStencilManager;
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
  MeshManager( const PropertyDatabase<dim>&, const FiniteElementManager&, const VSet<dim>& );
  
  /// MeshManager is not copy constructible
  MeshManager( const MeshManager& ) = delete;
  MeshManager&  operator=( const MeshManager& ) = delete;

  ~MeshManager();

  /// sets up distributed storage for variables, finite elements, and mesh connectivity, returns vectors of pointers remembering index-pointer mapping
  bool Initialize( const PropertyDatabase<dim>&,
                   const FiniteElementManager&,
                   const VSet<dim>& );

  /// assigns finite volume stencils to the FV pointers stored in each element
  void InitializeFiniteVolumeStencils( const PropertyDatabase<dim>&,
                                       const FiniteElementManager&,
                                       FiniteVolumeStencilManager<dim>& );
                                       
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

  
  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================

// TODO: check whether there are any use cases where the neighbor information can be applied during construction; else remove this optional argument

  /// by location only, no parent element  gets connected
  Node<dim>* const		 AddNodeAt( const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT );

  /// only if there is not already a node at this location, else a pointer to that node is returned, no parent element  gets connected
  Node<dim>* const		 AddNodeAtUniqueLocation( const Point<dim>&, size_t nearby_node,
                                                const LocalVariables&,
                                                BOX_BOUNDARY = NOT );

  /// method tries to find neighbors through the parent connectivity of the nodes
  Element<dim>*	const AddElement( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                  const LocalVariables&, const IntegrationPointVariables&,
                                  const std::vector<Node<dim>*>& nodes, int32 material_id );

  /// puts lower-dimensional element inside of an InterFace, connecting it to its base pointer; the neighbors are not connected yet
  Element<dim>*	const AddInterveningElement( csmp::InterFace<dim>* const,
                                             const LocalVariables&, const IntegrationPointVariables&,
                                             const std::vector<Node<dim>*>& nodes,
                                             int32 material_id );

  /// compatibility checks are performed
  Face<dim>* const ReplaceElementByFace( csmp::Element<dim>* eptr,
                                         csmp::Element<dim>* inner_eptr,
                                         csmp::Element<dim>* outer_eptr,
                                         const LocalVariables&,
                                         const IntegrationPointVariables& );
     
  /// the neighbor element pointers are not assigned; @note node pointers must be supplied in CCW order from outside looking in
  Face<dim>* const AddFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                            Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                            const LocalVariables&,
                            const IntegrationPointVariables&,
                            const std::vector<Node<dim>*>& nodes );

  /// adds Face that caps a higher-dimensional Element at the model boundary
  Face<dim>* const AddBoundaryFace( csmp::Element<dim>* const innerParent,
                                    size_t local_face_id,
 //                                   const FiniteElementManager& fe_manager, // put into MeshManager
                                    const LocalVariables&,
                                    const IntegrationPointVariables& ); ///< optional

  /// like AddFace, but with double the nodes (inside & outside) and neighbors; extra option to assign a precreated intervening element
  InterFace<dim>*	const	AddInterFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                      Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                                      Element<dim>* const intervening_elmt,
                                      const LocalVariables&,
                                      const IntegrationPointVariables& );
   /// compatibility checks are performed
  InterFace<dim>* const ReplaceFaceByInterFace( csmp::Face<dim>* eptr,
                                                const LocalVariables&,
                                                const IntegrationPointVariables& );
 
   /// duplicates Node, automatically creating a node manifold or adding it to an existing one.
  Node<dim>* const      Duplicate( Node<dim>* const nptr_inside,
                                   INTERFACE_SIDE new_node_side,
                                   ManifoldType geometry );

 
  /// updates all connectivity (elements, faces, interfaces, nodes to parents); call after  mesh modification
  void UpdateConnectivity();
  // TODO: create version of method that permits selective update of cells

  /// after disconnecting the nodes from potential manifolds, and parent elements, these are deleted
  size_t Delete( typename std::deque<Node<dim>*>::iterator first,
                 typename std::deque<Node<dim>*>::iterator last );

  /// deletes supplied sequence of elements returning their number; pointers are nulled for erase via  EraseNullPointerCells()
  // TODO: update connectivity of remaining mesh, in the meantime, call UpdateConnectivity()
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

  /// compacts deques, first filling in deleted cells with cells from the back; then erasing cells at the back
  size_t EraseNullPointerCells();

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence=false ) const;

  /// computes deques of numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>& ) const;

  /// adds distributed variables to the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, VSet<dim>& ) const;
  
  /// reads distributed variables from VSet
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

  /// re-establishes the neighbor connectivity between cells of the same dimensionality (Elements & Faces)
  /// @todo disambiguate connectivity between Face and InterFace object at manifolds
  template<template<size_t> class CELL>
  void RebuildConnectivity(  typename std::deque<CELL<dim>*>::iterator first,
                             typename std::deque<CELL<dim>*>::iterator last );

  /// same for node-to-node connectivity
  void RebuildConnectivity(  typename std::deque<Node<dim>*>::iterator first,
                             typename std::deque<Node<dim>*>::iterator last );

  /// Rebuild node-to-element parent relationships, for example after a region was removed
  void RebuildParentRelationships( typename std::vector<Node<dim>*>::iterator begin, typename std::vector<Node<dim>*>::iterator end );
  
  /// JCK's method to test the connectivity of a mesh after it had been read from binary file
  int32 CheckElementConnectivity( const MeshManager<dim>& );

  /// prints stored objects and their connectivity to screen
  void Out() const;
  
  
private:

  /// detecting and counting potentially empty cells or nodes in storage for prompting an update
  std::pair<std::array<size_t,4>,bool>  NullPointersInStorage() const;

private:
  /// access is via root node or element only    
  bool hybrid_element_mesh_;	///< true if the mesh consists of different FE types

  // root pointers to contiguous mesh patches
  std::deque<Node<dim>*>      nodes_;          ///<  nodes
  std::deque<Element<dim>*>   elements_;       ///<  pointers elements
  std::deque<Face<dim>*>      faces_;          ///<  pointers faces making up the boundaries
  std::deque<InterFace<dim>*> interfaces_;     ///<  pointers to interfaces making up the split boundaries
  // only used in models that contain node SplitBoundaries / IterFace objects
  NodeManifoldManager<dim>*   node_manifold_manager_ = nullptr; ///<  node manifolds of SplitBoundaries
};

} // end namespace csmp

#endif
