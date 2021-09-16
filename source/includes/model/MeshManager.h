#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
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
  
  // not copy constructible
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

  /// MeshManager swallows the mesh objects created by a user; const pointer cannot be modified
  Node<dim>*	    const	AddNode( const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT );
  
  Element<dim>*	  const	AddElement( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                    const LocalVariables&, const IntegrationPointVariables& );
  Element<dim>*	  const	AddElement( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                    const LocalVariables&, const IntegrationPointVariables&, size_t material_id );
                                    
  Face<dim>*	    const	AddFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                 const LocalVariables&, const IntegrationPointVariables& );
  Face<dim>*	    const	AddFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                 Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                                 const LocalVariables&,
                                 const IntegrationPointVariables& );
                                 
  // TODO: do we need to also supply the nodes?
  Face<dim>*	    const	AddFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                 Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                                 const std::vector<Node<dim>*>&  edge_nodes,
                                 const LocalVariables&,
                                 const IntegrationPointVariables& );

  InterFace<dim>*	const	AddInterFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                      const LocalVariables&, const IntegrationPointVariables& );
  InterFace<dim>*	const	AddInterFace( csmp::FiniteElement* const, const csmp::FiniteVolumeStencil<dim>* const,
                                      Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                                      Element<dim>* const intervening_elmt,
                                      const std::vector<Node<dim>*>&  edge_nodes, // do these need to be supplied?
                                      const LocalVariables&,
                                      const IntegrationPointVariables& );
  
   /// inserts new  object if it does not already exist in the tree, otherwise returns pointer to existing one.
  Node<dim>* const      Duplicate( const Node<dim>* const );
  Element<dim>* const   Duplicate( const Element<dim>* const );

  /// by location only, the parent element storage is not initialised; TODO: is this method needed
  Node<dim>* const			AddNodeAt( const Point<dim>&, const LocalVariables&,
                                   bool only_add_if_not_collocated, BOX_BOUNDARY = NOT );
 
  /// updates the connectivity of the mesh after its modification mesh
  void Update();

  /// after disconnecting the nodes from potential manifolds, and parent elements, these are deleted
  size_t Erase( typename std::vector<Node<dim>*>::iterator first,
                typename std::vector<Node<dim>*>::iterator last );

  /// erases the supplied sequence of elements returning the number of erasures, the pointers to the erased elements are nulled. The connectivity of the affected mesh neighborhood will get fixed.
  size_t Erase( typename std::vector<Element<dim>*>::iterator first,
                typename std::vector<Element<dim>*>::iterator last );

  size_t Erase( typename std::vector<InterFace<dim>*>::iterator first,
                typename std::vector<InterFace<dim>*>::iterator last );

  size_t Erase( typename std::vector<Face<dim>*>::iterator first,
                typename std::vector<Face<dim>*>::iterator last );

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence=false ) const;

  /// computes deques of numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>&, std::deque<const Node<dim>*>&, std::deque<const Element<dim>*>&,
                     std::deque<const Face<dim>*>&, std::deque<const InterFace<dim>*>&  ) const;

  /// adds distributed variables to the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, 
                                const std::deque<const Node<dim>*>&, 
                                const std::deque<const Element<dim>*>&,
                                const std::deque<const Face<dim>*>&,
                                const std::deque<const InterFace<dim>*>&,
                                VSet<dim>& ) const;
  
  /// reads distributed variables from VSet
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

  /// Rebuild parent relationships, for example after a region was removed
  void RebuildParentRelationships( typename std::vector<Node<dim>*>::iterator begin, typename std::vector<Node<dim>*>::iterator end );
  
  /// JCK's method to test the connectivity of a mesh after it had been read from binary file
  int32 CheckElementConnectivity( const MeshManager<dim>& );

  /// prints stored objects and their connectivity to screen
  void Out() const;
  
  
private:

  /// re-establishes the neighbor connectivity between cells of the same dimensionality (Elements, Faces, InterFaces)
  template<template<size_t> class CELL>
  void RebuildConnectivity(  typename std::deque<CELL<dim>*>::iterator first,
                             typename std::deque<CELL<dim>*>::iterator last );

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
  
  NodeManifoldManager<dim>    node_manifold_manager_; ///<  node manifolds of SplitBoundaries
};

// POTENTIALLY NEEDED METHODS
/*
  /// uses same nodes and data all existing connections and data; argument vector is only used if non-empty
  Face<dim>* const      FaceFromElement( const Element<dim>* const,
                                         const LocalVariables&, const IntegrationPointVariables&,
                                         Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                                         const std::vector<Face<dim>*>& nbor_faces );
               
  /// generates  extra nodes, if these are not already present as indicated by NodeManifolds
  InterFace<dim>* const InterFaceFromFace( const Face<dim>* const,
                                           const LocalVariables&, const IntegrationPointVariables&,
                                           Element<dim>* const inner_parent, Element<dim>* const outer_parent,
                                           const std::vector<InterFace<dim>*>& nbor_faces );

  /// replace lower-dimensional element with Face object, deleting the Elements and establishing the neighbor connectivity of the new Faces
  size_t          ReplaceElementsByFaces( const PropertyDatabase<dim>&,
                                          typename std::vector<Element<dim>*>::iterator first,
                                          typename std::vector<Element<dim>*>::iterator last );
*/

} // end namespace csmp

#endif
