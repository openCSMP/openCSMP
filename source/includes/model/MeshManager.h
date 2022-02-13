#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include "Node.h"
#include "NodeManifold.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "plf_colony.h"

namespace csmp {

class MeshManager_Test;
class ModelTopology;
template<size_t> struct IndexToPointerMapping;
template<size_t> class PropertyDatabase;
template<size_t> class VSet;
template<size_t> class NodeManifoldManager;
template<size_t,template<size_t> class> class ModelSubDomain;

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
  
  typename plf::colony<Node<dim> >::iterator      NodesBegin();
  typename plf::colony<Node<dim> >::iterator      NodesEnd();

  typename plf::colony<Element<dim> >::iterator   ElementsBegin();
  typename plf::colony<Element<dim> >::iterator   ElementsEnd();

  typename plf::colony<Face<dim> >::iterator      FacesBegin();
  typename plf::colony<Face<dim> >::iterator      FacesEnd();

  typename plf::colony<InterFace<dim> >::iterator InterFacesBegin();
  typename plf::colony<InterFace<dim> >::iterator InterFacesEnd();

  typename plf::colony<NodeManifold<dim> >::iterator NodeManifoldsBegin();
  typename plf::colony<NodeManifold<dim> >::iterator NodeManifoldsEnd();

  // const versions
  typename plf::colony<Node<dim> >::const_iterator      NodesBegin() const;
  typename plf::colony<Node<dim> >::const_iterator      NodesEnd() const;

  typename plf::colony<Element<dim> >::const_iterator   ElementsBegin() const;
  typename plf::colony<Element<dim> >::const_iterator   ElementsEnd() const;

  typename plf::colony<Face<dim> >::const_iterator      FacesBegin() const;
  typename plf::colony<Face<dim> >::const_iterator      FacesEnd() const;

  typename plf::colony<InterFace<dim> >::const_iterator InterFacesBegin() const;
  typename plf::colony<InterFace<dim> >::const_iterator InterFacesEnd() const;
  
  typename plf::colony<NodeManifold<dim> >::const_iterator NodeManifoldsBegin() const;
  typename plf::colony<NodeManifold<dim> >::const_iterator NodeManifoldsEnd() const;

/* not for colony
  Node<dim>* const      N( size_t ) const;
  Element<dim>* const   E( size_t ) const;
  Face<dim>* const      F( size_t ) const;
  InterFace<dim>* const I( size_t ) const;
*/
  /// direct access for backward compatibility
  const FiniteElementManager& FiniteElements() const { return fem_manager_; }

  /// direct access for backward compatibility
  const FiniteVolumeStencilManager<dim>& FiniteVolumes() const { return fvm_manager_; }

  
  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================
  
  /// for all outside-facing CELL neighbor perimeter face pointers that are not nullptr,/ set  neighbors of the corresponding cell to null, getting subdomain ready for deletion
  template<template<size_t> class CELL>
  size_t DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<dim,CELL>& );
  
  /// replaces supplied lower-dimensional elements with Face objects, establishing their connectivity; the Elements are deleted afterwards, setting input pointers to NULL
  std::vector<Face<dim>*>  ReplaceElementsByFaces( const PropertyDatabase<dim>&,
                                                   typename std::vector<Element<dim>*>::iterator first,
                                                   typename std::vector<Element<dim>*>::iterator last );

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
                                  const std::vector<Node<dim>*>& nodes, int32_t material_id );

  /// puts lower-dimensional element inside of an InterFace, connecting it to its base pointer; the neighbors are not connected yet
  Element<dim>*	const AddInterveningElement( csmp::InterFace<dim>* const,
                                             const LocalVariables&,
                                             const IntegrationPointVariables&,
                                             const std::vector<Node<dim>*>& nodes,
                                             int32_t material_id );

  /// creates a Face matching the current lower-dimensional element and deletes the element subsequently
  Face<dim>* const ReplaceElementByFace( csmp::Element<dim>* eptr,
                                         csmp::Element<dim>* inner_eptr,
                                         csmp::Element<dim>* outer_eptr,
                                         size_t adjacent_face_of_inner_element,
                                         size_t adjacent_face_of_outer_element,
                                         const LocalVariables& face_variables,
                                         const IntegrationPointVariables& face_integration_point_variables );

  /// creates Face matching the supplied lower-dimensional element but without deleting the underlying element 
  Face<dim>* const ConstructFaceFromElement( csmp::Element<dim>* eptr,
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

  /// updates all connectivity (elements, faces, interfaces, nodes to parents); however, node manifolds are not reconstructed
  void UpdateConnectivity();
  // TODO: create version of method that permits selective update of cells

  /// re-establishes the neighbor connectivity between cells of the same dimensionality (Elements & Faces)
  /// @todo disambiguate connectivity between Face and InterFace object at manifolds
  template<template<size_t> class CELL>
  void BuildConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                          typename std::vector<CELL<dim>*>::iterator last );

  template<template<size_t> class CELL>
  void BuildVolumeConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                                typename std::vector<CELL<dim>*>::iterator last );
                                       
  template<template<size_t> class CELL>
  void BuildSurfaceConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                                 typename std::vector<CELL<dim>*>::iterator last );
                                        
  template<template<size_t> class CELL>
  void BuildLineConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                              typename std::vector<CELL<dim>*>::iterator last );

  /// Starting with an existing node-to-parent element relationships, these are validated, removing excess connections, for example after a region was removed
  void RebuildNodeParentElementRelationships();


  // DELETIONS & MAINTANANCE OF MESH CONNECTIVITY
  // --------------------------------------------
  // NB: elements are responsible for the nodes, nodes for their manifolds

  /// after disconnecting the nodes from potential manifolds, the supplied range of nodes is deleted
  size_t Delete( typename std::vector<Node<dim>*>::iterator first,
                 typename std::vector<Node<dim>*>::iterator last );

  /// deletes the supplied range of elements, and orphaned nodes if any; the parent element storage of the nodes is rebuild; @note all input pointers are nulled
  size_t Delete( typename std::vector<Element<dim>*>::iterator first,
                 typename std::vector<Element<dim>*>::iterator last );

  size_t Delete( typename std::vector<InterFace<dim>*>::iterator first,
                 typename std::vector<InterFace<dim>*>::iterator last );

  size_t Delete( typename std::vector<Face<dim>*>::iterator first,
                 typename std::vector<Face<dim>*>::iterator last );

  /// JCK's method to test the connectivity of a mesh after it had been read from binary file
  int32_t CheckElementConnectivity() const;

  /// prints stored objects and their connectivity to screen
  void Out() const;
  
  
private:

  void Erase( typename plf::colony<Node<dim>>::const_iterator nit ) { nodes_.erase(nit); }
  void Erase( typename plf::colony<Element<dim>>::const_iterator it ) { elements_.erase(it); }
  void Erase( typename plf::colony<Face<dim>>::const_iterator it ) { faces_.erase(it); }
  void Erase( typename plf::colony<InterFace<dim>>::const_iterator it ) { interfaces_.erase(it); }

  size_t Erase( typename plf::colony<Node<dim>>::const_iterator first,
                typename plf::colony<Node<dim>>::const_iterator last ) { nodes_.erase(first,last); return nodes_.size(); }
              
  size_t Erase( typename plf::colony<Element<dim>>::const_iterator first,
                typename plf::colony<Element<dim>>::const_iterator last ) { elements_.erase(first,last); return elements_.size(); }
              
  size_t Erase( typename plf::colony<Face<dim>>::const_iterator first,
                typename plf::colony<Face<dim>>::const_iterator last ) { faces_.erase(first,last); return faces_.size(); }
              
  size_t Erase( typename plf::colony<InterFace<dim>>::const_iterator first,
                typename plf::colony<InterFace<dim>>::const_iterator last ) { interfaces_.erase(first,last); return interfaces_.size(); }

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence=false );
  
  /// puts nodes, elements, faces, and interfaces into the order given by Idx() variables; removes nullptr cells first
  void ReorderObjectsByIndexes();

  /// returns numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>&, bool get_indices_from_stored_variables=false );

  /// adds distributed variables to the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, VSet<dim>& ) const;
  
  /// reads distributed variables from VSet
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

private:

  FiniteElementManager             fem_manager_;
  FiniteVolumeStencilManager<dim>  fvm_manager_; ///< current finite volume specifications // TODO: make this a trait class because it needs no dynamic data!

  /// access is via root node or element only
  bool hybrid_element_mesh_;	///< true if the mesh consists of different FE types

  // root pointers to contiguous mesh patches; mutable to allow for behind scene updates
  plf::colony<Node<dim>>      nodes_;          ///<  nodes
  plf::colony<Element<dim>>   elements_;       ///<  pointers elements
  plf::colony<Face<dim>>      faces_;          ///<  pointers faces making up the boundaries
  plf::colony<InterFace<dim>> interfaces_;     ///<  pointers to interfaces making up the split boundaries
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
