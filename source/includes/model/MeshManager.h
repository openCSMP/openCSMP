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
#include "Region.h" 

namespace csmp {

class MeshManager_Test;
class ModelTopology;
template<uint32_t> struct IndexToPointerMapping;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class VSet;
template<uint32_t> class NodeManifoldManager;
template<uint32_t,template<uint32_t> class> class ModelSubDomain;
template<uint32_t> class FaceConstructionData;

/**
    @brief Manages the creation, deletion, storage and connectivity of Element, Face and InterFace objects.
    
    Contains and updates FiniteElementManager, FiniteVolumeStencilManager, and NodeManifoldManager.
    
    @pre nodes must be assigned to the elements as specified in CSMP_FEM_conventions.pdf following a right-hand coordinate system
    @pre the elements must be compatible with one-another
    @pre no hanging nodes
    @pre no degenerate elements with collocated nodes, non-manifold vertices, triangle boxes etc. must be contained in the input mesh
    @pre no interpenetrating elements
    
    The connectivity that is maintained at all times consists of:

    -  cells (Element, Face, InterFace) are connected to their equidimensional neighbors; Elements to Elements, Faces to Faces etc.

    -  there is only one such neighbor per cell; potential manifolds that can arise for lower-dimensional cells are disambiguated by choosing the most closely aligned line segments or surfaces with the nearest normal orientation

    - nodes are connected to their neighbors except for NodeManifolds connecting collocated but topologically separated nodes

    - nodes know their parent Element objects (not Faces nor Interfaces)

    - Face objects are connected to their higher dimensional Element neighbors; if located on the outside of the model, their inner neighbor will be connected so that their normal points outside of the model

    - InterFace objects can only occur on the inside of the model, being connected to two higher-dimensional Element neighbors
    
    When the MeshManager creates or deletes cells, the connectivity in the surroundings is updated
    
    @attention MESHMANAGER IS AGNOSTIC ABOUT REGIONS, BOUNDARIES, AND SPLIT BOUNDARIES; it only handles creation, deletion and connectivity of cells and nodes!
    
    @attention Changes in the mesh affect Region, Boundary and SplitBoundary objects meaning that they must updated, using
    corresponding functionality of ModelSubdomain and the Region-, Boundary- and SplitBoundary interfaces (policies) of Model.
    To support this, the MeshManager sets the supplied pointers to deleted items to NULLPOINTER.
    This way they can be deleted from the supplied cell ranges.
    If an object that should be deleted cannot, the pointer to it will not be nulled.
    
    @attention the pointers to items that shall be deleted from the colony are wrapped in iterators that allow setting the pointers to zero rather than just copies of them see method interfaces;
    For objects that are not going to be deleted, const-pointer access is granted.

    @author S.K. Matthai
    @date 2021 ( rewrite in 2024)

    @remark gain access to mesh using Mesh() public interface of Model.

    @attention the MeshManager takes care of the creation and destruction of Elements, Faces or Interfaces.
    Region or Boundary objects merely contain pointers to these.

*/
template<uint32_t dim>
class MeshManager {
public:
  MeshManager();
  MeshManager( const PropertyDatabase<dim>&, const VSet<dim>& );
  
  /// MeshManager is not copy constructible
  MeshManager( const MeshManager& ) = delete;
  MeshManager&  operator=( const MeshManager& ) = delete;

  ~MeshManager();

  /// sets up distributed storage for variables, finite elements, and mesh connectivity, returns vectors of pointers remembering index-pointer mapping
  bool Initialize( const PropertyDatabase<dim>&, const VSet<dim>&, bool initialise_FV_stencils );

  // ==============================================================
  //
  // MESH DIAGNOSTICS & ACCESS
  //
  // ==============================================================

  /// returns true if the mesh consists of multiple element types
  bool HybridElementMesh() const;
  
  /// uses a floodfill on the highest-dimensional elements in the mesh to identify whether the model consists  of disconnected mesh patches
  bool IsContiguous() const;
  
  /// has the NodeManifoldManager been initialised which is true if the original model contained any split boundaries
  bool HasNodeManifolds() const;

  /// returns number of nodes=vertices in the current mesh
  size_t Nodes() const;

  /// returns number of elements in the current mesh
  size_t Elements() const;

  /// returns number of Faces=lower-dimensional elements in current mesh
  size_t Faces() const;

  /// returns number of InterFaces=lower-dimensional elements in current mesh
  size_t Interfaces() const;
  
  /// number of Manifold objects existing in conjuction with the split boundaries
  size_t NodeManifolds() const;
  
  /// order of shape functions used by the elements in the mesh; @return set to record potentially different levels of P-refinement
  std::set<uint32_t> OrderOfShapeFunctions() const;
  
  typename plf::colony<Node<dim> >::iterator      NodesBegin();
  typename plf::colony<Node<dim> >::iterator      NodesEnd();

  typename plf::colony<Element<dim> >::iterator   ElementsBegin();
  typename plf::colony<Element<dim> >::iterator   ElementsEnd();

  typename plf::colony<Face<dim> >::iterator      FacesBegin();
  typename plf::colony<Face<dim> >::iterator      FacesEnd();

  typename plf::colony<InterFace<dim> >::iterator InterfacesBegin();
  typename plf::colony<InterFace<dim> >::iterator InterfacesEnd();

  typename plf::colony<NodeManifold<dim> >::iterator NodeManifoldsBegin();
  typename plf::colony<NodeManifold<dim> >::iterator NodeManifoldsEnd();

  // const versions
  typename plf::colony<Node<dim> >::const_iterator      NodesBegin() const;
  typename plf::colony<Node<dim> >::const_iterator      NodesEnd() const;

  typename plf::colony<Element<dim> >::const_iterator   ElementsBegin() const;
  typename plf::colony<Element<dim> >::const_iterator   ElementsEnd() const;

  typename plf::colony<Face<dim> >::const_iterator      FacesBegin() const;
  typename plf::colony<Face<dim> >::const_iterator      FacesEnd() const;

  typename plf::colony<InterFace<dim> >::const_iterator InterfacesBegin() const;
  typename plf::colony<InterFace<dim> >::const_iterator InterfacesEnd() const;
  
  typename plf::colony<NodeManifold<dim> >::const_iterator NodeManifoldsBegin() const;
  typename plf::colony<NodeManifold<dim> >::const_iterator NodeManifoldsEnd() const;

  /// direct access for backward compatibility
  const FiniteElementManager& FiniteElements() const { return fem_manager_; }

  /// direct access for backward compatibility
  const FiniteVolumeStencilManager<dim>* const FiniteVolumes() const { return fvm_manager_; }

  ///  assigns the finite volume stencils to the finite volume policies of the element, face, and interface so that this functionality can be used
  void InitializeFiniteVolumeStencils( const PropertyDatabase<dim>&, bool assign_stencils_to_elements );
  
 
  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================
  
  /// replaces supplied lower-dimensional elements inside of a model with Face objects, establishing their connectivity; the input Elements are deleted
  std::vector<Face<dim>*>  ReplaceInteriorElementsByFaces( const PropertyDatabase<dim>&,
                                                           typename std::vector<Element<dim>*>::iterator first,
                                                           typename std::vector<Element<dim>*>::iterator last );

   /// replaces supplied lower-dimensional elements at the external boundary of a model with Face objects, establishing their connectivity; the input Elements are deleted
  std::vector<Face<dim>*>  ReplaceBoundaryElementsByFaces( const PropertyDatabase<dim>&,
                                                           typename std::vector<Element<dim>*>::iterator first,
                                                           typename std::vector<Element<dim>*>::iterator last );

 /// creates Face objects between face/node sharing Elements updating the connectivity; inside elements are first in pair
  std::vector<Face<dim>*>  CreateFacesBetweenNodeSharingElements( const PropertyDatabase<dim>&,
                                          const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,std::pair<Element<dim>*,uint32_t> > >& );

  /// replaces supplied Face objects with InterFace ones adding  necessary nodes and node manifolds, establishing new connectivity; the input Faces are deleted
  std::vector<InterFace<dim>*>  ReplaceFacesByInterFaces( const PropertyDatabase<dim>&,
                                                          typename std::vector<Face<dim>*>::iterator first,
                                                          typename std::vector<Face<dim>*>::iterator first_at_boundary,
                                                          typename std::vector<Face<dim>*>::iterator last,
                                                          typename std::vector<Node<dim>*>::const_iterator perim_first,
                                                          typename std::vector<Node<dim>*>::const_iterator perim_last );

  /// replaces supplied Face objects with InterFace ones adding  necessary nodes and node manifolds, establishing new connectivity; the input Faces are deleted
  std::vector<InterFace<dim>*>  ReplaceElementsByInterFaces( const PropertyDatabase<dim>&,
                                                             typename std::vector<FaceConstructionData<dim>>::iterator first,
                                                             typename std::vector<FaceConstructionData<dim>>::iterator last,
                                                             typename std::vector<Node<dim>*>::const_iterator perim_first,
                                                             typename std::vector<Node<dim>*>::const_iterator perim_last,
                                                             std::set<Node<dim>*> & split_perimeter_nodes,
                                                             std::set<size_t>& region_material_ids,
                                                             bool convert_original_elements_to_intervening_elements );

  /// creates InterFace objects between face/node sharing Elements adding the necessary nodes, node manifolds, and InterFace connectivity, updating overall connectivity as well; inside elements are first in pair
  std::vector<InterFace<dim>*>  CreateInterfacesBetweenNodeSharingElements( const PropertyDatabase<dim>&,
                                           const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,std::pair<Element<dim>*,uint32_t> > >&,
                                           bool multiplicate_perimeter_nodes,
                                           std::unordered_set<Element<dim>*>& outside_elmts );

  /// creates InterFace objects between face/node sharing Elements adding the necessary node manifolds and InterFace connectivity; inside elements are first in pair
  std::vector<InterFace<dim>*>  CreateInterfacesBetweenNodeMatchingElements( const PropertyDatabase<dim>&,
                                           const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,std::pair<Element<dim>*,uint32_t> > >& );

  /// by location only, no parent element  gets connected
  Node<dim>* const		 AddNodeAt( const Point<dim>&, const LocalVariables&,
                                  BOX_BOUNDARY = NOT, TOPOTYPE = MESH_VERTEX );

   /// duplicates Node, automatically creating a node manifold or adding it to an existing one; manifold type is established
  Node<dim>* const     Duplicate( Node<dim>* const nptr_inside, const LocalVariables& lvars );

  /// method tries to find neighbors through the parent connectivity of the nodes
  Element<dim>*	const  AddElement( CSMP_FEM_TYPE,
                                   const LocalVariables& element_variables,
                                   const IntegrationPointVariables& element_integration_point_variables,
                                   const std::vector<Node<dim>*>& nodes, int32_t material_id );

  /// puts lower-dimensional element inside of an InterFace, connecting it to its base pointer; the neighbors are not connected yet
  Element<dim>*	const AddInterveningElement( csmp::InterFace<dim>* const,
                                             const LocalVariables&,
                                             const IntegrationPointVariables&,
                                             const std::vector<Node<dim>*>& nodes,
                                             int32_t material_id );

  /// for the creation of SplitBoundaries consisting of InterFaces containing lower-dimensional elements
  InterFace<dim>* const WrapInterFaceAroundElement( csmp::Element<dim>* const eptr,
                                                    csmp::Element<dim>* inner_eptr,
                                                    csmp::Element<dim>* outer_eptr,
                                                    uint32_t adjacent_face_of_inner_element,
                                                    uint32_t adjacent_face_of_outer_element,
                                                    const LocalVariables& interface_variables,
                                                    const IntegrationPointVariables& interface__ipoint_vars );

 /// creates a Face matching the current lower-dimensional element and deletes the element subsequently; setting element pointer to zero if deletion was successful
  Face<dim>* const ReplaceElementByFace( typename std::vector<Element<dim>*>::iterator,
                                         csmp::Element<dim>* inner_eptr,
                                         csmp::Element<dim>* outer_eptr,
                                         uint32_t adjacent_face_of_inner_element,
                                         uint32_t adjacent_face_of_outer_element,
                                         const LocalVariables& face_variables,
                                         const IntegrationPointVariables& face_integration_point_variables,
                                         bool delete_original_face=true );

/// multiplicates nodes and replaces (deletes)  lower-dimensional Element with InterFace object, setting element pointer to zero if deletion was successful
InterFace<dim>* const ReplaceElementByInterFace( typename std::vector<Element<dim>*>::iterator,
                                                 csmp::Element<dim>* inner_eptr,
                                                 csmp::Element<dim>* outer_eptr,
                                                 uint32_t adjacent_face_of_inner_element,
                                                 uint32_t adjacent_face_of_outer_element,
                                                 const LocalVariables& interface_variables,
                                                 const IntegrationPointVariables& interface__ipoint_vars );

  /// creates face and assigns nodes and higher-dimensional neighbors @note the face neighbor element pointers are not assigned;
  Face<dim>* const AddFace( Element<dim>* const inner_parent, uint32_t inner_parent_face_id,
                            Element<dim>* const outer_parent, uint32_t outer_parent_face_id,
                            const LocalVariables&,
                            const IntegrationPointVariables& );

  /// creates a lower-dimensional face with parents of adjacent higher-dimensional Face objects as parents; deduces element type
  // TODO: this has to be replaced by Edge else we have a permanent inconsistency
  Face<dim>* const AddEdgeFace( Face<dim>* const adjacent_face1, uint32_t parent_elmt1_segm_id,
                                Face<dim>* const adjacent_face2, uint32_t parent_elmt2_segm_id,
                                const LocalVariables&,
                                const IntegrationPointVariables&,
                                const std::vector<Node<dim>*>& nodes );

  /// adds Face that caps a higher-dimensional Element at the model boundary
  Face<dim>* const AddBoundaryFace( csmp::Element<dim>* const innerParent,
                                    uint32_t local_face_id,
                                    const LocalVariables&,
                                    const IntegrationPointVariables& ); ///< optional

   /// assuming that the nodes on either side of the interface are already there, the face gets replaced
  InterFace<dim>* const ReplaceFaceByInterFace( typename std::vector<Face<dim>*>::iterator,
                                                const LocalVariables&,
                                                const IntegrationPointVariables&,
                                                std::vector<Node<dim>*> outside_nodes );
 
  /// for reconstrunction / creation of interfaces when the supplied elements are already disconnected from one-another, having separate nodes
  InterFace<dim>*	const	AddInterFace( Element<dim>* const inner_parent, uint32_t inner_element_face_id,
                                      Element<dim>* const outer_parent, uint32_t outer_element_face_id,
                                      const LocalVariables& interface_variables,
                                      const IntegrationPointVariables& interface_integration_point_variables );

  /// For connecting node-matched, node-sharing mesh patches, after the creating of the necessary outside nodes and manifolds
  InterFace<dim>*	const	AddInterFace( Element<dim>* const inner_parent, uint32_t inner_element_face_id,
                                      Element<dim>* const outer_parent, uint32_t outer_element_face_id,
                                      const LocalVariables& interface_variables,
                                      const IntegrationPointVariables& interface_integration_point_variables,
                                      std::vector<Node<dim>*> outside_nodes );

  // DELETIONS & MAINTANANCE OF MESH CONNECTIVITY
  // --------------------------------------------
  // NB: elements are responsible for their nodes, nodes for their manifolds
  
  /// disconnects neighbors, removes Element from node-parent container and deletes element, returns iterator to next Element if deletion succeeded
  auto Delete( typename std::vector<Node<dim>*>::iterator ) -> typename plf::colony< Node<dim> >::iterator;

  /// disconnects neighbors, removes Element from node-parent container and deletes element, returns iterator to next Element if deletion succeeded
  auto Delete( typename std::vector<Element<dim>*>::iterator ) -> typename plf::colony< Element<dim> >::iterator;

  /// disconnects neighbors, deletes Face, , returns iterator to next Face if deletion succeeded
  auto Delete( typename std::vector<Face<dim>*>::iterator ) -> typename plf::colony< Face<dim> >::iterator;

  /// disconnects neighbors, deletes InterFace, , returns iterator to next InterFace if deletion succeeded
  auto Delete( typename std::vector<InterFace<dim>*>::iterator ) -> typename plf::colony< InterFace<dim> >::iterator;
  
  /// updates all connectivity (elements, faces, interfaces, nodes to parents); rebulding node manifolds if necessary
  void UpdateConnectivity();

  /// reconnects InterFaces, disconnects adjacent elements,  updates node-parents, node-neighbors and node-manifolds; similar to UpdateConnectivity() but only for the surroundings of interfaces
  void UpdateConnectivity(  typename std::vector<InterFace<dim>*>::const_iterator first,
                            typename std::vector<InterFace<dim>*>::const_iterator last );

  /// reconnects Faces and adjacent Elements,  rebuilds node-parents for Face nodes
  void UpdateConnectivity(  typename std::vector<Face<dim>*>::const_iterator first,
                            typename std::vector<Face<dim>*>::const_iterator last );

  /// re-establishes the neighbor connectivity between cells of the same dimensionality (Elements & Faces)
  /// @todo disambiguate connectivity between Face and InterFace object at manifolds
  template<template<uint32_t> class CELL>
  void BuildConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                          typename std::vector<CELL<dim>*>::const_iterator last );

  template<template<uint32_t> class CELL>
  void BuildVolumeConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                typename std::vector<CELL<dim>*>::const_iterator last );
                                       
  template<template<uint32_t> class CELL>
  void BuildSurfaceConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                 typename std::vector<CELL<dim>*>::const_iterator last );
                                        
  template<template<uint32_t> class CELL>
  void BuildLineConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                              typename std::vector<CELL<dim>*>::const_iterator last );
                              
  /// specialisation of the method that finds matching faces by using manifold-pointer based keys rather than node-pointer based ones
  void BuildInterFaceConnectivity( typename std::vector<InterFace<dim>*>::const_iterator first,
                                   typename std::vector<InterFace<dim>*>::const_iterator last );

  ///  for nodes attached to elements in the supplied element range, the parent and the neighbor connectivity is reconstructed from scratch
  void ConnectNodesToParentsAndNeighbors( typename std::vector<Element<dim>*>::iterator first,
                                          typename std::vector<Element<dim>*>::iterator last );

  /// deletes elements and potentially orphaned nodes if any;  parent element storage of the nodes is rebuild and connectivity repaired;  input pointers are nulled
  size_t DeleteElementsAfterDisconnectingRemainingOnes( typename std::vector<Element<dim>*>::iterator first,
                                                        typename std::vector<Element<dim>*>::iterator last );

  /// disconnects face patch from potential adjacent faces before deleting faces; input pointers are nulled
  // TODO: not used: only when Faces are converted into interfaces
 size_t DeleteFacesAndRepairConnnectivity( typename std::vector<Face<dim>*>::iterator first,
                                           typename std::vector<Face<dim>*>::iterator last );

  /// disconnectes interfaces from not-targeted neighbors before deleting them;  does not remove multiplicated nodes or manifolds;  input pointers are nulled
  size_t DeleteInterfacesAndRepairConnnectivity( typename std::vector<InterFace<dim>*>::iterator first,
                                                 typename std::vector<InterFace<dim>*>::iterator last );

  /// disconnects nodes from potential manifolds, deletes the nodes as well as manifolds that no longer are
  size_t DeleteNodesAndRepairNodeConnnectivity( typename std::vector<Node<dim>*>::iterator first,
                                                typename std::vector<Node<dim>*>::iterator last );

  ///  Reports nodes that do not belong to any parent elements, faces or interfaces
  size_t OrphanNodes() const;
  std::vector<const Node<dim>*> OrphanNodeVector() const;

  /// relying on the parent element information of its nodes,  finds equidimensional neighbors of element face and connects itself with them and vice versa; returns number of neighbors found
  uint32_t ConnectNeighborsUsingNodeParents( Element<dim>* const );

  /// method to test the connectivity of a mesh
  size_t CheckElementConnectivity() const;

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


  /// disconnected neighbor elements from cell and itself from its neighbors
  template<template<uint32_t> class CELL>
  void DetachNeighborsFrom( CELL<dim>* const );

  /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
  void AssignUniqueNumbers( bool in_a_single_sequence=false ) const;
  
  /// puts nodes, elements, faces, and interfaces into the order given by Idx() variables; removes nullptr cells first
  void ReorderObjectsByIndexes();
  
  /// returns numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>&, bool get_indices_from_stored_variables=false );

  /// adds distributed variables to the VSet
  void OutputStoredVariablesTo( const PropertyDatabase<dim>&, VSet<dim>& ) const;
  
  /// reads distributed variables from VSet
  void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

private:

  FiniteElementManager              fem_manager_;
  FiniteVolumeStencilManager<dim>*  fvm_manager_ = nullptr;           ///<  finite volume specifications
  NodeManifoldManager<dim>*         node_manifold_manager_ = nullptr; ///<  node manifolds in case there are InterFace objects making up SplitBoundaries

  /// access is via root node or element only
  bool hybrid_element_mesh_;	///< true if the mesh consists of different FE types

  plf::colony<Node<dim>>      nodes_;          ///<  nodes
  plf::colony<Element<dim>>   elements_;       ///<  pointers to elements in the model
  plf::colony<Face<dim>>      faces_;          ///<  pointers faces making up the boundaries
  plf::colony<InterFace<dim>> interfaces_;     ///<  pointers to interfaces making up the split boundaries

  friend class MeshManager_Test; ///< so that private methods can be tested
  friend class Model<dim>;       ///<  exclusive access to private member functions
};



// NON-MEMBER FUNCTIONS

/// Adjusts halo cells and perimeter node parent-element vectors in mesh domains which is going to be deleted; this domain is identified by supplied range of iterators
template<uint32_t dim>
void updateHaloElementConnectivity( typename std::vector<Element<dim>*>::iterator first,
                                    typename std::vector<Element<dim>*>::iterator last );

/// for Face and InterFace regions
template<uint32_t dim, template<uint32_t> class CELL>
void updateHaloCellConnectivity( typename std::vector<CELL<dim>*>::iterator first,
                                 typename std::vector<CELL<dim>*>::iterator last );


/// finds cells sharing the same nodes and reports their numbers; verbose reports the duplicate cells
template<uint32_t dim, template<uint32_t> class CELL>
size_t detectDuplicateCells( typename plf::colony<CELL<dim>>::const_iterator begin,
                             typename plf::colony<CELL<dim>>::const_iterator end,
                             bool verbose );

/// Rotates node numbers while keeping the separation of nodes into corner nodes, midside nodes (quadratic elements)
template<uint32_t dim>
void reverseOrderOfFaceNodes( uint32_t n_corner_nodes_of_face,
                              typename std::vector<Node<dim>*>::iterator first,
                              typename std::vector<Node<dim>*>::iterator last );


} // end namespace csmp

#endif
