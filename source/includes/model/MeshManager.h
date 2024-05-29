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
    @brief Helper class of the Model which takes care of the storage of Element, Face and InterFace objects;
    internal application is hidden and may vary between models (tree-storage is default).

    @author S.K. Matthai
    @date 2021 (complete rewrite)

    @remark gain access via Mesh() public interface of Model.

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
  size_t InterFaces() const;
  
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
                                                          typename std::vector<Node<dim>*>::const_iterator perim_last);

  /// replaces supplied Face objects with InterFace ones adding  necessary nodes and node manifolds, establishing new connectivity; the input Faces are deleted
  std::vector<InterFace<dim>*>  ReplaceElementsByInterFaces( const PropertyDatabase<dim>&,
                                                          typename std::vector<FaceConstructionData<dim>>::iterator first,
                                                          typename std::vector<FaceConstructionData<dim>>::iterator last,
                                                          typename std::vector<Node<dim>*>::const_iterator perim_first,
                                                          typename std::vector<Node<dim>*>::const_iterator perim_last,
                                                          std::set<Node<dim>*> & split_perimeter_nodes,
                                                          std::set<size_t>& region_material_ids);


  /// creates InterFace objects between face/node sharing Elements adding the necessary nodes, node manifolds, and InterFace connectivity, updating overall connectivity as well; inside elements are first in pair
  std::vector<InterFace<dim>*>  CreateInterfacesBetweenNodeSharingElements( const PropertyDatabase<dim>&,
                                           const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,std::pair<Element<dim>*,uint32_t> > >&,
                                           bool multiplicate_perimeter_nodes,
                                           Region<dim>& out_region );

  /// creates InterFace objects between face/node sharing Elements adding the necessary node manifolds and InterFace connectivity; inside elements are first in pair
  std::vector<InterFace<dim>*>  CreateInterfacesBetweenNodeMatchingElements( const PropertyDatabase<dim>&,
                                           const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,std::pair<Element<dim>*,uint32_t> > >& );

  /// by location only, no parent element  gets connected
  Node<dim>* const		 AddNodeAt( const Point<dim>&, const LocalVariables&,
                                  BOX_BOUNDARY = NOT, TOPOTYPE = MESH_VERTEX );

  /// only if there is not already a node at this location, else a pointer to that node is returned, no parent element  gets connected
  Node<dim>* const		 AddNodeAtUniqueLocation( const Point<dim>&, size_t nearby_node,
                                                const LocalVariables& node_variables,
                                                BOX_BOUNDARY = NOT,
                                                TOPOTYPE = MESH_VERTEX );

   /// duplicates Node, automatically creating a node manifold or adding it to an existing one; manifold type is established
  Node<dim>* const     Duplicate( Node<dim>* const nptr_inside, const LocalVariables& lvars );

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
                                         uint32_t adjacent_face_of_inner_element,
                                         uint32_t adjacent_face_of_outer_element,
                                         const LocalVariables& face_variables,
                                         const IntegrationPointVariables& face_integration_point_variables,
                                         bool delete_original_face=true );

/// multiplicates nodes and replaces (deletes)  lower-dimensional Element with InterFace object
InterFace<dim>* const ReplaceElementByInterFace( csmp::Element<dim>* eptr,
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
  InterFace<dim>* const ReplaceFaceByInterFace( csmp::Face<dim>* eptr,
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

  /// updates all connectivity (elements, faces, interfaces, nodes to parents); however, node manifolds are not reconstructed
  void UpdateConnectivity();

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
  size_t DeleteCellsAndRepairConnnectivity( typename std::vector<Element<dim>*>::iterator first,
                                            typename std::vector<Element<dim>*>::iterator last );

  /// disconnects face patch from potential adjacent faces before deleting faces; input pointers are nulled
  size_t DeleteCellsAndRepairConnnectivity( typename std::vector<Face<dim>*>::iterator first,
                                            typename std::vector<Face<dim>*>::iterator last );

  /// disconnectes interfaces from not-targeted neighbors before deleting them;  does not remove multiplicated nodes or manifolds;  input pointers are nulled
  size_t DeleteCellsAndRepairConnnectivity( typename std::vector<InterFace<dim>*>::iterator first,
                                            typename std::vector<InterFace<dim>*>::iterator last );

  /// disconnects nodes from potential manifolds and deletes the latter
  size_t DeleteCellsAndRepairConnnectivity( typename std::vector<Node<dim>*>::iterator first,
                                            typename std::vector<Node<dim>*>::iterator last );

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

// POTENTIAL METHODS?


   /// replaces face, constructing new nodes & manifolds where indicated by vector
//  InterFace<dim>* const ReplaceFaceByInterFace( csmp::Face<dim>* eptr,
//                                                const std::vector<bool>&  nodes_to_duplicate,
//                                                const LocalVariables&,
//                                                const IntegrationPointVariables& );



} // end namespace csmp

#endif
