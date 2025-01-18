#ifndef CSMP_FACE_H
#define CSMP_FACE_H

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Visitor;
class FiniteElementManager;
template<uint32_t> class FiniteVolumeStencilManager;

/**
    Lower dimensional surface (in 3D) or line (2D) elements  serve as boundaries (Face) or connect disconnected mesh domains (InterFace).
    The Face is used for material interfaces that are welded together by sharing nodes. Boundary objects of a CSMP model consist of Face objects.
    
    The Face is rather similar in its functionality to the Element with many member functions
    sharing their names. As an important difference, Face objects know their higher-dimensional neighbor Element objects.
    This means that Face objects are not only connected to each other, but have an extra set of pointers that connect them to Element objects.
    
    Face objects also have a unit normal that helps with the application of tractions etc.
    Since the Face is not inherited from Element, it constitutes a separate type and can store unique variables, especially ones like tractions,
    influxes or transfer coefficients, which would unnecessarily consume storage if discretised evenly over a domain.
    Element variables like permeability can still be accessed easily because the Face has pointers to its inner and outer Element neighbors.
    
    Face objects have been designed to support operations on internal or external model boundaries.
    If a Face is located at an external (Model) boundary its normal will point outward and only the inner higher-dimensional neighbor 0 will be initialised.
    Representing internal boundaries, Face objects provide access to both neighbors, also helping to distinguish inside from outside domains.
    This distinction is possible and indicated by to the Face normal that points to the outside.
    Thus, the Element objects that the Face normal points to are referred to as outside and the ones on the opposite side are the internal ones.
    This also means that the sense of node numbering of the Face matches that of
    the face of the Element on the inside.
    
    @note without a topology it is impossible to "compute" what the inside or outside of a surface in a boundary representation is. This info
    comes from CAD or geomodel into CSMP++ which can only perform consistency checks.
    
    @attention Since so much information is required to construct complete Face objects, this is done in 2 steps:
    1) the various constructors establish  Face connectivity with higher-dimensional neighbors, associated indices and policies.
    2) separate operations assign idx, nodes, neighbors, etc.
    
    @author Stephan Matthai
    @date 3/3/2016
*/
template<uint32_t dim>
class Face : public FiniteElementPolicy<dim,Face>,
             public FiniteVolumePolicy<dim,Face>,
             public LocalVariableStorage<dim,Face>
{
  public:
  
    Face() = delete;

    // ------------------------------------------------------------------------
    // Functionality used in Face construction process (in that order)
    // ------------------------------------------------------------------------

    /// constructs model-interior face as an exact copy of the supplied lower-dimensional element; no neighbor faces yet
    Face( const Element<dim>& dim_minus1_element, ///< supplies finite element policy & finite volume stencil information
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          uint32_t inner_parent_face_id,
          uint32_t outer_parent_face_id,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs model-boundary face as n-th (external) boundary face of the supplied higher dimensional parent element; no neighbor faces yet
    Face( Element<dim>& inner_parent,
          csmp::FiniteElement* FE_type_of_boundary_face,
          const FiniteVolumeStencilManager<dim>*,
          uint32_t n_boundary_face,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs face shared by the two volumetric elements inside of the model 
    Face( const FiniteElementManager&,
          const FiniteVolumeStencilManager<dim>*,
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          uint32_t inner_parent_face_id,
          uint32_t outer_parent_face_id,
          const LocalVariables&,
          const IntegrationPointVariables& );

   /// constructs face shared by the two volumetric elements inside of the model auto-detecting shared faces and nodes
    Face( const FiniteElementManager&,
          const FiniteVolumeStencilManager<dim>*,
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs model-edge line-element face connected with one or two Element objects that share their edge nodes with the Face on the model boundary
    Face( csmp::FiniteElement* FE_type_of_boundary_face,
          const FiniteVolumeStencilManager<dim>*,
          Element<dim>* const parent_of_face1,
          Element<dim>* const parent_of_face2,
          uint32_t parent_elmt1_segm_id,
          uint32_t parent_elmt2_segm_id,
          const std::vector<Node<dim>*>&  edge_nodes,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// for (RE)CONSTRUCTION of model from binary file; with storage but without connectivity
    Face( size_t index,
          csmp::FiniteElement*,
          const csmp::FiniteVolumeStencil<dim>*,
          const LocalVariables&,
          const IntegrationPointVariables& );
          
    Face( const Face<dim>& );
    Face<dim>& operator=( const Face<dim>& );
    // move semantics for colony
    Face( Face<dim>&& );
    Face<dim>& operator=( Face<dim>&& );
   
    ~Face() = default;

    /// connects face to the supplied node
    void Assign( uint32_t node, Node<dim>* const );
    
    /// disconnecting the Node without deleting it; its pointer is set to nullptr
    void Unassign( const csmp::Node<dim>* const );

    /// assign higher-dimensional neighbor elements to either side of face (outside is optional); needs nodes to be assigned first
    void Assign( Element<dim>* const innerElement, Element<dim>* const outerElement );
  
    /// tell face about its face neighbors
    void Assign( uint32_t nbor, Face<dim>* const );
      
    /// sets the pointer to given neighbor face to zero, reports whether removal was made
    bool Unassign( const Face<dim>* const );
    /// removal of anything attached to face neighbor
	  void UnassignNeighbor( uint32_t nbor );
  
    /// compares faces with one-another
    bool operator==( const Face<dim>& ) const;

    /// for use in range loops
    typename  std::vector<csmp::Face<dim>*>& NeighborElementVector();


    // ------------------------------------------------------------------------
    //  User interface of Face
    // ------------------------------------------------------------------------

    /// Local variable storage interface; required by LocalVariableStorage
    PLACEMENT Placement() const { return FACE; }

    uint32_t  Nodes() const;
    uint32_t  Neighbors() const;
	  uint32_t  ConnectedNeighbors() const;
    
    /// sides of Face object by analogy with Element
    uint32_t  Faces() const;

    /// only constant iterators are provided because the user is not supposed to change the node pr neighbor connectivity (done by MeshManager)
    typename std::vector<csmp::Node<dim>*>::const_iterator   NodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator   NodesEnd() const;

    typename std::vector<csmp::Node<dim>*>::const_iterator   CornerNodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator   CornerNodesEnd() const;

    typename std::vector<csmp::Face<dim>*>::const_iterator   NeighborsBegin() const;
    typename std::vector<csmp::Face<dim>*>::const_iterator   NeighborsEnd() const;

    /// to apply visitors whose application level is Boundary and target is Face
    void Accept( csmp::Visitor<dim>& );

    /// access the nodes that are connected to the Face
    csmp::Node<dim>* const N( uint32_t n_local ) const;
  
    /// access the neighbor faces of this face
    csmp::Face<dim>* const Neighbor( uint32_t ) const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void           Idx( size_t ) const;
    size_t         Idx() const;

    /// access the higher dimensional elements on either side of face; @attention returns nullptr if outside is not present
    Element<dim>* const Parent( INTERFACE_SIDE ) const;
  
    /// higher-dimensional element located on side opposite to where the unit normal points; will always be present
    Element<dim>* const InnerParent() const;
  
    /// higher-dimensional element located on the side of the face to which the unit normal points; @attention does not exist on model boundary
    Element<dim>* const OuterParent() const;
  
    /// returns which Face of the higher dimensional inner neighbor element this Face shares its nodes with
    /// local number of the face in the inner parent element, which borders against the interface
    void           ParentFaceID( INTERFACE_SIDE, uint32_t idx );
    uint32_t       InnerParentFaceID() const;
    uint32_t       OuterParentFaceID() const;
    uint32_t       ParentFaceID( INTERFACE_SIDE side ) const;

    // ------------------------------------------------------------------------
    //  Face geometry operations
    // ------------------------------------------------------------------------
  
    /// returns area of the face; method assumes same role as Volume() for the element
    double       Area() const;
  
    /// not a face-normal vector, but the shortest path between the barycenters of face and element
    void         VectorToInnerElementBaryCenter( VectorVariable<dim>& ) const;
  
    // unit normal computations for Face are handled by its FiniteElementPolicy the options are
    // Point<dim> UnitNormal() const;
    // void       UnitNormal( std::vector<double>& nrml ) const;
    // void       UnitNormal( VectorVariable<dim>& nrml ) const;

    /// returns a vector of the property of interest discretized on the node
    template<class Var>
    void        NodePropertyVector( const csmp::Index&, std::vector<Var>& ) const;

    /// inputs node coordinates into supplied matrix
    void        NodeCoordinateMatrix( DenseMatrix<DM_MIN>& ) const;

    /// the centre of gravity of the element
    Point<dim>  BaryCenter() const;

    /// projects node points onto line returning max distance between them; vec direction can have any length
    double    LengthInDirection( const VectorVariable<dim>& vecDirection ) const;

    /// prints state of this object
    void  Out() const;

  private:
    
    /// for exclusive use by MeshManager
    template<uint32_t> friend class MeshManager;
    void* operator new( size_t size );
    void operator delete( void* p );
  
    // TODO: SKM: deprecate this inefficient method
    void AssignFaceID( const std::vector<Node<dim>*>& faceNodes );
  
    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t           idx_;
    uint32_t                 inner_parent_face_id_ = std::numeric_limits<uint32_t>::max(); ///< face number of inside higher-dimensional parent element, segm id if Face is line element in 3D
    uint32_t                 outer_parent_face_id_ = std::numeric_limits<uint32_t>::max(); ///< face number of outside higher-dimensional parent element, segm id if Face is line element in 3D
    Element<dim>*            innerParent_;        ///< higher-dimensional neighbor in opposite direction of unit normal (always there)
    Element<dim>*            outerParent_;        ///< higher-dimensional neighbor element in direction of interface normal
    std::vector<Node<dim>*>  node_connector_;     ///< pointers to the nodes of the face
    std::vector<Face<dim>*>  face_connector_;     ///< the (equidimensional) neighbors of the face
};

} // csmp

#endif



