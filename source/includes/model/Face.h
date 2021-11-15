#ifndef CSMP_FACE_H
#define CSMP_FACE_H

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Visitor;
class FiniteElementManager;
template<size_t> class FiniteVolumeStencilManager;

/**
    Lower dimensional surface (3D) or line (2D) element that serves as interface (Face) or connector (InterFace)
    between higher dimensional mesh domains. The Face is used for material interfaces that are welded
    together. Boundary objects of a CSMP model consist of Face objects.
    
    The Face is rather similar in its functionality to the Element with many member functions
    sharing their names. However, Face objects know their higher-dimensional neighbors which are Elemen objects.
    This means that Face objects have an extra set of pointers that connect them to Element objects.
    
    Face objects also have a unit normal that helps with the application of tractions etc.
    
    Use Face objects for operations targeted on internal or external model boundaries.
    For external boundaries only the higher-dimensional neighbor 0 will be defined.
    For internal boundaries the Face provides access to both neighbors, distinguishing
    inside from outside neighbors. This distinction is made with regard to the Face normal.
    
    The Element objects that the Face normal points to are referred to as outside
    and the ones on the opposite side are the internal ones.
    This also means that the sense of node numbering of the Face matches that of
    the face of the Element on the inside.
    
    @author Stephan Matthai
    @date 3/3/2016
*/
template<size_t dim>
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
          size_t inner_parent_face_id,
          size_t outer_parent_face_id,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs model-boundary face as n-th (external) boundary face of the supplied higher dimensional parent element; no neighbor faces yet
    Face( Element<dim>& inner_parent,
          csmp::FiniteElement* FE_type_of_boundary_face,
          const FiniteVolumeStencilManager<dim>&,
          size_t n_boundary_face,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs face shared by the two volumetric elements inside of the model 
    Face( const FiniteElementManager&,
          const FiniteVolumeStencilManager<dim>&,
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          size_t inner_parent_face_id,
          size_t outer_parent_face_id,
          const LocalVariables&,
          const IntegrationPointVariables& );

   /// constructs face shared by the two volumetric elements inside of the model auto-detecting shared faces and nodes
    Face( const FiniteElementManager&,
          const FiniteVolumeStencilManager<dim>&,
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs model-edge line-element face connected with one or two Element objects that share their edge nodes with the Face on the model boundary
    Face( csmp::FiniteElement* FE_type_of_boundary_face,
          const FiniteVolumeStencilManager<dim>&,
          Element<dim>* const parent_of_face1,
          Element<dim>* const parent_of_face2,
          size_t parent_elmt1_segm_id,
          size_t parent_elmt2_segm_id,
          const std::vector<Node<dim>*>&  edge_nodes,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// for (RE)CONSTRUCTION of model from binary file; with storage but without connectivity
    Face( size_t index,
          csmp::FiniteElement*,
          const csmp::FiniteVolumeStencil<dim>*,
          const LocalVariables&,
          const IntegrationPointVariables& );

    Face( const Face& );
  
    /// hand-coded move constructor that is important since pointers need to be assigned
    Face( Face&& );

    ~Face();

    /// connects face to the supplied node
    void Assign( size_t node, Node<dim>* const );
    
    /// disconnecting the Node without deleting it; its pointer is set to nullptr
    void Unassign( csmp::Node<dim>* const );

    /// assign higher-dimensional neighbor elements to either side of face (outside is optional); needs nodes to be assigned first
    void Assign( Element<dim>* const innerElement, Element<dim>* const outerElement );
  
    /// tell face about its face neighbors
    void Assign( size_t nbor, Face<dim>* const );
      
    /// unassign its face neighbors
	  bool Unassign( Face<dim>* );

    /// @attention because of the pointers, this assignment makes sense only in the rarest cases
    Face& operator=( const Face<dim>& );

    /// hand-coded move assignment; important since pointers need to be assigned
    Face& operator=( Face<dim>&& );
  
    /// compares faces with one-another
    bool operator==( const Face<dim>& ) const;

    // TODO: poor design; rather tell Face what to do than taking over its functionality
    typename  std::vector<csmp::Face<dim>*>& NeighborElementVector();


    // ------------------------------------------------------------------------
    //  User interface of Face
    // ------------------------------------------------------------------------

    /// Local variable storage interface; required by LocalVariableStorage
    PLACEMENT Placement() const { return FACE; }

    size_t  Nodes() const;
    size_t  Neighbors() const;
	  size_t  ConnectedNeighbors() const;
    
    /// sides of Face object by analogy with Element
    size_t  Faces() const;

    /// to apply visitors whose application level is Boundary and target is Face
    void Accept( csmp::Visitor<dim>& );

    /// access the nodes that are connected to the Face
    csmp::Node<dim>*  N( size_t n_local ) const;
  
    /// access the neighbor faces of this face
    csmp::Face<dim>*  Neighbor( size_t ) const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void           Idx( size_t ) const;
    size_t         Idx() const;

    /// access the higher dimensional elements on either side of face; @attention returns nullptr if outside is not present
    Element<dim>*  Parent( INTERFACE_SIDE ) const;
  
    /// higher-dimensional element located on side opposite to where the unit normal points; will always be present
    Element<dim>*  InnerParent() const;
  
    /// higher-dimensional element located on the side of the face to which the unit normal points; @attention does not exist on model boundary
    Element<dim>*  OuterParent() const;
  
    /// returns which Face of the higher dimensional inner neighbor element this Face shares its nodes with
    /// local number of the face in the inner parent element, which borders against the interface
    void           ParentFaceID( INTERFACE_SIDE, size_t idx );
    size_t         InnerParentFaceID() const;
    size_t         OuterParentFaceID() const;
    size_t         ParentFaceID( INTERFACE_SIDE side ) const;
  
    /// returns the number of the desired node in the inner parent element of the Face
    size_t         ParentNodeNumber( size_t n_local ) const;

    // ------------------------------------------------------------------------
    //  Face geometry operations
    // ------------------------------------------------------------------------
  
    /// returns area of the face; method assumes same role as Volume() for the element
    double       Area() const;
  
    /// not a face-normal vector, but the shortest path between the barycenters of face and element
    void           VectorToInnerElementBaryCenter( VectorVariable<dim>& ) const;
  
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
    template<size_t> friend class MeshManager;
    void* operator new( size_t size );
    void operator delete( void* p );
  
    // TODO: SKM: deprecate this inefficient method
    void AssignFaceID( const std::vector<Node<dim>*>& faceNodes );

    // TODO: move this method to unit test: void CheckNodeOrderingAccordingToUnitNormalOrientation();
  
    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t           idx_;
    std::vector<Node<dim>*>  node_connector_;     ///< pointers to the nodes of the face
    std::vector<Face<dim>*>  face_connector_;     ///< the (equidimensional) neighbors of the face
    // not references or constant pointers because these may need to change during remeshing
    Element<dim>*            innerParent_;        ///< higher-dimensional neighbor in opposite direction of unit normal (always there)
    Element<dim>*            outerParent_;        ///< higher-dimensional neighbor element in direction of interface normal
    size_t inner_parent_face_id_ = UNSPECIFIED;   ///< face number of inside higher-dimensional parent element, segm id if Face is line element in 3D
    size_t outer_parent_face_id_ = UNSPECIFIED;   ///< face number of outside higher-dimensional parent element, segm id if Face is line element in 3D
};

} // csmp

#endif



