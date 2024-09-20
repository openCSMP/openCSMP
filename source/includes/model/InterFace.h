#ifndef CSMP_INTER_FACE_H
#define CSMP_INTER_FACE_H

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"


namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Face;


/**

@brief Surface (3D) or line (2D) element connector for elements which share faces but not nodes.

To connect disconnected patches of a mesh,
facilitating the application of tractions, boundary fluxes etc.
and for the representation of jump discontinuities of variables discretised on the nodes 
and interpolated with the element interpolation functions

@section conventions Conventions

For compatibility (usage conformance) with the Face and the Element, a range of methods are retained
which make sense only, however, if the side of the interface is specified. 
For this purpose,

current_side_

allows to set the state of the InterFace to INSIDE, OUTSIDE or MIDDLE; default is INSIDE.

If one wants to acces the nodes of an interface in the way one is used to from Element or Face,
the method 

InterFace::N( n_node, interface_side ) has to be used.

The nodes are ordered sequentially, INSIDE nodes first, then outside nodes.
The inside nodes follow the numbering of the nodes of the face of the higher-dimensional
element on the inside.
The outside nodes follow the face-node number of the OUTSIDE parent element.
This means that they are stored in the opposite order as those on the inside.

Inside and outside nodes start out to be matched (by point location).
Because of their different ordering (demanded to get the correct normals),
the first inside node matches the last outside node.
It follows that
@code
    auto n_nodes = FE()->Nodes() // of the finite element type of the interface.
    for ( auto i{0U}; i<n_nodes; i++ )
      assert( N(0,INSIDE)->Coordinate() == N(n_nodes-i-1U,OUTSIDE)->Coordinate() );
@endcode
should not fail.

The normal to the interface points from the INSIDE to the OUTSIDE element of the InterFace
and it is inherited from the lower dimensional element that the InterFace (or Face) was built from.
The boolean variable collocated_nodes, tells whether the nodes on opposite sides of the InterFace are
collocated, meaning that they have the same location (collocated_nodes=true) or whether they were
separated over the course of a simulation (collocated_nodes=false).
Default is true.

If the nodes are not collocated, the UnitNormal() method will return the normal to the mirror
symmetry plane between the 2 sides of the InterFace.

Nodes in the middle can only be accessed if the MIDDLE element pointer has been assigned
to a lower-dimensional intervening element inside of the InterFace.
If this is indeed the case, then access it via BaseElement().
The base (csmp::Element) of the InterFace is a nullptr by default, but it can be connected
to an intervening lower-dimensional element mesh once the InterFace has been created.

@note InterFaces are not listed as parents in their nodes.

@author SKM & Junchul Kim, refactored changing design to node-pointer based etc.
@author first version by P. Lang
@date 2011, 2014, 2016, 2019, 2022.

*/
template<uint32_t dim>
class InterFace : public FiniteElementPolicy<dim,InterFace>,
                  public FiniteVolumePolicy<dim,InterFace>, // TODO: what for? - flow parallel to interface on either side?
                  public LocalVariableStorage<dim,InterFace>
{
  public:

    // ------------------------------------------------------------------------
    // InterFace construction & connection with parents and neighbor Interfaces
    // ------------------------------------------------------------------------

    InterFace() = delete;

    /// New! SKM 29/7/2022: constructs  InterFace using the nodes and their numbering in the InterFace's higher-dimensional neighbors
    InterFace( csmp::Element<dim>&,
               csmp::Element<dim>* inner_parent,
               csmp::Element<dim>* outer_parent,
               uint32_t adjacent_face_of_inner_element,
               uint32_t adjacent_face_of_outer_element,
               const LocalVariables&  interface_props,
               const IntegrationPointVariables&  interface_integration_point_props );

    /// constructs complete InterFace using the Face nodes as inside nodes (in same order) and outside nodes supplied opposite order matching the inside nodes
    InterFace( csmp::Face<dim>*,
               const LocalVariables&  interface_props,
               const IntegrationPointVariables&  interface_integration_point_props,
               std::vector<Node<dim>*> outside_nodes );

    /// default: incomplete construction without connection to nodes
    InterFace( csmp::FiniteElement*, 
               const csmp::FiniteVolumeStencil<dim>*,
               const LocalVariables&  interface_props,
               const IntegrationPointVariables&  interface_integration_point_props );

    /// reconstruction of model from native binary file
    InterFace( size_t index,
               csmp::FiniteElement*,
               const csmp::FiniteVolumeStencil<dim>*,
               const LocalVariables&  interface_props,
               const IntegrationPointVariables&  interface_integration_point_props );

    InterFace( const InterFace<dim>& );
    InterFace( InterFace<dim>&& );
    ~InterFace() = default;
    InterFace& operator=( const InterFace<dim>& );
    InterFace& operator=( InterFace<dim>&& );
  
    /// self-detection in the interface construction process
    bool operator==( const InterFace<dim>& );

    /// reference to provide efficiency hack in MeshManager
    typename  std::vector<csmp::InterFace<dim>*>& NeighborElementVector();

    /// assigns nodes of  higher-dimensional neigbor elements to interface in the order they have in their corresponding faces
    void Assign( uint32_t n_local, Node<dim>*, INTERFACE_SIDE side );
  
    /// connects interFace to its higher-dimensional neighbors
    void Assign( Element<dim>* const inner_elmt, uint32_t inner_local_face_id,
                 Element<dim>* const outer_elmt, uint32_t outer_local_face_id );

    /// connects interFace to its higher-dimensional neighbors and nodes establishing connections by itself
    void AssignElementsAndNodes( Element<dim>* const inner_elmt,
                                 Element<dim>* const outer_elmt );

    /// connects InterFace to its higher-dimensional neighbor on the given side
    void Assign( Element<dim>* const parent, uint32_t faceId, INTERFACE_SIDE side );

    /// connects InterFace to its equidimensional neighbors (=number of InterFace finite-element faces)
    void Assign( uint32_t nbor, InterFace<dim>* const );
    
    /// sets neighbor pointer that was pointing to the argument object to 'nullptr'
    bool Unassign( const InterFace<dim>* const );
    /// removal of any InterFace attached to this neighbor slot
    void UnassignNeighbor( uint32_t nbor );
  
    /// connects interface to a lower-dimensional element with extra nodes situated inside the InterFace in a triple-layer mesh representation for fractures
    void Assign( Element<dim>* const intervening_elmt );


    // ------------------------------------------------------------------------
    //  Usage of InterFace objects in computations
    // ------------------------------------------------------------------------

    /// to apply visitors whose application level is Boundary and target is Face
    void Accept( csmp::Visitor<dim>& vis );

    /// of variables that are discretised on interfaces; required by LocalVariableStorage
    PLACEMENT Placement() const { return INTER_FACE; }

    /// node_connector_.size() = total nodes on both sides of InterFace
    uint32_t  Nodes() const;
    
    uint32_t  Faces() const;
  
    /// the InterFace object neighbors of the InterFace (one per face of interface)
    uint32_t  Neighbors() const;

    /// the InterFace object neighbors which are connected with the InterFace and not null.
    uint32_t  ConnectedNeighbors() const;

    /// only constant iterators are provided because the user is not supposed to change the node pr neighbor connectivity (done by MeshManager);
    typename std::vector<csmp::Node<dim>*>::const_iterator        NodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator        NodesEnd() const;

    typename std::vector<csmp::Node<dim>*>::const_iterator        CornerNodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator        CornerNodesEnd() const;

    typename std::vector<csmp::InterFace<dim>*>::const_iterator   NeighborsBegin() const;
    typename std::vector<csmp::InterFace<dim>*>::const_iterator   NeighborsEnd() const;

    /// access ONLY to the nodes on the Current_Side of the interface (determined by member current_side_)
    csmp::Node<dim>* const N( uint32_t n_local ) const;
    
    /// access the nodes that are connected to either, the inside or the outside of the Face
    csmp::Node<dim>* const N( uint32_t n_local, INTERFACE_SIDE side ) const;

    csmp::Node<dim>* const MatchingN( uint32_t n_local, INTERFACE_SIDE side ) const;

    /// switches internal state variable that sets interface side
    void            CurrentSide( INTERFACE_SIDE side );
    INTERFACE_SIDE  CurrentSide() const;
    
    /// END_POINT is a  classifier that applies on the perimeter of SplitBoundary objects terminating within models where INSIDE and OUTSIDE nodes are identical
    // TODO: review this functionality / adapt to manifolds
    bool IsEndPointNode( uint32_t n_local ) const;

    /// returns neighbor InterFace of interface
    csmp::InterFace<dim>* const Neighbor( uint32_t ) const;
    
    /// on-the-fly 0..n-1 numbering stored in mutable local variable and used for computations in interfaces (displacement gradients etc.)
    void           Idx( size_t ) const;
    size_t         Idx() const;

    /// access the higher dimensional elements on either side of interface; @attention returns nullptr if outside is not present
    Element<dim>* const Parent( INTERFACE_SIDE ) const;
  
    /// higher-dimensional element located on side opposite to where the unit normal points; will always be present
    Element<dim>* const InnerParent() const;
  
    /// higher-dimensional element located on the side of the interface to which the unit normal points; @attention does not exist on model boundary
    Element<dim>* const OuterParent() const;

    /// Local node numbers in higher-dimensional adjacent elements; costly to compute
    uint32_t       ParentNodeNumber( uint32_t n_local, INTERFACE_SIDE side ) const;

    /// access to the Element object from which the original face was created if it still is there (use HasBase()
    Element<dim>*  InterveningElement() const;

    /// is an equi-dimensional element connected to the MIDDLE element pointer of this InterFace
    bool           HasInterveningElement() const { return middleElement_!=nullptr; }
    
    /// tests whether the INSIDE nodes match the OUTSIDE nodes w.r.t. their position
    bool           AreNodesCollocated(double tolerance=std::numeric_limits<double>::epsilon()) const;
  
    /// local number of the face in the inner parent element, which borders against the interface
    void           ParentFaceID( INTERFACE_SIDE, uint32_t idx );
    uint32_t       InnerParentFaceID() const;
    uint32_t       OuterParentFaceID() const;
    uint32_t       ParentFaceID( INTERFACE_SIDE side ) const;

    // ------------------------------------------------------------------------
    //  InterFace geometric properties
    // ------------------------------------------------------------------------

    /// returns area of the interface; MIDDLE case is returned only if there is an intervening element
    double  Area( INTERFACE_SIDE=INSIDE ) const;
    
    /// returns the outward-pointing unit normal (from INSIDE to OUTSIDE)  using the middle element (if any) or inside node coordinates to construct face
    csmp::Point<dim>  UnitNormal( INTERFACE_SIDE side ) const;
  
    /// returns the outward-pointing unit normal (from INSIDE to OUTSIDE)  using the middle element (if any) or bisector node coordinates to construct face
    csmp::Point<dim>  UnitNormal() const;

    /// computes distance between corresponding pairs of nodes; @return false if nodes overlap, true if they are separated
    double NodeSpacing( uint32_t n_local ) const;

    /// returns a vector of the property of interest discretized on the node
    template<class Var>
    void    NodePropertyVector( const csmp::Index&, std::vector<Var>&, INTERFACE_SIDE=INSIDE ) const;

    /// returns a vector of the property of interest discretized on the node
    template<class Var>
    void    MatchingNodePropertyVector( const csmp::Index&, std::vector<Var>&, INTERFACE_SIDE ) const;


    /// the centre of gravity of the element (returns the mid-point of the 2-sides if detached)
    Point<dim>  BaryCenter() const;

    /// projects node points onto line returning max distance between them; vec direction can have any length
    double  LengthInDirection( const VectorVariable<dim>& vecDirection ) const;
    
    /// as needed by FiniteElementPolicy (uses current_side_ to retrieve matrix)
    void NodeCoordinateMatrix( DenseMatrix<DM_MIN>& ) const;

    /// as needed in the accumulation process in pde operators for a specifc side
    void NodeCoordinateMatrix( DenseMatrix<DM_MIN>&, INTERFACE_SIDE side ) const;


    // ------------------------------------------------------------------------
    // Screen Output
    // ------------------------------------------------------------------------

    void Out() const;
    
  protected:

    /// calculates the coordinate matrix from node coordinates representing the average of the inside and outside nodes of the interface
    void    BisectorCoordinateMatrix() const;

    /// calculates the coordinate matrix from node coordinates representing the average of the inside and outside nodes of the interface
    void    BisectorCoordinateMatrix(DenseMatrix<DM_MIN>& XY) const;

    /// finds the local numbers of the faces of the higher-dimensional element that will be connected by the interface; uses point coordinates that must be matched
    std::pair<uint32_t,uint32_t>  SharedElementFacesAndFaceIDs(); // TODO: replace with more generic functionality?
  
    /// connects the nodes of already connected higher dimensional neighbor elements to the InterFace
    void InitialiseNodeVector();

  private:
  
    /// for exclusive use by MeshManager
    template<uint32_t> friend class MeshManager;
    void* operator new( size_t size );
    void operator delete( void* p );

    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    Element<dim>*  innerParent_;                 ///< higher-dimensional parent element on side opposite to normal direction
    Element<dim>*  outerParent_;                 ///< higher-dimensional neighbor element in direction of interface normal
    Element<dim>*  middleElement_ = nullptr;     ///< pointer to lower-dimensional element that may have been inserted between the sides of InterFace.
    std::vector<Node<dim>*>       node_connector_;      ///< pointers to the nodes on inside followed by those on the outside
    std::vector<InterFace<dim>*>  interface_connector_; ///< neighbor interfaces; inside ones first, outside ones next, order determined by interface normal
    mutable size_t idx_;                                ///< unique identifier for indexing operations
    uint32_t       inner_parent_face_id_ = std::numeric_limits<uint32_t>::max(); ///< face number of inside higher-dimensional parent element
    uint32_t       outer_parent_face_id_ = std::numeric_limits<uint32_t>::max(); ///< face number of outside higher-dimensional parent element
    // used for compatibility with Element and Face methods (Neighbor etc.)
    mutable INTERFACE_SIDE current_side_;        ///< switch to return information from INSIDE, OUTSIDE or MIDDLE side of interface (default=INSIDE)

    friend class InterFace_Test; ///< friend declaration needed for the testing of private methods
};


} // csmp

#endif



