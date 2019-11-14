#ifndef CSMP_INTER_FACE_H
#define CSMP_INTER_FACE_H

#include "Element.h"

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"


namespace csmp {

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

Else the nodes are ordered sequentially, INSIDE nodes first, then outside nodes.
The nodes in the middle can only be accessed via BaseElement().

The normal to the interface points from the INSIDE to the OUTSIDE element of the InterFace
and it is inherited from the lower dimensional element that the InterFace (or Face) was built from.
The boolean variable collocated_nodes, tells whether the nodes on opposite sides of the InterFace are
collocated, meaning that they have the same location (collocated_nodes=true) or whether they were
separated over the course of a simulation (collocated_nodes=false).
Default is true.

If the nodes are not collocated, the UnitNormal() method will return the normal to the mirror
symmetry plane between the 2 sides of the InterFace.

The base (csmp::Element) of the InterFace is a nullptr by default, but it can be connected
to an intervening lower-dimensional element mesh once the InterFace has been created.

@note InterFaces are not listed as parents in their nodes.

@author SKM & Junchul Kim, refactored changing design to node-pointer based etc.
@author first version by P. Lang
@date 2011, 2014, 2016, 2019.

*/
template<size_t dim>
class InterFace : public FiniteElementPolicy<dim,InterFace>,
                  public FiniteVolumePolicy<dim,InterFace>,
                  public LocalVariableStorage<dim,InterFace<dim> >
{
  public:

    // ------------------------------------------------------------------------
    // InterFace construction & connection with parents and neighbor Interfaces
    // ------------------------------------------------------------------------

    InterFace() = delete;

    explicit InterFace( FiniteElement* );

    InterFace( FiniteElement*,
               FiniteVolumeStencil<dim>* );
    
    InterFace( csmp::FiniteElement*, 
               csmp::FiniteVolumeStencil<dim>*, 
               const LocalVariables& interface_props, 
               const IntegrationPointVariables& integration_point_props );

    /// for reconstruction of model from native binary file
    InterFace( size_t index,
               csmp::FiniteElement*,
               const LocalVariables& interface_props,
               const IntegrationPointVariables& integration_point_props );

    InterFace( const InterFace& );
  
    /// handcoded move contructor to deal with the pointers
    InterFace( InterFace&& );
  
    ~InterFace();

    /// connect InterFace to its equidimensional neighbors
    void Assign( size_t nbor, InterFace<dim>* const );
  
    /// connect interface to the element that it was created from
    void Assign( Element<dim>* const base_elmt );

    /// connect interface to its higher-dimensional neighbors
    void Assign( Element<dim>* const inner_elmt, Element<dim>* const outer_elmt, bool assign_nodes = true );
  
    /// as Assign, for the case that the shared faces are already known
    void Assign( Element<dim>* const inner_elmt, size_t inner_local_face_id,
                 Element<dim>* const outer_elmt, size_t outer_local_face_id,
                 bool assign_nodes );

    /// connect interface to its higher-dimensional neighbor on the given side
    void Assign( Element<dim>* const parent, size_t faceId, INTERFACE_SIDE side );

    void Assign( size_t n_local, size_t parent_node, INTERFACE_SIDE side );
    void Assign( size_t n_local, std::pair<size_t, size_t> parent_nodes );

    /// assign the corresponding node of the higher-dimensional neigbor element
    void Assign( size_t n_local, Node<dim>*, INTERFACE_SIDE side );
  
    /// TODO: just provide iterators to begin() and end(), or, still better, just tell Interface which neighbors to disconnect
    typename  std::vector<csmp::InterFace<dim>*>& NeighborElementVector();

	  /// unassigns interface neighbors
	  bool DisconnectNeighbor( InterFace<dim>* );

    InterFace& operator=( const InterFace& );

    /// hand-coded move assignment that deals with the pointers
    InterFace& operator=( InterFace&& );

    /// self-detection in the interface construction process
    bool operator==( const InterFace<dim>& );


    // ------------------------------------------------------------------------
    //  Usage of InterFace objects in computations
    // ------------------------------------------------------------------------

    /// to apply visitors whose application level is Boundary and target is Face
    void Accept( csmp::Visitor<dim>& vis );

    /// of variables that are discretised on interfaces
    PLACEMENT Placement() const { return INTER_FACE; }

    /// node_connector_.size() = total nodes on both sides of InterFace
    size_t  Nodes() const;
    
    size_t  Faces() const;
  
    /// the InterFace object neighbors of the InterFace (one per face of interface)
    size_t  Neighbors() const;

    /// the InterFace object neighbors which are connected with the InterFace and not null.
    size_t  ConnectedNeighbors() const;

    /// access to all nodes connected to the InterFace (inside nodes first)
    csmp::Node<dim>* N( size_t n_local ) const;
    
    /// access the nodes that are connected to either, the inside or the outside of the Face
    csmp::Node<dim>* N( size_t n_local, INTERFACE_SIDE side ) const;

    /// switches internal state variable that sets interface side
    void            CurrentSide( INTERFACE_SIDE side );
    INTERFACE_SIDE  CurrentSide() const;

    /// standard interface returns neighbor on current side of interface
    csmp::InterFace<dim>*  Neighbor( size_t ) const;
    
    /// neighbor on choosen side of interface is returned
    csmp::InterFace<dim>*  Neighbor( size_t, INTERFACE_SIDE ) const;

    /// on-the-fly 0..n-1 numbering stored in mutable local variable and used for computations in interfaces (displacement gradients etc.)
    void           Idx( size_t ) const;
    size_t         Idx() const;

    /// access the higher dimensional elements on either side of interface; @attention returns nullptr if outside is not present
    Element<dim>*  Parent( INTERFACE_SIDE ) const;
  
    /// higher-dimensional element located on side opposite to where the unit normal points; will always be present
    Element<dim>*  InnerParent() const;
  
    /// higher-dimensional element located on the side of the interface to which the unit normal points; @attention does not exist on model boundary
    Element<dim>*  OuterParent() const;

    /// access to the Element object from which the original face was created if it still is there (use HasBase()
    Element<dim>*  BaseElement() const;

    /// Local node numbers in higher-dimensional adjacent elements; costly to compute
    size_t         ParentNodeNumber( size_t n_local, INTERFACE_SIDE side ) const;

    /// is an equi-dimensional element connected to the MIDDLE element pointer of this InterFace
    bool           HasBase() const { return baseElement_!=nullptr; }
  
    // TODO: review this functionality, included SharedElementFaces()
    /// local number of the face in the inner parent element, which borders against the interface
    size_t         InnerParentFaceID() const;
    size_t         OuterParentFaceID() const;
    size_t         ParentFaceID( INTERFACE_SIDE side ) const;

    // ------------------------------------------------------------------------
    //  InterFace geometric properties
    // ------------------------------------------------------------------------

    /// returns area of the face; method assumes same role as Volume() for the element
    double64       Area( INTERFACE_SIDE=MIDDLE ) const;
    
    /// unit normals on either side point from INSIDE to OUTSIDE, but have different orientation when nodes are spatially separated 
    void           UnitNormal( VectorVariable<dim>&, INTERFACE_SIDE side ) const;
  
    /// returns normal to 'current' side of interface
    void           UnitNormal( VectorVariable<dim>& ) const;
  
    /// returns the normal pointing from the inside to the outside higher-dimensional Element of the InterFace, calculated for bisector plane.
    csmp::Point<dim>  UnitNormal() const;

    /// computes distance between corresponding pairs of nodes; @return false if nodes overlap, true if they are separated
    // TODO: review this functionality
    bool           NodeSpacing( size_t n_local, VectorVariable<dim>& ) const;

    /// returns a vector of the property of interest discretized on the node
    template<class Var>
    void           NodePropertyVector( const csmp::Index&, std::vector<Var>&, INTERFACE_SIDE=INSIDE ) const;

    /// inputs node coordinates into supplied matrix; for MIDDLE the midpoints between the nodes on either side are used
    void           NodeCoordinateMatrix( DenseMatrix<DM_MIN>&, INTERFACE_SIDE=MIDDLE ) const;

    /// the centre of gravity of the element (returns the mid-point of the 2-sides if detached)
    Point<dim>     BaryCenter() const;

    /// projects node points onto line returning max distance between them; vec direction can have any length
    double64       LengthInDirection( const VectorVariable<dim>& vecDirection ) const;

    // ------------------------------------------------------------------------
    // Screen Output
    // ------------------------------------------------------------------------

    void Out() const;

  private:
    /// finds the local numbers of the faces of the higher-dimensional element that will be connected by the interface; uses point coordinates that must be matched
    std::pair<size_t,size_t>  SharedElementFaces();
  
    /// connect the nodes of the higher dimensional neighbor elements to the InterFace; @note can also be done individually with Assign
    void InitializeNodeVector();

    /// as above when the local indices of the shared faces of the higher-dimensional elements adjacent to the face are already known
    void InitializeNodeVector( size_t inner_face_ID, size_t outer_face_ID );

    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t                idx_;                 ///< unique identifier for indexing operations
    std::vector<Node<dim>*>       node_connector_;      ///< pointers to the nodes on inside followed by those on the outside
    std::vector<InterFace<dim>*>  interface_connector_; ///< neighbor interfaces on inside followed by those on outside
    bool                          collocated_nodes_;    ///< nodes on both sides of InterFace are co-located = default
  
    // used for compatibility with Element and Face methods (Neighbor etc.)
    INTERFACE_SIDE  current_side_;  ///< switch to return information from INSIDE, OUTSIDE or MIDDLE side of interface (default=INSIDE)

    // pointers to the higher-dimensional elements this face sits in between
    // the (same dimensional) neighbors are stored by the base class
    Element<dim>* innerParent_;            ///< higher-dimensional parent element on side opposite to normal direction
    size_t        inner_parent_face_id_;   ///< face number of inside higher-dimensional parent element

    Element<dim>* outerParent_;            ///< higher-dimensional neighbor element in direction of interface normal
    size_t        outer_parent_face_id_;   ///< face number of outside higher-dimensional parent element

    Element<dim>* baseElement_ = nullptr;  ///< pointer to lower-dimensional element that sits between sides of InterFace.
};


} // csmp

#endif



