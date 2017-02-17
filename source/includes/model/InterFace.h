#ifndef CSMP_INTER_FACE_H
#define CSMP_INTER_FACE_H

#include "Element.h"

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"
#include "InterFaceRemeshingTraits.h"


namespace csmp {

/// surface (3D) or line (2D) element connector for elements which share faces but not nodes
template<size_t dim>
class InterFace : public InterFaceRemeshingTraits<dim,InterFace>,
                  public FiniteElementPolicy<dim,InterFace>,
                  public FiniteVolumePolicy<dim,InterFace>,
                  public LocalVariableStorage<dim,InterFace<dim> >
{

  public:

    // ------------------------------------------------------------------------
    // Functionality used in Face construction process
    // ------------------------------------------------------------------------

    InterFace();

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
//    InterFace( InterFace&& );
    ~InterFace();

    /// connect InterFace to its equidimensional neighbors
    void Assign( size_t nbor, InterFace<dim>* const );
  
    /// connect interface to the element that it was created from
    void Assign( Element<dim>* const base_elmt );

    /// connect interface to its higher-dimensional neighbors
    void Assign( Element<dim>* const inner_elmt, Element<dim>* const outer_elmt, bool assign_nodes = true );

    /// connect interface to its higher-dimensional neighbors
    void Assign( Element<dim>* const parent, size_t faceId, INTERFACE_SIDE side );

    void Assign( size_t n_local, size_t parent_node, INTERFACE_SIDE side );
    void Assign( size_t n_local, std::pair<size_t, size_t> parent_nodes );

    // TODO: SKM FIX - assign the actual node
    void Assign( size_t n_local, Node<dim>*, INTERFACE_SIDE side );
  
    // ------------------------------------------------------------------------
    // Basic information
    // ------------------------------------------------------------------------

    InterFace& operator=( const InterFace& );

    /// relation operators
    bool operator==( const InterFace<dim>& );
  
    size_t  Nodes() const     { return node_connector_.size() / 2U; }; // since there is a duplicate set
    size_t  Neighbors() const { return interface_connector_.size(); };
    
    /// for element face, there can be a neighbor
    size_t  Faces() const { return interface_connector_.size(); };

    /// Local variable storage interface
    PLACEMENT Placement() const { return INTER_FACE; }
  
    // ------------------------------------------------------------------------
    // Member access
    // ------------------------------------------------------------------------

    /// to apply visitors whose application level is Boundary and target is Face
    void Accept( csmp::Visitor<dim>& vis );

    /// helper method for remeshing purposes; @todo move to remeshing policy
    typename  std::vector<csmp::InterFace<dim>*>& NeighborElementVector();

    /// access the nodes that are connected to inside of the Face
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

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
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

    /// Node numbers correspondance
    size_t         ParentNodeNumber( size_t n_local, INTERFACE_SIDE side ) const;
    size_t         ParentNodeNumberOppositeTo( size_t n_from_parent, INTERFACE_SIDE side ) const;


    // TODO: check and capitalize
    bool           HasBase() const { return baseElement_!=nullptr; }
    // TODO: review this functionality
    size_t         InnerParentFaceID() const;
    size_t         OuterParentFaceID() const;
    size_t         ParentFaceID( INTERFACE_SIDE side ) const;

    // ------------------------------------------------------------------------
    // Geometry
    // ------------------------------------------------------------------------

    /// returns area of the face; method assumes same role as Volume() for the element
    double64       Area() const;
    void           UnitNormal( VectorVariable<dim>&, INTERFACE_SIDE side ) const;
    void           UnitNormal( VectorVariable<dim>& )         const;

    /// computes spacing between pairs of nodes; @return false if overlap, true if separated
    bool           NodeSpacing( size_t n_local, VectorVariable<dim>& ) const;

    // ------------------------------------------------------------------------
    // Functionality
    // ------------------------------------------------------------------------

    /// returns a vector of the property of interest discretized on the node
    template<class Var>
    void           NodePropertyVector( const csmp::Index&, std::vector<Var>&, INTERFACE_SIDE=INSIDE ) const;

    /// inputs node coordinates into supplied matrix
    void           NodeCoordinateMatrix( DenseMatrix<DM_MIN>&, INTERFACE_SIDE=INSIDE ) const;

    /// the centre of gravity of the element (returns the mid-point of the 2-sides if detached)
    Point<dim>     BaryCenter() const;

    /// projects node points onto line returning max distance between them; vec direction can have any length
    double64       LengthInDirection( const VectorVariable<dim>& vecDirection ) const;

    // ------------------------------------------------------------------------
    // Screen Output
    // ------------------------------------------------------------------------

    void Out() const;

  private:

    void InitializeNodeCorrespondanceVector();
  
    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t                idx_;       ///< unique identifier for indexing operations
    std::vector<Node<dim>*>       node_connector_;      ///< pointers to the nodes on inside followed by those on the outsie
    std::vector<InterFace<dim>*>  interface_connector_; ///< nbor interfaces on inside followed by those on outside

    std::vector<std::pair< size_t, size_t > > 
      parent_elements_node_connector_; ///< parent element local node-ids on inside followed by those on outside
    
    Element<dim>* baseElement_; ///< Element object from which interface was constructed
    
    // used for compatibility with Element and Face methods (Neighbor etc.)
    INTERFACE_SIDE  current_side_;  ///< switch to return information from INNER or OUTER side of interface

    // pointers to the higher-dimensional elements this face sits in between
    // the (same dimensional) neighbors are stored by the base class
    Element<dim>* innerParent_;           ///< higher-dimensional parent element on side opposite to normal direction
    size_t        inner_parent_face_id_;  ///< face number of inside higher-dimensional parent element

    Element<dim>* outerParent_;           ///< higher-dimensional neighbor element in direction of interface normal
    size_t        outer_parent_face_id_;  ///< face number of outside higher-dimensional parent element
};



/**

@class InterFace InterFace "main_library/InterFace.h"

@author SKM refactored and changed design to node-pointer based etc.
@author P. Lang
@author Roman
@date 2011,2014, 2016

@note InterFaces are not listed as parents in their nodes.

@section motivation Motivation

To represent surfaces / interfaces in an mesh.  

@section conventions Conventions

The base (csmp::Element) class has in its node container the nodes of the inner element only.
If one wants to acces nodes of an interface, InterFace::N( i, j ) is to be used.

*/


} // csmp

#endif



