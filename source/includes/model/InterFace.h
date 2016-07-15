#ifndef CSMP_INTER_FACE_H
#define CSMP_INTER_FACE_H

#include "Element.h"

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementTraits.h"
#include "FiniteVolumeTraits.h"
#include "InterFaceRemeshingTraits.h"


namespace csmp {

/// surface (3D) or line (2D) element connector for elements which share faces but not nodes
template<size_t dim>
class InterFace : public InterFaceRemeshingTraits<dim, InterFace>,
                  public FiniteElementTraits<dim,InterFace >,
                  public FiniteVolumeTraits<dim,InterFace >,
                  public LocalVariableStorage<dim,InterFace<dim> >
{

  public:

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
    InterFace( InterFace&& );
    ~InterFace();
    
    InterFace& operator=( const InterFace& );

    /// relation operators
    bool operator==( const InterFace<dim>& );

    /// Local variable storage interface
    PLACEMENT Placement() const { return INTER_FACE; }

    /// Visitor
    void Accept( csmp::Visitor<dim>& vis );
    
    // ------------------------------------------------------------------------
    // Functionality of construction process
    // ------------------------------------------------------------------------

    /// Assigns parents and initializes interface's node vector
    void Assign( const FiniteVolumeStencil<dim>* const );
    void Assign( size_t nbor, InterFace<dim>* const );
    void Assign( Element<dim>* const base_elmt );
    void Assign( Element<dim>* const inner_elmt, Element<dim>* const outer_elmt, bool assign_nodes = true );
    void Assign( Element<dim>* const parent, size_t faceId, INTERFACE_SIDE side );
    void Assign( size_t n_local, size_t parent_node, INTERFACE_SIDE side );
    void Assign( size_t n_local, std::pair<size_t, size_t> parent_nodes );
    void Assign( FiniteElement*); /// this method was added to aid OpenMP usage.

    // TODO: SKM FIX - assign the actual node
    void Assign( size_t n_local, Node<dim>*, INTERFACE_SIDE side );

    // ------------------------------------------------------------------------
    // Member access
    // ------------------------------------------------------------------------

    // for remeshing purposes
    typename  std::vector<csmp::InterFace<dim>*>& NeighborElementVector();

    csmp::Node<dim>* N( size_t n_local, INTERFACE_SIDE side ) const;
    csmp::Node<dim>* N( size_t n_local ) const;
    
    /// switches internal state variable that sets interface side
    void            CurrentSide( INTERFACE_SIDE side );
    INTERFACE_SIDE  CurrentSide() const;

    /// standard interface returns neighbor on current side of interface
    csmp::InterFace<dim>*  Neighbor( size_t ) const;
    
    /// neighbor on choosen side of interface is returned
    csmp::InterFace<dim>*  Neighbor( size_t, INTERFACE_SIDE ) const;

    FiniteElement*                  FE() const;
    const FiniteVolumeStencil<dim>* FV_Stencil() const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void           Idx( size_t ) const;
    size_t         Idx() const;

    Element<dim>*  Parent( INTERFACE_SIDE side ) const;
    Element<dim>*  InnerParent() const;
    Element<dim>*  OuterParent() const;
    Element<dim>*  BaseElement() const;
    // TODO: check and capitalize
    bool           hasBase()      const;

    size_t         InnerParentFaceID() const;
    size_t         OuterParentFaceID() const;
    size_t         ParentFaceID( INTERFACE_SIDE side ) const;

    /// Node numbers correspondance
    size_t         ParentNodeNumber( size_t n_local, INTERFACE_SIDE side ) const;
    size_t         ParentNodeNumberOppositeTo( size_t n_from_parent, INTERFACE_SIDE side ) const;

    // ------------------------------------------------------------------------
    // Geometry
    // ------------------------------------------------------------------------

    double64       Area() const;
    void           UnitNormal( VectorVariable<dim>&, INTERFACE_SIDE side ) const;
    void           UnitNormal( VectorVariable<dim>& )         const;
    // computations for spacing between node pairs
    bool           SpacingNode( size_t n_local, VectorVariable<dim>& ) const;

    // ------------------------------------------------------------------------
    // Screen Output
    // ------------------------------------------------------------------------

    void Out() const;

  private:

    void InitializeNodeCorrespondanceVector();
  
    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t                          idx_;       ///< unique identifier for indexing operations
    csmp::FiniteElement*                    fptr_;      ///< bridge-(pattern) connector to FiniteElement
    const csmp::FiniteVolumeStencil<dim>*   fvptr_;     ///< connector to finite-volume functionality
    
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



