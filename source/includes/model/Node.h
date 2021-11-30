#ifndef CSMP_NODE_H
#define CSMP_NODE_H

#include "Box.h"
#include "Point.h"
#include "LocalVariableStorage.h"

namespace csmp {

template<size_t dim> class Visitor;
template<size_t dim> class Element;
template<size_t dim> class NodeManifold;

/**
 
@brief Node = lowest class of the model hierarchy.

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1999

@section motivation Motivation

Node class in the finite element mesh hierarchy
stores the coordinate of the node point and the variables
associated with the node.
 
@section design Design Intent

Mesh vertex with coordinates and associated flags.  
 
 
@section participants Participants

Uses Point<> template to represent the coordinate.  
 
 
@section collaborations Collaborations

Elements are registered as parents, Faces and InterFaces are not.
 
*/
template<size_t dim>
class Node : public LocalVariableStorage<dim,Node> {
  public:
    Node();
    /// custom constructor used when model is reconstructed from binary file
    Node( size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT );
    ~Node();
    /// constructs Node with same idx, position, property values, and pointer connections as the argument Node
    Node( const Node& );
    Node( Node&& );
    Node& operator=( const Node& );
    Node& operator=( Node&& );

    /// relation operators
    bool operator==( const Node<dim>& );

    /// Local variable storage interface
    PLACEMENT Placement() const { return NODE; }

    // maintaining node to parent element connectivity
    /// assign new parent element where there is a  NOT_INITIALISED  slot in the parent element storage
    void             Assign( size_t parent_elmt_node_number, Element<dim>* parent_elmt );
    /// if found, sets matching parent element pointer to nullptr and the corresponding node number to NOT_INITIALIZED
    bool             Unassign( Element<dim>* parent_elmt );
    /// changes parent element related containers to new size
    void             ResizeParentStorage( size_t parent_elements );
    /// removing parent elements that were previously assigned a nullptr
    void             EraseNullPointerParents();
    /// remove all current parent elements
    void             EraseParents();

    /// connects the node to other topologically collocated nodes if any
    void             Assign( NodeManifold<dim>* const );
    
    /// access to manifold if any; returns nullptr if the node is not a manifold
    bool             IsManifold() const;

    /// access to other topologically collocarted Node objects through manifold if any; returns nullptr if the node is not a manifold
    NodeManifold<dim>* const Manifold() const;

    /// support of the vistor design pattern
    void             Accept( csmp::Visitor<dim>& );

    /// returns how many elements share this node
    size_t           Parents() const;
    /// access to the (0..n-1) parent element 
    Element<dim>*    Parent( size_t ) const;
    /// the local number of this node within the node-numbering scheme of parent element (and equal to sector number)
    size_t           ParentNodeNumber( size_t parent_element ) const;
    /// the number of nodes that this node is connected with
    size_t           Neighbors() const;
    /// access to any of these nodes
    Node<dim>*       Neighbor( size_t ) const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void             Idx( size_t id_0_to_n_minus_1 ) const; // since idx is mutable
    size_t           Idx() const;
    Point<dim>       Coordinate() const;
    void             Coordinate( const Point<dim>& );
    /// flagging for box-shaped models: NOT, LEFT, BOTTOM, RIGHT, TOP, BACK, FRONT etc.
    void             AtBoundary( BOX_BOUNDARY );
    BOX_BOUNDARY     AtBoundary() const;

    /// accessors/mutators for specific node coordinates x=0, y=1, z=2 (z exists only in 3D)
    double         operator[]( size_t i ) const;
    double&        operator[]( size_t i );
    double&        operator()( size_t i );
    
    void             x( double );
    void             y( double );
    void             z( double );
    
    double         x() const;
    double         y() const;
    double         z() const;

    /// output current state of class Node
    void             Out() const;

  private:
    /// private because these operators are owned by the MeshManager
    template<size_t> friend class MeshManager;
    void* operator new( size_t size );
    void  operator delete( void* );
  
  private:
    mutable size_t                 idx_;                      ///< 0..n-1
    BOX_BOUNDARY                   at_boundary_;              ///< which model boundary the Node is on
    Point<dim>                     xyz_;                      ///< coordinate array
    std::vector<ONE_BYTE_NUMBER>   parent_node_indexes_;      ///< local parent node number (0...nodes-1)
    std::vector<Element<dim>*>     parent_element_pointers_;  ///< parent element pointers
    NodeManifold<dim>*             manifold_ = nullptr;       ///< node manifold pointer
    
    friend class FiniteElement_TestData; // for testing 
};


// FUNCTIONS INVOLVING NODES

/// returns parent elements shared by face, inner side is reported first; outer next else application: give nodes of lower-dimensional face to find element on either side
template<size_t dim>
std::pair<Element<dim>*,Element<dim>*>  parentElementsSharedByFace( const std::vector<Node<dim>*>& face_nodes_in_correct_order );

} // csmp

#endif






