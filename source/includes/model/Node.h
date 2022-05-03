#ifndef CSMP_NODE_H
#define CSMP_NODE_H

#include "Box.h"
#include "Point.h"
#include "LocalVariableStorage.h"

namespace csmp {

template<uint32_t dim> class Visitor;
template<uint32_t dim> class Element;
template<uint32_t dim> class NodeManifold;
class FiniteElement_TestData;

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
template<uint32_t dim>
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

    /// comparitor (Roman, 2014), very costly. @todo rethink logic and rewrite
    bool operator==( const Node<dim>& );

    /// Local variable storage interface
    PLACEMENT Placement() const { return NODE; }


    // node to parent element connectivity (sorted vector that is searchable)
    
    /// assign new parent element where there is a  NOT_INITIALISED  slot in the parent element storage
    void Assign( uint32_t parent_elmt_node_number, Element<dim>* parent_elmt );
    /// if found, sets matching parent element pointer to nullptr and the corresponding node number to NOT_INITIALIZED
    bool Unassign( Element<dim>* parent_elmt );
    /// changes parent element related containers to new size
    void ResizeParentStorage( uint32_t parent_elements );
    /// sorts parent vectors for searching
    void SortParents();
    /// removing parent elements that were previously assigned a nullptr
    void EraseNullPointerParents();
    /// remove all current parent elements
    void EraseParents();

    /// returns how many elements share this node
    uint32_t Parents() const;
    /// access to the (0..n-1) parent element
    Element<dim>* Parent( uint32_t ) const;
    /// the local number of this node within the node-numbering scheme of parent element (and equal to sector number)
    uint32_t ParentNodeNumber( uint32_t parent_element ) const;
    /// checks whether Element is a parent of the node
    bool IsParent( const Element<dim>* const ) const;
    
    
    // node manifolds (where nodes have been multiplicated at material interfaces)
    
    /// connects the node to other topologically collocated nodes if any
    void Assign( NodeManifold<dim>& );
    
    /// access to manifold if any; returns nullptr if the node is not a manifold
    bool IsManifold() const;

    /// access to other topologically collocarted Node objects through manifold if any; returns nullptr if the node is not a manifold
    NodeManifold<dim>* const Manifold() const;


    // node neighbors (on the opposite side of the finite element segments that the node is on)
    // (sorted vector that is searchable)
    
    /// initialises the corner-node to neighbor corner node pointer vector
    void Assign( std::set<Node<dim>*>& neighbor_nodes );

    /// expects a sorted vector without duplicates
    void Assign( std::vector<Node<dim>*>& neighbor_nodes, bool sort_neighbors=false );
    
    /// copies property values from the argument node to the current node
    void CopyPropertyValuesFrom( Node<dim>& );

    /// rebuilds the neighbor connectivity working through higher-dimensional parent element edges that the node is part of; returns new number of neighbors
    uint32_t ReassignNeighbors();

    /// removes null pointers and potential duplicates returning the resulting number of neighbors
    uint32_t UpdateNeighbors();
    
    bool IsNeighbor( const Node<dim>* const ) const;
    void AddNeighbor( Node<dim>* neighbor_node );
    void RemoveNeighbor( const Node<dim>* const neighbor_node );
    
    /// the number of corner nodes that this node is directly connected with via segments
    uint32_t Neighbors() const;
    
    /// access to any of the neighbor nodes
    Node<dim>* Neighbor( uint32_t ) const;


    // basic functionality of the Node
    
    /// support of the vistor design pattern
    void Accept( csmp::Visitor<dim>& );

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void Idx( size_t id_0_to_n_minus_1 ) const; // since idx is mutable
    size_t Idx() const;
    void Coordinate( const Point<dim>& );
    Point<dim> Coordinate() const;
    /// flagging for box-shaped models: NOT, LEFT, BOTTOM, RIGHT, TOP, BACK, FRONT etc.
    void AtBoundary( BOX_BOUNDARY );
    BOX_BOUNDARY AtBoundary() const;

    /// accessors/mutators for specific node coordinates x=0, y=1, z=2 (z exists only in 3D)
    double  operator[]( uint32_t i ) const;
    double& operator[]( uint32_t i );
    double& operator()( uint32_t i );
    
    void x( double );
    void y( double );
    void z( double );
    
    double x() const;
    double y() const;
    double z() const;

    /// output current state of class Node
    void Out() const;

  private:
    /// private because these operators are owned by the MeshManager
    template<uint32_t> friend class MeshManager;
    void* operator new( size_t size );
    void  operator delete( void* );
  
  private:
    Point<dim>                     xyz_;                      ///< coordinate array
    mutable size_t                 idx_;                      ///< 0..n-1
    std::vector<Element<dim>*>     parent_element_pointers_;  ///< parent element pointers
    std::vector<Node<dim>*>        neighbor_node_pointers_;   ///< corner node to corner node on opposite end of the segment pointer
    NodeManifold<dim>*             manifold_ = nullptr;       ///< node manifold pointer
    std::vector<ONE_BYTE_NUMBER>   parent_node_indexes_;      ///< local parent node number (0...nodes-1)
    BOX_BOUNDARY                   at_boundary_;              ///< which model boundary the Node is on
    
    friend class FiniteElement_TestData; // for testing
    friend std::istream& operator >> ( std::istream&, FiniteElement_TestData& );
};


// FUNCTIONS INVOLVING NODES

/// returns  elements that share face, inner side is reported first; outer next; face-nodes must be in correct order. Application: from nodes of lower-dimensional face find elements on in- and outside
template<uint32_t dim>
std::pair<Element<dim>*,Element<dim>*>  parentElementsSharedByFace( typename std::vector<Node<dim>*>::const_iterator first,
                                                                    typename std::vector<Node<dim>*>::const_iterator last );

/// if there is only one parent element expected, then use this method instead of 'parentElementsSharedByFace'
template<uint32_t dim>
std::pair<Element<dim>*,size_t>  parentElement( typename std::vector<Node<dim>*>::const_iterator first,
                                                typename std::vector<Node<dim>*>::const_iterator last );

/// prints Idx and boundary flag values of the nodes that this node is connected with
template<uint32_t dim>
void printNeighbors( const Node<dim>* const );

/// prints current parent information and checks for duplicate parents
template<uint32_t dim>
void printParents( const Node<dim>* const );

/// calculates the size of the Node excluding the dynamic contribution to the stored variables
template<uint32_t dim>
size_t sizeOf( const Node<dim>* const );

} // csmp

#endif






