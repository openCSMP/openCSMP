#ifndef CSMP_NODE_H
#define CSMP_NODE_H

#include "CSMP_number_types.h"
#include "Box.h"
#include "Point.h"
#include "LocalVariableStorage.h"

namespace csmp {

template<size_t dim> class Visitor;
template<size_t dim> class Element;

// NODE class for Model
template<size_t dim>
class Node : public LocalVariableStorage<dim,Node<dim> >
  {
  public:
    Node();
    /// custom constructor used when model is reconstructed from binary file
    Node( size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY=NOT );
    ~Node();
    Node( const Node& nd );
    Node( Node&& nd );
    Node& operator=( const Node& );

    /// relation operators
    bool operator==( const Node<dim>& );

    /// Local variable storage interface
    PLACEMENT Placement() const { return NODE; }

    void             Assign( size_t parent_elmt_node_number, Element<dim>* parent_elmt );
    bool             Unassign( Element<dim>* parent_elmt );
    void             ResizeParentStorage( size_t parent_elements );
    void             EraseParents();
    
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
    /// accessors/mutators for specific node coordinates x=0, y=1, z=2 (z exists only in 3D)
    double64         operator[]( size_t i ) const;
    double64&        operator[]( size_t i );
    double64&        operator()( size_t i );
    void             x( double64 );
    void             y( double64 );
    void             z( double64 );
    double64         x() const;
    double64         y() const;
    double64         z() const;

    void             Out() const { Out(std::cout); }
    void             Out(std::ostream& os) const;

    /// an flagging to be deprecated in the future
    void             AtBoundary( BOX_BOUNDARY );
    BOX_BOUNDARY     AtBoundary() const;

  private:
    mutable size_t                 idx_;                      ///< 0..n-1
    BOX_BOUNDARY                   at_boundary_;              ///< which model boundary the Node is on
    Point<dim>                     xyz_;                      ///< coordinate array
    std::vector<ONE_BYTE_NUMBER>   parent_node_indexes_;      ///< local parent node number (0...nodes-1)
    std::vector<Element<dim>*>     parent_element_pointers_;  ///< parent element pointers
};



/**
 
@class Node Node "main_library/Node.h"

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1999

@section motivation Motivation

Node class in the finite element mesh hierarchy.  
 
 
@section design Design Intent

Mesh vertex with coordinates and associated flags.  
 
 
@section participants Participants

Uses Point<> template to represent the coordinate.  
 
 
@section collaborations Collaborations

Elements are registered as parents, Faces and InterFaces are not.
 
*/


} // csmp

#endif






