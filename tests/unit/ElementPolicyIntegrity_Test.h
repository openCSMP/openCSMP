//
//  ElementPolicyIntegrity_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 9/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef ELEMENT_POLICY_INTEGRITY_TEST_H
#define ELEMENT_POLICY_INTEGRITY_TEST_H

/**

Can the Element class be used without finite volumes? - does it contain too much in this case.
Can we build a purely geometry element without finite element functionality.

This test distangles the cross dependencies of the policies of the Element class.

*/

#include "Element.h"
#include "Test.h"
#include "Box.h"

namespace csmp {

/**
    Building an Element only element
    with only the functionality that is stand-alone.
*/
template<size_t dim>
class ElementOnly {
  public:
    ElementOnly( size_t id, size_t nodes, size_t nbors )
     : idx_(0), at_boundary_(NOT),
       elmt_connector_(nbors,nullptr), ///< a triangle
       node_connector_(nodes,nullptr)
     {}
  // which member functions make sense here
    ElementOnly&  operator=( const ElementOnly& );
    void Accept( csmp::Visitor<dim>& );
    void Assign( size_t nbor, ElementOnly<dim>* const e_ptr ) { elmt_connector_[nbor] = e_ptr; }
    void Assign( size_t node, Node<dim>* const n_ptr ) { node_connector_[node] = n_ptr; }
    // iterators
    // accessors
    csmp::Node<dim>*         N( size_t n_local )    const { return node_connector_[n_local]; }
    csmp::ElementOnly<dim>*  Neighbor( size_t nbor ) const { return elmt_connector_[nbor]; };

    void         Idx( size_t idx ) const { idx_=idx; }
    size_t       Idx() const { return idx_; }
    void         AtBoundary( BOX_BOUNDARY b ) { at_boundary_=b; }
    BOX_BOUNDARY AtBoundary() const { return at_boundary_; }
  
    // can we do these?
    size_t     Nodes() const { return node_connector_.size(); }
    size_t     Neighbors() const { return elmt_connector_.size(); }
    size_t     Faces() const { return elmt_connector_.size(); }
    void       CoordinateMatrix() const;
    void       CoordinateMatrix( DenseMatrix<DM_MIN>& XY )  const;
  
    // maybe?    double64   AspectRatio() const;
    double64   LengthInDirection( const VectorVariable<dim>& vecDirection ) const;
    Point<dim> BaryCenter() const;
    size_t     Segments() const { return 3U; /* type dependent number, static FEM_specifications class? */ }
    Point<dim> SegmentMidPoint( size_t segm ) const;
    double64   SegmentLength( size_t segm ) const;
    void       SegmentLengths( std::vector<double64>& ) const;
    double64   AspectRatio() const;
    double64   InnerRadius() const;
  
  private:
    mutable size_t                  idx_;
    std::vector<ElementOnly<dim>*>  elmt_connector_; ///< neighbors
    std::vector<csmp::Node<dim>*>   node_connector_; ///< nodes
    BOX_BOUNDARY                    at_boundary_;
};


template<size_t dim>
class ElementWithVariableStorage : public LocalVariableStorage<dim,Element> {
  public:
    // as above?
    /// returns property values at the nodes
    template<class Var>
    void       NodePropertyVector( const csmp::Index&, std::vector<Var>& ) const;
};

template<size_t dim>
class ElementWithStorageAndFEM : public LocalVariableStorage<dim,Element>,
                                 public FiniteElementPolicy<dim,Element> {
};

/** Testing
*/
class ElementPolicyIntegrity_Test : public Test {
  public:
    virtual void run();

};


} // end csmp

#endif /* ELEMENT_POLICY_INTEGRITY_TEST_H */
