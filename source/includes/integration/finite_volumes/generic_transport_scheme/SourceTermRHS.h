//
//  SourceTermRHS.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_SOURCE_TERM_RHS_H
#define CSMP_SOURCE_TERM_RHS_H

#include "VectorOperator.h"

namespace csmp {

template<size_t dim>
class SourceTermRHS : public VectorOperator<dim> {
  public:
    /// note that this source variable may be an element or a nodal one
    SourceTermRHS( const Model<dim>&, const char* source_variable );
  
    /// accumulates integrated (single) values of source term onto the finite volume centered on the Node
    virtual void AccumulateFV( const Node<dim>*, std::vector<double64>& rhs ) const;
  
    /// accumulates distributed values of the source term on the finite volume
    virtual void AccumulateStencil( const Element<dim>*, std::vector<double64>& rhs ) const;
  
  private:
    const csmp::Index src_key_;

};

} // end csmp

#endif /* CSMP_SOURCE_TERM_RHS_H */
