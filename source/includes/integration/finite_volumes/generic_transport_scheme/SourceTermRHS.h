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

/**
   Either an areal/volumetric source specified on the Element or Face or
   a nodal (point) source term to be accumulated into the righthand vector.
*/
template<uint32_t dim>
class SourceTermRHS : public VectorOperator<dim> {
  public:
    /// note that this source variable may be an element or a nodal one
    SourceTermRHS( const csmp::INDEX<SCALAR,ELEMENT>& source_variable_key,
                   const csmp::INDEX<SCALAR,NODE>& nodal_source_variable_key,
                   const csmp::INDEX<SCALAR,NODE>& transported_variable_key );
                   
    virtual ~SourceTermRHS() {}
  
    /// accumulates integrated (single) values of source term onto the finite volume centered on the Node
    virtual void AccumulateFiniteVolume( const Node<dim>&, std::vector<double>& rhs ) const;
  
    /// accumulates distributed values of the source term on the finite volume
    virtual void AccumulateStencil( const Element<dim>&, std::vector<double>& rhs ) const;
    //virtual void AccumulateStencil( const Face<dim>&, std::vector<double>& rhs ) const;
    //virtual void AccumulateStencil( const InterFace<dim>&, std::vector<double>& rhs ) const;
  
  private:
    const csmp::Index& esrc_key_;   ///<  area/volume source  placed on the element 
    const csmp::Index& nsrc_key_;   ///< point source
    const csmp::Index& adv_key_;    ///<  transport variable
};

} // end csmp

#endif /* CSMP_SOURCE_TERM_RHS_H */
