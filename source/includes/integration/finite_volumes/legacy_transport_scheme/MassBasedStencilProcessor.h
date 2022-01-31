#ifndef MASS_BASED_STENCIL_PROCESSOR_H
#define MASS_BASED_STENCIL_PROCESSOR_H

#include "ExplicitStencilProcessor.h"
#include "Index.h"

namespace csmp {


template<size_t dim>
class MassBasedStencilProcessor : public ExplicitStencilProcessor<dim>{
public:

  MassBasedStencilProcessor( const csmp::Index& adv_lhs_key, 
                             const csmp::Index& adv_rhs_key, 
                             const csmp::Index& velo_key );

  virtual void InitializeFirstOrder( const FV_Parameter& param, const Element<dim>& e, const VARIABLE_TYPE &vt=SCALAR, const size_t var_comp_nr=0 );
  
  csmp::Index adv_rhs_key_;

  std::vector<double> psi2_;  // rhs advection variable
};


}// end csmp  
#endif //MASS_BASED_STENCIL_PROCESSOR_H
