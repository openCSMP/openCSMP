// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MASS_BASED_STENCIL_PROCESSOR_H
#define MASS_BASED_STENCIL_PROCESSOR_H

#include "ExplicitStencilProcessor.h"
#include "Index.h"

namespace csmp {


template<uint32_t dim>
class MassBasedStencilProcessor : public ExplicitStencilProcessor<dim>{
public:

  MassBasedStencilProcessor( const csmp::Index& adv_lhs_key, 
                             const csmp::Index& adv_rhs_key, 
                             const csmp::Index& velo_key );

  virtual void InitializeFirstOrder( const FV_Parameter& param, const Element<dim>& e, const VARIABLE_TYPE &vt=SCALAR, uint32_t var_comp_nr=0 );
  
  csmp::Index adv_rhs_key_;

  std::vector<double> psi2_;  // rhs advection variable
};


}// end csmp  
#endif //MASS_BASED_STENCIL_PROCESSOR_H
