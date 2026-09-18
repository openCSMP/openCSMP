// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef EXTRACT_TENSOR_VARIABLE_EIGENVALUE_H
#define EXTRACT_TENSOR_VARIABLE_EIGENVALUE_H

#include "Interrelation.h"


namespace csmp {

  /** Interrelation which maps the user-specified eigenvalue of a tensor variable
  into the user-supplied scalar variable.

  @author A Paluszny
  @date 6/2006
  */

template<uint32_t dim>
class ExtractTensorVariableEigenvalue : public Interrelation<dim> {
    Operand<dim>&        t;
    Operand<dim>&        s;
    TensorVariable<dim>  ts;
    const uint32_t       eval;
    
  public:
    ExtractTensorVariableEigenvalue( const PropertyDatabase<dim>& p,
                                     const char* tens_var, const char* scalar_var,
                                     uint32_t eval );
                                    
    void Calculate() override final;     
};




template<uint32_t dim>
inline void ExtractTensorVariableEigenvalue<dim>::Calculate()
 {
    t.AssignTo( ts );
    VectorVariable<dim> evals;
    ts.EigenValues(evals);
    s = evals[eval];
 } 

} // csmp

#endif

