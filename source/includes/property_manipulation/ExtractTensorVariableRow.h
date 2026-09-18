// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef EXTRACT_TENSOR_VARIABLE_ROW_H
#define EXTRACT_TENSOR_VARIABLE_ROW_H

#include "Interrelation.h"


namespace csmp {

  /** Interrelation which maps the user-specified row of a tensor variable
  into the user-supplied vector variable.

  @author A Paluszny
  @date 6/2006
  */

template<uint32_t dim>
class ExtractTensorVariableRow : public Interrelation<dim> {
    Operand<dim>&        t;
    Operand<dim>&        v;
    TensorVariable<dim>  ts;
    const uint32_t       row;
    
  public:
    ExtractTensorVariableRow( const PropertyDatabase<dim>& p,
                              const char* tens_var, const char* vector_var, 
                              uint32_t r );
                                    
    void Calculate() override final;
};




template<uint32_t dim>
inline void ExtractTensorVariableRow<dim>::Calculate()
 {
    t.AssignTo( ts );
    v = ts.Row(row);
 } 

} // csmp

#endif

