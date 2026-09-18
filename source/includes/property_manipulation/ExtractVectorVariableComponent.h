// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef EXTRACT_VECTOR_VARIABLE_COMPONENT_H
#define EXTRACT_VECTOR_VARIABLE_COMPONENT_H

#include "Interrelation.h"


namespace csmp {

  /**

  Interrelation which maps the user-specified component of a vector variable
  into the user-supplied scalar variable.
  This may be useful for instance if one wants to output and display the
  variable using a color scheme rather than vectors.

  @author S.K. Matthaei
  @date 6/2000
  */

template<uint32_t dim>
class ExtractVectorVariableComponent : public Interrelation<dim> {
    Operand<dim>&        V;
    Operand<dim>&        S;
    VectorVariable<dim>  vc;
    const uint32_t       component;
    
  public:
    ExtractVectorVariableComponent( const PropertyDatabase<dim>& p,
                                    const char* vec_var, const char* scalar_var, 
                                    uint32_t comp );
    
    void Calculate() override final;
};




template<uint32_t dim>
inline void ExtractVectorVariableComponent<dim>::Calculate()
 {
    V.AssignTo( vc );
    S = vc(component);
 } 

} // csmp

#endif

