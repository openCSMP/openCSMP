#ifndef EXTRACT_TENSOR_VARIABLE_COMPONENT_H
#define EXTRACT_TENSOR_VARIABLE_COMPONENT_H

#include "Interrelation.h"


namespace csmp {

  /**

  Interrelation which maps the user-specified component of a tensor variable
  into the user-supplied scalar variable.
  This may be useful for instance if one wants to output and display the
  variable using a color scheme rather than ellipsoids.

  @author S.K. Matthaei
  @date 6/2000
  */

template<uint32_t dim>
class ExtractTensorVariableComponent : public Interrelation<dim> {
    Operand<dim>&        T;
    Operand<dim>&        S;
    TensorVariable<dim>  ts;
    const size_t         comp_i, comp_j;
    
  public:
    ExtractTensorVariableComponent( const PropertyDatabase<dim>& p,
                                    const char* tens_var, const char* scalar_var, 
                                    size_t i, size_t j );
                                    
    ~ExtractTensorVariableComponent();
    
    void Calculate();
};




template<uint32_t dim>
inline void ExtractTensorVariableComponent<dim>::Calculate()
 {
    T.AssignTo( ts );
    S = ts(comp_i,comp_j);
 } 

} // csmp

#endif

