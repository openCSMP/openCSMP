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

template<size_t dim>
class ExtractVectorVariableComponent : public Interrelation<dim> {
    Operand<dim>&        V;
    Operand<dim>&        S;
    VectorVariable<dim>  vc;
    const size_t         component;
    
  public:
    ExtractVectorVariableComponent( const PropertyDatabase<dim>& p,
                                    const char* vec_var, const char* scalar_var, 
                                    size_t comp );
                                    
    ~ExtractVectorVariableComponent();
    
    void Calculate();
};




template<size_t dim>
inline void ExtractVectorVariableComponent<dim>::Calculate()
 {
    V.AssignTo( vc );
    S = vc(component);
 } 

} // csmp

#endif

