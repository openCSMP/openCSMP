#ifndef EXTRACT_VECTOR_VARIABLE_LENGTH_H
#define EXTRACT_VECTOR_VARIABLE_LENGTH_H

#include "Interrelation.h"



namespace csmp {

  /**

  Interrelation to map the length of a vector variable into a scalar variable.
  This may be useful for instance if one wants to output and display the
  variable using a color scheme rather than vectors.

  @author S.K. Matthaei
  @date 6/2000
  */

template<size_t dim>
class ExtractVectorVariableLength : public Interrelation<dim> {
    Operand<dim>&        V;
    Operand<dim>&        S;
    VectorVariable<dim>  vc;
    
  public:
    ExtractVectorVariableLength( const PropertyDatabase<dim>& p,
                                 const char* vec_var, const char* to_scalar_var );
                                    
    ~ExtractVectorVariableLength();
    
    void Calculate();
};





template<size_t dim>
inline void ExtractVectorVariableLength<dim>::Calculate()
 {
    V.AssignTo( vc );
    S = vc.Length();
 } 

} // csmp

#endif

