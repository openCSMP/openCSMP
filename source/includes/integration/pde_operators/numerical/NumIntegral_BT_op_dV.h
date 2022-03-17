#ifndef NumIntegral_BT_op_dV_h
#define NumIntegral_BT_op_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

/** 

@brief classical stiffness matrix, i.e., interpolation function derivative matrix
for a vector solution variable like displacement which has u, v, w components.

@author S.K. Matthai
@author S. Roberts
@date 1999 

*/
template<uint32_t dim,class CELL=Element<dim> >
class NumIntegral_BT_op_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_BT_op_dV( const PropertyDatabase<dim>&,
                          const char* oper, const char* test );

    virtual void ComputeContribution( const CELL& );
    
    virtual NumIntegral_BT_op_dV<dim,CELL>* clone() const { return new NumIntegral_BT_op_dV<dim,CELL> (*this); }
    
  private:  
    DenseMatrix<DM_MIN>  B, BT, STR; 
};

} // csmp

#endif
