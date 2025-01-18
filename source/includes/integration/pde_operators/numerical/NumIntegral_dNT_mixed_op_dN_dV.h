#ifndef NumIntegral_dNT_mixed_op_dN_dV_h
#define NumIntegral_dNT_mixed_op_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// PDE operator:  K div^2 T
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_mixed_op_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_dNT_mixed_op_dN_dV( const PropertyDatabase<dim>&,
                                    const char* oper, 
                                    const char* nodal_oper_multiplier,
                                    const char* basic, 
                                    const char* test );
    
    virtual void GetOperands( const CELL<dim>& );
    virtual void ComputeContribution( const CELL<dim>& );

  private:
    DenseMatrix<DM_MIN>     B, BT, NVAL;
    Index                   nkey; ///< nodal multiplier for operand
    std::vector<double>   ip_nmult;
};

} // csmp

#endif
