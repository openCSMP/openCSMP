#ifndef NumIntegral_dNT_mixed_op_dN_dV_h
#define NumIntegral_dNT_mixed_op_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// PDE operator:  K div^2 T
template<size_t dim,class CELL>
class NumIntegral_dNT_mixed_op_dN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_dNT_mixed_op_dN_dV( const PropertyDatabase<dim>& pref, 
                        const char* oper, 
                        const char* nodal_oper_multiplier,
                        const char* basic, 
                        const char* test );
    
    virtual void GetOperands( CELL& e );
    virtual void ComputeContribution( CELL& e );

  private:
    DenseMatrix<DM_MIN>     B, BT, NVAL;
    Index                   nkey; ///< nodal multiplier for operand
    std::vector<double64>   ip_nmult;
};

} // csmp

#endif
