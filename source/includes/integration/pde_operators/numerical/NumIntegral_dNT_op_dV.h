#ifndef NUM_INTEGRAL_DNT_OP_DV_H
#define NUM_INTEGRAL_DNT_OP_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/**
    To integrate over a gradient represented by a vector property.
    
    @example Compute   gravityTerm = rho_w * gravityVector
    and project it onto the dip-vector of a lower dimensional element.
    The resulting vector goes into the righthandside integrated numerically via this integral:
    
    NumIntegral_dNT_op_dV(  model.Database(), "gravity term", "fluid pressure" );,
    
    @author Shaho Bazr-Afkan
    @date 2011
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_op_dV final : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_dNT_op_dV( const PropertyDatabase<dim>&,
                           const char* oper,    ///< (vector) gradient property, e.g., rho g grad z
                           const char* test );  ///< scalar, for instance fluid pressure
    
    void ComputeContribution( const CELL<dim>& ) override final;
    
    NumIntegral_dNT_op_dV<dim,CELL>* clone() const override final { return new NumIntegral_dNT_op_dV<dim,CELL> (*this); }
    
  private:
    DenseMatrix<DM_MIN>  B_, BT_; 
};

} // csmp

#endif
