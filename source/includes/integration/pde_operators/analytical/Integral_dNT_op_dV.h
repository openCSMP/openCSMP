// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INTEGRAL_DNT_OP_DV_H
#define INTEGRAL_DNT_OP_DV_H

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
    
    @attention modelled on the numeric version by Shaho Bazr-Afkan
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_op_dV final : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_dNT_op_dV( const PropertyDatabase<dim>&,
                        const char* oper,    ///< (vector) gradient property, e.g., rho g grad z
                        const char* test );  ///< scalar, for instance fluid pressure
    
    void GetOperands( const CELL<dim>& ) noexcept override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    Integral_dNT_op_dV<dim,CELL>* clone() const override final { return new Integral_dNT_op_dV<dim,CELL> (*this); }
    
  private:
    DenseMatrix<DM_MIN>  DN_, DNT_;
    VectorVariable<dim>  grad_prop_;
};

// INLINE FUNCTIONS

/**
    Only needs to read a single VectorVariable representing the operand.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline void Integral_dNT_op_dV<dim,CELL>::GetOperands( const CELL<dim>& e ) noexcept
 {
    e.Read( this->MaterialOperandKey(), grad_prop_ );
 }




} // csmp

#endif
