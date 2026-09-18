// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NumIntegral_NT_lhsop_N_dV_h
#define NumIntegral_NT_lhsop_N_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
     Mass or capacitance matrix
     
     C_jk = ∫_Ω N_j [σ] N_k dV
     
     To represent heat capacity, storativity (specific storage).
     Used as transient storage term in pressure diffusion, compressibility matrix in poromechanics.
     
     @note Supported via LumpedFormulation(true) — produces a diagonal matrix whose entries equal the row sums of the consistent matrix.
 */
template<uint32_t dim, template<uint32_t> class CELL=csmp::Element>
class NumIntegral_NT_lhsop_N_dV final : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>&,
                               const char* oper,
                               const char* basic,
                               const char* test );

    void GetOperands(         const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

    NumIntegral_NT_lhsop_N_dV<dim,CELL>* clone() const override final
      { return new NumIntegral_NT_lhsop_N_dV<dim,CELL>(*this); }

  private:
    std::vector<double> op_;    ///< scalar material operand values —
                                ///< size 1 for element/face/region placement,
                                ///< size n_integration_points otherwise
    DenseMatrix<DM_MIN> NT_;    ///< interpolation function column vector (nodes x 1)
    DenseMatrix<DM_MIN> N_;     ///< interpolation function row vector    (1 x nodes)
    DenseMatrix<DM_MIN> TEMP_;  ///< accumulated mass matrix              (nodes x nodes)
};

/* copyright (c) 1997 by Dr. Stephan K. Matthaei & Stephen G. Roberts */

} // csmp

#endif
















