// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_DUDN_OP_U_DS_IFACE_DIM1_T_H
#define NUM_INTEGRAL_DUDN_OP_U_DS_IFACE_DIM1_T_H

#include "MathOperatorLHS.h"
#include "Index.h"
#include "InterFace.h"

namespace csmp {

class FiniteElement;
template<uint32_t> class Model;
template<uint32_t> class InterFace;


template<uint32_t dim>
class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature : public MathOperatorLHS<dim,InterFace> {
  public:
    NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature( const Model<dim>&,
                                                       const char* oper,
                                                       const char* basic,
                                                       const char* test);
    /**
       Computes difference between the variable values at the topologically collocated nodes, including that in the intervening
       lower-dimensional region, using it to compute   the transfer coefficient (=coupling coefficient) at the time level t + delta t.
    */
    virtual void GetOperands( const InterFace<dim>& ) override;
    
    virtual void ComputeContribution( const InterFace<dim>& ) override;
    
    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions, Luat Khoa Tran)
    virtual void AssignToGlobal( const InterFace<dim>&, SparseMatrix&, std::vector<double>&, const std::vector<size_t>& ) override;
    
    
    virtual NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>* clone() const override
      { return new NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region_temperature<dim>(*this); }
      

  private:
      csmp::INDEX<SCALAR,ELEMENT>  thickness_key_;  //Benoit add
      std::vector<double>    transfer_coefficients_;  ///< as computed from delta t and delta test-function operand value at the nodes
      const bool  order1_FD_approximation_ = true;
};

} // csmp


#endif
