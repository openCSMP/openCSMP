// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_DUDN_RHSOP_U_DV_DS_H
#define NUM_INTEGRAL_DUDN_RHSOP_U_DV_DS_H

#include "MathOperatorRHS.h"
#include "Index.h"
#include "InterFace.h"

namespace csmp {

class FiniteElement;
template<uint32_t> class Model;

/**
       Use to apply a Robin type (flux) boundary condition coupling the sides of the Interface with a lower-dimensional
       region in the middle between them.
       
       @note use for transient computations involving fractures
       
       @author SKM
       @date 18/8/22
*/
template<uint32_t dim>
class NumIntegral_dudn_rhsop_u_dS : public MathOperatorRHS<dim,InterFace> {
  public:
    /// interface variable has to be flagged Robin and elmt variable will be read as a matching variable placed on the element
    NumIntegral_dudn_rhsop_u_dS( const Model<dim>&, const char* interface_oper, const char* elmt_oper, const char* test );
    /**
       Computes difference between the variable values at the topologically collocated nodes, including that in the intervening
       lower-dimensional region, using it to compute   the transfer coefficient (=coupling coefficient) at the time level t + delta t.
    */
    void GetOperands( const InterFace<dim>& ) override final;
    
    void ComputeContribution( const InterFace<dim>& ) override final;
    
    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions, Luat Khoa Tran)
 //   void AssignToGlobal( const InterFace<dim>&,
 //                                std::vector<double>& rhs,
 //                                const std::vector<size_t>& DOF_indexes ) override final;
    
    /// adjusting the time increment in case it changes during the transient calculation
    void UpdateTimeIncrement( double dt ) { delta_t_ = dt; }
    
    NumIntegral_dudn_rhsop_u_dS<dim>* clone() const override final;
      
    void ResetMinMaxTransferTerms() {
         transfer_min_ =  1e30;
         transfer_max_ = -1e30;
      }

  private:
      // for the computation of fluxes from the 1D analytic solution
      double Diffusion1D( double val_farfield, double diffusivity, double x, double t );
      double Grad_Var_AtX0( double val_farfield, double diffusivity, double t );
      double Flux1DAtX0( double val_farfield, double diffusivity, double t, double conductivity );

  private:
      csmp::INDEX<SCALAR,ELEMENT>  diffusivity_key_;  ///< extra key for diffusion in higher-dim elements adjacent to interface
      csmp::INDEX<SCALAR,ELEMENT>  conductivity_key_; ///< extra key for conduction in higher-dim elements adjacent to interface
      std::vector<double>    transfer_coefficients_;  ///< as computed from delta t and delta test-function operand value at the nodes
      double  delta_t_;                               ///< current time increment (between t0 and t + delta t)
      double  transfer_min_                =  1e30;
      double  transfer_max_                = -1e30;
      const bool  order1_FD_approximation_ = true;
};

} // csmp


#endif /* NUM_INTEGRAL_DUDN_OP_U_DV_H */
