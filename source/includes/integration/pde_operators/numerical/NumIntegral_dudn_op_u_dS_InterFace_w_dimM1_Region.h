#ifndef NUM_INTEGRAL_DUDN_OP_U_DV_H
#define NUM_INTEGRAL_DUDN_OP_U_DV_H

#include "MathOperatorLHS.h"
#include "Index.h"

namespace csmp {

class FiniteElement;
template<uint32_t> class Model;

double diffusion1D( double val_farfield, double diffusivity, double x, double t );
double flux1DAtX0( double val_farfield, double diffusivity, double t, double transfer_coefficient );


/**
       Use to apply a Robin type (flux) boundary condition coupling the sides of the Interface with a lower-dimensional
       region in the middle between them.
       
       @note use for transient computations involving fractures
       
       @author SKM
       @date 18/8/22
*/
template<uint32_t dim>
class NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region : public MathOperatorLHS<dim,InterFace> {
  public:
    NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region( const Model<dim>&,
                                                       const char* oper,
                                                       const char* basic,
                                                       const char* test,
                                                       double delta_t );
    /**
       Computes difference between the variable values at the topologically collocated nodes, including that in the intervening
       lower-dimensional region, using it to compute   the transfer coefficient (=coupling coefficient) at the time level t + delta t.
    */
    virtual void GetOperands( const InterFace<dim>& ) override;
    
    virtual void ComputeContribution( const InterFace<dim>& ) override;
    
    /// used by PDE_IntegratorUoM for assembly of a pre-eliminated solution matrix and RH vector (scalar versions, Luat Khoa Tran)
    virtual void AssignToGlobal( const InterFace<dim>&, SparseMatrix&, std::vector<double>&, const std::vector<size_t>& ) override;
    
    /// adjusting the time increment in case it changes during the transient calculation
    void UpdateTimeIncrement( double dt ) { delta_t_ = dt; }
    
    virtual NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>* clone() const override
      { return new NumIntegral_dudn_op_u_dS_InterFace_w_dimM1_Region<dim>(*this); }

  private:
      csmp::INDEX<SCALAR,ELEMENT>  diffusivity_key_;  ///< test function key (in base class)
      std::vector<double>    transfer_coefficients_;  ///< as computed from delta t and delta test-function operand value at the nodes
      double  delta_t_;                               ///< current time increment (between t0 and t + delta t)
};

} // csmp


#endif /* NUM_INTEGRAL_DUDN_OP_U_DV_H */
