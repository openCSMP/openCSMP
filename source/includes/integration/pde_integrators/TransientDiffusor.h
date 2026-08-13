#ifndef TRANSIENT_DIFFUSOR_H
#define TRANSIENT_DIFFUSOR_H

#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_dNT_op_dV.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "SAMG_Exception.h"
#else
#include "LinearSolver.h"
#endif

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;

/** transient diffusion solver with Backward-Euler time-stepping

  [C/dt + K]{p}_t+dt = [C/dt]{p}_t + q   
   ------------------------------------
   specific to the CSP reservoir simulator
   only a single equation is solved for total pressure
   result variable 'fluid pressure'
   initialises hydrostatic pressure to zero
   
   @attention gradient_multiplier currently is not used.
   
*/
template<uint32_t dim,template<uint32_t> class CELLTYPE=Element>
class TransientDiffusor : public PDE_Integrator<dim,CELLTYPE> {
  public:
    TransientDiffusor( Model<dim>&,
                        const char* diffusivity,
                        const char* diffusing_variable,
                        const char* storage_variable,
                        const char* element_source_variable );
                                     
    TransientDiffusor( Model<dim>&,
                        const char* diffusivity,
                        const char* diffusing_variable,
                        const char* storage_variable,
                        const char* element_source_variable,
                        const char* point_source_variable );
                              
    TransientDiffusor( Model<dim>&, 
                       const char* lhs_diffusivity,
                       const char* rhs_diffusivity,
                       const char* diffusing_variable,
                       const char* storage_variable,
                       const char* spatial_source_variable,
                       const char* point_source_variable );
    
    // including gravity
    TransientDiffusor( Model<dim>&, 
                        const char* diffusivity,
                        const char* diffusing_variable,
                        const char* storage_variable,
                        const char* spatial_source_variable,
                        const char* gradient_variable, double gradient_multiplier );

    TransientDiffusor( Model<dim>&, 
                        const char* lhs_diffusivity,
                        const char* rhs_diffusivity,
                        const char* diffusing_variable,
                        const char* storage_variable,
                        const char* spatial_source_variable,
                        const char* gradient_variable, double gradient_multiplier );
                              
    TransientDiffusor( Model<dim>&, 
                      const char* lhs_diffusivity,
                      const char* rhs_diffusivity,
                      const char* diffusing_variable,
                      const char* storage_variable,
                      const char* spatial_source_variable,
                      const char* point_source_variable,
                      const char* gradient_variable, double gradient_multiplier );

    virtual ~TransientDiffusor();
  
    /// [C/dt + K]{p}_t+dt = [C/dt]{p}_t + q + g [d(var)/dy]
    void ComputeTransientStateFullyImplicit( Model<dim>& sg, double time_increment, bool verbose=false ); 

    virtual void AdjustSolverSettings();

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings&  GetSolverSettings();
#else
    /// add extra functionality for alternative solver if needed
#endif

  protected:

    NumIntegral_dNT_lhsop_dN_dV<dim,CELLTYPE>   conductance_;
    NumIntegral_NT_rhsop_N_dV<dim,CELLTYPE>* source_;
    NumIntegral_NT_lhsop_N_dV<dim,CELLTYPE>  capacitance_lhs_;
    NumIntegral_NT_rhsop_N_dV<dim,CELLTYPE>  capacitance_rhs_;
    PointSource_rhsop<dim,CELLTYPE>*         nodal_source_;
    NumIntegral_dNT_op_dV<dim,CELLTYPE>*     gravity_;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings_;
    SAMG_Solver    solver_;
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver_;
#endif
    const double  grad_multiplier_;
    std::string   dep_var_name_;
    
  private:
    TransientDiffusor();
};

} // end csmp

#endif
