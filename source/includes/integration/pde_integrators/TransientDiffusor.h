#ifndef TRANSIENT_DIFFUSOR_H
#define TRANSIENT_DIFFUSOR_H

#include "CSMP_definitions.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "PDE_Integrator.h"

namespace csmp {

template<uint32_t> class Model;


/** transient diffusion solver with Backward-Euler time-stepping

  [C/dt + K]{p}_t+dt = [C/dt]{p}_t + q   
   ------------------------------------
   specific to the CSP reservoir simulator
   only a single equation is solved for total pressure
   result variable 'fluid pressure'
   initialises hydrostatic pressure to zero
   
   @attention gradient_multiplier currently is not used.
   
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
class TransientDiffusor : public PDE_Integrator<dim,COMPUTATION_DOMAIN> {
  public:
    typedef typename COMPUTATION_DOMAIN<dim>::CellType  ComputationCell;

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

    NumIntegral_dNT_op_dN_dV<dim,ComputationCell>   conductance_;
    NumIntegral_NT_op_N_dV<dim,ComputationCell>*    source_;
    NumIntegral_NT_lhsop_N_dV<dim,ComputationCell>  capacitance_lhs_;
    NumIntegral_NT_op_N_dV<dim,ComputationCell>     capacitance_rhs_;
    PointSource_rhsop<dim,ComputationCell>*         nodal_source_;
    NumIntegral_dNT_op_dV<dim,ComputationCell>*     gravity_;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings                   settings_;
#else
    /// add extra functionality for alternative solver if needed
#endif

    const double                  grad_multiplier_;
    std::string                     dep_var_name_;
    
  private:
    TransientDiffusor();
};

} // end csmp

#endif
