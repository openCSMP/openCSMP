#ifndef STEADY_STATE_DIFFUSOR_H
#define STEADY_STATE_DIFFUSOR_H

#include "CSMP_definitions.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "PointSource_rhsop.h"
#include "PDE_Integrator_CRM.h"

namespace csmp {

template<size_t> class Model;


/** Spatial integration of steady-state partial-differential diffusion eqns: k div^2=q+..., 

    @attention Gradient multiplier is currently not used and has no influence on results.
    @todo (3) Make generic so that it can take both numerically and analytically integrated FEs
*/
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
class SteadyStateDiffusor : public PDE_Integrator_CRM<dim,COMPUTATION_DOMAIN> {
  public:
    typedef typename COMPUTATION_DOMAIN<dim>::CellType  ComputationCell;

	  SteadyStateDiffusor(  Model<dim>&,
                          const char* diffusivity,
                          const char* diffusing_variable,
                          const char* spatial_source_variable,
                          bool LumpedRHS);

    SteadyStateDiffusor( Model<dim>&, 
                                const char* diffusivity,
                                const char* diffusing_variable,
                                const char* spatial_source_variable );
                              
    SteadyStateDiffusor( Model<dim>&, 
                                const char* diffusivity,
                                const char* diffusing_variable,
                                const char* spatial_source_variable,  // on element or integration points
                                const char* point_source_variable );  // on the nodes

    /// integrating gradient terms applied to the righthand sides; no sources
    SteadyStateDiffusor( Model<dim>&, 
                                const char* diffusivity,
                                const char* diffusing_variable,
                                const char* gradient_variable, double64 gradient_multiplier );

    /// as above, but with FE source terms
    SteadyStateDiffusor( Model<dim>&, 
                                const char* diffusivity,
                                const char* diffusing_variable,
                                const char* spatial_source_variable,
                                const char* gradient_variable, double64 gradient_multiplier );

    SteadyStateDiffusor( Model<dim>&, 
                                const char* lhs_diffusivity,
                                const char* rhs_diffusivity,
                                const char* diffusing_variable,
                                const char* spatial_source_variable,
                                const char* gradient_variable, double64 gradient_multiplier );

    virtual ~SteadyStateDiffusor();

    /// K div^2 = q
    void ComputeSteadyState( COMPUTATION_DOMAIN<dim>& subDomain, bool verbose=false );

    virtual void AdjustSolverSettings();

#ifdef CSMP_WITH_SAMG_SOLVER
    /// Access to SAMG solver settings owned by SteadyStateDiffusor object
    SAMG_Settings& GetSolverSettings();
#else
    /// add extra functionality for alternative solver if needed
#endif

    
  protected:

    NumIntegral_dNT_op_dN_dV<dim,ComputationCell>   conductance_;
    NumIntegral_NT_op_N_dV<dim,ComputationCell>*    source_;
    PointSource_rhsop<dim,ComputationCell>*         nodal_source_;
    NumIntegral_dNT_op_dV<dim,ComputationCell>*     gravity_;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings                   settings_; // only this derived class knows about SAMG
#else
    /// add extra functionality for alternative solver if needed
#endif

    const double64                  grad_multiplier_;
    std::string                     dep_var_name_;

    
  private:
    SteadyStateDiffusor();
    SteadyStateDiffusor( const SteadyStateDiffusor& );
#ifdef CSMP_WITH_SAMG_SOLVER
    void Adjust_SAMG_ForSubsequentSolves();
#else
    /// add extra functionality for alternative solver if needed
#endif

  private:
    bool   firstCall_;
};

} // end csmp

#endif 
