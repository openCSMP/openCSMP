#include "TwoPhaseImplicitNodeCenteredFVTransProfile.h"
#include "PL_Utilities.h"

#ifdef SAMG_WITH_CSMP_SOLVER
namespace csmp{


/// First-order constructor

template<size_t dim>
TwoPhaseImplicitNodeCenteredFVTransProfile<dim>::
              TwoPhaseImplicitNodeCenteredFVTransProfile( const char* group_name,
                                                          Model<dim>& model,
                                                          bool with_capillary_spreading,
                                                          bool with_gravitational_forces,
                                                          const char* porosity,
                                                          const char* diffusivity,
                                                          const char* viscosity_n,
                                                          const char* viscosity_w,
                                                          const char* density_n,
                                                          const char* density_w,
                                                          const char* phase1_to_update,
                                                          const char* phase2_to_advect,
                                                          const char* transp_velocity,
                                                          const char* nodal_source,
                                                          const char* reference_variable_to_no_flow_bc,
                                                          bool nonlinear_scheme)
: model_( model ),
  advector_( group_name, model, porosity, diffusivity, viscosity_n, viscosity_w, density_n, density_w,
             phase1_to_update, phase2_to_advect, transp_velocity, nodal_source, with_capillary_spreading, with_gravitational_forces ,reference_variable_to_no_flow_bc, nonlinear_scheme),
  dumpFileName_( "SAMG_Advection_" ),
  firstCall_(true)

{
  /// Adjust SAMG settings
  this->Solver().AdjustSolverSettings();
  AdjustSolverSettings();
}

/// Second-order constructor

template<size_t dim>
TwoPhaseImplicitNodeCenteredFVTransProfile<dim>::
              TwoPhaseImplicitNodeCenteredFVTransProfile( const char* group_name,
                                                          Model<dim>& model,
                                                          bool second_order_in_space,
                                                          bool with_capillary_spreading,
                                                          bool with_gravitational_forces,
                                                          const char* porosity,
                                                          const char* diffusivity,
                                                          const char* viscosity_n,
                                                          const char* viscosity_w,
                                                          const char* density_n,
                                                          const char* density_w,
                                                          const char* phase1_to_update,
                                                          const char* phase2_to_advect,
                                                          const char* transp_velocity,
                                                          const char* nodal_source,
                                                          const char* reference_variable_to_no_flow_bc,
                                                          bool nonlinear_scheme)
: model_( model ),
  advector_( group_name, model, second_order_in_space, false, porosity, diffusivity, viscosity_n, viscosity_w, density_n, density_w,
             phase1_to_update, phase2_to_advect, transp_velocity, nodal_source, with_capillary_spreading, with_gravitational_forces , reference_variable_to_no_flow_bc, nonlinear_scheme),
  dumpFileName_( "SAMG_Advection_" ),
  firstCall_(true)

{
  /// Adjust SAMG settings
  this->Solver().AdjustSolverSettings();
  AdjustSolverSettings();
}


template<size_t dim>
bool TwoPhaseImplicitNodeCenteredFVTransProfile<dim>::AdjustSolverSettings() {

  std::cout <<"\n\n*** TwoPhaseImplicitNodeCenteredFVTransProfile::AdjustSolverSettings ***\n\n";

  this->Solver().GetSolverSettings().Set_iout1( 1 );
  this->Solver().GetSolverSettings().Set_iout2( 0 );

  /// SAMG solution criterion
  this->Solver().GetSolverSettings().Set_eps(0.);
  this->Solver().GetSolverSettings().Set_rel_eps(1.E-10);

  /// SAMG output to file
  #ifdef SAMG_OUTPUT_TO_FILE
      this->Solver().GetSolverSettings().Set_idmp( 8 );                    // define SAMG command and file output
      this->Solver().GetSolverSettings().Set_ioform( "f" );                // define SAMG file output format for reduced file size, idmp > 1 is required
      this->Solver().GetSolverSettings().Set_filnam_dump( dumpFileName_ ); // set filename for SAMG file output other than default "level", idmp > 1 is required
  #endif

  #ifndef SAMG_MULTIPLE_INSTANCES
     /// IMP-IMS without SAMG Multiple Instances: setup information from previous solver call has to be deleted (purge memory)
     this->Solver().GetSolverSettings().Set_iswit(5);
  #endif

  /// Renounce coarsening
  /**
  @section arguments Input Arguments

  Gauss-Seidel preconditioned BI-CGstab without coarsening is permitted only for
  nealy diagonal matrices
  */
  #ifdef RENOUNCE_COARSENING
     this->Solver().GetSolverSettings().Set_levelx(1);
  #endif

  #ifdef SAMG_MULTIPLE_INSTANCES

     /// Select SAMG instance
     this->Solver().GetSolverSettings().SetSolverInstance(1);

     /// IMP-IMS with SAMG Multiple Instances: iswit(4) required for first call (without purging memory)
     this->Solver().GetSolverSettings().Set_iswit(4);

  #endif

  return true;
} // AdjustSolverSettings

/** Reconfiguration of SAMG solver settings after first call

@section arguments Input Arguments

Under provision of a separate SAMG instance (SAMG_MULTIPLE_INSTANCES) or applicaton of the IMPES solution
scheme the below defined SAMG input parameters are re-adjusted to gain additional solver performance.

- SAMG input parameters reconfigured after first SAMG solver call:
1. Provide solution of the previous time step as initial guess -> SAMG_Settings::Set_itypu
2. Switch to relative convergence criterion -> SAMG_Settings::Set_eps()

The gained functionality concerns SAMG input arguments that are required to be set in a different way at
the initial solver call. The if-loop "firstCall_" has been implemented to cater for changed SAMG solver
settings after the solver routine has been called for the first time. Changed settings take effect with
the second call of the the SAMG solver instance.

For faster numeric convergence SAMG settings can be changed in order to provide the result of the previous
solver call as initial guess. The setting (Set_itypu(0) has to be called succeeding an initial solver call.

Re-use of the coarsening setup created at previous solver calls is not feasible since matrix entries are
not constant. The SAMG solver input parameter iswit therefore remains at its initial setting.

After the first solution has converged with an absolut convergence set at the round-off error, a less strict
convergence criterion, either absolute or residual, is subsequently used as stopping criterion. Absolute convergence
is defined by "res <= eps". Relative convergence is defined by "res <= eps.res0" (res0 = starting residual).
*/
template<size_t dim>
bool TwoPhaseImplicitNodeCenteredFVTransProfile<dim>::Solve( double modelTime )
{
  throw;
}

template<size_t dim>
bool TwoPhaseImplicitNodeCenteredFVTransProfile<dim>::Solve( double modelTime, TwoPhaseModel<dim>& saturationFunctions, double timeInterval )
{
    /// SAMG output to file
    #ifdef SAMG_OUTPUT_TO_FILE
        std::string currentDumpFileName( dumpFileName_ );
        currentDumpFileName.append( numberToString( static_cast<size_t>( modelTime ) ) );
        this->Solver().GetSolverSettings().Set_filnam_dump( currentDumpFileName );
    #endif

    /// Solve
    std::cout <<"\n\n*** TwoPhaseImplicitNodeCenteredFVTransProfile::Solve ***\n\n";
    advector_.TransportPhase( saturationFunctions, timeInterval );


    #ifdef SAMG_MULTIPLE_INSTANCES
        if ( firstCall_ ) {

          /// Use solution of previous timestep as an initial guess
          this->Solver().GetSolverSettings().Set_itypu(0); // use solution of previous timestep as an initial guess
          std::cout << "\n\n*** TwoPhaseImplicitNodeCenteredFVTransProfile::Solve 'First Call: Set itypu = 0 ";

        #ifdef RELATIVE_CONVERGENCE
            /// Relative convergence is used as stopping criterion "res <= eps.res0" (res0 = starting residual) - for IMPES without SAMG Multiple Instances only
            this->Solver().GetSolverSettings().Set_eps( this->Solver().GetSolverSettings().Get_rel_eps() );
            std::cout << " and relative solution criterion eps = " << this->Solver().GetSolverSettings().Get_eps() << " after first solver call' ***\n\n";
        #else
            /// Absolute convergence is used as stopping criterion "res <= eps" (res0 = starting residual)
            //this->Solver().GetSolverSettings().Set_eps( 1.E-10 );
            //std::cout << " and absolute solution criterion eps = " << -this->Solver().GetSolverSettings().Get_eps() << " after first solver call' ***\n\n";
        #endif

          firstCall_ = false;
        }
    #endif

    return true;
}

template class TwoPhaseImplicitNodeCenteredFVTransProfile<1U>;
template class TwoPhaseImplicitNodeCenteredFVTransProfile<2U>;
template class TwoPhaseImplicitNodeCenteredFVTransProfile<3U>;

} // csmp

#endif // CSMP_WITH_SAMG_SOLVER
