#include "SteadyStateDiffusionProfile.h"
#include "Model.h"
#include "PL_Utilities.h"

namespace csmp {

template<size_t dim>
SteadyStateDiffusionProfile<dim>::SteadyStateDiffusionProfile( Model<dim>& model,
                                                               const char* diffusivity,
                                                               const char* diffusingVariable,
                                                               const char* spatialSourceVariable )
: model_( model ),
  ssds_( model, diffusivity, diffusingVariable, spatialSourceVariable ),
  dumpFileName_( "SAMG_Diffusion_" ),
  firstCall_(true)
{
  /// Adjust SAMG Settings
  ssds_.AdjustSolverSettings();
  AdjustSolverSettings();
} // end constructor


template<size_t dim>
SteadyStateDiffusionProfile<dim>::SteadyStateDiffusionProfile( Model<dim>& model,
                                                               const char* diffusivity,
                                                               const char* diffusingVariable,
                                                               const char* spatialSourceVariable,
                                                               const char* gradientVariable,
                                                               double64 gradient_multiplier )
: model_( model ),
  ssds_( model, diffusivity, diffusingVariable, spatialSourceVariable, gradientVariable, gradient_multiplier),
  dumpFileName_( "SAMG_Diffusion_" ),
  firstCall_(true)
{
  /// Adjust SAMG Settings
  ssds_.AdjustSolverSettings();
  AdjustSolverSettings();
} // end constructor


template<size_t dim>
bool SteadyStateDiffusionProfile<dim>::AdjustSolverSettings()
{
#ifdef CSMP_WITH_SAMG_SOLVER
  std::cout <<"\n\n*** SteadyStateDiffusionProfile::AdjustSolverSettings ***\n\n";

  this->Solver().GetSolverSettings().Set_iout1( 1 );
  this->Solver().GetSolverSettings().Set_iout2( 0 );

  /// SAMG solution criterion
  this->Solver().GetSolverSettings().Set_eps(0.);
  this->Solver().GetSolverSettings().Set_rel_eps(1.E-10);

  /// SAMG output to file settings
  #ifdef SAMG_OUTPUT_TO_FILE
      this->Solver().GetSolverSettings().Set_idmp( 8 );                       // define SAMG command and file output
      this->Solver().GetSolverSettings().Set_ioform( "f" );                   // define SAMG file output format for reduced file size, idmp > 1 is required
      this->Solver().GetSolverSettings().Set_filnam_dump( dumpFileName_ );    // set filename for SAMG file output other than default "level", idmp > 1 is required
  #endif

  #ifdef SAMG_MULTIPLE_INSTANCES
        /// Select SAMG instance
        this->Solver().GetSolverSettings().SetSolverInstance(0);

        /// IMP-IMS with SAMG Multiple Instances: Reusing setup is possible
        /// Re-use solver setup from previous timestep, internal checks force setup when required
        #ifndef NO_PRIMARY_SOLVER_CONTROL
            this->Solver().GetSolverSettings().Set_iswit(7);
        #else
            /// Re-use solver setup from previous timestep, without internal checks forcing a new setup when required
            /// iswit(4) is required for first call (without purging memory) to store previous SAMG setup information
            this->Solver().GetSolverSettings().Set_iswit(4);
        #endif

        /// Agressive first level coarsening nredlev(1) for decreased setup time and reduced no. of cycles
        //this->Solver().GetSolverSettings().Set_nred(1);

        /// Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
        this->Solver().GetSolverSettings().Set_a_cmplx(2);

        /// Pre-adjust SAMG mesh complexity, based on solver output
        this->Solver().GetSolverSettings().Set_g_cmplx(1.5);
        this->Solver().GetSolverSettings().Set_w_avrge(2);

  #endif

  #ifndef SAMG_MULTIPLE_INSTANCES
     #ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
         #ifndef NO_PRIMARY_SOLVER_CONTROL
            /// Re-use solver setup from previous timestep, internal checks force setup when required
            this->Solver().GetSolverSettings().Set_iswit(7);
         #else
            /// Re-use solver setup from previous timestep, without internal checks forcing a new setup when required
            /// iswit(4) is required for first call (without purging memory) to store previous SAMG setup information
            this->Solver().GetSolverSettings().Set_iswit(4);
         #endif
     #else
     /// IMP-IMS without SAMG Multiple Instances: setup information from previous solver call has to be deleted (purge memory)
        this->Solver().GetSolverSettings().Set_iswit(5);
     #endif

  /// Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
  this->Solver().GetSolverSettings().Set_a_cmplx(2);

  /// Pre-adjust SAMG mesh complexity, based on solver output
  this->Solver().GetSolverSettings().Set_g_cmplx(1.5);
  this->Solver().GetSolverSettings().Set_w_avrge(2);

  #endif
#endif

  return true;
} // AdjustSolverSettings


/** Reconfiguration of SAMG solver settings after completing first solver call

    @section arguments Input Arguments

    Under provision of a separate SAMG instance (SAMG_MULTIPLE_INSTANCES) or applicaton of the IMPES solution
    scheme the below defined SAMG input parameters are re-adjusted to gain additional solver performance.

    - SAMG input parameters reconfigured after first SAMG solver call:
    1. Provide solution of the previous time step as initial guess -> Set_itypu
    2. Reuse SAMG solver setup if not initialized in automatic mode Set_iswit(7)
    3. Change absolute convergence criterion -> Set_eps() or
    4. Switch to relative convergence criterion -> Set_eps()

    The gained functionality concerns SAMG input arguments that are required to be set in a different way at
    the initial solver call. The if-loop "firstCall_" has been implemented to cater for changed SAMG solver
    settings after the solver routine has been called for the first time. Changed settings take effect with
    the second call of the the SAMG solver instance.

    ad 1. For faster numeric convergence SAMG settings can be changed in order to provide the result of the previous
    solver call as initial guess. The setting (Set_itypu(0) has to be called succeeding an initial solver call.

    ad 2. Handling the SAMG solver memory is done in a similar manner: The coarsening setup created at the initial
    solver call is kept in the memory with any setting other than Set_iswit(5). Subsequent calls shall reuse the setup.

    Manual (no primary solver control) setup refresh check, iswit(4) is adjusted immediately after the first SAMG solver call as seen
    below to iswit(3). In the event of exceeding a defined benchmark of iteration cycles in SAMG_Solver::CheckCycleCriterion a new setup
    will be forced in the upcoming solver call. Since the paramter iswit remains unchanged, the new setup is reused from the following
    solver calls until the criterion is failed again. The benchmark case is set within the firstCall_ sequence as well.

    Automatic primary solver control performs setup refresh check (convegence control, ...) with iswit(7) which corresponds
    iswit(3) functionality. Whenever possible and reasonable, SAMG runs with partial setup. Otherwise, SAMG makes a full setup.
    In contrast to this, the automatic control does not assume a full hierarchy to be available from a previous run. If no setup
    is available, SAMG will automatically perform a fresh setup.

    ad 3. After the first solution has converged with an absolut convergence set to the round-off error, a less strict
    convergence criterion, either absolute or residual, is subsequently used as stopping criterion. Absolute convergence
    is defined by "res <= eps".

    ad.4 Relative convergence is defined by "res <= eps.res0" (res0 = starting residual).
    */
template<size_t dim>
bool SteadyStateDiffusionProfile<dim>::Solve( double64 modelTime )
{ 
    /// SAMG output to file
    #ifdef SAMG_OUTPUT_TO_FILE
        std::string currentDumpFileName( dumpFileName_ );
        currentDumpFileName.append( numberToString( static_cast<size_t>( modelTime ) ) );
        this->Solver().GetSolverSettings().Set_filnam_dump( currentDumpFileName );
    #endif

    /// Solve
    std::cout <<"\n\n*** SteadyStateDiffusionProfile::Solve ***\n\n";
    ssds_.ComputeSteadyState( model_ );

#ifdef CSMP_WITH_SAMG_SOLVER
    #ifdef SAMG_MULTIPLE_INSTANCES
        if ( firstCall_ ) {
           /// Use solution of previous timestep as an initial guess
          this->Solver().GetSolverSettings().Set_itypu(0);
          std::cout << "\n\n*** SteadyStateDiffusionProfile::Solve: 'First Call: Set SAMG input parameter itypu = 0";

          /// Reuse SAMG solver setup, this step is only necessary if initial setting was iswit(4)
          #ifdef NO_PRIMARY_SOLVER_CONTROL
              this->Solver().GetSolverSettings().Set_iswit(3); // re-use solver setup from previous timestep
              std::cout << ", iswit = 3,";
          #endif

          #ifdef RELATIVE_CONVERGENCE
              /// Relative convergence is used as stopping criterion "res <= eps.res0" (res0 = starting residual) - for IMPES without SAMG Multiple Instances only
              this->Solver().GetSolverSettings().Set_eps( this->Solver().GetSolverSettings().Get_rel_eps() );
              std::cout << " and relative solution criterion eps = " << this->Solver().GetSolverSettings().Get_eps() << " after first solver call' ***\n\n";
          #else
              /// Absolute convergence is used as stopping criterion "res <= eps" (res0 = starting residual)
              this->Solver().GetSolverSettings().Set_eps( 1.E-14 );
              std::cout << " and absolute solution criterion eps = " << -this->Solver().GetSolverSettings().Get_eps() << " after first solver call' ***\n\n";
          #endif
          firstCall_ = false;
        }
    #endif

    #ifndef SAMG_MULTIPLE_INSTANCES
         if ( firstCall_ ) {
         #ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES

              /// Use solution of previous timestep as an initial guess - for IMPES without SAMG Multiple Instances only
              this->Solver().GetSolverSettings().Set_itypu(0);
              std::cout << "\n\n*** SteadyStateDiffusionProfile::Solve 'First Call: Set SAMG input parameter itypu = 0' ***\n";

              /// Reuse SAMG solver setup - for IMPES without SAMG Multiple Instances only - this step is only necessary if initial setting was iswit(4)
              #ifdef NO_PRIMARY_SOLVER_CONTROL
                  this->Solver().GetSolverSettings().Set_iswit(3); // re-use solver setup from previous timestep
                  std::cout << "\n\n*** SteadyStateDiffusionProfile::Solve 'First Call: Set SAMG input parameter iswit = 3' ***\n";
              #endif

          #endif

          #ifdef RELATIVE_CONVERGENCE
              /// Relative convergence is used as stopping criterion "res <= eps.res0" (res0 = starting residual) - for IMPES without SAMG Multiple Instances only
              this->Solver().GetSolverSettings().Set_eps( this->Solver().GetSolverSettings().Get_rel_eps() );
              std::cout << "\n\n*** SteadyStateDiffusionProfile::Solve 'First Call: Set SAMG relative solution criterion eps = " << this->Solver().GetSolverSettings().Get_eps() << "' ***\n\n";
          #else
              /// Absolute convergence is used as stopping criterion "res <= eps" (res0 = starting residual)
              this->Solver().GetSolverSettings().Set_eps( 1.E-14 );
              std::cout << "\n\n*** SteadyStateDiffusionProfile::Solve 'First Call: Set SAMG absolute solution criterion eps = " << -this->Solver().GetSolverSettings().Get_eps() << "' ***\n\n";
          #endif

          firstCall_ = false;
          }
    #endif
#endif

    return true;
} // Solve

template class SteadyStateDiffusionProfile<1U>;
template class SteadyStateDiffusionProfile<2U>;
template class SteadyStateDiffusionProfile<3U>;

} // csmp
