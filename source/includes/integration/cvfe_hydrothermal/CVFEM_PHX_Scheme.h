#ifndef CVFEM_PHX_SCHEME_H
#define CVFEM_PHX_SCHEME_H

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

#include "CVFEM_PHX_VariableNames.h"
#include "CVFEM_PHX_VariableNames.h"
#include "TwoPhaseTransportPHX.h"
#include "CVFEM_PressureGradientVisitor.h"
#include "PoreVolumeVisitor.h"
#include "NaClH2OPropertiesVisitorPHX.h"
#include "CVFEM_Visitor.h"
#include "UpwindControlVisitor.h"

#include "PropertyHandle.h"

#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "CVFEM_NumIntegral_dNT_op_dN_dV.h"
#include "CVFEM_PointSource_rhsop.h"
#include "NumIntegral_NT_lhs_nodal_op_N_dV.h"
#include "NumIntegral_NT_rhs_nodal_op_N_dV.h"
#include "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.h"
#include "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.h"

#include "ExplicitFiniteVolumeTransportPHX.h"
#include "Limiter.h"
#include "PDE_Integrator.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;
template<uint32_t> class PropertyDatabase;
//template<uint32_t, template> class PDE_Integrator;

  /**
     @class CVFEM_PHX_Scheme CVFEM_PHX_Scheme.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes
  
     @section motivation Motivation
      Class to perform calculations of CVFEM PHX scheme.

     @section usage Usage
      Uses CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
          
     @endcode
     
     @section dependencies Dependencies
      CVFEM_PHX_VariableNames
      TwoPhaseTransportPHX
      CVFEM_PressureGradientVisitor
      PoreVolumeVisitor
      NumIntegral_NT_lhsop_N_dV
      NumIntegral_NT_op_N_dV
      CVFEM_NumIntegral_dNT_op_dN_dV
      CVFEM_PointSource_rhsop
      NumIntegral_NT_lhs_nodal_op_N_dV
      NumIntegral_NT_rhs_nodal_op_N_dV
      CVFEM_Upwind_NumIntegral_dNT_op_dN_dV
      CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV
      NaClH2OPropertiesVisitorPHX
      CVFEM_Visitor
      Limiter
      PDE_Integrator
      Brine
      ConvertConcentrationUnitsNaCl
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
template<uint32_t dim>
class CVFEM_PHX_Scheme {
  public:

    explicit CVFEM_PHX_Scheme( Model<dim>&, bool with_gravity = true );
    ~CVFEM_PHX_Scheme( );
    
    void SetLargestTimeStep( double timestep ); // modifying maximum size of time step
    void ChangeTimeStepTo( double timestep ); // adjusting timestep
    void InitialFluidPropertiesFromPTX(); // initialize fluid properties from current PTX conditions
    void PrepareTransientCalculations(); // preparation before transient calculations, calculating pressure gradient and updwind nodes from current status
    double Apply(); // main function to apply CVFEM scheme in transeint calculations, returns time step used for calculations

    void TemperatureDependentHeatCapacityRock( double cpr_min_ext, double t_min_ext,
                                               double cpr_max_ext, double t_max_ext ); // modify calculations of temperature-dependent heat capacity of the rock
    void OpenBoundaries( double T_gradC, double p_Pa, double wt ); // switch on open top, specifying temperature, pressure and salinity of inflowing fluid
    void OpenBoundaries( double wt ); // switch on open top, specifying salinity of inflowing fluid
    void CheckConsistency( bool check ); // switch consostency check on or off - currently disabled

    void SetBrickWallLimiterTo( bool limit ); // switch brick wall limiter for fluid pressure on or off

    void AddAdvectionVariable( const char* balanced_variable,
                               const char* new_lhs_liquid, const char* new_rhs_liquid,
                               const char* new_lhs_vapor,  const char* new_rhs_vapor  ); // add further variables for FV calculations
    void Adjust_CFL_Criterion( double scale_factor, bool take_pore_velocity); // modifying cfl criterion

    void GetFacetFluxFromInsideNodeToOutsideNode( Element<dim>& e, unsigned int facet_idx,
                                                  double& flux_liquid, double& flux_vapor ); // access function to transient fluxes


private:

// sequence matters  
    Model<dim>& model;
    const PropertyDatabase<dim>& p_ref;
    
    bool verbose, move_on, pressure_loop;
    bool check_consistency, open_top;
    bool time_tracking;
    bool brick_wall_limiter;

    CVFEM_PHX_VariableNames names;

    CVFEM_PressureGradientVisitor<dim>    pres_grad;
    PoreVolumeVisitor<dim>                pore_visitor;
    ExplicitFiniteVolumeTransportPHX<dim> fv_transport_vapor, fv_transport_liquid;
    UpwindControlVisitor<dim>             upwind_control;
    TwoPhaseTransportPHX<dim>             transport;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings p_settings, T_settings;
    SAMG_Solver p_SAMG_solver, T_SAMG_solver;
#else
    /// add extra functionality for alternative solver if needed
    CSMP_DEFAULT_LINEAR_SOLVER p_LINEAR_solver, T_LINEAR_solver;
#endif
    PDE_Integrator<dim,Element> p_FE_SAMG;
    PDE_Integrator<dim,Element> T_FE_SAMG;
    PDE_Integrator<dim,Element> p_FE_Gauss;

    NumIntegral_NT_lhs_nodal_op_N_dV<dim> capacitance_lhs;
    CVFEM_NumIntegral_dNT_op_dN_dV<dim>   conductance;
    NumIntegral_NT_rhs_nodal_op_N_dV<dim> capacitance_rhs;
    CVFEM_PointSource_rhsop<dim>          heat_bottom;

    NumIntegral_NT_lhs_nodal_op_N_dV<dim>        capacitance_lhs_p;
    CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>   conductance_p_upwind_liquid,conductance_p_upwind_vapor;
    NumIntegral_NT_rhs_nodal_op_N_dV<dim>        capacitance_rhs_p;
    CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim> grav_liq, grav_vap;
    CVFEM_PointSource_rhsop<dim>                 source_p, source_p2;

    NaClH2OPropertiesVisitorPHX<dim>  equilibrator_properties;

    PropertyHandle<dim> diff_mass, mass, fluid_density,
                           diff_enthalpy, enthalpy, total_enthalpy,
                           bfm, bfe, bfs, dhc;
                           
    std::vector<PropertyHandle<dim>* > reset_properties;

//    CVFEM_Visitor<dim> mass_visitor, enthalpy_visitor, conduction_visitor;

    Limiter<dim>  pressure_limiter_transport, pressure_limiter_fluid;

    double dt, cfl_dt, largest_timestep;
    double current_dt, control_dt, old_dt;
    double min_value, max_value;
    double cfl_min_l, cfl_min_v, cfl_max_l, cfl_max_v;
    double temp;

    int timestep;

    void AdvanceTransientVariables(); // book keeping of variables for transient calculations
    void ResetVariables(); // reset variables for transient pressure calculations
    void FullReset(); // reset variables for transient calculations after thermal equilibration
    void AdvectionDiffusionLoops(); // outer loop including advection and pressure diffion until mass-based time step criterion is met
    void PressureLoop(); // inner loop including pressure diffion until cfl-based time step criterion is met
//    void ApplyCVFEM_Visitors(); // application of CVFEM_visitors for consistency check - currently not in use
    void FluidRockEquilibration(); // thermal quilibration between fluid and rock
//    void CheckForConsistency(); // performing consistency check between CVFEM_visitors and FV calculations - currently not in use

};

} // end namespace csmp



#endif
