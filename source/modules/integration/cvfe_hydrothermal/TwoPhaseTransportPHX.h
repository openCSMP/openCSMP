#ifndef TWO_PHASE_TRANSPORT_PHX_H
#define TWO_PHASE_TRANSPORT_PHX_H

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

#include "ExplicitFiniteVolumeTransportPHX.h"
#include "UpwindControlVisitor.h"

namespace csmp {

/// class to handle two-phase flow caclulations of FV part within the CVFEM scheme (Weis et al., Geofluids, 2014).

template<size_t dim>
class TwoPhaseTransportPHX {
  public:

    TwoPhaseTransportPHX( Model<dim>& model,
                          UpwindControlVisitor<dim>& upwind_visitor, // visitor defining upwind nodes
                          ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor, // class for transport with vapor phase
                          ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid); // class for transport with liquid phase

    ~TwoPhaseTransportPHX( );
    
    double64  AdvectMassConserved( const double64& time_increment); // main function to coordinate two-phase flow
    void SetLargestTimeStep( const double64& max_time_step ); // modify maximum size of time step
    void WithGravityComponentLiquidAndVapor(); // active gravity component
    void UpdateProjection( ); // update projection of velocity onto facet normal
    void AddAdvectionVariable( const char* new_lhs_liquid, const char* new_rhs_liquid,
                               const char* new_lhs_vapor,  const char* new_rhs_vapor ); // adding further variables for FV calculations

  private:
  
    Model<dim>& model_ref;
    UpwindControlVisitor<dim>& UpwindVisitor;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_vapor;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_liquid;

    void  DetermineFacetFluxes( ); // calculating facet fluxes for both phases
    void  CheckDryFiniteVolumesAndAdjustTimestep( ); // checking mass-based time step criterion and adjusting time step if necessary
    void  PerformFacetFluxes( ); // finalize FV calculations
    void  TrackFluxes( ); // special function to track fluxes of magmatic fluids

    double64    cfl_dt, time_step_factor;
    double64    upper_shell_T, lower_shell_T;

    // tracking fluxes
    csmp::Index          tffi_key, tffo_key,
                        mlfi_key, mlfo_key,
                        mvfi_key, mvfo_key,
                        slfi_key, slfo_key,
                        svfi_key, svfo_key,
                        T_key;

	ScalarVariable  tffi, tffo, mlfi, mlfo, mvfi, mvfo;
    ScalarVariable  slfi, slfo, svfi, svfo, temperature;

};

/**
     @class TwoPhaseTransportPHX TwoPhaseTransportPHX.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Class for two-phase flow calculations within CVFEM scheme.

     @section usage Usage
      To be used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Controls two-phase flow calculations by coordating flow of vapor and liquid phase
          
     @endcode
     
     @section dependencies Dependencies
	 Needs ExplicitFiniteVolumeTransportPHX and UpwindControlVisitor
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */


} // end namespace csmp



#endif
