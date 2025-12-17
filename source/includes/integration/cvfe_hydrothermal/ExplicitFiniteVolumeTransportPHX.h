#ifndef EXPLICIT_FINITE_VOLUME_TRANSPORT_PHX_H
#define EXPLICIT_FINITE_VOLUME_TRANSPORT_PHX_H

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

#include "NodeCenteredFiniteVolumeTransport.h"
#include "StencilProcessorPHX.h"

namespace csmp {

template<uint32_t> class Model;

  /**
     @class ExplicitFiniteVolumeTransportPHX ExplicitFiniteVolumeTransportPHX.h

     Finitive volume calculations for one phase within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes
  
     @section motivation Motivation
      Specialized version of explicit finite volume transport class for PHX scheme.

     @section usage Usage
      To be used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Performs mass balanced finite volume tansport
          
     @endcode
     
     @section dependencies Dependencies
	 The functions are tailored for use in TwoPhaseTransportPHX.
	 Needs StencilProcessorPHX
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
template<uint32_t dim>
class ExplicitFiniteVolumeTransportPHX : public NodeCenteredFiniteVolumeTransport<dim> {

  public:

    ExplicitFiniteVolumeTransportPHX( Model<dim>& model,
                                      std::vector<std::string>& lhs_property, // vector with left-hand side variables for finite volume transport
                                      std::vector<std::string>& rhs_property, // vector with right-hand side variables for finite volume transport - matching lhs_property ordering
                                      const char* velocity, // "velocity" in this particular case should only be K * grad( P )
                                      const char* source, // currently not used but needed for constructor of base class
                                      const char* density); // density of the fluid phase

                                       
    virtual ~ExplicitFiniteVolumeTransportPHX();

    double  AdvectMassConservedSinglePhase( const double& time_increment ); // For use with single-phase fluid
	                                                                           // the function takes a time step suggestions
	                                                                           // performs transport with this time step
	                                                                           // if flux with this time step violates mass-based time step criterium, flux is performed with a smaller time step
	                                                                           // the time step used for transport calculations is returned
    void  AddAdvectionVariable( const char* new_lhs, const char* new_rhs); // Adding further variables for finite volume calculations
	void WithGravityComponent( ); // activate gravity component
    void  SetMaximumTimeStep( const double& new_max_time_step ); // manually set the maximum time step size

	// the remaining function are used by TwoPhaseTransportPHX
    csmp::Index GetDensityKey( ); 
    void UpdateProjection( );     // updating projection of "velocity", i.e. pressure gradient, onto facte normals
    
	void  DetermineFacetFlux( const double& time_increment,
		                      std::vector<DenseMatrix<DM_MIN> >& upwind_visitor );
    void  DetermineFacetFluxSinglePhase( );
    void  AdjustAndPerformFacetFlux( const double& time_factor );

    double GetPropertyValue( size_t idx );
    double GetFluxOut( size_t idx ); // return the main property
    double GetFluxIn(  size_t idx ); // returns the main property
    double GetFluxOut( size_t idx, uint32_t property ); // returns the indicated property of the property-vector
    double GetFluxIn(  size_t idx, uint32_t property ); // returns the indicated property of the property-vector
    double GetMainPropertyLHS( size_t idx );
    double GetFacetFlux( Element<dim>& e, uint32_t facet_idx, uint32_t property = 0U);
    // TODO: not defined: std::vector<double> GetUpwindCoefficients( size_t eidx );
    double GetProjectedVelocities( Element<dim>& e, uint32_t i );

    double GetFacetNormalVelocity( size_t element, uint32_t facet );
    double GetFacetNormalComponent( size_t element, uint32_t facet, uint32_t x_or_y_or_z );
    double GetFacetArea( size_t element, uint32_t facet );

  protected:

// TODO: Why store these in here?
    Model<dim>& model_ref;
    Region<dim>& region_ref;
    
    StencilProcessorPHX<dim>  stencilPHX;

    void  ComposeAdvection( );
    void  WriteResults( );
    void  CalculateOutflowPerPoreVolume();
    void  CalculateInflowPerPoreVolume();
    void  AdjustFluxOut( const double& time_factor );
    void  CalculateFluxIn();
    void  GetUpwindMatrix( const Element<dim>& e, std::vector<DenseMatrix<DM_MIN> >& upwind_visitor );
    void  CheckDryFiniteVolumesAndAdjustTimestepForSinglePhase( );

    double  internal_time_step, max_time_step, single_phase_time_step_factor;
    bool with_gravity;
    DenseMatrix<DM_MIN>  upwind;
    
    csmp::Index               k_key, rho_key, pv_key;
    std::vector<csmp::Index>  pl_key, pr_key;

    // flux vectors
    std::vector<std::vector<double> > flux_in_vectors, flux_out_vectors,
                                      property_vectors, mass_balance_vectors;
                                        
    std::vector<std::vector<std::vector<double> > > facet_flux_vectors;
};

} // end namespace csmp

#endif
