#ifndef STENCIL_PROCESSOR_H
#define STENCIL_PROCESSOR_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"

namespace csmp {

class FV_Parameter;
template<size_t>  class Element;
template<size_t>  class TwoPhaseModel;

template<size_t dim>
struct StencilProcessor {
    ~StencilProcessor();
    /// advection only
    StencilProcessor( const csmp::Index& adv1_key,    // the advected nodal property
                      const csmp::Index& velo_key,    // transport velocity
                      const csmp::Index& src_key,     // nodal source  
                      const csmp::Index& velo_mult_key); // velocity multiplier

    /// advection and diffusion
    StencilProcessor( const csmp::Index& adv1_key,          // the advected nodal property
                      const csmp::Index& velo_key,          // the Darcian transport velocity
                      const csmp::Index& diffusivity_key,   // diffusivity of adv. prop.
                      const csmp::Index& src_key,           // nodal source  
                      const csmp::Index& velo_mult_key);    // velocity multiplier

    StencilProcessor( const StencilProcessor& tfs );

    /// Advection Diffusion Equation

    /// first-order scheme for advection-diffusion equation:  sector volumes, advected variable values, and facet fluxes
    void InitializeFirstOrder(const FV_Parameter& param,
                              const Element<dim>& e ,
                              const VARIABLE_TYPE &vt=SCALAR, const size_t var_comp_nr=0);
    
    /// Computes Flux Mismatch for FVs at the Boundaries
    void ComputeBoundaryFluxMismatch(const FV_Parameter& param, const Element<dim>& e);
     
    /// second-order scheme for advection-diffusion equation: sector volumes, advected variable values, facet fluxes, and facet advected variable values
    void InitializeSecondOrder( const FV_Parameter& param,
                                const Element<dim>& e );

    /// second-order in time scheme for advection-diffusion equation: advected variable values at nodes and interpolated to facets
    void InitializeAdvectedVariableValues( const Element<dim>& e );


    /// for advection (isat1)
    void IsotropicallyLimitTransportProperties( const Element<dim>& e, const std::vector<std::pair<double,double> >& SMINMAX );

    void ApplyLeastSquareMethodToLimitTransportProperties( const Element<dim>& e,
                                                           const std::vector<std::pair<double,double> >& SMINMAX,
                                                           const csmp::Index& mass_center_key,
                                                           const csmp::Index& grad_sn_key,
                                                           const csmp::Index& grad_sn_limiter_key );

    /// second-order version
    double   LimitExplicitSourceTerm( const Element<dim>& e,
                                  const std::vector<std::pair<double,double> >& SMINMAX,
                                  double inside_var_value,
                                  double outside_var_value,
                                  double psi_dash_f,        // average value at facet
                                  double flux );

    /// theta-limiting (in time)
    void EvaluateThetaValues( const Element<dim>& e,
                              const std::vector<double>&  FVPOREVOL, // vols of FV's
                              const std::vector<double>&  SAT0, // advected variable value at t0
                              const std::vector<std::vector<double> >& FACETFLUXES0, // old timestep
                              const std::vector<std::vector<double> >& LTDSATS0,     // old timestep
                              double time_increment,
                              bool  use_max_theta );

    /// Incompressible Two Phase Flow Equation


    /// first-order scheme for two phase flow ( linearazed form ): as in reference documentation
    void AccumulateImplicitTwoPhaseSolution1( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm );

    /// second-order scheme for two phase flow ( linearized form ): as in InitializeSecondOrder for single phase
    void AccumulateImplicitTwoPhaseSolution2( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm );

    /// first-order scheme for two-phase flow ( nonlinear Newton Raphson approach ):
    void AccumulateImplicitTwoPhaseSolution1_NonlinearNewtonRaphson( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading);

    /// second-order scheme for two-phase flow ( nonlinear Newton Raphson approach, standart limiter ):
    void AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson( const std::vector<std::pair<double,double> >& SMINMAX,
                                              const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading);


    /// second-order scheme for two-phase flow ( nonlinear Newton Raphson approach, lsm gradient limiter ):
    void AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson( const std::vector<std::pair<double,double> >& SMINMAX,
                                              const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading,
                                              const csmp::Index& mass_center_key,
                                              const csmp::Index& grad_sn_key,
                                              const csmp::Index& grad_sn_limiter_key);

    /// inflow & outflow boundary correction for two-phase flow
    void CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson(const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              size_t pnid,
                                              double& flux,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading);


     void Out() const;
    
    // public variables   
    size_t                 eidx_;               ///< local array index for the stencil
    double               diff_coeff_, velo_mult_;
    std::vector<double>  sector_pore_volume_; ///< the volume of each finite-volume sector multiplied by its porosity
    std::vector<double>  facet_flux_;         ///< scalar volume flux across finite-volume facets
    std::vector<double>  facet_flux_rhs_;
    std::vector<size_t>    upstream_node_;      ///< inside or outside node local index

    // interpolated values of advected quantity at facet integration points
    std::vector<double>  psi1_,      ///< transported variable at t+dt
                           ipsi1_,     ///< interpolated transported variable (t+dt) (can be limited)
                           theta_;     ///< limiter values for each FV facet
                           
    // fluid sources (+) or sinks (-) due to deviations from conservative fluxes
    std::vector<double>  nodal_src_; ///< nodal source/sink term
    std::vector<double>  src_;       ///< sector / divergece related sources for lhs
    std::vector<double>  rhs_src_;   ///< sector / divergece related sources for rhs
    /// node position w.r.t. the finite volume facet normal
    mutable size_t         inside_node_, outside_node_;
     
    // variable keys 
    csmp::Index            adv1_key_;      ///<  advected nodal quantity at t + dt
    csmp::Index            vel_key_;       ///< transport velocity at t
    csmp::Index            diff_key_;      ///<  diffusivity of advected nodal quantity
    csmp::Index            src_key_;       ///<   nodal source   
    csmp::Index            velo_mult_key_; ///< velocity multiplier key
    DenseMatrix<DM_MIN>    DN_,DNT_, dpcds_grad_;
    std::vector<double>  dpcdsn_, dsdn_; ///< capillary pressure and saturation gradients

  private:
    StencilProcessor();
    StencilProcessor& operator=( const StencilProcessor& tfs );
};





/**
 
@class StencilProcessor  StencilProcessor "generic_node_centered_finite_volumes\StencilProcessor.h"
@author S.K. Matthaei
@date 2004

@section motivation Motivation
 

When visiting the edges of each finite element, this object deals with
the accumulation of the fluxes, saturations etc. for the transport
schemes. 
*/


} // namespace csmp

#endif




























