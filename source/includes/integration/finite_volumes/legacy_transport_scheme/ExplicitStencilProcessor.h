#ifndef EXPLICIT_STENCIL_PROCESSOR_H
#define EXPLICIT_STENCIL_PROCESSOR_H

#include "DenseMatrix.h"
#include "Index.h"

namespace csmp {

class FV_Parameter; 
template<size_t> class Element;
template<size_t> class TwoPhaseModel;

/// for the calculation of FV equations to  be solved explicitly
template<size_t dim>
struct ExplicitStencilProcessor {
    /// for advection only
    ExplicitStencilProcessor( const csmp::Index& adv_key, 
                              const csmp::Index& velo_key );
    /// for advection diffusion                          
    ExplicitStencilProcessor( const csmp::Index& adv_key, 
                              const csmp::Index& velo_key,
                              const csmp::Index& diffusivity_key );
                              
    virtual ~ExplicitStencilProcessor();


    /// Advection Diffusion Equation

    /// first-order accurate solution for advection equation
    virtual void InitializeFirstOrder( const FV_Parameter& param, const Element<dim>& e ,const VARIABLE_TYPE &vt=SCALAR, const size_t var_comp_nr=0);

    /// limited second-order explicit advection equation ( standart limiter, without diffusion )
    void AccumulateExplicitAdvectionSolution2( 
                                     const std::vector<std::pair<double64,double64> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     std::vector<double64>& res);

    /// limited second-order explicit advection equation ( lsm gradient limiter, without diffusion )
    void AccumulateExplicitAdvectionSolution2(
                                     const std::vector<std::pair<double64,double64> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     std::vector<double64>& res,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_psi_key,
                                     const csmp::Index& grad_psi_limiter_key);

    /// limited second-order explicit advection-diffusion equation ( standart limiter )
    void AccumulateExplicitAdvectionDiffusionSolution1(
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     std::vector<double64>& res, 
                                     const VARIABLE_TYPE &vt=SCALAR, 
                                     const size_t var_comp_nr=0);

    /// limited second-order explicit advection-diffusion equation ( standart limiter )
    void AccumulateExplicitAdvectionDiffusionSolution2( 
                                     const std::vector<std::pair<double64,double64> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     std::vector<double64>& res);

    /// limited second-order explicit advection-diffusion equation ( lsm gradient limiter )
    void AccumulateExplicitAdvectionDiffusionSolution2(
                                     const std::vector<std::pair<double64,double64> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     std::vector<double64>& res,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_psi_key,
                                     const csmp::Index& grad_psi_limiter_key);

    /// Incompressible Two Phase Flow Equation

    /// first-order accurate solution for two phase flow ( only viscous flow )
    void AccumulateExplicitTwoPhaseSolution1_Visc( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              std::vector<double64>& res);

    /// second-order accurate solution for two phase flow ( standart limiter, only viscous flow )
    void AccumulateExplicitTwoPhaseSolution2_Visc( const std::vector<std::pair<double64,double64> >& SMINMAX,
                                                   const FV_Parameter& param,
                                                   const Element<dim>& e,
                                                   TwoPhaseModel<dim>& relperm,
                                                   std::vector<double64>& res);


    /// second-order accurate solution for two phase flow ( lsm gradient limiter, only viscous flow )
    void AccumulateExplicitTwoPhaseSolution2_Visc( const std::vector<std::pair<double64,double64> >& SMINMAX,
                                                   const FV_Parameter& param,
                                                   const Element<dim>& e,
                                                   TwoPhaseModel<dim>& relperm,
                                                   std::vector<double64>& res,
                                                   const csmp::Index& mass_center_key,
                                                   const csmp::Index& grad_sn_key,
                                                   const csmp::Index& grad_sn_limiter_key);

    /// first-order accurate solution for two-phase flow with gravitational and capillary effects
    void AccumulateExplicitTwoPhaseSolution1( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              std::vector<double64>& res,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading);

    /// second-order accurate solution for two phase flow ( standart limiter )
    void AccumulateExplicitTwoPhaseSolution2( const std::vector<std::pair<double64,double64> >& SMINMAX,
                                              const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              std::vector<double64>& res);

    /// second-order accurate solution for two phase flow ( lsm gradient limiter )
    void AccumulateExplicitTwoPhaseSolution2( const std::vector<std::pair<double64,double64> >& SMINMAX,
                                              const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              std::vector<double64>& res,
                                              const csmp::Index& mass_center_key,
                                              const csmp::Index& grad_sn_key,
                                              const csmp::Index& grad_sn_limiter_key);


    /// second-order accurate solution for two phase flow with standart limiter ( with gravity forces )
    void AccumulateExplicitTwoPhaseSolution2 ( const std::vector<std::pair<double64,double64> >& SMINMAX,
                                              const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              std::vector<double64>& res,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading);

    /// second-order accurate solution for two phase flow with lsm gradient limiter ( with gravity forces )
    void AccumulateExplicitTwoPhaseSolution2 ( const std::vector<std::pair<double64,double64> >& SMINMAX,
                                              const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              std::vector<double64>& res,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading,
                                              const csmp::Index& mass_center_key,
                                              const csmp::Index& grad_sn_key,
                                              const csmp::Index& grad_sn_limiter_key);


    /// inflow & outflow boundary correction for two phase flow
    void AccumulateExplicitTwoPhaseSolutionAtBoundary(const FV_Parameter& param,
                                              const Element<dim>& e,
                                              TwoPhaseModel<dim>& relperm,
                                              size_t pnid,
                                              double64& flux,
                                              bool with_gravity_forces,
                                              bool with_capillary_spreading);
                                    
    void Out() const { Out(std::cout); }
    void Out(std::ostream& os) const;

    // local array index for the stencil
    size_t               eidx_;
    // the volume of each finite-volume sector multiplied by its porosity
    std::vector<double64>      sector_pore_volume_; // pore volumes, index corresponds to nodes
    // scalar volume flux across finite-volume facets
    std::vector<double64>      facet_flux_;
    // interpolated values of advected quantity at facet integration points
    std::vector<double64>      psi1_;  // transported variable at t, and t+dt
    // fluid sources (+) or sinks (-) due to deviations from conservative fluxes
    std::vector<double64>      src_, pc_; // divergence for each facet
    // node position w.r.t. the finite volume facet normal
    mutable size_t       inside_node_, outside_node_;
     
    csmp::Index  adv1_key_;   // advected nodal quantity at t + dt
    csmp::Index  vel_key_;   // transport velocity at t
    csmp::Index  diff_key_;   // diffusivity of advected nodal quantity

    std::vector<double64>  grad_, pc_grad_, dsdn_;
    DenseMatrix<DM_MIN>    DN_;

  protected:
    ExplicitStencilProcessor();
  
  private:
    // the usual stuff
    ExplicitStencilProcessor( const ExplicitStencilProcessor& tfs );
    ExplicitStencilProcessor& operator=( const ExplicitStencilProcessor& tfs );
};


/**
 
@class ExplicitStencilProcessor  ExplicitStencilProcessor "generic_node_centered_finite_volumes\ExplicitStencilProcessor.h"
@author S.K. Matthaei
@date 2004

 
@section motivation Motivation

When visiting the edges of each finite element, this object deals with
the accumulation of the fluxes, saturations etc. for the transport
schemes. 
*/


} // namespace csmp

#endif




























