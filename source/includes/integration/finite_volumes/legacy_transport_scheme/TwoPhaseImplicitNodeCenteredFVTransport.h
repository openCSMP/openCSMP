// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TWO_PHASE_IMPLICIT_NODE_CENTERED_FINITE_VOLUME_TRANSPORT_H
#define TWO_PHASE_IMPLICIT_NODE_CENTERED_FINITE_VOLUME_TRANSPORT_H

#include "CSMP_definitions.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "NodeCenteredFiniteVolumeAlgorithm.h"

namespace csmp {

template<uint32_t dim, template<uint32_t> class STP>
class TwoPhaseImplicitNodeCenteredFVTransport : public NodeCenteredFiniteVolumeTransport<dim> {
  public:

    /// First-order two-phase flow constructor, no capillary spreading, wetting phase=phase1
    TwoPhaseImplicitNodeCenteredFVTransport( const char* group_name,
                                             Model<dim>&,
                                             const char* porosity="porosity",
                                             const char* viscosity_n="viscosity carbonic phase",
                                             const char* viscosity_w="viscosity aqueous phase",
                                             const char* density_n="density carbonic phase",
                                             const char* density_w="density aqueous phase",
                                             const char* phase1="saturation aqueous phase",
                                             const char* phase2="saturation carbonic phase",
                                             const char* transp_velocity="velocity",
                                             const char* nodal_source="nodal fluid volume source",
                                             bool with_capillary_spreading  = false,
                                             bool with_gravitational_forces = false,
                                             const char* reference_variable_to_no_flow_bc="fluid pressure",
                                             bool nonlinear_scheme = false);

    /// Second-order two-phase flow constructor, no capillary spreading, wetting phase=phase1
    TwoPhaseImplicitNodeCenteredFVTransport( const char* group_name,
                                             Model<dim>&,
                                             bool second_order_in_space,
                                             bool second_order_in_time,
                                             const char* porosity="porosity",
                                             const char* viscosity_n="viscosity carbonic phase",
                                             const char* viscosity_w="viscosity aqueous phase",
                                             const char* density_n="density carbonic phase",
                                             const char* density_w="density aqueous phase",
                                             const char* phase1="saturation aqueous phase",
                                             const char* phase2="saturation carbonic phase",
                                             const char* transp_velocity="velocity",
                                             const char* nodal_source="nodal fluid volume source",
                                             bool with_capillary_spreading  = false,
                                             bool with_gravitational_forces = false,
                                             const char* reference_variable_to_no_flow_bc="fluid pressure",
                                             bool nonlinear_scheme = false);

    /// First-order two-phase flow constructor, with capillary spreading modelled as non-linear diffusion
    TwoPhaseImplicitNodeCenteredFVTransport( const char* group_name,
                                             Model<dim>&,
                                             const char* porosity,
                                             const char* diffusivity, // <- difference to first constructor
                                             const char* viscosity_n,
                                             const char* viscosity_w,
                                             const char* density_n,
                                             const char* density_w,
                                             const char* phase1,
                                             const char* phase2,
                                             const char* transp_velocity,
                                             const char* nodal_source,
                                             bool with_capillary_spreading  = false,
                                             bool with_gravitational_forces = false,
                                             const char* reference_variable_to_no_flow_bc="fluid pressure",
                                             bool nonlinear_scheme = false);

    /// Second-order two-phase flow constructor, with capillary spreading modelled as non-linear diffusion
    TwoPhaseImplicitNodeCenteredFVTransport( const char* group_name,
                                             Model<dim>&,
                                             bool second_order_in_space,
                                             bool second_order_in_time,
                                             const char* porosity,
                                             const char* diffusivity, // <- difference to first constructor
                                             const char* viscosity_n,
                                             const char* viscosity_w,
                                             const char* density_n,
                                             const char* density_w,
                                             const char* phase1,
                                             const char* phase2,
                                             const char* transp_velocity,
                                             const char* nodal_source,
                                             bool with_capillary_spreading  = false,
                                             bool with_gravitational_forces = false,
                                             const char* reference_variable_to_no_flow_bc="fluid pressure",
                                             bool nonlinear_scheme = false);


    /// for field-scale simulation: with standard variables, gravity and thickness attribute for lower-dim elements
    TwoPhaseImplicitNodeCenteredFVTransport( const char* group_name,
                                             Model<dim>&,
                                             bool second_order_in_space,
                                             bool second_order_in_time,
                                             const char* transp_velocity,
                                             bool with_capillary_forces,
                                             const char* thickness_low_dim_elmt ="thickness",
                                             bool nonlinear_scheme = false );


    virtual ~TwoPhaseImplicitNodeCenteredFVTransport();
    
    /// transport of a non-wetting phase in 2-phase flow
    virtual double TransportPhase( TwoPhaseModel<dim>&, double time_interval );
                               
    /// two-phase, takes into account element shape in 2D and 3D unless one deals with 1D domain
    virtual double AnisotropicCourantIncrement( TwoPhaseModel<dim>&, double max_time_increment ); 

    virtual void AdjustSolverSettings();

#ifdef CSMP_WITH_SAMG_SOLVER
    /// Access to advector's samg solver settings
    virtual SAMG_Settings& GetSolverSettings();
#else
    /// add extra functionality for alternative solver if needed
#endif

    void  DisableCapillarySpreading();
    void  EnableCapillarySpreading();
    void  DisableGravitationalForces();
    void  EnableGravitationalForces();

    void  SetNoFlowBoundaryConditionKey(Model<dim>&, const char* no_flow_bc="fluid pressure");
    void  ApplyLinearScheme();
    void  ApplyNonlinearScheme();

    void  MaxNewtonRahsonIterations(size_t max_newton_raphson_iterations);
    void  MaxLineSearchIterations(size_t max_line_search_iterations);
    void  TargetNewtonRaphsonResidual(double target_newton_raphson_residual);

  private:

    /// test for correct variable type and placement
    virtual void CheckTransportVariables() const; 
     
    bool IsInteriorStencil( const Element<dim>* const ) const;

    void MinMaxAdvectedProperty();

    void MinMaxAdvectedPropertyIncludingTheCurrentNode();

    void MinMaxAdvectedPropertyExceptTheCurrentNode();

    void ComputePiecewiseConstantNodalSource();

    void AddSourceTerm();


    /// advects non-wetting phase in 2-phase flow (1st-order accurate)
    void SolveTransportEquation1stOrder( TwoPhaseModel<dim>&,
                                         double time_increment,
                                         bool account_for_nodal_sources = true );

    /// advects non-wetting phase in 2-phase flow without gravity (1st-order accurate)
    void SolveTransportEquation1stOrder_NonlinearNewtonRaphson( TwoPhaseModel<dim>& relperm,
                                                                double time_increment);

    /// advects non-wetting phase in 2-phase flow without gravity (2nd-order accurate)
    void SolveTransportEquation2ndOrderInSpace_NonlinearNewtonRaphson( Model<dim>&, TwoPhaseModel<dim>& relperm,
                                                                       double time_increment);
 
  private:
  
    mutable std::vector<double>  LHSSRC, RHSSRC, /// < source terms for limiting
                                   OMEGA;          /// < balancing parameter for LHS and RHS

    double target_newton_raphson_residual_;
    size_t max_newton_raphson_iterations_, max_line_search_iterations_;

    const csmp::Index  rhn_key,
                       rhw_key,
                       mun_key,
                       muw_key,
                       ph1_key; ///< the non-wetting phase that is not advected but updated

    bool with_capillary_spreading_,
         with_gravitational_forces_;

    csmp::Index  reference_variable_to_no_flow_bc_key_;

    bool nonlinear_scheme_;

    // vectors used in line-search algorithm
    std::vector<double>  SN_, DS_;

    NodeCenteredFiniteVolumeAlgorithm<dim>  advector_;
    std::set<csmp::Element<dim>*>  halo_stencils_;
};



/**

@class TwoPhaseImplicitNodeCenteredFVTransport TwoPhaseImplicitNodeCenteredFVTransport.h
"generic_node_centered_finite_volumes/TwoPhaseImplicitNodeCenteredFVTransport.h"

@author S.K. Matthaei
@date 2002
 

@section motivation Motivation

The TwoPhaseImplicitNodeCenteredFVTransport module is provided to allow users
to perform passive advection and advection of a single non-wetting 
phase in two-phase flow simulations. It represents the fastest but
CFL dependent transport scheme for arbitrary FE meshes including 
different FE types and mixes of volume and surface elements.  

The scheme is also very tolerant in that no specific boundary conditions 
need to be applied. The value of the advected variable at inflow 
boundaries is not modified and governs what gets transported into
the model domain.  

 
@section design Design Intent

The scheme has been designed to capitalize on the new generic polytype
CVFE technology that is part of CSMP.  
 

@section applicability Applicability

Use the scheme if the transport velocity field is rapidly changing or
invalidated by the transport itself. In this case the scheme is faster
that the implicit scheme embodied in the FiniteVolumeTransport module
which starts to pay off at about CFL=20 or greater.  

Do not use this scheme when the CFL condition based transport increment
leads global advection steps that are smaller than about 1/1,000 of the 
model length in the direction of the flow. In this case the computations
would be prohibitively slow.
 

@section structure Structure

The user will only interact with the TwoPhaseImplicitNodeCenteredFVTransport module,
but the FiniteVolumeAlgorithm is used to manage the accumulation
and output of the solution in each transport step, and the 
FiniteVolumeProcessor class facilitates the generation of the 
finite-volume stencils from each finite element, so that the accumulation
can proceed on an element by element basis and no FV mesh needs to be 
stored.
 
 
@section collaboration Collaboration

The scheme collaborates with the finite elements of the current model.
 
 
@section consequences Consequences

This implementation of the scheme cannot be restricted to Regions, but 
always needs to be applied to the entire model, i.e. Model.  

Currently it also only works for 3D models. Surface computations can
still be done if the mesh consists entirely of surface elements in
3D space. 

*/
 
} // end namespace csmp

#endif
