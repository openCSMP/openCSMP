#ifndef TWO_PHASE_EXPLICIT_NODE_CENTERED_FV_TRANSPORT_H
#define TWO_PHASE_EXPLICIT_NODE_CENTERED_FV_TRANSPORT_H

#include "NodeCenteredFiniteVolumeTransport.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class TwoPhaseModel;

/// two-phase explicit transport scheme
template<size_t dim, template<size_t> class STP>
class TwoPhaseExplicitNodeCenteredFVTransport : public NodeCenteredFiniteVolumeTransport<dim> {
  public:
    TwoPhaseExplicitNodeCenteredFVTransport( const char* region_name,
                                             Model<dim>&,
                                             const char* porosity="porosity",
                                             const char* diffusivity="capillary diffusivity",
                                             const char* advected_prop1="saturation aqueous phase",
                                             const char* advected_prop2="saturation carbonic phase",
                                             const char* transp_velocity="velocity",
                                             const char* nodal_source="nodal fluid volume source",
                                             bool second_order_accuracy=false,
                                             bool with_capillary_spreading=true,
                                             bool with_gravitational_forces=false,
                                             const char* reference_variable_to_no_flow_bc="fluid pressure",
                                             const char* thickness_attribute=NULL );

    virtual ~TwoPhaseExplicitNodeCenteredFVTransport();
    
    /// 2-phase flow, transport of the non-wetting phase
    virtual double TransportPhase( TwoPhaseModel<dim>&, double time_interval );
                               
    /// two-phase, takes into account element shape in 2D and 3D unless one deals with 1D domain
    virtual double AnisotropicCourantIncrement( TwoPhaseModel<dim>&, double max_time_increment );

    void  DisableCapillarySpreading();
    void  EnableCapillarySpreading();
    void  DisableGravitationalForces();
    void  EnableGravitationalForces();

    void  SetNoFlowBoundaryConditionKey(Model<dim>&, const char* no_flow_bc="fluid pressure");

private:

    bool IsInteriorStencil( const Element<dim>* const ) const;

    /// writes results back to model
    double OutputResults( const csmp::Index& adv_key, bool show_range ) const;

    double OutputResults( const csmp::Index& adv1_key,
                            const csmp::Index& adv2_key, bool show_range, bool do_range_check ) const;


    void MinMaxAdvectedProperty();

    void MinMaxAdvectedPropertyIncludingTheCurrentNode();

    void MinMaxAdvectedPropertyExceptTheCurrentNode();

    /// builds solution taking into account nodal sources and sinks
    void Compose2PhaseSolution( TwoPhaseModel<dim>&, double time_interval, bool divergence_free_correction = false );


    /// compensate inflow & outflow
    bool FractionalFlowThroughBoundaryFiniteVolume( const Node<dim>*,
                                                    TwoPhaseModel<dim>&,
                                                    double& inflow,
                                                    double& flux_balance ) const;

    void AssignFractionalFlowBoundaryConditions( TwoPhaseModel<dim>& );

    void AssignGenericFlowBoundaryConditions( TwoPhaseModel<dim>& );

    /// Makes divirgence free correction
    void DivergenceFreeCorrection( TwoPhaseModel<dim>&, double time_interval );


    /// advects non-wetting phase in 2-phase flow (1st-order accurate)
    void AdvectVariable1stOrder( TwoPhaseModel<dim>&,
	                               double time_increment );

    /// advects non-wetting phase in 2-phase flow (2nd-order accurate)
    void AdvectVariable2ndOrder( TwoPhaseModel<dim>&,
	                               double time_increment );

 
  private:

    std::set<csmp::Element<dim>*>  halo_stencils_;
    STP<dim>                       stencil_;
    std::vector<double>          RESULT;

    bool                           with_capillary_spreading_,
                                   with_gravitational_forces_;

    csmp::Index                    sw_key_,
                                   no_flow_bc_key_;
};


/**
 
@class TwoPhaseExplicitNodeCenteredFVTransport TwoPhaseExplicitNodeCenteredFVTransport.h
"generic_node_centered_finite_volumes/TwoPhaseExplicitNodeCenteredFVTransport.h" 

@author S.K. Matthaei
@author S. Geiger
@date 2002
 
@section motivation Motivation 

The TwoPhaseExplicitNodeCenteredFVTransport module performs 
advection a single non-wetting phase in two-phase flow simulations. 
It is a CFL dependent transport scheme for arbitrary FE meshes including 
different FE types and mixes of volume and surface elements.  

The scheme is tolerant in that no specific boundary conditions 
need to be applied. The value of the advected variable at inflow 
boundaries is not modified and governs what gets transported into
the model domain.  
 
@section applicability Applicability

Use the explicit scheme if the transport velocity field is rapidly changing /
invalidated by the transport itself. In this case the scheme is faster
than the implicit scheme. The use of the latter will
start to pay off at about CFL=20 or greater.  

Do not use this scheme when the CFL condition based transport increment
leads global advection steps that are smaller than about 1/1,000 of the 
model length in the direction of the flow. In this case the computations
would be prohibitively slow.
 
@section structure Structure

The explicit advection scheme consists out of the three modules which are
documented in a power-point documentation in the documentation folder.  
 
TwoPhaseExplicitNodeCenteredFVTransport
ExplicitStencilProcessor
  
The user will only interact with the TwoPhaseExplicitNodeCenteredFVTransport module.
The ExplicitStenProcessor class facilitates the generation of the 
finite-volume stencil contributions from each finite element, so that the accumulation
can proceed on an element by element basis. No FV mesh needs to be stored.
 
 
@section implementation Implementation

Flux-balances are computed for all finite volumes including those 
which are located on the boundaries of the Region or Model to 
which the transport scheme is applied. They are compensated for also
on all domain boundaries.  


@section examples Application Examples

This is an example of how the second order accurate scheme is applied:  
 
@code

TwoPhaseExplicitNodeCenteredFVTransport<1U,ExplicitStencilProcessor>( "Model", 
                                                                       reservoir_model, 
                                                                      "porosity", 
                                                                      "capillary diffusivity", 
                                                                      "saturation water", "saturation oil", 
                                                                      "velocity", 
                                                                      "nodal fluid volume source",
                                                                       second_order_accurate );

 cout <<"\n\nTransportPhase: Configuring TRANSPORT simulation: ";
 if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
 else                cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
 cout <<"\nThe grid Courant number is "<< explicit_advector.AnisotropicCourantIncrement( reservoir_model ) << endl;
 cout <<"\nEnter advection time: ";
 double time_interval;
 cin >> time_interval;
 cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
 clock_t ticks = clock();
 explicit_advector.TransportPhase( relperm_model, time_interval );
 ticks = clock() - ticks;
 cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

@endcode 

@todo (3) Remove NCFVTAlgorithm depencency (B)
 
*/


} // end namespace csmp

#endif

