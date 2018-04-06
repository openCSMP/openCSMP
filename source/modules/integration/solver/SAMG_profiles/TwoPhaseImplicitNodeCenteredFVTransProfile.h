#ifndef TWOPHASEIMPLICITNODECENTEREDFVTRANSPROFILE_H
#define TWOPHASEIMPLICITNODECENTEREDFVTRANSPROFILE_H

#include "SAMG_Profile.h"
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"

namespace csmp{

template<size_t dim>
class TwoPhaseImplicitNodeCenteredFVTransProfile : public SAMG_Profile
{
public:

    /// First-order two-phase flow constructor
    TwoPhaseImplicitNodeCenteredFVTransProfile( const char* group_name,
                                                Model<dim>& model,
                                                bool with_capillary_spreading,
                                                bool with_gravitational_forces,
                                                const char* porosity  = "porosity",
                                                const char* diffusivity = "capillary diffusivity",
                                                const char* viscosity_n = "viscosity carbonic phase",
                                                const char* viscosity_w = "viscosity aqueous phase",
                                                const char* density_n = "density carbonic phase",
                                                const char* density_w = "density aqueous phase",
                                                const char* phase1_to_update = "saturation aqueous phase",
                                                const char* phase2_to_advect = "saturation carbonic phase",
                                                const char* transp_velocity = "velocity",
                                                const char* nodal_source= "nodal fluid volume source",
                                                const char* reference_variable_to_no_flow_bc="fluid pressure",
                                                bool nonlinear_scheme = false);


    /// Second-order two-phase flow constructor
    TwoPhaseImplicitNodeCenteredFVTransProfile( const char* group_name,
                                                Model<dim>& model,
                                                bool second_order_in_space,
                                                bool with_capillary_spreading,
                                                bool with_gravitational_forces,
                                                const char* porosity  = "porosity",
                                                const char* diffusivity = "capillary diffusivity",
                                                const char* viscosity_n = "viscosity carbonic phase",
                                                const char* viscosity_w = "viscosity aqueous phase",
                                                const char* density_n = "density carbonic phase",
                                                const char* density_w = "density aqueous phase",
                                                const char* phase1_to_update = "saturation aqueous phase",
                                                const char* phase2_to_advect = "saturation carbonic phase",
                                                const char* transp_velocity = "velocity",
                                                const char* nodal_source= "nodal fluid volume source",
                                                const char* reference_variable_to_no_flow_bc="fluid pressure",
                                                bool nonlinear_scheme = false);


    /// Solve for transport equation
    virtual bool Solve( double64 modelTime );

    /// Overloaded for transport
    bool Solve( double64 modelTime, TwoPhaseModel<dim>& saturationFunctions, double64 timeInterval );

    TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>& Solver() { return advector_; }

private:
    TwoPhaseImplicitNodeCenteredFVTransProfile();

    /// Configure SAMG solver settings
    virtual bool AdjustSolverSettings();

private:
    Model<dim>& model_;
    TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor> advector_;

    bool firstCall_;
    std::string   dumpFileName_;
};


} // csmp

#endif // TWOPHASEIMPLICITNODECENTEREDFVTRANSPROFILE_H
