// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PermeabilityVisitor_h
#define PermeabilityVisitor_h

#include "Visitor.h"
#include "Model.h"
#include "Element.h"
#include "Exception.h"

#include "compareFloats.h"

#include <vector>

using namespace std;
using namespace csmp;

namespace csmp
{

template<uint32_t dim>
class PermeabilityVisitor: public Visitor<dim>
{

public:
    PermeabilityVisitor( Model<dim> &model);

    ~PermeabilityVisitor();

    virtual void Visit(Region<dim> *r);
    virtual void Visit(Element<dim> *e);
    virtual void Visit(Model<dim> *m);

    //! Setter functions
    // depth-dependent
    void DepthDependent( bool d_dep );
    void SetReferenceDepthTo( bool ref_dep, double refernce_depth_from_bottom );

    // temperature-dependent
    void TemperatureDependent( bool T_dep);
    void ChangeBrittleDuctileTransitionTemperature( bool change_values,
                                                   double  transition_start,
                                                   double transition_ductile, double  transition_end, double  log_k_start,
                                                   double  log_k_duct, double log_k_end );

    // pore-fluid-factor / hydrofracturing
    void PoreFluidFactorDependent( bool Pff_dep, bool average_Pff_ );
    void Hydrofracturing( bool hydro_frac, double log_perm_max, double log_perm_min);
    void SetImmediateClosureTo( bool closure );
    void NoFracturingBelowIntrusion( bool no_frac_below_intrusion_, std::string region_name );

    // mineral precipitation/dissolution
    void MineralDependent( bool min_dep );

    // anisotropy
    void AnisotropicPermeabilityTensor( bool anisotropic_k_ );

    // intrusion region tracking
    void SetIntrusionRegionID( uint32_t intrusion_ID );

    //! Standalone batch operations (not related to Visit function)
    void AssignPermeabilityTensor(Model<dim> &model);
    void ComputeSimpleFracturableFlag();

private:

    PermeabilityVisitor();

    Model<dim> &model;

    Index
        // core / bookkeeping
        permeability_key_,
        vertical_permeability_key_,
        horizontal_permeability_key_,
        permeability_tensor_key_,
        permeability_ID_key_,
        log_k_increase_key_,
        fracturing_ref_key_,
        fracturing_events_key_,
        k_start_key_,
        region_ID_key_,

        // temperature-dependent
        temperature_e_key_,

        // pore-fluid-factor / hydrofracturing
        fluid_pressure_key_,
        failure_pressure_key_,

        // mineral precipitation/dissolution
        pore_volume_element_key_,
        bulk_volume_element_key_,
        porosity_key_,

        // anisotropy
        anisotropy_factor_key_,

        // grad-P scaling
        gp_sc_key_;


    //! Functions called within the visitor
    double CalculateDepthDependentPermeability( Element<dim> *e );
    double CalculateTemperatureDependentPermeability( Element<dim> *e );
    double CalculatePoroAndPermChangeByMineralPrecipitationAndDissolution(Element<dim> *e);
    void   CalculatePoreFluidFactor( Element<dim> *e, bool average_pff );
    double CalculatePoreFluidFactorDependentPermeability( Element<dim> *e );
    double CalculateHydrofracturedPermeability( Element<dim> *e );
    void   CalculateGradPScaling(Element<dim> *e);


    bool
        // intrusion region tracking
        with_intrusion_,

        // depth-dependent
        depth_dependent_,
        with_reference_depth_,

        // temperature-dependent
        temperature_dependent_,

        // pore-fluid-factor / hydrofracturing
        pore_fluid_factor_dependent_,
        average_Pff_,
        overpressured_,
        hydro_fracturing_,
        immediate_closure_,
        no_frac_below_intrusion_,

        // mineral precipitation/dissolution
        mineral_dependent_,

        // anisotropy
        anisotropic_k_,

        // grad-P scaling
        avoid_shallow_and_reverse_gradients_;


    double
        // core / bookkeeping
        log_k_start_,
        log_k_,
        log_k_previous_,
        log_k_min_,

        // depth-dependent
        log_k_D_,
        reference_depth_,

        // temperature-dependent
        log_k_T_,
        log_k_back_,
        log_k_ductile_,
        hydro_failure_temperature_,
        ductile_failure_temperature_,
        lithos_failure_temperature_,

        // pore-fluid-factor / hydrofracturing
        log_k_Pff_,
        log_k_max_,
        log_k_Hydrofracture_,
        pff_,
        pore_fluid_factor_,
        chamber_y_min_,

        // mineral precipitation/dissolution
        log_k_mineral_;


    uint32_t
        // intrusion region tracking
        intrusion_ID_;


    std::vector<ScalarVariable>
        // pore-fluid-factor / hydrofracturing
        fluid_pressure_vec_,
        failure_pressure_vec_;


    ScalarVariable
        // core / bookkeeping
        k_,
        kV_, kH_,
        log_k_increase_,
        k_ID_,
        fracturing_ref_,

        // temperature-dependent
        T_,

        // mineral precipitation/dissolution
        porosity_,
        pore_volume_element_,
        bulk_volume_element_,
        deltaQuartzPrecipitated_element_,
        deltaGoldPrecipitated_element_,

        // anisotropy
        anisotropy_factor_,

        // grad-P scaling
        gradP_sc_;


    //densities (specific masses)
    double au_density_          = 19320.;      //kg/m³
    double quartz_density_      = 2650.;       //kg/m³
    double bulk_rock_density_   = 2500.;       //kg/m³


};
} //csp
#endif
