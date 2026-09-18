// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAGMA_MODEL_H
#define MAGMA_MODEL_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"
#include "compareFloats.h"
#include "H2ONaClFluidProperties.h"
#include "ConvertConcentrationUnitsNaCl.h"


namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/**
@author BLC
@date 2016
*/

template<size_t dim>
class MagmaModel : public Visitor<dim> {

public:
    MagmaModel(Model<dim>& model,
                               double initial_wtpercent_volatile_in_magma,
                               double crystal_density,
                               double crystallization_curve_exponent,
                               double initial_volatile_vol_frac_chamber,
                               double compressibility_rock,
                               double compressibility_magma,
                               double sigma1,
                               int saturation_model,
                               double saturation,
                               double salinity_magmatic_fluids);

    virtual ~MagmaModel();

    void ReadVariables(Node<dim>* node);
    void StoreVariables(Node<dim>* node);
    bool ConsistencyChecks1(bool pause_at_end, double density_deformation, double deformation_overpressure);
    bool ConsistencyChecks2(bool pause_at_end);
    bool FractionSumChecks(bool pause_at_end);
    void InitializeMagma(Model<dim>* node, std::string regionname);
    void CalculateAverageCrystallinity(Model<dim>* mdl, std::string regionname);
    void WithHaliteTrend(double fraction_magma);
    void SetTimeIncrement(double &time_step );
    virtual void Visit(Node<dim>* node);
    virtual void Visit(Model<dim>* node);
    virtual void Visit(Region<dim>* node);

private:
    std::string                group_name;  // if visitor is restricted to Region
    const PropertyDatabase<dim>&  prop_ref_;
    Model< dim>&                  model_ref_;
    csmp::Index

    liquidus_temperature_key,                    // volume deformation when fluid pressure is above threshold pressure
    solidus_temperature_key,                     // volume deformation when fluid pressure is above threshold pressure

    volume_deformation_key,           // volume deformation when fluid pressure is above threshold pressure
    total_volume_deformation_key,
    rock_density_scaling_key,         // rock density scaling to keep constant rock mass during volume deformation
    // the problem is that volume deformation translates into a change in 1-phi,
    // the latter being used as a proxy for rock volume later on (thermal equilibrator)

    nodal_density_rock_key,           // density rock
    threshold_pressure_key,           // threshold pressure
    lithostatic_pressure_key,         // lithostatic pressure
    bulk_volume_key,                  // nodal volume (control volume used in FV scheme)
    minimum_pore_volume_key,             // baseline pore volume - crystallization evolution with potential contraction of magma at ex<0.4
    minimum_mass_fluid_key,              //
    melt_volume_key,                  // melt volume
    crystal_volume_key,               // crystal volume
    delta_melt_volume_key,            // change in melt volume over a crystallization step (over a timestep)
    delta_crystal_volume_key,         // change in crystal volume over a crystallization step (over a timestep)
    crystallinity_key,                // crystallinity (change with temperature according to a law)
    mass_volatile_produced_key,       // mass of fluid exsolved (positive) or dissolved (negative) given a change in crystallinity
    total_mass_volatile_produced_key, // the total mass of fluid produced over the course of the simulation
    melt_mass_key,                    // melt mass
    crystal_mass_key,                 // crystal mass
    pressure_key,                     // fluid pressure
    magma_density_key,                // magma density
    melt_density_key,                 // melt density
    temperature_key,                  // temperature
    previous_temperature_key,         // previous temperature
    temperature_difference_key,       // change in temperature over one timestep
    nodal_porosity_key,               // nodal porosity
    volatile_dissolved_fraction_key,  // volatile dissolved fraction in melt in wt% (can potentially vary with pressure)
    volatile_dissolved_mass_key,       // mass of volatile dissolved in melt
    volatile_dissolved_eq_key,         // volatile dissolved fraction in melt at equilibrium occordind to Henry's law (pressure dependence)
    pore_volume_key,                  // pore volume
    previous_pore_volume_key,         // previous pore volume
    volatile_saturation_key,          // volatile saturation
    relative_saturation_key,          // relative saturation
    volume_change_factor_key,         // scaling parameter for "volumic" variables derived from the change in pore volume
    crystal_fraction_key,             // crystal volume fraction
    melt_fraction_key,                // melt volume fraction
    volatile_fraction_key,            // volatile volume fraction = nodal porosity, better for visualization of the chamber region
    volumetric_mass_key,                            // fluid density
    Htp_key,                           // previous total enthalpy
    msp_key,                           // previous mass salt
    m_key,                             // mass
    H_key,                             // enthalpy
    hCl_key, //enthalpy content liquid
    hCv_key, //enthalpy content vapor
    mm_key,                            // magmatic fluid mass
    mf_key,                            // fluid mass
    ml_key,                            // fluid mass liquid
    mv_key,                            // fluid mass vapor
    mms_key,                           // magmatic mass salt
    cCf_key,                            // copper content fluid
    znCf_key,                            // zinc content fluid
    kcl_ms_key,                        // mass KCl
    mpr_key,                           // magmatic production rate
    nodal_compressibility_rock_key,   // nodal_compressibility_rock
    minimum_nodal_density_rock_key,
    wt_key,                            // salinity
    average_crystallinity_key,
    delta_average_crystallinity_key,
    ddt_magmatic_ratio_key,
    depth_key,                            // salinity
    porous_flag_key;                // nodal porous flag (flow region);

    ScalarVariable volume_deformation;
    ScalarVariable total_volume_deformation;
    ScalarVariable rock_density_scaling;
    ScalarVariable minimum_rock_density_scaling;
    ScalarVariable nodal_density_rock;
    ScalarVariable nodal_density_rock_before;
    ScalarVariable threshold_pressure;
    ScalarVariable lithostatic_pressure;
    ScalarVariable nodal_compressibility_rock;

    ScalarVariable bulk_volume;
    ScalarVariable minimum_pore_volume;
    ScalarVariable minimum_mass_fluid;
    ScalarVariable melt_volume;
    ScalarVariable melt_volume_before;
    ScalarVariable crystal_volume;
    ScalarVariable crystal_volume_before;
    ScalarVariable delta_melt_volume;
    ScalarVariable delta_melt_volume_density;
    ScalarVariable delta_crystal_volume;
    ScalarVariable crystallinity;
    ScalarVariable mass_volatile_produced;
    ScalarVariable total_mass_volatile_produced;
    ScalarVariable melt_mass;
    ScalarVariable melt_mass_before;
    ScalarVariable crystal_mass;
    ScalarVariable crystal_mass_before;
    ScalarVariable pressure;
    ScalarVariable magma_density;
    ScalarVariable melt_density;
    ScalarVariable temperature;
    ScalarVariable previous_temperature;
    ScalarVariable temperature_difference;
    ScalarVariable fluid_density;
    ScalarVariable nodal_porosity;
    ScalarVariable volatile_dissolved_fraction;
    ScalarVariable volatile_dissolved_mass_before;
    ScalarVariable volatile_dissolved_mass;
    ScalarVariable volatile_dissolved_eq;
    ScalarVariable pore_volume;
    ScalarVariable volatile_saturation;
    ScalarVariable relative_saturation;
    ScalarVariable volume_change_factor;
    ScalarVariable crystal_fraction;
    ScalarVariable melt_fraction;
    ScalarVariable volatile_fraction;
    ScalarVariable ddt_magmatic_ratio;
    ScalarVariable average_crystallinity;
    ScalarVariable delta_average_crystallinity;


    ScalarVariable crystallinity_before;
    ScalarVariable nodal_porosity_before;
    ScalarVariable volatile_dissolved_fraction_before;
    ScalarVariable pore_volume_before;
    ScalarVariable previous_pore_volume;

    ScalarVariable volumetric_mass;
    ScalarVariable ml;                            // fluid mass liquid
    ScalarVariable mv;                          // fluid mass vapor
    ScalarVariable Htp;
    ScalarVariable msp;
    ScalarVariable phi;
    ScalarVariable mass;
    ScalarVariable enthalpy;
    ScalarVariable hCl;
    ScalarVariable hCv;
    ScalarVariable magmatic_mass;
    ScalarVariable fluid_mass;
    ScalarVariable magmatic_salt;
    ScalarVariable nmf;
    ScalarVariable mpr;
    ScalarVariable tmf;
    ScalarVariable mass_KCl;
    ScalarVariable cCf, znCf;
    ScalarVariable injected_fluid;
    ScalarVariable minimum_nodal_density_rock;
    ScalarVariable wt;
    ScalarVariable depth;
    ScalarVariable porous_flag;

    ScalarVariable liquidus_temperature;
    ScalarVariable solidus_temperature;


    double       initial_wtpercent_volatile_in_magma_;
    double       melt_density_;
    double       crystal_density_;
    double       crystallization_curve_exponent_;
    double       initial_volatile_vol_frac_chamber_;
    double       liquidus_temperature_;
    double       solidus_temperature_;
    double       salinity_magmatic_fluids;
    double       dt;
    double       fraction_KCl_magma;
    double       pre_b, tem_b, sal_b, x_NaCl_b;
    double       sigma1_;
    double       compressibility_rock_;
    double       compressibility_magma_;
    double       brine_enthalpy;
    double        saturation;
    int             saturation_model;

    bool           with_halite_trend, power_law_,sill_emplaced;
    bool           with_pressure_volume_change, with_copper, with_zinc;

    H2ONaClFluidProperties fluid;
};
} // end namespace csmp

#endif
