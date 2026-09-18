// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef Quartz_MODEL_H
#define Quartz_MODEL_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"
#include "H2OLookup.h"
#include "compareFloats.h"
#include "LinearSolver.h"
#include <cmath>
#include <iostream>


namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/** @file QuartzModel.h
 *  @author Nicolas Krattiger
 *  @brief Model to calculate quartz precipitation and dissolution in hydrothermal system
 *  @date 2024
 */

template<uint32_t dim>
class QuartzModel : public Visitor<dim> {

public:
    QuartzModel(Model<dim>& model);

    virtual ~QuartzModel();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //functions related to transport of solutes
    void CalculateSpeciesFractionsAndMobilities(Model<dim>* model);
    void ApplyBoundaryFlowCorrections(Model<dim>* model);
    void Partitioning(Model<dim>* model);
    void InitializeQuartz(Model<dim>* model);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    virtual void Visit(Node<dim>* node);
    virtual void Visit(Model<dim>* node);
    virtual void Visit(Region<dim>* node);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //functions for calculation of geochemical equilibria
    void SolveForQuartzEquilibriumState();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //functions related to transport of solutes
    void GetFluidMasses();
    void ReadVariables( Node<dim>* node );
    void GetSpeciesContentsFluid();
    void InitializeConcentrationsAndGasMoleNumbers();
    void InitializeMineralMoleNumbers();
    void CalculateTotalMoleNumbersOfBasisSpecies();
    void ReadEquilibriumConstantsFromLiterature();
    void CheckNonnegativityOfBasisSpecies();

    void CalculateMobilities_aux (
        csmp::Index vapormobilitykey,csmp::Index liquidmobilitykey,csmp::Index vaporcontentkey,csmp::Index liquidcontentkey,
        ScalarVariable vapormobility,ScalarVariable liquidmobility,ScalarVariable vaporcontent,ScalarVariable liquidcontent,
        typename std::vector<Node<dim>*>::const_iterator nit);

    void CalculateSpeciesFractions_aux (
        csmp::Index vaporcontentkey,csmp::Index liquidcontentkey,csmp::Index vaporfractionkey,csmp::Index liquidfractionkey,
        ScalarVariable vaporcontent,ScalarVariable liquidcontent,ScalarVariable vaporfraction,ScalarVariable liquidfraction,
        typename std::vector<Node<dim>*>::const_iterator nit);

    void ApplyBoundaryFlowCorrections_aux (
        csmp::Index vaporcontentkey,csmp::Index liquidcontentkey,csmp::Index fluidcontentkey,
        ScalarVariable vaporcontent,ScalarVariable liquidcontent,ScalarVariable fluidcontent,
        typename std::vector<Node<dim>*>::const_iterator nit);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //functions related to storage of variables
    void CalculateAndStoreConcentrations( Node<dim>* node );
    void CalculateAndStoreQuartzMassesFromMoleNumbers( Node<dim>* node );

    std::string                     group_name; // if visitor is restricted to Region
    const PropertyDatabase<dim>&    prop_ref_;
    Model< dim>&                    model_ref_;

    // temperature ranges for equilibrium calculations in degC
    //  quartz
    double lowercutofftemperature_quartzcalculations = 0.;
    double uppercutofftemperature_quartzcalculations = 1000.;
    //

    csmp::Index
    temperature_key,                            // temperature
    pressure_key,                               // pressure
    xnacl_l_key,                                // salt fraction liquid
    xnacl_v_key,                                // salt fraction vapor
    pore_volume_key,                            // pore volume
    rho_l_key,                                  // liquid density
    rho_v_key,                                  // vapor density
    mt_key,                                     // fluid density (total fluid mass)
    ml_key,                                     // fluid mass liquid
    mv_key,                                     // fluid mass vapor
    sl_key,                                     // liquid saturation
    sv_key,                                     // vapor saturation
    mml_key,                                    // liquid mass mobility
    mmv_key,                                    // vapor mass mobility
    bfm_key,                                    // boundary mass flow
    meltmass_key,                               // mass of the melt
    crystalmass_key,                            // mass of the crystals

    //SiO2 species
    sio2_content_fluid_key, sio2_content_liquid_key, sio2_content_vapor_key, // sio2 content in kg sio2/pore volume
    sio2_concentration_fluid_key, sio2_concentration_liquid_key, sio2_concentration_vapor_key, sio2_concentration_melt_key, sio2_concentration_crystals_key, // sio2 concentration in mol sio2/kg water
    sio2_fraction_fluid_key, sio2_fraction_liquid_key, sio2_fraction_vapor_key,                      // sio2 mass in kg sio2/(kg sio2 + kg water)
    sio2_mass_total_key, sio2_mass_fluid_key, sio2_mass_liquid_key, sio2_mass_vapor_key, sio2_mass_melt_key, sio2_mass_crystals_key, // sio2 mass in kg
    sio2_mobility_liquid_key, sio2_mobility_vapor_key,

    //quartz
    vein_quartz_solid_key;              // proxy for quartz repricitated as vein filling

    ScalarVariable
    temperature,                            // temperature
    pressure,                               // pressure
    xnacl_l,                                // salt fraction liquid
    xnacl_v,                                // salt fraction vapor
    pore_volume,                            // pore volume
    rho_l,                                  // liquid density
    rho_v,                                  // vapor density
    mt,                                     // fluid density (total fluid mass)
    ml,                                     // fluid mass liquid
    mv,                                     // fluid mass vapor
    sl,                                     // liquid saturation
    sv,                                     // vapor saturation
    mml,                                    // liquid mass mobility
    mmv,                                    // vapor mass mobility
    bfm,                                    // boundary mass flow
    meltmass,                               // mass of the melt
    crystalmass,                            // mass of the crystals

    //SiO2 species
    sio2_content_fluid, sio2_content_liquid, sio2_content_vapor,  // sio2 content in kg sio2/pore volume
    sio2_concentration_fluid, sio2_concentration_liquid, sio2_concentration_vapor, sio2_concentration_melt, sio2_concentration_crystals, // sio2 concentration in mol sio2/kg water
    sio2_fraction_fluid, sio2_fraction_liquid, sio2_fraction_vapor,                      // sio2 mass in kg sio2/(kg sio2 + kg water)
    sio2_mass_total, sio2_mass_fluid, sio2_mass_liquid, sio2_mass_vapor, sio2_mass_melt, sio2_mass_crystals, // sio2 mass in kg
    sio2_mobility_liquid, sio2_mobility_vapor,                                                // sio2 mobility

    //quartz
    vein_quartz_solid;              // proxy for quartz repricitated as vein filling

    // declaration of concentrations, masses and mole numbers
    double liquid_mass, vapor_mass;
    double massh2oliq, massh2ogas;
    double vein_quartz_solid_mass, vein_quartz_solid_mol;

    double
            sio2_tot_mol=0, sio2_tot_concentration=0.;

    double sio2_concentration=0, sio2_gas_mol=0;

    bool is_sio2_tot_mol_zero=false;

    // equilibrium constants
    // aqueous species
    double Ksio2_lit = 0., logKsio2_lit = 0.,

    // minerals
            Kquartzaq_lit = 0., logKquartzaq_lit = 0.,
            Kquartzgas_lit = 0., logKquartzgas_lit = 0.,

    // gas/liquid partition coefficients
            D_sio2 = 0.;

    //species molar masses
    double sio2_molar_mass = 60.08e-3;   //kg/mol

};
} // end namespace csmp

#endif
