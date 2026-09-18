// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "QuartzModel.h"
#include "Model.h"
#include "Region.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;
namespace csmp {

template<uint32_t dim>
QuartzModel<dim>::QuartzModel(Model<dim>& model)

    : Visitor<dim>(MODEL, NODE),

      model_ref_(model),
      prop_ref_(model.Database()),

      //general variables keys
      temperature_key                (model.Database().StorageKey("temperature")),
      pressure_key                   (model.Database().StorageKey("fluid pressure")),
      pore_volume_key                (model.Database().StorageKey("pore volume")),
      rho_l_key                      (model.Database().StorageKey("density liquid")),
      rho_v_key                      (model.Database().StorageKey("density vapor")),
      sl_key                         (model.Database().StorageKey("saturation liquid")),
      sv_key                         (model.Database().StorageKey("saturation vapor")),
      ml_key                         (model.Database().StorageKey("fluid mass liquid")),
      mv_key                         (model.Database().StorageKey("fluid mass vapor")),
      mml_key                        (model.Database().StorageKey("liquid mass mobility")),
      mmv_key                        (model.Database().StorageKey("vapor mass mobility")),
      bfm_key                        (model.Database().StorageKey("boundary flow mass")),
      meltmass_key                   (model.Database().StorageKey("melt mass")),
      crystalmass_key                (model.Database().StorageKey("crystal mass")),

      sio2_fraction_liquid_key        (model.Database().StorageKey("sio2 fraction liquid")),
      sio2_fraction_vapor_key         (model.Database().StorageKey("sio2 fraction vapor")),
      sio2_fraction_fluid_key         (model.Database().StorageKey("sio2 fraction fluid")),
      sio2_content_liquid_key         (model.Database().StorageKey("sio2 content liquid")),
      sio2_content_vapor_key          (model.Database().StorageKey("sio2 content vapor")),
      sio2_content_fluid_key          (model.Database().StorageKey("sio2 content fluid")),
      sio2_mobility_liquid_key        (model.Database().StorageKey("liquid sio2 mobility")),
      sio2_mobility_vapor_key         (model.Database().StorageKey("vapor sio2 mobility")),
      sio2_concentration_liquid_key   (model.Database().StorageKey("sio2 concentration liquid")),
      sio2_concentration_vapor_key    (model.Database().StorageKey("sio2 concentration vapor")),
      sio2_concentration_fluid_key    (model.Database().StorageKey("sio2 concentration fluid")),
      sio2_mass_liquid_key            (model.Database().StorageKey("sio2 mass liquid")),
      sio2_mass_vapor_key             (model.Database().StorageKey("sio2 mass vapor")),
      sio2_mass_fluid_key             (model.Database().StorageKey("sio2 mass fluid")),
      sio2_mass_total_key             (model.Database().StorageKey("sio2 mass total")),

      vein_quartz_solid_key          (model.Database().StorageKey("vein quartz solid"))


{} // end QuartzModel

template<uint32_t dim>
QuartzModel<dim>::~QuartzModel()
{} // end ~QuartzModel

template<uint32_t dim>
void QuartzModel<dim>::Visit(Model<dim>* node)
{}

template<uint32_t dim>
void QuartzModel<dim>::Visit(Region<dim>* node)
{}

template<uint32_t dim>
void QuartzModel<dim>::Visit(Node<dim>* node)
{    

    ReadVariables( node );
    GetFluidMasses();
    GetSpeciesContentsFluid();

    //initialization and calculation of equilibrium constants and vapor-liquid partition coefficients
    ReadEquilibriumConstantsFromLiterature();

    //calculation and conversion of aqueous species quantities
    InitializeConcentrationsAndGasMoleNumbers();

    //calculation of masses and mole numbers of minerals
    InitializeMineralMoleNumbers();

    CalculateTotalMoleNumbersOfBasisSpecies();

    CheckNonnegativityOfBasisSpecies();

    //quartz
    if (temperature()>lowercutofftemperature_quartzcalculations
            && temperature()<uppercutofftemperature_quartzcalculations) {

        SolveForQuartzEquilibriumState();

        //calculate mineral masses from their mole numbers
        CalculateAndStoreQuartzMassesFromMoleNumbers( node );
    }

    CalculateAndStoreConcentrations( node );
}
// end Visit

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template<uint32_t dim>
void QuartzModel<dim>::SolveForQuartzEquilibriumState()
{

    double saturation_sio2_aq=Kquartzaq_lit;
    double saturation_sio2_gas=Kquartzgas_lit;

    vein_quartz_solid_mol += std::max(sio2_tot_mol
                                      -massh2oliq*saturation_sio2_aq
                                      -massh2ogas*saturation_sio2_gas,
                                      0.);

    if (massh2oliq!=0.) {
        sio2_concentration = saturation_sio2_aq;}
    else {
        sio2_concentration = 0.;}

    sio2_gas_mol=massh2ogas*saturation_sio2_gas;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::ReadEquilibriumConstantsFromLiterature()
{
    double T(temperature());

    // SiO2(solid, Quartz) = SiO2(aq)
    // from C. Manning, GCA, 1994; in mol/kg H2O

    if (!essentiallyEqual(rho_l(),0.,numeric_limits<double>::epsilon()))
    {
        logKquartzaq_lit = 4.2620-5764.2/T + 1.7513e6/T/T - 2.2869e8/T/T/T
                +(2.8454-1006.9/T+3.5689e5/T/T)*log10(rho_l()/1000.);
        Kquartzaq_lit = pow(10.,logKquartzaq_lit);
    }
    else
    {
        Kquartzaq_lit = 0.;
    }

    if (!essentiallyEqual(rho_v(),0.,numeric_limits<double>::epsilon()))
    {
        logKquartzgas_lit = 4.2620-5764.2/T + 1.7513e6/T/T - 2.2869e8/T/T/T
                +(2.8454-1006.9/T+3.5689e5/T/T)*log10(rho_v()/1000.);
        Kquartzgas_lit = pow(10.,logKquartzgas_lit);
    }
    else
    {
        Kquartzgas_lit = 0.;
    }

    if (Kquartzaq_lit==0.)
    {
        D_sio2=0.;
    }
    else
    {
        D_sio2 = Kquartzgas_lit/Kquartzaq_lit;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::ApplyBoundaryFlowCorrections(Model<dim>* model)
{    
    const Region<dim>&   ref = model_ref_.Region("Model");
    for( typename std::vector<Node<dim>*>::const_iterator
         nit = ref.NodesBegin();
         nit != ref.NodesEnd();
         ++nit )
    {
        ApplyBoundaryFlowCorrections_aux(sio2_content_vapor_key,sio2_content_liquid_key,sio2_content_fluid_key,
                                         sio2_content_vapor,sio2_content_liquid,sio2_content_fluid,
                                         nit);
    }
}


template<uint32_t dim>
void QuartzModel<dim>::ApplyBoundaryFlowCorrections_aux

(csmp::Index vaporcontentkey,csmp::Index liquidcontentkey,csmp::Index fluidcontentkey,
 ScalarVariable vaporcontent,ScalarVariable liquidcontent,ScalarVariable fluidcontent,
 typename std::vector<Node<dim>*>::const_iterator nit)

{
    (*nit)->Read(bfm_key,bfm);

    (*nit)->Read(ml_key, ml);
    (*nit)->Read(mv_key, mv);

    (*nit)->Read(liquidcontentkey, liquidcontent);
    (*nit)->Read(vaporcontentkey, vaporcontent);

    double scalingfactor(0.);
    double tr_frac_fluid_inflow(0.); //must be adapted if recharge water has tracer content


    if (bfm()>0. && !essentiallyEqual(mml()+mmv(),0.,numeric_limits<double>::epsilon()))
    { //outflow

        scalingfactor = 1./(1.+((mml()/(mmv()+mml()))*bfm()/(mv()+ml())));

        //NK may 2024, : the relative mass proportions of vapor and liquid contributing to the total boundary flow should be approximately proportional
        //to the relative proportions of mobilities. Since the previous mobilities are not accessible, we take the current mobilities
        //as an approximation. Therefore the contribution of the liquid to the total boundary flow is (mmv/(mmv+mml))*bfm,
        //analogous for vapor

        liquidcontent() *= scalingfactor;

        scalingfactor = 1./(1.+((mmv()/(mmv()+mml()))*bfm()/(mv()+ml())));
        vaporcontent()  *= scalingfactor;

    } else if (bfm()<=0.) { //inflow

        liquidcontent() -= bfm()*tr_frac_fluid_inflow;
        vaporcontent()  -= bfm()*tr_frac_fluid_inflow;
    }

    fluidcontent = vaporcontent + liquidcontent;

    (*nit)->Store(liquidcontentkey, liquidcontent);
    (*nit)->Store(vaporcontentkey, vaporcontent);
    (*nit)->Store(fluidcontentkey, fluidcontent);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::CalculateSpeciesFractionsAndMobilities(Model<dim>* model)
{
    const Region<dim>&   ref = model_ref_.Region("Model");
    for( typename std::vector<Node<dim>*>::const_iterator
         nit = ref.NodesBegin();
         nit != ref.NodesEnd();
         ++nit )
    {
        CalculateSpeciesFractions_aux(sio2_content_vapor_key,sio2_content_liquid_key,
                                      sio2_fraction_vapor_key,sio2_fraction_liquid_key,
                                      sio2_content_vapor,sio2_content_liquid,
                                      sio2_fraction_vapor,sio2_fraction_liquid,
                                      nit);

        CalculateMobilities_aux(sio2_mobility_vapor_key,sio2_mobility_liquid_key,
                                sio2_fraction_vapor_key,sio2_fraction_liquid_key,
                                sio2_mobility_vapor,sio2_mobility_liquid,
                                sio2_fraction_vapor,sio2_fraction_liquid,
                                nit);
    }
}

template<uint32_t dim>
void QuartzModel<dim>::CalculateSpeciesFractions_aux

(csmp::Index vaporcontentkey,csmp::Index liquidcontentkey,
 csmp::Index vaporfractionkey,csmp::Index liquidfractionkey,
 ScalarVariable vaporcontent,ScalarVariable liquidcontent,
 ScalarVariable vaporfraction,ScalarVariable liquidfraction,
 typename std::vector<Node<dim>*>::const_iterator nit)

{
    (*nit)->Read(ml_key, ml);
    (*nit)->Read(mv_key, mv);

    (*nit)->Read(liquidcontentkey, liquidcontent);
    (*nit)->Read(vaporcontentkey, vaporcontent);

    if (ml()!=0.) {liquidfraction = liquidcontent / ml;}
    else {liquidfraction = 0.;}

    if (mv()!=0.) {vaporfraction = vaporcontent / mv;}
    else {vaporfraction = 0.;}

    (*nit)->Store(liquidfractionkey, liquidfraction);
    (*nit)->Store(vaporfractionkey, vaporfraction);

}

template<uint32_t dim>
void QuartzModel<dim>::CalculateMobilities_aux

(csmp::Index vapormobilitykey,csmp::Index liquidmobilitykey,
 csmp::Index vaporfractionkey,csmp::Index liquidfractionkey,
 ScalarVariable vapormobility,ScalarVariable liquidmobility,
 ScalarVariable vaporfraction,ScalarVariable liquidfraction,
 typename std::vector<Node<dim>*>::const_iterator nit)

{
    (*nit)->Read(mml_key,mml);
    (*nit)->Read(mmv_key,mmv);

    (*nit)->Read(liquidfractionkey, liquidfraction);
    (*nit)->Read(vaporfractionkey, vaporfraction);

    liquidmobility() = mml()*liquidfraction();
    vapormobility() = mmv()*vaporfraction();

    (*nit)->Store(liquidmobilitykey, liquidmobility);
    (*nit)->Store(vapormobilitykey, vapormobility);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::InitializeQuartz(Model<dim>* model)
{
    const Region<dim>&   ref = model_ref_.Region("Model");
    for( typename std::vector<Node<dim>*>::const_iterator
         nit = ref.NodesBegin();
         nit != ref.NodesEnd();
         ++nit )
    {
        (*nit)->Read(temperature_key, temperature);
        (*nit)->Read(rho_l_key, rho_l);
        (*nit)->Read(rho_v_key, rho_v);
        (*nit)->Read(sl_key, sl);
        (*nit)->Read(sv_key, sv);
        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);
        (*nit)->Read(pore_volume_key, pore_volume);

        (*nit)->Read(sio2_content_liquid_key, sio2_content_liquid);
        (*nit)->Read(sio2_content_vapor_key, sio2_content_vapor);

        ReadEquilibriumConstantsFromLiterature();

        massh2ogas = mv()*pore_volume();
        massh2oliq = ml()*pore_volume();

        double sio2_saturation_aq = Kquartzaq_lit;
        double sio2_saturation_gas = Kquartzgas_lit;

        sio2_content_liquid()=massh2oliq/pore_volume()*sio2_saturation_aq*sio2_molar_mass;
        sio2_content_vapor()=massh2ogas/pore_volume()*sio2_saturation_gas*sio2_molar_mass;
        vein_quartz_solid()=0.;

        (*nit)->Store(sio2_content_liquid_key, sio2_content_liquid);
        (*nit)->Store(sio2_content_vapor_key, sio2_content_vapor);
        (*nit)->Store(vein_quartz_solid_key, vein_quartz_solid);

    }
}

template<uint32_t dim>
void QuartzModel<dim>::Partitioning(Model<dim>* model)
{
    const Region<dim>&   ref = model_ref_.Region("Model");
    for( typename std::vector<Node<dim>*>::const_iterator
         nit = ref.NodesBegin();
         nit != ref.NodesEnd();
         ++nit )
    {
        (*nit)->Read(sio2_content_liquid_key, sio2_content_liquid);
        (*nit)->Read(sio2_content_vapor_key, sio2_content_vapor);
        (*nit)->Read(sio2_content_fluid_key, sio2_content_fluid);

        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(rho_l_key, rho_l);
        (*nit)->Read(rho_v_key, rho_v);

        (*nit)->Read(meltmass_key, meltmass );
        (*nit)->Read(crystalmass_key, crystalmass );
        (*nit)->Read(pore_volume_key, pore_volume );

        double sio2_mass_tot(0.);

        double massfractionwaterliquid = ml()/(ml()+mv());

        sio2_mass_liquid() = sio2_content_liquid()*pore_volume();
        sio2_mass_vapor() = sio2_content_vapor()*pore_volume();
        sio2_mass_fluid() = sio2_mass_liquid()+sio2_mass_vapor();

        //partitioning between the aqueous phases
        if (mv()==0.)
        {
            sio2_mass_liquid() = sio2_mass_fluid();
            sio2_mass_vapor() = 0.;
        }
        else if (ml()==0.)
        {
            sio2_mass_liquid() = 0.;
            sio2_mass_vapor() = sio2_mass_fluid();
        }
        else
        {
            sio2_mass_liquid() = sio2_mass_fluid()/(1+D_sio2*mv()/ml());
            sio2_mass_vapor() = sio2_mass_fluid() - sio2_mass_liquid();
        }

        sio2_content_liquid() = sio2_mass_liquid()/pore_volume();
        sio2_content_vapor() = sio2_mass_vapor()/pore_volume();
        sio2_content_fluid() = sio2_content_liquid()+sio2_content_vapor();

        //calculate concentrations
        //all concentrations in ppm (mg Li/kg fluid or melt)
        /////
        if (ml()!=0.)
        {
            sio2_concentration_liquid()=1.e6*sio2_content_liquid()/ml();
        }
        else
        {
            sio2_concentration_liquid()=0.;
        }

        /////
        if (mv()!=0.)
        {
            sio2_concentration_vapor()=1.e6*sio2_content_vapor()/mv();
        }
        else
        {
            sio2_concentration_vapor()=0.;
        }


        /////
        sio2_concentration_fluid()=1.e6*(sio2_content_liquid()+sio2_content_vapor())/(mv()+ml());

        (*nit)->Store(sio2_concentration_liquid_key, sio2_concentration_liquid);
        (*nit)->Store(sio2_concentration_vapor_key, sio2_concentration_vapor);
        (*nit)->Store(sio2_concentration_fluid_key, sio2_concentration_fluid);

        (*nit)->Store(sio2_content_liquid_key, sio2_content_liquid);
        (*nit)->Store(sio2_content_vapor_key, sio2_content_vapor);
        (*nit)->Store(sio2_content_fluid_key, sio2_content_fluid);

        (*nit)->Store(sio2_mass_liquid_key, sio2_mass_liquid);
        (*nit)->Store(sio2_mass_vapor_key, sio2_mass_vapor);
        (*nit)->Store(sio2_mass_fluid_key, sio2_mass_fluid);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::CalculateAndStoreConcentrations( Node<dim>* node )
{

    sio2_concentration_liquid=1.e6*sio2_concentration*sio2_molar_mass;
    if (massh2ogas!=0.)
    {
        sio2_concentration_vapor=1.e6*sio2_gas_mol*sio2_molar_mass/(massh2ogas);
    }
    else
    {
        sio2_concentration_vapor=0.;
    }
    sio2_concentration_fluid=1.e6*sio2_molar_mass * (sio2_concentration*massh2oliq+sio2_gas_mol)/(massh2ogas+massh2oliq);

    node->Store(sio2_concentration_vapor_key,sio2_concentration_vapor);
    node->Store(sio2_concentration_liquid_key,sio2_concentration_liquid);
    node->Store(sio2_concentration_fluid_key,sio2_concentration_fluid);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::ReadVariables(Node<dim>* node)
{
    //Read general variables
    node->Read(temperature_key, temperature);                           //°C
    node->Read(pressure_key, pressure);                                 //Pa
    node->Read(pore_volume_key, pore_volume);                           //m³
    node->Read(rho_l_key, rho_l);                                       //kg/m³
    node->Read(rho_v_key, rho_v);                                       //kg/m³
    node->Read(ml_key, ml);                                             //kg/m³
    node->Read(mv_key, mv);                                             //kg/m³
    node->Read(sl_key, sl);                                             //-/-
    node->Read(sv_key, sv);                                             //-/-
    node->Read(mml_key, mml);                                           //m²/s
    node->Read(mmv_key, mmv);                                           //m²/s
    node->Read(bfm_key, bfm);                                           //kg/m³

    node->Read(sio2_content_liquid_key, sio2_content_liquid);           //kg/m³
    node->Read(sio2_content_vapor_key, sio2_content_vapor);             //kg/m³

    //Read other variables
    node->Read(vein_quartz_solid_key, vein_quartz_solid);               //kg/m³
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::GetSpeciesContentsFluid()
{
    sio2_content_fluid() =          sio2_content_liquid()+sio2_content_vapor();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::GetFluidMasses()
{
    massh2oliq =                    ml()*pore_volume(); //kg
    massh2ogas =                    mv()*pore_volume();  //kg
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::InitializeConcentrationsAndGasMoleNumbers()
{
    if (definitelyGreaterThan(massh2oliq,0.,numeric_limits<double>::epsilon())) //massh2oliq nonzero
    {
        sio2_concentration = sio2_content_liquid () * pore_volume ()/(sio2_molar_mass*massh2oliq );             //mol/kg
        //mol/kg
    }
    else //if massh2oliq=0
    {
        sio2_concentration = 0.;
    }

    //calculation and conversion of gaseous species quantities
    sio2_gas_mol = sio2_content_vapor () * pore_volume ()/sio2_molar_mass;
    //mol
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::InitializeMineralMoleNumbers()
{
    vein_quartz_solid_mass   = vein_quartz_solid();
    vein_quartz_solid_mol    = vein_quartz_solid_mass/sio2_molar_mass;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::CalculateTotalMoleNumbersOfBasisSpecies()
{
    //These are the total mole numbers of the basis species used in the SolveForEquilibriumState function.
    sio2_tot_mol = sio2_content_fluid()*pore_volume()/sio2_molar_mass;
    //mol
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::CalculateAndStoreQuartzMassesFromMoleNumbers( Node<dim>* node )
{    
    vein_quartz_solid_mass = vein_quartz_solid_mol * sio2_molar_mass;
    vein_quartz_solid() = vein_quartz_solid_mass;

    node->Store(vein_quartz_solid_key, vein_quartz_solid);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<uint32_t dim>
void QuartzModel<dim>::CheckNonnegativityOfBasisSpecies()
{
    is_sio2_tot_mol_zero=false;

    //check for nonnegativity/non-zero values of species concentrations
    if (!definitelyGreaterThan(sio2_tot_mol,0.,numeric_limits<double>::epsilon()))
        is_sio2_tot_mol_zero=true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class QuartzModel<1U>;
template class QuartzModel<2U>;
template class QuartzModel<3U>;

} // end namespace csmp
