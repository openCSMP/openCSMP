// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "LithiumModel.h"
#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "H2OLookup.h"

using namespace std;
namespace csmp {

template<size_t dim>
LithiumModel<dim>::LithiumModel(Model<dim> &model, bool with_kinetics)

    : Visitor<dim>(MODEL, NODE),

    prop_ref_(model.Database()),
    model_ref_(model),

    xCl_key(model.Database().StorageKey("salt content liquid")),
    p_xCl_key(model.Database().StorageKey("previous salt content liquid")),
    pore_volume_key(model.Database().StorageKey("pore volume")),
    bulk_volume_key(model.Database().StorageKey("bulk volume")),
    ml_key(model.Database().StorageKey("fluid mass liquid")),
    p_ml_key(model.Database().StorageKey("previous fluid mass liquid")),
    mml_key(model.Database().StorageKey("liquid mass mobility")),
    mmv_key(model.Database().StorageKey("vapor mass mobility")),
    liml_key(model.Database().StorageKey("liquid lithium mobility")),
    limv_key(model.Database().StorageKey("vapor lithium mobility")),
    lithium_content_fluid_key(model.Database().StorageKey("lithium content fluid")),
    lithium_content_liquid_key(model.Database().StorageKey("lithium content liquid")),
    lithium_content_vapor_key(model.Database().StorageKey("lithium content vapor")),
    lithium_concentration_liquid_key(model.Database().StorageKey("lithium concentration liquid")),
    p_lithium_content_liquid_key(model.Database().StorageKey("previous lithium content liquid")),
    lithium_solubility_liquid_key(model.Database().StorageKey("lithium solubility liquid")),
    lithium_solid_key(model.Database().StorageKey("lithium solid")),
    lithium_bulk_volumic_mass_key(model.Database().StorageKey("lithium bulk volumic mass")),
    bfm_key(model.Database().StorageKey("boundary flow mass")),

    with_kinetics(with_kinetics),

    timestep(0.)

{} // end LithiumModel

template<size_t dim>
LithiumModel<dim>::~LithiumModel()
{} // end ~LithiumModel

template<size_t dim>
void LithiumModel<dim>::Visit(Model<dim> *node)
{}

template<size_t dim>
void LithiumModel<dim>::Visit(Region<dim> *node)
{}

template<size_t dim>
void LithiumModel<dim>::Visit(Node<dim> *node) {
    //TODO: check if lithium variables scaling with volume is done properly
    //Careful! here liquid = water + NaCl dissolved
    //Careful! here vapor = water + NaCl dissolved
    //Careful! liquid_water and vapor_water refers to H2O only

    //There is no lithium input at boundaries

    //Read
    node->Read(xCl_key, xCl);
    node->Read(p_xCl_key, p_xCl);
    node->Read(pore_volume_key, pore_volume);
    node->Read(bulk_volume_key, bulk_volume);
    node->Read(ml_key, ml);
    node->Read(p_ml_key, p_ml);
    node->Read(mml_key, mml);
    node->Read(mmv_key, mmv);
    node->Read(lithium_content_liquid_key, lithium_content_liquid);
    node->Read(p_lithium_content_liquid_key, p_lithium_content_liquid);
    node->Read(lithium_content_fluid_key, lithium_content_fluid);
    node->Read(lithium_solid_key, lithium_solid);
    node->Read(bfm_key, bfm);
    node->Read(lithium_solubility_liquid_key, lithium_solubility_liquid);

    if (lithium_content_liquid() < 0.) cerr << endl << "START Visit LiCl NEGATIVE!!!!";

    //double li_molar_mass = 6.941e-3; //kg/mol

    double dissolution_rate(1.e-20); //kg Li per second per m3 (bulk) do -20, -10, -11, -12

    //Compute starting lithium and fluid phases mass values, from transport step
    double liquid_water_mass(0.),
        prev_liquid_water_mass(0.),
        lithium_solid_mass(0.),
        lithium_dissolved_mass_liquid(0.),
        prev_lithium_dissolved_mass_liquid(0.);

    prev_liquid_water_mass              = (p_ml() - p_xCl()) * pore_volume(); //kg WATER
    liquid_water_mass                   = (ml() - xCl()) * pore_volume(); //kg WATER
    lithium_solid_mass                  = lithium_solid() * bulk_volume();
    lithium_dissolved_mass_liquid       = lithium_content_liquid() * pore_volume();//kg li, Careful! here liquid = water + NaCl dissolved, transported variable


    if (p_lithium_content_liquid() > 0.)
        prev_lithium_dissolved_mass_liquid  = p_lithium_content_liquid() * pore_volume();//kg li, Careful! here liquid = water + NaCl dissolved, transported variable


    //BOUNDARY FLOW
    double li_mass_out_l(0.);

    if (bfm() > 0.) {
        // (Benoit 18/08/2026) BASIS FIX + zero-guard.
        // BASIS: the carrier bfm*mml/(mml+mmv) is LIQUID mass (water + dissolved
        // NaCl), but this was multiplied by a concentration expressed per WATER
        // mass (prev_lithium_dissolved / prev_liquid_water_mass). Mismatched
        // bases over-removed lithium at outflow boundaries by 1/(1-X), with X
        // the liquid salt mass fraction — a few % for dilute brine, >40% near
        // halite saturation. The well path already pairs these correctly
        // (Coupled_reservoir_well_..._Scheme.cpp:2182: it subtracts
        // salt_mass_transfer_rate to get a WATER carrier before applying the
        // per-water concentration).
        // Fixed by expressing the concentration per LIQUID mass instead — the
        // salt cancels exactly, since
        //   bfm*f*PV * (prev_water/prev_liquid) * (prev_Li/prev_water)
        //     = bfm*f*PV * prev_Li/prev_liquid_mass
        // which also makes this line structurally identical to TracerModel's.
        // ("boundary flow salt" bfs is NOT the right vehicle here: it lumps in
        //  vapor salt and halite, so (bfm-bfs)*f would be less accurate. bfs is
        //  therefore genuinely unused by this model and has been removed.)
        //
        // ZERO-GUARD: mml+mmv is zero where neither phase is mobile (kr = 0 in
        // both, or a fully dried node); with bfm > 0 that produced a silent NaN.
        const double mm_tot          = mml() + mmv();
        const double prev_liquid_mass = p_ml() * pore_volume(); //kg liquid (water + dissolved NaCl)
        if (prev_liquid_mass > 0. && mm_tot > 0.)
            li_mass_out_l = std::max(0.0, bfm() * mml() / mm_tot * (prev_lithium_dissolved_mass_liquid / prev_liquid_mass) * pore_volume());
    }

    //Careful! here bfm = water + NaCl dissolved. We use previous as boundary flow was calculated with previous values of density, etc...

    if (li_mass_out_l > prev_lithium_dissolved_mass_liquid) {
        cerr << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        cerr << endl << "li_mass_out_l > prev_lithium_dissolved_mass_liquid: " << li_mass_out_l << " vs prev: " << prev_lithium_dissolved_mass_liquid
             << " vs current: " << lithium_dissolved_mass_liquid;
        li_mass_out_l = prev_lithium_dissolved_mass_liquid;
    }

    lithium_dissolved_mass_liquid -= li_mass_out_l;

    //END BOUNDARY FLOW

    // Potential phase change occured during thermal equilibration!
    if (liquid_water_mass <= 0. && lithium_dissolved_mass_liquid > 0.)
    //Precipitates all lithium that was left in liquid in previous timestep if there is no more liquid
    {
        //cerr<<endl<<"ml: "<<ml()<<", xCl: "<<xCl();
        //cerr<<endl<<"No more liquid after thermal equilibration, precipitating the lithium that was there previously";
        lithium_solid_mass += lithium_dissolved_mass_liquid;
        lithium_dissolved_mass_liquid = 0.;
    }

    // //////////////////////////////////////
    // Calculating concentration (in mass of lithium per mass of liquid phase, the latter excludes dissolved NaCl
    lithium_concentration_liquid() = 0.;

    if (liquid_water_mass > 0.) lithium_concentration_liquid() = lithium_dissolved_mass_liquid / liquid_water_mass; //kg li per kg WATER

    //END Calculate concentration


    // Calculate lithium solubility
    // NO CALCULATION, FIXED AT INITIALZATION, EQUAL TO INITIAL CONCENTRATION
    // END Calculate lithium solubility

    // Calculate dissolution/precipitation
    double lithium_mass_dissolving(0.);
    //double lithium_mass_precipitating (0.);

    //Dissolution
    if (liquid_water_mass > 0. && lithium_concentration_liquid() < lithium_solubility_liquid()) {
        //lithium_mass_dissolving = -(lithium_concentration_liquid() - lithium_solubility_liquid()) * liquid_water_mass;

        if (with_kinetics) {
            lithium_mass_dissolving = dissolution_rate * timestep * bulk_volume();

            if ((lithium_dissolved_mass_liquid + lithium_mass_dissolving) > lithium_solubility_liquid() * liquid_water_mass) {
                lithium_mass_dissolving = lithium_solubility_liquid() * liquid_water_mass - lithium_dissolved_mass_liquid; //solubility limit
            }
        }

        if (lithium_mass_dissolving >= lithium_solid_mass)
        //all lithium dissolves
        {
            lithium_mass_dissolving = lithium_solid_mass;

            lithium_dissolved_mass_liquid += lithium_mass_dissolving;
            lithium_solid_mass = 0.;
        }

        else {
            lithium_dissolved_mass_liquid += lithium_mass_dissolving;
            lithium_solid_mass -= lithium_mass_dissolving;
        }
    }

    // //Precipitation  — kept commented, but the ORDER BUG is fixed so this is
    // //safe to re-enable as written. (Benoit 18/08/2026)
    // //Previously the "all precipitates" branch read:
    // //    lithium_dissolved_mass_liquid = 0.;
    // //    lithium_solid_mass += lithium_dissolved_mass_liquid;   // adds ZERO
    // //which zeroed the dissolved load before transferring it, silently
    // //destroying all of it instead of moving it to the solid phase.
    // if( liquid_water_mass > 0. && lithium_concentration_liquid() > lithium_solubility_liquid() )
    // {
    //     lithium_mass_precipitating = (lithium_concentration_liquid() - lithium_solubility_liquid()) * liquid_water_mass;


    //     if( lithium_mass_precipitating >= lithium_dissolved_mass_liquid )
    //     //all lithium precipitates
    //     {
    //         lithium_mass_precipitating = lithium_dissolved_mass_liquid;

    //         lithium_solid_mass += lithium_mass_precipitating;   // transfer FIRST
    //         lithium_dissolved_mass_liquid = 0.;                 // then zero
    //     }

    //     else
    //     {
    //         lithium_dissolved_mass_liquid -= lithium_mass_precipitating;
    //         lithium_solid_mass += lithium_mass_precipitating;
    //     }
    // }

    //Update
    lithium_content_liquid() = lithium_dissolved_mass_liquid / pore_volume();
    lithium_content_fluid()  = lithium_dissolved_mass_liquid / pore_volume();

    if (liquid_water_mass > 0.) lithium_concentration_liquid() = lithium_dissolved_mass_liquid / liquid_water_mass; //kg li per kg WATER

    // if((lithium_concentration_liquid()-lithium_solubility_liquid()) > 1.e-6) cerr<<"CONCENTRATION ABOVE SOLUBILITY!!! Conc: "<<lithium_concentration_liquid()
    //          <<", Solub: "<<lithium_solubility_liquid()<<endl;


    if (lithium_content_liquid() < 0.) cerr << endl << "END Visit LiCl NEGATIVE!!!!";
    if (lithium_concentration_liquid() < 0.) cerr << endl << "END Visit Li Conc NEGATIVE!!!!";

    // (Benoit 18/08/2026) Clamping a negative content CREATES mass. Report how
    // much, so a slow leak is distinguishable from a clean run instead of
    // looking identical to one. Concentration is deliberately left unclamped,
    // matching previous behaviour.
    if (lithium_content_fluid() < 0.) {
        cerr << endl << "END Visit LiCf clamping to 0, mass created: "
             << -lithium_content_fluid() * pore_volume() << " kg";
        lithium_content_fluid() = 0.;
    }
    if (lithium_content_liquid() < 0.) {
        cerr << endl << "END Visit LiCl clamping to 0, mass created: "
             << -lithium_content_liquid() * pore_volume() << " kg";
        lithium_content_liquid() = 0.;
    }
    //NEED TO CHECK EXPRESSIONS ABOVE!!!!

    lithium_solid() = lithium_solid_mass / bulk_volume();
    lithium_bulk_volumic_mass() = (lithium_dissolved_mass_liquid + lithium_solid_mass) / bulk_volume();

    // (Benoit 18/08/2026) Denominator mt() -> ml(). Same fix as TracerModel.
    // mobility = mml * C / ml, not / mt: dividing by the TOTAL fluid mass
    // mt = ml + mv (+ma) let vapour that carries no lithium dilute the
    // lithium mobility, a spurious retardation of mt/ml. The guard already
    // tested ml() > 0. while the denominator used mt() — the mismatch was
    // the tell. Invisible while mt == ml (single-phase liquid).
    //
    // NOTE: DIAGNOSTIC ONLY — transport rebuilds the lithium flux from the
    // mass flux at the donor's lithium/ml ratio and never reads this value.
    // See ExplicitFiniteVolumeTransportPHX::DetermineFacetFluxRatesOnce.
    liml() = 0.;
    if (ml() > 0.) {
        liml() =  mml() * lithium_content_liquid() / ml(); //kg Li / kg liquid
    }

    //////////////////////////////////////

    // Store results
    node->Store(lithium_content_liquid_key, lithium_content_liquid);
    node->Store(lithium_content_vapor_key, makeScalar(ANY, 0.));
    node->Store(lithium_concentration_liquid_key, lithium_concentration_liquid);
    node->Store(lithium_content_fluid_key, lithium_content_fluid);


    node->Store(liml_key, liml);
    node->Store(limv_key, makeScalar(ANY, 0.));

    node->Store(lithium_solid_key, lithium_solid);
    node->Store(lithium_bulk_volumic_mass_key, lithium_bulk_volumic_mass);
    //////////////////////////////////////

} // end Visit

template<size_t dim>
void LithiumModel<dim>::InitializeFromConcentration(Model<dim> *node, double concentration) {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for (typename std::vector<Node<dim>*>::const_iterator it = ref.NodesBegin(); it != ref.NodesEnd(); ++it) {
        (*it)->Read(pore_volume_key, pore_volume);
        (*it)->Read(ml_key, ml);
        (*it)->Read(xCl_key, xCl);

        double liquid_water_mass = (ml() - xCl()) * pore_volume(); //kg WATER
        double lithium_dissolved_mass_liquid = liquid_water_mass * concentration;

        lithium_content_fluid() = lithium_content_liquid() = lithium_dissolved_mass_liquid / pore_volume();
        lithium_concentration_liquid() = concentration;
        lithium_solubility_liquid() = concentration;

        (*it)->Store(lithium_content_liquid_key, lithium_content_liquid);
        (*it)->Store(lithium_concentration_liquid_key, lithium_concentration_liquid);
        (*it)->Store(lithium_content_fluid_key, lithium_content_fluid);
        (*it)->Store(lithium_solubility_liquid_key, lithium_solubility_liquid);
    }
}

template<size_t dim>
void LithiumModel<dim>::SetTimeStep(double timestep_ext) {
    timestep = timestep_ext;
}

template class LithiumModel<1U>;
template class LithiumModel<2U>;
template class LithiumModel<3U>;

} // end namespace csmp