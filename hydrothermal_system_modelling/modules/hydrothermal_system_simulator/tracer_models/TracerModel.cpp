// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "TracerModel.h"
#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"

using namespace std;
namespace csmp {

template<size_t dim>
TracerModel<dim>::TracerModel(Model<dim> &model)

    : Visitor<dim>(MODEL, NODE),
    prop_ref_(model.Database()),
    model_ref_(model),

    pore_volume_key(model.Database().StorageKey("pore volume")),
    ml_key(model.Database().StorageKey("fluid mass liquid")),
    p_ml_key(model.Database().StorageKey("previous fluid mass liquid")),
    mml_key(model.Database().StorageKey("liquid mass mobility")),
    mmv_key(model.Database().StorageKey("vapor mass mobility")),
    trml_key(model.Database().StorageKey("liquid tracer mobility")),
    trmv_key(model.Database().StorageKey("vapor tracer mobility")),

    tracer_content_fluid_key(model.Database().StorageKey("tracer content fluid")),
    tracer_content_liquid_key(model.Database().StorageKey("tracer content liquid")),
    tracer_content_vapor_key(model.Database().StorageKey("tracer content vapor")),
    tracer_concentration_liquid_key(model.Database().StorageKey("tracer concentration")),
    p_tracer_content_liquid_key(model.Database().StorageKey("previous tracer content liquid")),
    bfm_key(model.Database().StorageKey("boundary flow mass"))


{} // end TracerModel

template<size_t dim>
TracerModel<dim>::~TracerModel()
{} // end ~TracerModel

template<size_t dim>
void TracerModel<dim>::Visit(Model<dim> *node)
{}

template<size_t dim>
void TracerModel<dim>::Visit(Region<dim> *node)
{}

template<size_t dim>
void TracerModel<dim>::Visit(Node<dim> *node) {

    //Careful! here liquid = water + NaCl dissolved
    //Careful! here vapor = water + NaCl dissolved
    //Careful! liquid and vapor_water refers to H2O only

    //Read
    node->Read(pore_volume_key, pore_volume);
    node->Read(ml_key, ml);
    node->Read(p_ml_key, p_ml);
    node->Read(mml_key, mml);
    node->Read(mmv_key, mmv);
    node->Read(tracer_content_liquid_key, tracer_content_liquid);
    node->Read(p_tracer_content_liquid_key, p_tracer_content_liquid);
    node->Read(tracer_content_fluid_key, tracer_content_fluid);
    node->Read(bfm_key, bfm);

    if (tracer_content_liquid() < 0.) cerr << endl << "START Visit tracer content NEGATIVE!!!!";

    //Compute starting tracer and fluid phases mass values, from transport step
    double liquid_mass(0.),
        prev_liquid_mass(0.),
        tracer_dissolved_mass_liquid(0.),
        prev_tracer_dissolved_mass_liquid(0.);

    prev_liquid_mass              = p_ml() * pore_volume(); //kg liquid
    liquid_mass                   = ml() * pore_volume(); //kg liquid
    tracer_dissolved_mass_liquid  = tracer_content_liquid() * pore_volume();//kg tracer

    if (p_tracer_content_liquid() > 0.)
        prev_tracer_dissolved_mass_liquid  = p_tracer_content_liquid() * pore_volume();//kg tracer


    //BOUNDARY FLOW
    double tr_mass_out_l(0.);

    if (bfm() > 0.)
    {
        const double mm_tot = mml() + mmv();
        if (prev_liquid_mass > 0. && mm_tot > 0.)
            tr_mass_out_l = std::max(0.0, bfm() * mml() / mm_tot * (prev_tracer_dissolved_mass_liquid / prev_liquid_mass) * pore_volume());
    }

    if (tr_mass_out_l > prev_tracer_dissolved_mass_liquid) {
        cerr << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        cerr << endl << "tr_mass_out_l > prev_tracer_dissolved_mass_liquid: " << tr_mass_out_l << " vs prev: " << prev_tracer_dissolved_mass_liquid
             << " vs current: " << tracer_dissolved_mass_liquid;
        tr_mass_out_l = prev_tracer_dissolved_mass_liquid;
    }

    tracer_dissolved_mass_liquid -= tr_mass_out_l;

    //END BOUNDARY FLOW


    // DRY-OUT (Benoit 18/08/2026)
    // Explicit statement of what happens when the cell holds no liquid.
    // LithiumModel precipitates its dissolved load to "lithium solid" here;
    // the tracer has no solid store, so it stays as "tracer content liquid"
    // in a cell with ml == 0. Transport then computes its donor ratio as
    // TrCl/ml = 0, so the tracer is IMMOBILE until the cell rewets, at which
    // point it re-mobilises. Mass is conserved throughout — nothing is lost
    // or created — but the tracer does not follow the vapour out of the cell.
    // That is correct for a liquid-only conservative tracer; it is written
    // out here so it reads as a decision rather than an omission.
    // NO BEHAVIOUR CHANGE: this branch only zeroes the concentration, which
    // the guard below already did.
    if (liquid_mass <= 0.) tracer_concentration_liquid() = 0.;

    //Update
    tracer_content_liquid() = tracer_dissolved_mass_liquid / pore_volume();
    tracer_content_fluid()  = tracer_dissolved_mass_liquid / pore_volume();

    if (liquid_mass > 0.) tracer_concentration_liquid() = tracer_dissolved_mass_liquid / liquid_mass; //kg tracer per kg liquid


    // (Benoit 18/08/2026) Clamping a negative content CREATES mass. Report how
    // much, so a slow leak is distinguishable from a clean run instead of
    // looking identical to one.
    if (tracer_content_fluid() < 0.) {
        cerr << endl << "END Visit TrCf NEGATIVE, clamping to 0, mass created: "
             << -tracer_content_fluid() * pore_volume() << " kg";
        tracer_content_fluid() = 0.;
    }
    if (tracer_content_liquid() < 0.) {
        cerr << endl << "END Visit TrCl NEGATIVE, clamping to 0, mass created: "
             << -tracer_content_liquid() * pore_volume() << " kg";
        tracer_content_liquid() = 0.;
    }
    if (tracer_concentration_liquid() < 0.) {
        // warn only — deliberately NOT clamped, matching previous behaviour
        cerr << endl << "END Visit tracerConc NEGATIVE: " << tracer_concentration_liquid();
    }


    // (Benoit 18/08/2026) Denominator mt() -> ml();
    // The mobility convention in NaClH2OPropertiesVisitorPHX is
    // mobility = (krl/mul) * [content per unit LIQUID VOLUME], e.g.
    // xml = krl*xVl/mul with xVl = xCl/sl. Converting a per-pore-volume
    // content C gives mobility = mml * C / ml, NOT / mt. Dividing by the
    // TOTAL fluid mass mt = ml + mv (+ma) scaled the tracer down by the
    // liquid mass fraction — vapour carrying no tracer still diluted it, a
    // spurious retardation of exactly mt/ml. Invisible in single-phase
    // liquid runs (mt == ml), wrong as soon as vapour appears.
    //
    // NOTE: this field is DIAGNOSTIC ONLY. Transport rebuilds the tracer flux
    // from the mass flux at the donor's tracer/ml ratio and never reads this
    // value.
    // See ExplicitFiniteVolumeTransportPHX::DetermineFacetFluxRatesOnce.
    trml() = 0.;
    if (ml() > 0.) {
        trml() =  mml() * tracer_content_liquid() / ml(); //kg tracer / kg liquid
    }

    // Store results
    node->Store(tracer_content_liquid_key, tracer_content_liquid);
    node->Store(tracer_content_vapor_key, makeScalar(ANY, 0.));
    node->Store(tracer_concentration_liquid_key, tracer_concentration_liquid);
    node->Store(tracer_content_fluid_key, tracer_content_fluid);


    node->Store(trml_key, trml);
    node->Store(trmv_key, makeScalar(ANY, 0.));
    //////////////////////////////////////

} // end Visit

template class TracerModel<1U>;
template class TracerModel<2U>;
template class TracerModel<3U>;

} // end namespace csmp