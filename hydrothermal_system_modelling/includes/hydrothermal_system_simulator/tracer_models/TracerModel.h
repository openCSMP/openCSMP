// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef tracer_MODEL_H
#define tracer_MODEL_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"
#include "H2OLookup.h"


namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/**
        @author BLC
        @date 2023
    */

template<size_t dim>
class TracerModel : public Visitor<dim> {

public:
    TracerModel(Model<dim> &model);

    virtual ~TracerModel();

    virtual void Visit(Node<dim> *node);
    virtual void Visit(Model<dim> *node);
    virtual void Visit(Region<dim> *node);

private:
    // (Benoit 18/08/2026) Removed unused members:
    //   keys   : bulk_volume_key, mt_key, p_mt_key, sl_key, rvl_key, rl_key
    //   scalars: mt, p_mt, mv, sl, trmv, rvl, rl
    //   other  : group_name, timestep, SetTimeStep
    // mt/p_mt fell out when the mobility denominator moved to ml(); the rest
    // were never read. trmv_key is KEPT (the vapor stores use it) but the
    // trmv ScalarVariable is not — those stores pass makeScalar(ANY, 0.).
    // SetTimeStep is gone because this class has no kinetics; the two driver
    // calls (Geothermal_lithium_example.cpp) were removed in the same commit.
    const PropertyDatabase<dim>  &prop_ref_;
    Model< dim>                  &model_ref_;

    csmp::Index
        pore_volume_key,                  // pore volume
        ml_key,                           // fluid mass liquid
        p_ml_key,                           // prev fluid mass liquid

        mml_key, mmv_key,
        trml_key,trmv_key,
        tracer_content_fluid_key,           // tracer mass kg/m3 in pore fluid
        tracer_content_liquid_key,          // tracer mass kg/m3 in pore liquid
        tracer_content_vapor_key,          // tracer mass kg/m3 in pore vapor
        tracer_concentration_liquid_key,          // tracer mass kg per mass liquid
        p_tracer_content_liquid_key,          // tracer mass kg/m3 in pore liquid
        bfm_key;                          // boundary mass flow

    ScalarVariable
        pore_volume,                  // pore volume
        ml,
        p_ml,
        mml, mmv,
        trml,
        tracer_content_fluid,          // tracer mass kg/m3 in pore fluid
        tracer_content_liquid,          // tracer mass kg/m3 in pore liquid
        tracer_concentration_liquid,
        p_tracer_content_liquid,          // tracer mass kg/m3 in pore liquid
        bfm;

};
} // end namespace csmp

#endif