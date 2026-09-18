// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef Lithium_MODEL_H
#define Lithium_MODEL_H

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
class LithiumModel : public Visitor<dim> {

public:
    LithiumModel(Model<dim> &model, bool with_kinetics = true);

    virtual ~LithiumModel();

    void InitializeFromConcentration(Model<dim> *node, double concentration);
    void SetTimeStep(double timestep_ext);

    virtual void Visit(Node<dim> *node);
    virtual void Visit(Model<dim> *node);
    virtual void Visit(Region<dim> *node);

private:
    // (Benoit 18/08/2026) Removed unused members:
    //   keys   : mt_key, p_mt_key, sl_key, nodal_depth_key, bfs_key, rvl_key, rl_key
    //   scalars: temperature, pressure, mt, p_mt, mv, sl, bfs, rvl, rl
    //   other  : group_name, total_lithium_initial, H2OLookup water
    // mt/p_mt fell out when the mobility denominator moved to ml(). bfs
    // ("boundary flow salt") was read every node and never used — see the
    // basis fix in the boundary block of the .cpp: once carrier and
    // concentration are on the same basis the salt cancels, so bfs is not
    // needed. limv_key is KEPT (the vapor stores use it); there is no limv
    // ScalarVariable, those stores pass makeScalar(ANY, 0.).
    // `water` was a full H2OLookup instance per visitor, never referenced.
    const PropertyDatabase<dim>  &prop_ref_;
    Model< dim>                  &model_ref_;

    csmp::Index
        xCl_key,                           // salt content liquid
        p_xCl_key,                           // salt content liquid
        pore_volume_key,                  // pore volume
        bulk_volume_key,                  // bulk volume
        ml_key,                           // fluid mass liquid
        p_ml_key,                           // prev fluid mass liquid

        mml_key, mmv_key,
        liml_key,limv_key,
        lithium_content_fluid_key,           // lithium mass kg/m3 in pore fluid
        lithium_content_liquid_key,          // lithium mass kg/m3 in pore liquid
        lithium_content_vapor_key,          // lithium mass kg/m3 in pore vapor
        lithium_concentration_liquid_key,          // lithium mass kg per mass liquid WATER
        p_lithium_content_liquid_key,          // lithium mass kg/m3 in pore liquid
        lithium_solubility_liquid_key,       // lithium solubility kg/kg in liquid
        lithium_solid_key,                   // lithium mass kg/m3 in solid form, in the bulk volume
        lithium_bulk_volumic_mass_key,       // lithium mass kg/m3 in all forms, in the bulk volume
        bfm_key;                          // boundary mass flow

    ScalarVariable
        xCl, p_xCl,                           // salt content liquid
        pore_volume,                  // pore volume
        bulk_volume,                  // bulk volume
        ml,
        p_ml,
        mml, mmv,
        liml,
        lithium_content_fluid,          // lithium mass kg/m3 in pore fluid
        lithium_content_liquid,          // lithium mass kg/m3 in pore liquid
        lithium_concentration_liquid,
        p_lithium_content_liquid,          // lithium mass kg/m3 in pore liquid
        lithium_solubility_liquid,       // lithium solubility kg/kg in liquid
        lithium_solid,                   // lithium mass kg/m3 in solid form, in the bulk volume
        lithium_bulk_volumic_mass,
        bfm;

    bool with_kinetics;

    double timestep;

};
} // end namespace csmp

#endif