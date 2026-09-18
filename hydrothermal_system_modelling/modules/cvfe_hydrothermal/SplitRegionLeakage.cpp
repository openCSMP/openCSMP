// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SplitRegionLeakage.h"
#include "MathOperatorLHS.h"
#include "Model.h"
#include "SplitBoundary.h"
#include "compareFloats.h"
#include "Exception.h"
#include <set>

using namespace std;
namespace csmp {

// ── SB-leakage CFL tuning knobs (all in one place) ───────────────────────────
// Drainage & stuffing are checked per-phase; overdraining (exact mass, margin-
// guarded) is the always-on backstop. Adjust these for threshold sweeps without
// touching the logic below.
static constexpr double sb_drain_margin    = 0.0; // overdrain fires within this fraction of empty (leave >= margin of holdup)
static constexpr double sb_drain_threshold = 2000.;  // per-phase drainage CFL fraction that triggers a cut
static constexpr double sb_stuff_threshold = 2000.;  // per-phase stuffing CFL fraction that triggers a cut

// NOTE: cut-vs-clamp is NOT a file-scope knob. It is set by the caller via
// beat_mode_ (see header): ExchangeSubstep -> clamp; one-shot legacy pass -> cut.

template<uint32_t dim>
SplitRegionLeakage<dim>::SplitRegionLeakage(Model<dim> &model,
                                            bool open_space,
                                            double fault_perm_anisotropy)
    :
    // 1D leakage computation
    // Only works with split boundaries

    model_ref_(model),

    prop_ref_(model.Database()),

    open_space_(open_space),
    fault_perm_anisotropy_(fault_perm_anisotropy),

    p_key(model.Database().StorageKey("fluid pressure")),
    pp_key(model.Database().StorageKey("previous fluid pressure")),

    nfvs_key(model.Database().StorageKey("nodal fluid volume source")),
    qM_key(model.Database().StorageKey("split region mass source")),
    qE_key(model.Database().StorageKey("split region energy source")),
    qS_key(model.Database().StorageKey("split region salt source")),
    k_key(model.Database().StorageKey("permeability")),

    viscosity_l_key(model.Database().StorageKey("viscosity liquid")),
    viscosity_v_key(model.Database().StorageKey("viscosity vapor")),
    sat_l_key(model.Database().StorageKey("saturation liquid")),
    sat_v_key(model.Database().StorageKey("saturation vapor")),
    sat_halite_key(model.Database().StorageKey("saturation halite")),
    density_l_key(model.Database().StorageKey("density liquid")),
    density_v_key(model.Database().StorageKey("density vapor")),

    pore_volume_key(model.Database().StorageKey("pore volume")),

    //Advection variables:
    ml_key(model.Database().StorageKey("fluid mass liquid")),
    mv_key(model.Database().StorageKey("fluid mass vapor")),
    hCl_key(model.Database().StorageKey("enthalpy content liquid")),
    hCv_key(model.Database().StorageKey("enthalpy content vapor")),
    hv_key(model.Database().StorageKey("enthalpy vapor")),
    xCl_key(model.Database().StorageKey("salt content liquid")),
    xCv_key(model.Database().StorageKey("salt content vapor")),

    thickness_key(model.Database().StorageKey("thickness")),

    //Legacy mass source for visualization, not used in P equation anymore
    grav_mass_source_key(model.Database().StorageKey("split region gravity mass source")),

    split_grad_pot_l_key(model.Database().StorageKey("split region potential gradient l")),
    split_grad_pot_v_key(model.Database().StorageKey("split region potential gradient v")),

    split_relperm_l_key(model.Database().StorageKey("split region relperm l")),
    split_relperm_v_key(model.Database().StorageKey("split region relperm v")),

    with_SB_leakage_report(false),

    dt(0.)

{

} // end SplitRegionLeakage

template<uint32_t dim>
SplitRegionLeakage<dim>::~SplitRegionLeakage()
{} // end ~SplitRegionLeakage


template<uint32_t dim>
void SplitRegionLeakage<dim>::SetTimeStep(double dt_ext) {
    dt = dt_ext;
}

template<uint32_t dim>
void SplitRegionLeakage<dim>::ReadVariables()
{
    nfvs_in  = in_node->Read(nfvs_key);
    nfvs_out  = out_node->Read(nfvs_key);
    nfvs_mid  = mid_node->Read(nfvs_key);

    mul_in = in_node->Read(viscosity_l_key);
    mul_out = out_node->Read(viscosity_l_key);
    mul_mid = mid_node->Read(viscosity_l_key);

    muv_in = in_node->Read(viscosity_v_key);
    muv_out = out_node->Read(viscosity_v_key);
    muv_mid = mid_node->Read(viscosity_v_key);

    Sl_in = in_node->Read(sat_l_key);
    Sl_out = out_node->Read(sat_l_key);
    Sl_mid = mid_node->Read(sat_l_key);

    Sv_in = in_node->Read(sat_v_key);
    Sv_out = out_node->Read(sat_v_key);
    Sv_mid = mid_node->Read(sat_v_key);

    Sh_in = in_node->Read(sat_halite_key);
    Sh_out = out_node->Read(sat_halite_key);
    Sh_mid = mid_node->Read(sat_halite_key);

    rhol_in = in_node->Read(density_l_key);
    rhol_out = out_node->Read(density_l_key);
    rhol_mid = mid_node->Read(density_l_key);

    rhov_in = in_node->Read(density_v_key);
    rhov_out = out_node->Read(density_v_key);
    rhov_mid = mid_node->Read(density_v_key);

    pore_volume_in = in_node->Read(pore_volume_key);
    pore_volume_out = out_node->Read(pore_volume_key);
    pore_volume_mid = mid_node->Read(pore_volume_key);

    grav_mass_source_in_tot = in_node->Read(grav_mass_source_key);
    grav_mass_source_out_tot = out_node->Read(grav_mass_source_key);
    grav_mass_source_mid_tot = mid_node->Read(grav_mass_source_key);

    qM_in = in_node->Read(qM_key);
    qM_out = out_node->Read(qM_key);
    qM_mid = mid_node->Read(qM_key);

    qE_in = in_node->Read(qE_key);
    qE_out = out_node->Read(qE_key);
    qE_mid = mid_node->Read(qE_key);

    qS_in = in_node->Read(qS_key);
    qS_out = out_node->Read(qS_key);
    qS_mid = mid_node->Read(qS_key);

    if(with_gold_)
    {
        qSO2_in  = in_node ->Read(qSO2_key);
        qSO2_out = out_node->Read(qSO2_key);
        qSO2_mid = mid_node->Read(qSO2_key);

        qAuHS0_in  = in_node ->Read(qAuHS0_key);
        qAuHS0_out = out_node->Read(qAuHS0_key);
        qAuHS0_mid = mid_node->Read(qAuHS0_key);
    }

    // Transport variables
    pressure_in  = in_node->Read(p_key);
    pressure_out  = out_node->Read(p_key);
    pressure_mid  = mid_node->Read(p_key);

    ppressure_in  = in_node->Read(pp_key);
    ppressure_out  = out_node->Read(pp_key);
    ppressure_mid  = mid_node->Read(pp_key);

    ml_in  = in_node->Read(ml_key);
    ml_out = out_node->Read(ml_key);
    ml_mid = mid_node->Read(ml_key);

    mv_in  = in_node->Read(mv_key);
    mv_out = out_node->Read(mv_key);
    mv_mid = mid_node->Read(mv_key);

    hCl_in = in_node->Read(hCl_key);
    hCl_out = out_node->Read(hCl_key);
    hCl_mid = mid_node->Read(hCl_key);

    hCv_in = in_node->Read(hCv_key);
    hCv_out = out_node->Read(hCv_key);
    hCv_mid = mid_node->Read(hCv_key);

    xCl_in = in_node->Read(xCl_key);
    xCl_out = out_node->Read(xCl_key);
    xCl_mid = mid_node->Read(xCl_key);

    xCv_in = in_node->Read(xCv_key);
    xCv_out = out_node->Read(xCv_key);
    xCv_mid = mid_node->Read(xCv_key);

    if (with_tracer_)
    {
        TrCl_in = in_node->Read(TrCl_key);
        TrCl_out = out_node->Read(TrCl_key);
        TrCl_mid = mid_node->Read(TrCl_key);

        TrCv_in = in_node->Read(TrCv_key);
        TrCv_out = out_node->Read(TrCv_key);
        TrCv_mid = mid_node->Read(TrCv_key);
    }

    if (with_lithium_)
    {
        LiCl_in = in_node->Read(LiCl_key);
        LiCl_out = out_node->Read(LiCl_key);
        LiCl_mid = mid_node->Read(LiCl_key);

        LiCv_in = in_node->Read(LiCv_key);
        LiCv_out = out_node->Read(LiCv_key);
        LiCv_mid = mid_node->Read(LiCv_key);
    }

    if(with_air_)
    {
        mua_in  = in_node->Read(viscosity_a_key);
        mua_out = out_node->Read(viscosity_a_key);
        mua_mid = mid_node->Read(viscosity_a_key);

        Sa_in  = in_node->Read(sat_a_key);
        Sa_out = out_node->Read(sat_a_key);
        Sa_mid = mid_node->Read(sat_a_key);

        rhoa_in  = in_node->Read(density_a_key);
        rhoa_out = out_node->Read(density_a_key);
        rhoa_mid = mid_node->Read(density_a_key);

        ma_in  = in_node->Read(ma_key);
        ma_out = out_node->Read(ma_key);
        ma_mid = mid_node->Read(ma_key);

        hCa_in  = in_node->Read(hCa_key);
        hCa_out = out_node->Read(hCa_key);
        hCa_mid = mid_node->Read(hCa_key);
    }

    if (with_gold_)
    {
        so2Cl_in = in_node->Read(so2Cl_key);
        so2Cl_out = out_node->Read(so2Cl_key);
        so2Cl_mid = mid_node->Read(so2Cl_key);

        so2Cv_in = in_node->Read(so2Cv_key);
        so2Cv_out = out_node->Read(so2Cv_key);
        so2Cv_mid = mid_node->Read(so2Cv_key);

        auhs0Cl_in = in_node->Read(auhs0Cl_key);
        auhs0Cl_out = out_node->Read(auhs0Cl_key);
        auhs0Cl_mid = mid_node->Read(auhs0Cl_key);

        auhs0Cv_in = in_node->Read(auhs0Cv_key);
        auhs0Cv_out = out_node->Read(auhs0Cv_key);
        auhs0Cv_mid = mid_node->Read(auhs0Cv_key);
    }


}

template<uint32_t dim>
void SplitRegionLeakage<dim>::ComputeRelPerm() {
    //effective permeability calculatiom

    if(Sh_in()<1.)
    {
        effective_sat_in_l  = Sl_in()  / (1. - Sh_in());
        effective_sat_in_v  = Sv_in()  / (1. - Sh_in());
        if(with_air_) effective_sat_in_a  = Sa_in()  / (1. - Sh_in());
    }
    if(Sh_mid()<1.)
    {
        effective_sat_mid_l = Sl_mid() / (1. - Sh_mid());
        effective_sat_mid_v = Sv_mid() / (1. - Sh_mid());
        if(with_air_) effective_sat_mid_a = Sa_mid() / (1. - Sh_mid());
    }
    if(Sh_out()<1.)
    {
        effective_sat_out_l = Sl_out() / (1. - Sh_out());
        effective_sat_out_v = Sv_out() / (1. - Sh_out());
        if(with_air_) effective_sat_out_a = Sa_out() / (1. - Sh_out());
    }

    // Hard coded residual saturation l of 30%!!!!
    relative_perm_in_l()  = 1. / 0.7 * effective_sat_in_l - 0.3 / 0.7;
    relative_perm_mid_l() = 1. / 0.7 * effective_sat_mid_l - 0.3 / 0.7;
    relative_perm_out_l() = 1. / 0.7 * effective_sat_out_l - 0.3 / 0.7;

    if (effective_sat_in_l < 0.3)  relative_perm_in_l()  = 0.;
    if (effective_sat_mid_l < 0.3) relative_perm_mid_l() = 0.;
    if (effective_sat_out_l < 0.3) relative_perm_out_l() = 0.;

    if(!with_air_)
    {
        // Original vapor kr formula (no air)
        relative_perm_in_v()  = 1. / 0.7 * effective_sat_in_v;
        relative_perm_mid_v() = 1. / 0.7 * effective_sat_mid_v;
        relative_perm_out_v() = 1. / 0.7 * effective_sat_out_v;

        if (effective_sat_in_l < 0.3)  relative_perm_in_v()  = 1.;
        if (effective_sat_mid_l < 0.3) relative_perm_mid_v() = 1.;
        if (effective_sat_out_l < 0.3) relative_perm_out_v() = 1.;
    }
    else
    {

        // Gas kr = 1 - kr_l, partitioned between vapor and air by their effective saturations
        double kr_gas_in  = 1. - relative_perm_in_l();
        double kr_gas_mid = 1. - relative_perm_mid_l();
        double kr_gas_out = 1. - relative_perm_out_l();

        double s_gas_in  = effective_sat_in_v  + effective_sat_in_a;
        double s_gas_mid = effective_sat_mid_v + effective_sat_mid_a;
        double s_gas_out = effective_sat_out_v + effective_sat_out_a;

        if (s_gas_in > 0.)
        {
            relative_perm_in_v()  = kr_gas_in  * effective_sat_in_v  / s_gas_in;
            relative_perm_in_a()  = kr_gas_in  * effective_sat_in_a  / s_gas_in;
        }
        else { relative_perm_in_v()  = 0.; relative_perm_in_a()  = 0.; }

        if (s_gas_mid > 0.)
        {
            relative_perm_mid_v() = kr_gas_mid * effective_sat_mid_v / s_gas_mid;
            relative_perm_mid_a() = kr_gas_mid * effective_sat_mid_a / s_gas_mid;
        }
        else { relative_perm_mid_v() = 0.; relative_perm_mid_a() = 0.; }

        if (s_gas_out > 0.)
        {
            relative_perm_out_v() = kr_gas_out * effective_sat_out_v / s_gas_out;
            relative_perm_out_a() = kr_gas_out * effective_sat_out_a / s_gas_out;
        }
        else { relative_perm_out_v() = 0.; relative_perm_out_a() = 0.; }
    }

    // Case of saturation halite = 1, will probably break elswhere anyway if that's the case...
    if(essentiallyEqual(Sh_in(), 1., numeric_limits<double>::epsilon()))
    {
        relative_perm_in_l() = 0.;
        relative_perm_in_v() = 0.;
        if(with_air_) relative_perm_in_a() = 0.;
    }
    if(essentiallyEqual(Sh_out(), 1., numeric_limits<double>::epsilon()))
    {
        relative_perm_out_l() = 0.;
        relative_perm_out_v() = 0.;
        if(with_air_) relative_perm_out_a() = 0.;
    }
    if(essentiallyEqual(Sh_mid(), 1., numeric_limits<double>::epsilon()))
    {
        relative_perm_mid_l() = 0.;
        relative_perm_mid_v() = 0.;
        if(with_air_) relative_perm_mid_a() = 0.;
    }

    in_node->Store(split_relperm_l_key, relative_perm_in_l);
    mid_node->Store(split_relperm_l_key, relative_perm_mid_l);
    out_node->Store(split_relperm_l_key, relative_perm_out_l);
    in_node->Store(split_relperm_v_key, relative_perm_in_v);
    mid_node->Store(split_relperm_v_key, relative_perm_mid_v);
    out_node->Store(split_relperm_v_key, relative_perm_out_v);
    if(with_air_)
    {
        in_node->Store(split_relperm_a_key, relative_perm_in_a);
        mid_node->Store(split_relperm_a_key, relative_perm_mid_a);
        out_node->Store(split_relperm_a_key, relative_perm_out_a);
    }
    // Here we assumed Kr = S linear.

}

template<uint32_t dim>
void SplitRegionLeakage<dim>::ComputeHydraulicConductivity()
{
    hydraulic_conductivity_in_l = hydraulic_conductivity_in_v =
        hydraulic_conductivity_out_l = hydraulic_conductivity_out_v = 0.;
    if(with_air_) hydraulic_conductivity_in_a = hydraulic_conductivity_out_a = 0.;

    // - If we are modeling the fault as an open space
    //  should we use external permeabilities?

    // - We always upstream rel_perm and viscosity!

    //liquid
    if (grad_pot_in_l()<=0.)
    {
        //flow l in to mid
        if ( !open_space_ )
        {
            if (mul_in() > 0.) hydraulic_conductivity_in_l = permeability_m() * relative_perm_in_l() / mul_in();
        }
        else
            if (mul_in() > 0.) hydraulic_conductivity_in_l = permeability_i() * relative_perm_in_l() / mul_in();
    }
    else
    {
        //flow l mid to in
        if ( !open_space_ )
        {
            if (mul_mid() > 0.) hydraulic_conductivity_in_l = permeability_m() * relative_perm_mid_l() / mul_mid();
        }
        else
            if (mul_mid() > 0.) hydraulic_conductivity_in_l = permeability_i() * relative_perm_mid_l() / mul_mid();
    }

    //vapor
    if (grad_pot_in_v()<=0.)
    {
        //flow v in to mid
        if ( !open_space_ )
        {
            if (muv_in() > 0.) hydraulic_conductivity_in_v = permeability_m() * relative_perm_in_v() / muv_in();
        }
        else
            if (muv_in() > 0.) hydraulic_conductivity_in_v = permeability_i() * relative_perm_in_v() / muv_in();
    }
    else
    {
        //flow v mid to in
        if ( !open_space_ )
        {
            if (muv_mid() > 0.) hydraulic_conductivity_in_v = permeability_m() * relative_perm_mid_v() / muv_mid();
        }
        else
            if (muv_mid() > 0.) hydraulic_conductivity_in_v = permeability_i() * relative_perm_mid_v() / muv_mid();
    }

    //liquid
    if (grad_pot_out_l()<=0.)
    {
        //flow l out to mid
        if ( !open_space_ )
        {
            if (mul_out() > 0.) hydraulic_conductivity_out_l = permeability_m() * relative_perm_out_l() / mul_out();
        }
        else
            if (mul_out() > 0.) hydraulic_conductivity_out_l = permeability_o() * relative_perm_out_l() / mul_out();
    }
    else
    {
        //flow l mid to out
        if ( !open_space_ )
        {
            if (mul_mid() > 0.) hydraulic_conductivity_out_l = permeability_m() * relative_perm_mid_l() / mul_mid();
        }
        else
            if (mul_mid() > 0.) hydraulic_conductivity_out_l = permeability_o() * relative_perm_mid_l() / mul_mid();
    }

    //vapor
    if (grad_pot_out_v()<=0.)
    {
        //flow v out to mid
        if ( !open_space_ )
        {
            if (muv_out() > 0.) hydraulic_conductivity_out_v = permeability_m() * relative_perm_out_v() / muv_out();
        }
        else
            if (muv_out() > 0.) hydraulic_conductivity_out_v = permeability_o() * relative_perm_out_v() / muv_out();
    }
    else
    {
        //flow l mid to out
        if ( !open_space_ )
        {
            if (muv_mid() > 0.) hydraulic_conductivity_out_v = permeability_m() * relative_perm_mid_v() / muv_mid();
        }
        else
            if (muv_mid() > 0.) hydraulic_conductivity_out_v = permeability_o() * relative_perm_mid_v() / muv_mid();
    }

    //air
    if(with_air_)
    {
        if (grad_pot_in_a()<=0.)
        {
            //flow a in to mid
            if ( !open_space_ )
            {
                if (mua_in() > 0.) hydraulic_conductivity_in_a = permeability_m() * relative_perm_in_a() / mua_in();
            }
            else
                if (mua_in() > 0.) hydraulic_conductivity_in_a = permeability_i() * relative_perm_in_a() / mua_in();
        }
        else
        {
            //flow a mid to in
            if ( !open_space_ )
            {
                if (mua_mid() > 0.) hydraulic_conductivity_in_a = permeability_m() * relative_perm_mid_a() / mua_mid();
            }
            else
                if (mua_mid() > 0.) hydraulic_conductivity_in_a = permeability_i() * relative_perm_mid_a() / mua_mid();
        }

        if (grad_pot_out_a()<=0.)
        {
            //flow a out to mid
            if ( !open_space_ )
            {
                if (mua_out() > 0.) hydraulic_conductivity_out_a = permeability_m() * relative_perm_out_a() / mua_out();
            }
            else
                if (mua_out() > 0.) hydraulic_conductivity_out_a = permeability_o() * relative_perm_out_a() / mua_out();
        }
        else
        {
            //flow a mid to out
            if ( !open_space_ )
            {
                if (mua_mid() > 0.) hydraulic_conductivity_out_a = permeability_m() * relative_perm_mid_a() / mua_mid();
            }
            else
                if (mua_mid() > 0.) hydraulic_conductivity_out_a = permeability_o() * relative_perm_mid_a() / mua_mid();
        }
    }
}


template<uint32_t dim>
void SplitRegionLeakage<dim>::ComputeFlowPotential()
{

    grad_pot_in_l() = grad_pot_out_l() = grad_pot_in_v() = grad_pot_out_v() = 0.;
    if(with_air_) { grad_pot_in_a() = grad_pot_out_a() = 0.; }
    grav_component_in_l = grav_component_in_v = grav_component_out_l = grav_component_out_v = 0.;
    if(with_air_) { grav_component_in_a = grav_component_out_a = 0.; }
    gradP_in = gradP_out = 0.;

    Point<dim> in_mid    = in_node->Coordinate()  - mid_node->Coordinate();
    Point<dim> out_mid   = out_node->Coordinate() - mid_node->Coordinate();

    // dx
    distance_in_mid = distance_out_mid = 0.;

    if(in_node == out_node)
    {
        distance_in_mid  = thickness()/2;
        distance_out_mid  = thickness()/2;
    }
    else
    {
        distance_in_mid  = in_mid.Length();
        distance_out_mid = out_mid.Length();
        // We need to pull apart the SB by thickness for this to work!!!!!! Alternative: use thickness/2 instead
    }

    in_mid  = unrml * (distance_in_mid);
    out_mid = -1 * unrml * (distance_out_mid);

    // dy
    double dy_in_mid  = 0.;
    double dy_out_mid = 0.;

    if constexpr(dim == 3U) dy_in_mid = (dotProduct(in_mid, Point<3U>(0., 1., 0.)));
    if constexpr(dim == 2U) dy_in_mid = (dotProduct(in_mid, Point<2U>(0., 1.)));
    if constexpr(dim == 3U) dy_out_mid = (dotProduct(out_mid, Point<3U>(0., 1., 0.)));
    if constexpr(dim == 2U) dy_out_mid = (dotProduct(out_mid, Point<2U>(0., 1.)));

    //rho*g*dy/dx, here we use rho and not mass concentration!
    if( distance_in_mid > 0. && (in_node != out_node))//no calculation at periphery
    {
        // grav_component_in_l = 9.80665 * rhol_mid() * dy_in_mid / (distance_in_mid);
        // grav_component_in_v = 9.80665 * rhov_mid() * dy_in_mid / (distance_in_mid);
        // if(with_air_)
        //     grav_component_in_a = 9.80665 * rhoa_mid() * dy_in_mid / (distance_in_mid);

        //24-04-2026 test improve:
        double avg_rhol_in_mid = (rhol_in() + rhol_mid()) / 2.0;
        double avg_rhov_in_mid = (rhov_in() + rhov_mid()) / 2.0;
        grav_component_in_l = 9.80665 * avg_rhol_in_mid * dy_in_mid / (distance_in_mid);
        grav_component_in_v = 9.80665 * avg_rhov_in_mid * dy_in_mid / (distance_in_mid);
        if(with_air_)
        {
            double avg_rhoa_in_mid = (rhoa_in() + rhoa_mid()) / 2.0;
            grav_component_in_a = 9.80665 * avg_rhoa_in_mid * dy_in_mid / (distance_in_mid);
        }

        gradP_in  = (pressure_mid() - pressure_in())  / (distance_in_mid);
    }

    if( distance_out_mid > 0.)
    {
        if(in_node != out_node)//no gravity component at periphery
        {
            // grav_component_out_l = 9.80665 * rhol_mid() * dy_out_mid / (distance_out_mid);
            // grav_component_out_v = 9.80665 * rhov_mid() * dy_out_mid / (distance_out_mid);
            // if(with_air_)
            //     grav_component_out_a = 9.80665 * rhoa_mid() * dy_out_mid / (distance_out_mid);

            //24-04-2026 test improve:
            double avg_rhol_out_mid = (rhol_out() + rhol_mid()) / 2.0;
            double avg_rhov_out_mid = (rhov_out() + rhov_mid()) / 2.0;
            grav_component_out_l = 9.80665 * avg_rhol_out_mid * dy_out_mid / (distance_out_mid);
            grav_component_out_v = 9.80665 * avg_rhov_out_mid * dy_out_mid / (distance_out_mid);
            if(with_air_)
            {
                double avg_rhoa_out_mid = (rhoa_out() + rhoa_mid()) / 2.0;
                grav_component_out_a = 9.80665 * avg_rhoa_out_mid * dy_out_mid / (distance_out_mid);
            }
        }
        gradP_out = (pressure_mid() - pressure_out()) / (distance_out_mid);
    }

    grad_pot_in_l()  = gradP_in  + grav_component_in_l;
    grad_pot_in_v()  = gradP_in  + grav_component_in_v;
    grad_pot_out_l() = gradP_out + grav_component_out_l;
    grad_pot_out_v() = gradP_out + grav_component_out_v;
    if(with_air_)
    {
        grad_pot_in_a()  = gradP_in  + grav_component_in_a;
        grad_pot_out_a() = gradP_out + grav_component_out_a;
    }
}


//==============================================================================
//  Leakage-in-the-loop: frozen-rate split (see header).
//==============================================================================
template<uint32_t dim>
void SplitRegionLeakage<dim>::PrepareStageRates()
{
    // Run the full pass with dt=1 in capture mode: the flux branches then
    // produce RATES, which the capture gate records per triple instead of
    // clamping/applying. Geometry, upwinding and viz stores happen exactly as
    // in a normal pass (viz values are per-second in this mode — cosmetic).
    const double dt_saved = dt;
    dt = 1.0;
    capture_rates_ = true;
    stage_rates_.clear();
    node_drain_rates_perpv_.clear();
    stale_reported_.clear();
    ComputeGravity_PressureSourceTerm();
    capture_rates_ = false;
    dt = dt_saved;
}

template<uint32_t dim>
void SplitRegionLeakage<dim>::ExchangeSubstep(double sub_dt)
{
    // One relay beat of the interface: per triple, legs = frozen rate * sub_dt,
    // applied against FRESH node state (transport wrote the model between
    // beats) under per-BEAT drain budgets. grad_pot is restored from the
    // stage record for donor selection (frozen contract). Positivity is owned
    // entirely by the budgets here — the per-leg caps of the one-shot path
    // are rate-inert and not re-applied.
    drain_budgets_.clear();     // budgets are per beat
    reported_nodes.clear();
    cut_dt = false;             // clamp mode never requests cuts

    // A silent no-op here starves the fault of matrix recharge with no other
    // symptom (found the hard way) — never fail quietly.
    if ( stage_rates_.empty() )
    {
        static bool warned = false;
        if ( !warned )
        {
            cerr << "\n[LeakLoop][WARNING] ExchangeSubstep called with NO stage"
                    " rates — PrepareStageRates() missing or captured nothing."
                    " Interface exchange is a NO-OP this stage.";
            warned = true;
        }
        return;
    }

    beat_mode_ = true;      // per-beat clamp, never cut
    ledger_active_ = true;
    beat_withheld_ = 0.;    // this beat's standing queue (see ledger block below)
    for ( auto& r : stage_rates_ )
    {
        in_node = r.in_n;  out_node = r.out_n;  mid_node = r.mid_n;
        ReadVariables();        // fresh masses / pore volumes / concentrations

        // ORIGIN PROBE: is a node ALREADY vapor-mass-with-nonpositive-enthalpy
        // BEFORE this beat touches it? If yes, the stale state came from
        // upstream (equilibration / dust reset), not from leakage — leakage
        // conserves E and M exactly (symmetric UTV, single-factor clamp).
        // One-shot per node per stage to avoid spam.
        {
            auto stale = [&]( Node<dim>* n ) {
                if ( stale_reported_.count(n) ) return;
                const double mvn = n->Read(mv_key), hcvn = n->Read(hCv_key);
                const double hvn = n->Read(hv_key);   // specific vapor enthalpy
                if ( mvn > 1.e-6 && hcvn <= 0. ) {
                    cerr << "\n[LeakLoop][STALE-ENTRY] node=" << n->Idx()
                    << " arrives mv=" << mvn << " hCv=" << hcvn << " hv(spec)=" << hvn
                    << ( hvn <= 0. ? "  => specific enthalpy zero (equilibration/thermo origin)"
                                  : "  => content zero but hv>0 (mass/energy split origin)" )
                    << " — upstream of leakage.";
                    stale_reported_.insert(n);
                }
            };
            stale(in_node); stale(mid_node); stale(out_node);
        }

        // ── FLUX RECONSTRUCTION (velocity freeze) ───────────────────────────
        // Rebuild every leg from the frozen VELOCITY and the CURRENT donor
        // state, exactly mirroring the flux-assembly branches (net_Mflux =
        // vel * m_donor * sub_dt; intensive legs = net_Mflux * content_donor /
        // m_donor = vel * content_donor * sub_dt). Donor side is frozen (sign
        // of grad_pot at capture); donor MASS and CONTENT are read fresh this
        // beat, so the flux shrinks as a donor drains and E/M cannot desync.
        // Removal accumulators (mass_*_removed_from_*) are rebuilt here too,
        // since the budget clamp and the CheckDry drain map both consume them.
        auto donor_m = [&]( DonorSide side, double m_matrix, double m_mid )
        { return ( side == DONOR_MATRIX ) ? m_matrix : m_mid; };

        // reset removal accumulators for this triple/beat
        mass_l_removed_from_mid = mass_l_removed_from_in = mass_l_removed_from_out = 0.;
        mass_v_removed_from_mid = mass_v_removed_from_in = mass_v_removed_from_out = 0.;
        mass_a_removed_from_mid = mass_a_removed_from_in = mass_a_removed_from_out = 0.;

        // helper: build one leg (mass + intensive), tally the donor's removal.
        // in_leg==true means the "in" side leg; side picks matrix(in/out) vs mid.
        auto build_leg =
            [&]( double vel, DonorSide side,
                double m_matrix, double m_mid,
                double hC_matrix, double hC_mid,
                double xC_matrix, double xC_mid,
                double TrC_matrix, double TrC_mid,   // liquid only (else pass 0)
                double LiC_matrix, double LiC_mid,   // liquid only
                double& net_M, double& net_E, double& net_X,
                double* net_Tr, double* net_Li,
                double& rem_matrix, double& rem_mid,
                bool has_TrLi ) -> void
        {
            const double m  = donor_m( side, m_matrix, m_mid );
            const double hC = ( side == DONOR_MATRIX ) ? hC_matrix : hC_mid;
            const double xC = ( side == DONOR_MATRIX ) ? xC_matrix : xC_mid;
            net_M = vel * m * sub_dt;
            if ( m > 0. )
            {
                net_E = net_M * hC / m;   // = vel * hC * sub_dt (mass cancels)
                net_X = net_M * xC / m;
                if ( has_TrLi )
                {
                    const double TrC = ( side == DONOR_MATRIX ) ? TrC_matrix : TrC_mid;
                    const double LiC = ( side == DONOR_MATRIX ) ? LiC_matrix : LiC_mid;
                    if ( net_Tr ) *net_Tr = net_M * TrC / m;
                    if ( net_Li ) *net_Li = net_M * LiC / m;
                }
            }
            else { net_E = net_X = 0.; if ( net_Tr ) *net_Tr = 0.; if ( net_Li ) *net_Li = 0.; }
            // removal is tallied on the DONOR side (matrix->in/out or mid)
            ( side == DONOR_MATRIX ? rem_matrix : rem_mid ) += fabs( net_M );
        };

        // liquid legs (carry tracer + lithium)
        build_leg( r.vel_in_l,  r.don_in_l,
                  ml_in(),  ml_mid(),  hCl_in(),  hCl_mid(),  xCl_in(),  xCl_mid(),
                  TrCl_in(), TrCl_mid(), LiCl_in(), LiCl_mid(),
                  net_Mflux_in_l,  net_Eflux_in_l,  net_Xflux_in_l,
                  &net_Trflux_in_l, &net_Liflux_in_l,
                  mass_l_removed_from_in, mass_l_removed_from_mid, true );
        build_leg( r.vel_out_l, r.don_out_l,
                  ml_out(), ml_mid(),  hCl_out(), hCl_mid(),  xCl_out(), xCl_mid(),
                  TrCl_out(), TrCl_mid(), LiCl_out(), LiCl_mid(),
                  net_Mflux_out_l, net_Eflux_out_l, net_Xflux_out_l,
                  &net_Trflux_out_l, &net_Liflux_out_l,
                  mass_l_removed_from_out, mass_l_removed_from_mid, true );

        // vapor legs (no tracer/lithium)
        build_leg( r.vel_in_v,  r.don_in_v,
                  mv_in(),  mv_mid(),  hCv_in(),  hCv_mid(),  xCv_in(),  xCv_mid(),
                  0.,0.,0.,0.,
                  net_Mflux_in_v,  net_Eflux_in_v,  net_Xflux_in_v,
                  nullptr, nullptr,
                  mass_v_removed_from_in, mass_v_removed_from_mid, false );
        build_leg( r.vel_out_v, r.don_out_v,
                  mv_out(), mv_mid(),  hCv_out(), hCv_mid(),  xCv_out(), xCv_mid(),
                  0.,0.,0.,0.,
                  net_Mflux_out_v, net_Eflux_out_v, net_Xflux_out_v,
                  nullptr, nullptr,
                  mass_v_removed_from_out, mass_v_removed_from_mid, false );

        grad_pot_in_l()  = r.gp_in_l;   grad_pot_out_l() = r.gp_out_l;
        grad_pot_in_v()  = r.gp_in_v;   grad_pot_out_v() = r.gp_out_v;

        if ( with_air_ )
        {
            double dummyX_in = 0., dummyX_out = 0.;   // air carries no salt in the assembly
            build_leg( r.vel_in_a,  r.don_in_a,
                      ma_in(),  ma_mid(),  hCa_in(),  hCa_mid(),  0.,0.,
                      0.,0.,0.,0.,
                      net_Mflux_in_a,  net_Eflux_in_a,  dummyX_in,
                      nullptr, nullptr,
                      mass_a_removed_from_in, mass_a_removed_from_mid, false );
            build_leg( r.vel_out_a, r.don_out_a,
                      ma_out(), ma_mid(),  hCa_out(), hCa_mid(),  0.,0.,
                      0.,0.,0.,0.,
                      net_Mflux_out_a, net_Eflux_out_a, dummyX_out,
                      nullptr, nullptr,
                      mass_a_removed_from_out, mass_a_removed_from_mid, false );
            grad_pot_in_a() = r.gp_in_a;  grad_pot_out_a() = r.gp_out_a;
            available_air_mid = ma_mid() * pore_volume_mid();
            available_air_in  = ma_in()  * pore_volume_in();
            available_air_out = ma_out() * pore_volume_out();
        }

        available_liquid_mid = ml_mid() * pore_volume_mid();
        available_liquid_in  = ml_in()  * pore_volume_in();
        available_liquid_out = ml_out() * pore_volume_out();
        available_vapor_mid  = mv_mid() * pore_volume_mid();
        available_vapor_in   = mv_in()  * pore_volume_in();
        available_vapor_out  = mv_out() * pore_volume_out();

        overdraining      = false;
        overdrain_detail_ = "";
        CFL_drain_max = CFL_stuff_max = 0.;   // CFL machinery inert under the clamp

        // (Per-beat DRIFT diagnostic removed: after the velocity freeze the
        //  reconstructed E/M equals fresh specific enthalpy by construction, so
        //  it only reported "a mid drained" — expected and harmless. The E/M
        //  CORRUPT tripwire below stays as cheap insurance; the stage ledger
        //  carries the accuracy story, now split transport vs interface.)

        ApplyDrainBudgets();
        UpdateTransportVariables();
        // PATCH 16/09/2026 (D2.1): the one-shot path stores the viz sources, but in
        // leakage-in-the-loop mode nothing did (capture pass skips it, and it zeroes
        // them first), so "split region mass/energy/salt source" were always 0.
        // Accumulate the clamped legs of every beat: total exchanged over the stage.
        TotalNetFluxes_for_visualization();
    }

    // ── standing-queue update (per beat) ────────────────────────────────────
    // This beat's total withheld IS the queue standing at the fault right now.
    // Its final value (last beat) = stage carryover; its max = worst transient.
    iface_ledger_.parked_now = beat_withheld_;
    if ( beat_withheld_ > iface_ledger_.parked_max )
        iface_ledger_.parked_max = beat_withheld_;

    ledger_active_ = false;
    beat_mode_ = false;
}

//==============================================================================
//  ApplyDrainBudgets — the per-beat/per-pass budget clamp, extracted so the
//  one-shot apply path and ExchangeSubstep share one implementation.
//  Requires: in/out/mid nodes positioned, ReadVariables() done, net_*flux and
//  removal accumulators set (freshly computed OR restored from frozen rates),
//  available_* computed, grad_pot members valid for donor selection (current
//  visit, or restored from the stage-rate record).
//==============================================================================
template<uint32_t dim>
void SplitRegionLeakage<dim>::ApplyDrainBudgets()
{
    const double dm = 1. - sb_drain_margin;
    if ( !beat_mode_ )
        throw csmp::Exception( ERROR, "SplitRegionLeakage::ApplyDrainBudgets",
                              "budget clamp invoked outside a sub-step beat — the one-shot legacy "
                              "path must CUT on overdrain, not clamp. Caller wiring is wrong." );
    {
        // budget accessor: initialize from start-of-pass holdup on first touch
        auto budget_of = [&](Node<dim>* node, double avail_l, double avail_v, double avail_a) -> DrainBudget&
        {
            DrainBudget& b = drain_budgets_[node];
            if ( !b.init )
            {
                b.l = dm * avail_l;
                b.v = dm * avail_v;
                b.a = dm * avail_a;
                b.init = true;
            }
            return b;
        };
        DrainBudget& bud_mid = budget_of(mid_node, available_liquid_mid, available_vapor_mid, with_air_ ? available_air_mid : 0.);
        DrainBudget& bud_in  = budget_of(in_node,  available_liquid_in,  available_vapor_in,  with_air_ ? available_air_in  : 0.);
        DrainBudget& bud_out = budget_of(out_node, available_liquid_out, available_vapor_out, with_air_ ? available_air_out : 0.);
        // periphery: in_node == out_node -> bud_in and bud_out alias the SAME entry

        // scale for this visit's requested removal against remaining budget;
        // decrements the budget by the applied amount
        auto draw = [this](double requested, double& budget) -> double
        {
            requested = fabs(requested);
            if ( requested <= 0. ) return 1.;
            const double allowed = std::min(requested, std::max(0., budget));
            budget -= allowed;
            // interface queue: budget-withheld mass is standing at this node
            // this beat. Tallied into the per-beat accumulator (see the block
            // at the end of ExchangeSubstep). ledger_active_ gates to beats.
            if ( ledger_active_ )
            {
                iface_ledger_.applied += allowed;
                const double w = requested - allowed;
                if ( w > 0. )
                {
                    beat_withheld_ += w;
                    if ( w > iface_ledger_.worst )
                    { iface_ledger_.worst = w; iface_ledger_.worst_idx = mid_node->Idx(); }
                }
            }
            return allowed / requested;
        };

        const double s_l_mid = draw(mass_l_removed_from_mid, bud_mid.l);
        const double s_v_mid = draw(mass_v_removed_from_mid, bud_mid.v);
        const double s_l_in  = draw(mass_l_removed_from_in,  bud_in.l);
        const double s_v_in  = draw(mass_v_removed_from_in,  bud_in.v);
        const double s_l_out = draw(mass_l_removed_from_out, bud_out.l);
        const double s_v_out = draw(mass_v_removed_from_out, bud_out.v);

        // scale one leg: mass + everything riding on it
        auto clamp_leg = [](double s, double& M, double& E, double& X, double& Tr, double& Li)
        {
            if ( s < 1. ) { M *= s; E *= s; X *= s; Tr *= s; Li *= s; }
        };
        double s_min = 1.;

        // donor of each leg follows the flux-assembly branch condition:
        // grad_pot <= 0 : matrix side is donor;  > 0 : mid is donor
        {
            const double s = ( grad_pot_in_l()  <= 0. ) ? s_l_in  : s_l_mid;
            clamp_leg(s, net_Mflux_in_l,  net_Eflux_in_l,  net_Xflux_in_l,  net_Trflux_in_l,  net_Liflux_in_l);
            s_min = std::min(s_min, s);
        }
        {
            const double s = ( grad_pot_out_l() <= 0. ) ? s_l_out : s_l_mid;
            clamp_leg(s, net_Mflux_out_l, net_Eflux_out_l, net_Xflux_out_l, net_Trflux_out_l, net_Liflux_out_l);
            s_min = std::min(s_min, s);
        }
        {
            double dummyTr = 0., dummyLi = 0.;
            const double s_iv = ( grad_pot_in_v()  <= 0. ) ? s_v_in  : s_v_mid;
            clamp_leg(s_iv, net_Mflux_in_v,  net_Eflux_in_v,  net_Xflux_in_v,  dummyTr, dummyLi);
            const double s_ov = ( grad_pot_out_v() <= 0. ) ? s_v_out : s_v_mid;
            clamp_leg(s_ov, net_Mflux_out_v, net_Eflux_out_v, net_Xflux_out_v, dummyTr, dummyLi);
            s_min = std::min({s_min, s_iv, s_ov});
        }
        if ( with_air_ )
        {
            const double s_a_mid = draw(mass_a_removed_from_mid, bud_mid.a);
            const double s_a_in  = draw(mass_a_removed_from_in,  bud_in.a);
            const double s_a_out = draw(mass_a_removed_from_out, bud_out.a);
            double dummyX = 0., dummyTr = 0., dummyLi = 0.;
            const double s_ia = ( grad_pot_in_a()  <= 0. ) ? s_a_in  : s_a_mid;
            clamp_leg(s_ia, net_Mflux_in_a,  net_Eflux_in_a,  dummyX, dummyTr, dummyLi);
            dummyX = 0.;
            const double s_oa = ( grad_pot_out_a() <= 0. ) ? s_a_out : s_a_mid;
            clamp_leg(s_oa, net_Mflux_out_a, net_Eflux_out_a, dummyX, dummyTr, dummyLi);
            s_min = std::min({s_min, s_ia, s_oa});
        }

        if ( s_min < 1. )
        {
            cerr << endl << " OVERDRAIN CLAMP: node=" << mid_node->Idx()
            << (in_node == out_node ? " [PERIPHERY]" : "")
            << "  s_min=" << s_min
            << "  [" << (overdrain_detail_.empty() ? "budget" : overdrain_detail_) << "]";
        }
        overdraining = false;   // resolved by budget clamp — never request a cut
    }
}

template<uint32_t dim>
bool SplitRegionLeakage<dim>::ComputeGravity_PressureSourceTerm()
{ //The main task

    // ════════════════════════════════════════════════════════════════════════
    //  CONTROL FLOW (read this before the loops below)
    // ────────────────────────────────────────────────────────────────────────
    //  We iterate split boundaries -> cells -> nodes. For each (in, mid, out)
    //  node triple we:
    //     1. read state, form flow potential, upwind relperm/conductivity;
    //     2. compute the in<->mid and out<->mid mass/energy/salt fluxes;
    //     3. compute the per-node drainage & stuffing CFL fractions and the
    //        exact `overdraining` mass check;
    //     4. call UpdateTransportVariables(), which — and THIS is the non-obvious
    //        part — both:
    //          (a) tests the CFL/overdraining guards and may set cut_dt = true, and
    //          (b) APPLIES the fluxes to the model via Store(), but only inside its
    //              own `if (!cut_dt)` block. So on a cut, the node's fluxes are
    //              COMPUTED but NOT APPLIED.
    //
    //  cut_dt is the single latch that ties it together:
    //     * It starts false and is only ever set true inside UpdateTransportVariables.
    //     * The three nested `if (!cut_dt)` guards below (boundary/cell/node) then
    //       short-circuit ALL remaining iterations — once any node requests a cut,
    //       no further node is processed or applied.
    //     * We return cut_dt to the scheme; if true, the scheme resets and retries
    //       at a smaller dt (discarding any fluxes applied before the cutting node).
    //
    //  Net: a cut cleanly aborts the whole leakage pass with nothing partially
    //  committed that matters — applied-then-discarded earlier nodes are wiped by
    //  the scheme reset, and the cutting node onward was never applied.
    // ════════════════════════════════════════════════════════════════════════

    cut_dt = false;
    sb_cut_reason_ = "";   // reason captured when a cut is requested (for logging)

    // ── CFL TUNING ────────────────────────────────────────────────────────────
    // Drainage & stuffing are always checked per-phase; overdraining (exact mass,
    // margin-guarded) is the always-on backstop. The thresholds/margin are file-
    // scope constexpr at the top of this file (sb_drain_threshold, sb_stuff_threshold,
    // sb_drain_margin) so this function and UpdateTransportVariables() share them.

    //Legacy gravity mass source for visualization, not used in P equation anymore
    model_ref_.InputPropertyValue("split region gravity mass source", ScalarVariable(ANY, 0.));
    // No need to call model? can just store on SB interface?

    model_ref_.InputPropertyValue("split region mass source", ScalarVariable(ANY, 0.));
    model_ref_.InputPropertyValue("split region energy source", ScalarVariable(ANY, 0.));
    model_ref_.InputPropertyValue("split region salt source", ScalarVariable(ANY, 0.));
    // No need to call model? can just store on SB interface?

    if (with_SB_leakage_report) 
        cerr<<endl<<"open_space: "<<open_space_<<", fault perm anisotropy: "<<fault_perm_anisotropy_;

    // Track nodes that have already triggered the negative-variable warning
    // so each unique node reports only once, regardless of how many parent cells share it.
    reported_nodes.clear();
    drain_budgets_.clear();   // per-pass drain budgets (see header / clamp block)

    for (auto sb = model_ref_.SplitBoundariesBegin(); sb != model_ref_.SplitBoundariesEnd(); ++sb)
    {
        //if (cut_dt) cerr << " SKIP SB ";

        // Short-circuit guard: once UpdateTransportVariables (called per node, far
        // below) sets cut_dt, this and the two inner `if (!cut_dt)` guards (cell,
        // node) skip all remaining work so we fall straight through to `return
        // cut_dt`. The three guards are the same mechanism at three loop levels;
        // they are NOT commented individually below.
        if (!cut_dt)
        {
            for (auto sb_c = sb->second.CellsBegin(); sb_c != sb->second.CellsEnd(); ++sb_c)
            {
                //if (cut_dt) cerr << ", SKIP CELL";

                if (!cut_dt)
                {
                    thickness       = (*sb_c)->InterveningElement()->Read(thickness_key);
                    permeability_m  = (*sb_c)->InterveningElement()->Read(k_key);
                    permeability_i  = (*sb_c)->InnerParent()->Read(k_key);
                    permeability_o  = (*sb_c)->OuterParent()->Read(k_key);

                    uint32_t n_nodes = (*sb_c)->FE()->Nodes();
                    area             = (*sb_c)->InterveningElement()->Volume() / static_cast<double>((*sb_c)->InterveningElement()->Nodes());

                    //cerr<<endl<<"(*sb_c)->FE()->Nodes()="<<(*sb_c)->FE()->Nodes()<<", (*sb_c)->InterveningElement()->Nodes()="<<(*sb_c)->InterveningElement()->Nodes();

                    permeability_m *= fault_perm_anisotropy_;

                    unrml = (*sb_c)->UnitNormal(MIDDLE); //From inside to outside

                    // PATCH 16/09/2026 (D2.1): these two fields were declared but never
                    // written. Conductivity: k*kr/mu of the liquid legs, mean over nodes.
                    const auto viz_hc_key = model_ref_.Database().StorageKey("split region hydraulic conductivity");
                    const auto viz_nA_key = model_ref_.Database().StorageKey("split region nodal area");
                    double viz_hc_sum = 0.;


                    for (uint32_t n{0U}; n < n_nodes; ++n)
                    {
                        //if (cut_dt) cerr << ", SKIP NODE";

                        if (!cut_dt)
                        {
                            in_node   = (*sb_c)->MatchingN(n, INSIDE);
                            out_node  = (*sb_c)->MatchingN(n, OUTSIDE);
                            mid_node  = (*sb_c)->MatchingN(n, MIDDLE);

                            double tot_area = 0;
                            for ( auto p = 0; p<mid_node->Parents(); p++ )
                            {
                                tot_area += mid_node->Parent(p)->Volume() /(mid_node->Parent(p)->Nodes());
                                //cerr<<endl<<"Parent: "<<p<<" has volume: "<<mid_node->Parent(p)->Volume()<<" and nodes: "<<mid_node->Parent(p)->Nodes()<<", tot_area: "<<tot_area;
                            }

                            double area_scale = tot_area/area;
                            //cerr<<endl<<"area_scale:"<<area_scale;

                            ReadVariables();

                            ComputeFlowPotential();

                            ComputeRelPerm();

                            ComputeHydraulicConductivity();

                            // PATCH 16/09/2026 (D2.1): visualization only.
                            viz_hc_sum += 0.5 * (hydraulic_conductivity_in_l + hydraulic_conductivity_out_l);
                            {
                                ScalarVariable nA(ANY, tot_area);
                                mid_node->Store(viz_nA_key, nA);
                                in_node->Store(viz_nA_key, nA);
                                out_node->Store(viz_nA_key, nA);
                            }

                            // Initializing to 0

                            volume_l_removed_from_mid = volume_v_removed_from_mid = 0.;
                            volume_l_removed_from_in = volume_v_removed_from_in = 0.;
                            volume_l_removed_from_out = volume_v_removed_from_out = 0.;

                            volume_l_added_to_mid = volume_v_added_to_mid = 0.;
                            volume_l_added_to_in = volume_v_added_to_in = 0.;
                            volume_l_added_to_out = volume_v_added_to_out = 0.;

                            mass_l_removed_from_mid = mass_v_removed_from_mid = 0.;
                            mass_l_removed_from_in = mass_v_removed_from_in = 0.;
                            mass_l_removed_from_out = mass_v_removed_from_out = 0.;

                            grav_mass_source_in_l = grav_mass_source_out_l = grav_mass_source_in_v = grav_mass_source_out_v = 0.;

                            net_Mflux_in_l = net_Mflux_out_l = net_Eflux_in_l = net_Eflux_out_l = net_Xflux_in_l = net_Xflux_out_l = 0.;
                            net_Mflux_in_v = net_Mflux_out_v = net_Eflux_in_v = net_Eflux_out_v = net_Xflux_in_v = net_Xflux_out_v = 0.;

                            if(with_air_)
                            {
                                volume_a_removed_from_mid = 0.;
                                volume_a_removed_from_in  = 0.;
                                volume_a_removed_from_out = 0.;

                                volume_a_added_to_mid = 0.;
                                volume_a_added_to_in  = 0.;
                                volume_a_added_to_out = 0.;

                                mass_a_removed_from_mid   = 0.;
                                mass_a_removed_from_in    = 0.;
                                mass_a_removed_from_out   = 0.;
                                grav_mass_source_in_a = grav_mass_source_out_a = 0.;
                                net_Mflux_in_a = net_Mflux_out_a = net_Eflux_in_a = net_Eflux_out_a = 0.;
                            }

                            if (with_tracer_)
                                net_Trflux_in_l = net_Trflux_out_l = net_Trflux_in_v = net_Trflux_out_v = 0.;

                            if (with_lithium_)
                                net_Liflux_in_l = net_Liflux_out_l = net_Liflux_in_v = net_Liflux_out_v = 0.;

                            if(with_gold_)
                            {
                                net_so2flux_in_l = net_so2flux_out_l = net_so2flux_in_v = net_so2flux_out_v = 0.;
                                net_auhs0flux_in_l = net_auhs0flux_out_l = net_auhs0flux_in_v = net_auhs0flux_out_v = 0.;
                            }


                            // Fluxes calculation, we consistently use "concentration/content" variables and not thermodynamic variables
                            // to compute fluxes.

                            // liquid
                            // flow in to mid
                            if (grad_pot_in_l() <= 0.)
                            {
                                grav_mass_source_in_l  = grav_component_in_l * hydraulic_conductivity_in_l * area * ml_in() * dt;
                                net_Mflux_in_l  = grad_pot_in_l() * hydraulic_conductivity_in_l * area * ml_in() * dt;

                                //Cap mass flux
                                net_Mflux_in_l  = std::max(-ml_in()*pore_volume_in(), net_Mflux_in_l);
                                mass_l_removed_from_in += fabs(net_Mflux_in_l);

                                if (ml_in() > 0.)
                                {
                                    volume_l_removed_from_in += fabs(net_Mflux_in_l)/ml_in();
                                    volume_l_added_to_mid += fabs(net_Mflux_in_l)/ml_in();

                                    net_Eflux_in_l  = net_Mflux_in_l * hCl_in()  / ml_in();
                                    net_Xflux_in_l  = net_Mflux_in_l * xCl_in()  / ml_in();

                                    if (with_tracer_)
                                        net_Trflux_in_l = net_Mflux_in_l * TrCl_in() / ml_in();

                                    if (with_lithium_)
                                        net_Liflux_in_l = net_Mflux_in_l * LiCl_in() / ml_in();

                                    if(with_gold_)
                                    {
                                        net_so2flux_in_l   = net_Mflux_in_l * so2Cl_in() / ml_in();
                                        net_auhs0flux_in_l = net_Mflux_in_l * auhs0Cl_in() / ml_in();
                                    }
                                }
                            }

                            // flow mid to in
                            else
                            {
                                grav_mass_source_in_l  = grav_component_in_l * hydraulic_conductivity_in_l * area * ml_mid()* dt;
                                net_Mflux_in_l  = grad_pot_in_l() * hydraulic_conductivity_in_l * area * ml_mid() * dt;

                                //Cap mass flux
                                net_Mflux_in_l  = std::min(ml_mid()*pore_volume_mid(), net_Mflux_in_l);
                                mass_l_removed_from_mid += fabs(net_Mflux_in_l);

                                if (ml_mid() > 0.)
                                {
                                    volume_l_removed_from_mid += fabs(net_Mflux_in_l)/ml_mid();
                                    volume_l_added_to_in += fabs(net_Mflux_in_l)/ml_mid();

                                    net_Eflux_in_l  = net_Mflux_in_l * hCl_mid()  / ml_mid();
                                    net_Xflux_in_l  = net_Mflux_in_l * xCl_mid()  / ml_mid();

                                    if (with_tracer_)
                                        net_Trflux_in_l = net_Mflux_in_l * TrCl_mid() / ml_mid();

                                    if (with_lithium_)
                                        net_Liflux_in_l = net_Mflux_in_l * LiCl_mid() / ml_mid();

                                    if(with_gold_)
                                    {
                                        net_so2flux_in_l   = net_Mflux_in_l * so2Cl_mid() / ml_mid();
                                        net_auhs0flux_in_l = net_Mflux_in_l * auhs0Cl_mid() / ml_mid();
                                    }
                                }
                            }

                            // flow out to mid
                            if (grad_pot_out_l() <= 0.)
                            {
                                grav_mass_source_out_l  = grav_component_out_l * hydraulic_conductivity_out_l * area * ml_out()* dt;
                                net_Mflux_out_l  = grad_pot_out_l() * hydraulic_conductivity_out_l * area * ml_out() * dt;

                                //Cap mass flux
                                net_Mflux_out_l  = std::max(-ml_out()*pore_volume_out(), net_Mflux_out_l);
                                mass_l_removed_from_out += fabs(net_Mflux_out_l);

                                if (ml_out() > 0.)
                                {
                                    volume_l_removed_from_out += fabs(net_Mflux_out_l)/ml_out();
                                    volume_l_added_to_mid += fabs(net_Mflux_out_l)/ml_out();

                                    net_Eflux_out_l  = net_Mflux_out_l * hCl_out()  / ml_out();
                                    net_Xflux_out_l  = net_Mflux_out_l * xCl_out()  / ml_out();

                                    if (with_tracer_)
                                        net_Trflux_out_l = net_Mflux_out_l * TrCl_out() / ml_out();

                                    if (with_lithium_)
                                        net_Liflux_out_l = net_Mflux_out_l * LiCl_out() / ml_out();

                                    if(with_gold_)
                                    {
                                        net_so2flux_out_l   = net_Mflux_out_l * so2Cl_out() / ml_out();
                                        net_auhs0flux_out_l = net_Mflux_out_l * auhs0Cl_out() / ml_out();
                                    }
                                }
                            }
                            // flow mid to out
                            else
                            {
                                grav_mass_source_out_l  = grav_component_out_l * hydraulic_conductivity_out_l * area * ml_mid()* dt;
                                net_Mflux_out_l  = grad_pot_out_l() * hydraulic_conductivity_out_l * area * ml_mid() * dt;

                                //Cap mass flux
                                net_Mflux_out_l  = std::min(ml_mid()*pore_volume_mid(), net_Mflux_out_l);
                                mass_l_removed_from_mid += fabs(net_Mflux_out_l);

                                if (ml_mid() > 0.)
                                {
                                    volume_l_removed_from_mid += fabs(net_Mflux_out_l)/ml_mid();
                                    volume_l_added_to_out += fabs(net_Mflux_out_l)/ml_mid();

                                    net_Eflux_out_l = net_Mflux_out_l  * hCl_mid()  / ml_mid();
                                    net_Xflux_out_l = net_Mflux_out_l  * xCl_mid()  / ml_mid();

                                    if (with_tracer_)
                                        net_Trflux_out_l = net_Mflux_out_l * TrCl_mid() / ml_mid();

                                    if (with_lithium_)
                                        net_Liflux_out_l = net_Mflux_out_l * LiCl_mid() / ml_mid();

                                    if(with_gold_)
                                    {
                                        net_so2flux_out_l   = net_Mflux_out_l * so2Cl_mid() / ml_mid();
                                        net_auhs0flux_out_l = net_Mflux_out_l * auhs0Cl_mid() / ml_mid();
                                    }
                                }
                            }

                            // same for vapor
                            // flow in to mid
                            if (grad_pot_in_v() <= 0.)
                            {
                                grav_mass_source_in_v  = grav_component_in_v * hydraulic_conductivity_in_v * area * mv_in() * dt;
                                net_Mflux_in_v  = grad_pot_in_v() * hydraulic_conductivity_in_v * area * mv_in() * dt;

                                //Cap mass flux
                                net_Mflux_in_v  = std::max(-mv_in()*pore_volume_in(), net_Mflux_in_v);
                                mass_v_removed_from_in += fabs(net_Mflux_in_v);

                                if (mv_in() > 0.)
                                {
                                    volume_v_removed_from_in += fabs(net_Mflux_in_v)/mv_in();
                                    volume_v_added_to_mid += fabs(net_Mflux_in_v)/mv_in();

                                    net_Eflux_in_v  = net_Mflux_in_v * hCv_in()  / mv_in();
                                    net_Xflux_in_v  = net_Mflux_in_v * xCv_in()  / mv_in();

                                    if (with_tracer_)
                                        net_Trflux_in_v = net_Mflux_in_v * TrCv_in() / mv_in();

                                    if (with_lithium_)
                                        net_Liflux_in_v = net_Mflux_in_v * LiCv_in() / mv_in();

                                    if(with_gold_)
                                    {
                                        net_so2flux_in_v   = net_Mflux_in_v * so2Cv_in() / mv_in();
                                        net_auhs0flux_in_v = net_Mflux_in_v * auhs0Cv_in() / mv_in();
                                    }
                                }
                            }

                            // flow mid to in
                            else
                            {
                                grav_mass_source_in_v  = grav_component_in_v * hydraulic_conductivity_in_v * area * mv_mid()* dt;
                                net_Mflux_in_v  = grad_pot_in_v() * hydraulic_conductivity_in_v * area * mv_mid() * dt;

                                //Cap mass flux
                                net_Mflux_in_v  = std::min(mv_mid()*pore_volume_mid(), net_Mflux_in_v);
                                mass_v_removed_from_mid += fabs(net_Mflux_in_v);

                                if (mv_mid() > 0.)
                                {
                                    volume_v_removed_from_mid += fabs(net_Mflux_in_v)/mv_mid();
                                    volume_v_added_to_in += fabs(net_Mflux_in_v)/mv_mid();

                                    net_Eflux_in_v  = net_Mflux_in_v * hCv_mid()  / mv_mid();
                                    net_Xflux_in_v  = net_Mflux_in_v * xCv_mid()  / mv_mid();

                                    if (with_tracer_)
                                        net_Trflux_in_v = net_Mflux_in_v * TrCv_mid() / mv_mid();

                                    if (with_lithium_)
                                        net_Liflux_in_v = net_Mflux_in_v * LiCv_mid() / mv_mid();

                                    if(with_gold_)
                                    {
                                        net_so2flux_in_v   = net_Mflux_in_v * so2Cv_mid() / mv_mid();
                                        net_auhs0flux_in_v = net_Mflux_in_v * auhs0Cv_mid() / mv_mid();
                                    }
                                }
                            }

                            // flow out to mid
                            if (grad_pot_out_v() <= 0.)
                            {
                                grav_mass_source_out_v  = grav_component_out_v * hydraulic_conductivity_out_v * area * mv_out()* dt;
                                net_Mflux_out_v  = grad_pot_out_v() * hydraulic_conductivity_out_v * area * mv_out() * dt;

                                //Cap mass flux
                                net_Mflux_out_v  = std::max(-mv_out()*pore_volume_out(), net_Mflux_out_v);
                                mass_v_removed_from_out += fabs(net_Mflux_out_v);

                                if (mv_out() > 0.)
                                {
                                    volume_v_removed_from_out += fabs(net_Mflux_out_v)/mv_out();
                                    volume_v_added_to_mid += fabs(net_Mflux_out_v)/mv_out();

                                    net_Eflux_out_v  = net_Mflux_out_v * hCv_out()  / mv_out();
                                    net_Xflux_out_v  = net_Mflux_out_v * xCv_out()  / mv_out();

                                    if (with_tracer_)
                                        net_Trflux_out_v = net_Mflux_out_v * TrCv_out() / mv_out();

                                    if (with_lithium_)
                                        net_Liflux_out_v = net_Mflux_out_v * LiCv_out() / mv_out();

                                    if(with_gold_)
                                    {
                                        net_so2flux_out_v   = net_Mflux_out_v * so2Cv_out() / mv_out();
                                        net_auhs0flux_out_v = net_Mflux_out_v * auhs0Cv_out() / mv_out();
                                    }
                                }
                            }
                            // flow mid to out
                            else
                            {
                                grav_mass_source_out_v  = grav_component_out_v * hydraulic_conductivity_out_v * area * mv_mid()* dt;
                                net_Mflux_out_v  = grad_pot_out_v() * hydraulic_conductivity_out_v * area * mv_mid() * dt;

                                //Cap mass flux
                                net_Mflux_out_v  = std::min(mv_mid()*pore_volume_mid(), net_Mflux_out_v);
                                mass_v_removed_from_mid += fabs(net_Mflux_out_v);

                                if (mv_mid() > 0.)
                                {
                                    volume_v_removed_from_mid += fabs(net_Mflux_out_v)/mv_mid();
                                    volume_v_added_to_out += fabs(net_Mflux_out_v)/mv_mid();

                                    net_Eflux_out_v = net_Mflux_out_v  * hCv_mid()  / mv_mid();
                                    net_Xflux_out_v = net_Mflux_out_v  * xCv_mid()  / mv_mid();

                                    if (with_tracer_)
                                        net_Trflux_out_v = net_Mflux_out_v * TrCv_mid() / mv_mid();

                                    if (with_lithium_)
                                        net_Liflux_out_v = net_Mflux_out_v * LiCv_mid() / mv_mid();

                                    if(with_gold_)
                                    {
                                        net_so2flux_out_v   = net_Mflux_out_v * so2Cv_mid() / mv_mid();
                                        net_auhs0flux_out_v = net_Mflux_out_v * auhs0Cv_mid() / mv_mid();
                                    }

                                }
                            }

                            // air phase
                            if(with_air_)
                            {
                                // flow in to mid
                                if (grad_pot_in_a() <= 0.)
                                {
                                    grav_mass_source_in_a  = grav_component_in_a * hydraulic_conductivity_in_a * area * ma_in() * dt;
                                    net_Mflux_in_a  = grad_pot_in_a() * hydraulic_conductivity_in_a * area * ma_in() * dt;

                                    //Cap mass flux
                                    net_Mflux_in_a  = std::max(-ma_in()*pore_volume_in(), net_Mflux_in_a);
                                    mass_a_removed_from_in += fabs(net_Mflux_in_a);

                                    if (ma_in() > 0.)
                                    {
                                        volume_a_removed_from_in += fabs(net_Mflux_in_a)/ma_in();
                                        volume_a_added_to_mid += fabs(net_Mflux_in_a)/ma_in();

                                        net_Eflux_in_a  = net_Mflux_in_a * hCa_in()  / ma_in();
                                    }
                                }
                                // flow mid to in
                                else
                                {
                                    grav_mass_source_in_a  = grav_component_in_a * hydraulic_conductivity_in_a * area * ma_mid() * dt;
                                    net_Mflux_in_a  = grad_pot_in_a() * hydraulic_conductivity_in_a * area * ma_mid() * dt;

                                    //Cap mass flux
                                    net_Mflux_in_a  = std::min(ma_mid()*pore_volume_mid(), net_Mflux_in_a);
                                    mass_a_removed_from_mid += fabs(net_Mflux_in_a);

                                    if (ma_mid() > 0.)
                                    {
                                        volume_a_removed_from_mid += fabs(net_Mflux_in_a)/ma_mid();
                                        volume_a_added_to_in += fabs(net_Mflux_in_a)/ma_mid();

                                        net_Eflux_in_a  = net_Mflux_in_a * hCa_mid()  / ma_mid();
                                    }
                                }

                                // flow out to mid
                                if (grad_pot_out_a() <= 0.)
                                {
                                    grav_mass_source_out_a  = grav_component_out_a * hydraulic_conductivity_out_a * area * ma_out() * dt;
                                    net_Mflux_out_a  = grad_pot_out_a() * hydraulic_conductivity_out_a * area * ma_out() * dt;

                                    //Cap mass flux
                                    net_Mflux_out_a  = std::max(-ma_out()*pore_volume_out(), net_Mflux_out_a);
                                    mass_a_removed_from_out += fabs(net_Mflux_out_a);

                                    if (ma_out() > 0.)
                                    {
                                        volume_a_removed_from_out += fabs(net_Mflux_out_a)/ma_out();
                                        volume_a_added_to_mid += fabs(net_Mflux_out_a)/ma_out();

                                        net_Eflux_out_a  = net_Mflux_out_a * hCa_out()  / ma_out();
                                    }
                                }
                                // flow mid to out
                                else
                                {
                                    grav_mass_source_out_a  = grav_component_out_a * hydraulic_conductivity_out_a * area * ma_mid() * dt;
                                    net_Mflux_out_a  = grad_pot_out_a() * hydraulic_conductivity_out_a * area * ma_mid() * dt;

                                    //Cap mass flux
                                    net_Mflux_out_a  = std::min(ma_mid()*pore_volume_mid(), net_Mflux_out_a);
                                    mass_a_removed_from_mid += fabs(net_Mflux_out_a);

                                    if (ma_mid() > 0.)
                                    {
                                        volume_a_removed_from_mid += fabs(net_Mflux_out_a)/ma_mid();
                                        volume_a_added_to_out += fabs(net_Mflux_out_a)/ma_mid();

                                        net_Eflux_out_a = net_Mflux_out_a  * hCa_mid()  / ma_mid();
                                    }
                                }
                            }

                            //----------------------------------------------------------------------------------
                            // Draining and CFl checks

                            //==========================================================
                            //  CAPTURE MODE (PrepareStageRates): dt==1 above, so the
                            //  net_*flux members and removal accumulators now hold
                            //  RATES (per second). Record them per triple, accumulate
                            //  the per-PV drain rates for CheckDry voting, and skip
                            //  clamping/application entirely — nothing moves at Prepare.
                            //  NOTE: the per-leg std::max caps in the branches above are
                            //  inert at dt==1 (rate << absolute available); in per-substep
                            //  mode positivity is owned solely by the per-beat budgets.
                            //----------------------------------------------------------
                            if ( capture_rates_ )
                            {
                                // ── CAPTURE (velocity freeze) ─────────────────────
                                // net_Mflux = grad_pot*K*area*m_donor*dt (dt==1 here).
                                // Store vel = net_Mflux / m_donor : the frozen
                                // velocity coefficient. Donor picked by grad_pot sign
                                // (matrix side when <=0, mid when >0), matching the
                                // flux-assembly branches exactly. Guard against a
                                // zero-mass donor (flux is then 0 anyway -> vel 0).
                                TripleRates r;
                                r.in_n = in_node;  r.out_n = out_node;  r.mid_n = mid_node;

                                auto vel_of = [](double net_Mflux, double m_donor) -> double
                                { return ( m_donor > 0. ) ? net_Mflux / m_donor : 0.; };

                                // liquid in-leg: donor = in (gp<=0) or mid (gp>0)
                                r.don_in_l  = ( grad_pot_in_l()  <= 0. ) ? DONOR_MATRIX : DONOR_MID;
                                r.vel_in_l  = vel_of( net_Mflux_in_l,
                                                    r.don_in_l  == DONOR_MATRIX ? ml_in()  : ml_mid() );
                                r.don_out_l = ( grad_pot_out_l() <= 0. ) ? DONOR_MATRIX : DONOR_MID;
                                r.vel_out_l = vel_of( net_Mflux_out_l,
                                                     r.don_out_l == DONOR_MATRIX ? ml_out() : ml_mid() );
                                r.don_in_v  = ( grad_pot_in_v()  <= 0. ) ? DONOR_MATRIX : DONOR_MID;
                                r.vel_in_v  = vel_of( net_Mflux_in_v,
                                                    r.don_in_v  == DONOR_MATRIX ? mv_in()  : mv_mid() );
                                r.don_out_v = ( grad_pot_out_v() <= 0. ) ? DONOR_MATRIX : DONOR_MID;
                                r.vel_out_v = vel_of( net_Mflux_out_v,
                                                     r.don_out_v == DONOR_MATRIX ? mv_out() : mv_mid() );
                                if ( with_air_ )
                                {
                                    r.don_in_a  = ( grad_pot_in_a()  <= 0. ) ? DONOR_MATRIX : DONOR_MID;
                                    r.vel_in_a  = vel_of( net_Mflux_in_a,
                                                        r.don_in_a  == DONOR_MATRIX ? ma_in()  : ma_mid() );
                                    r.don_out_a = ( grad_pot_out_a() <= 0. ) ? DONOR_MATRIX : DONOR_MID;
                                    r.vel_out_a = vel_of( net_Mflux_out_a,
                                                         r.don_out_a == DONOR_MATRIX ? ma_out() : ma_mid() );
                                }
                                else
                                {
                                    r.don_in_a = r.don_out_a = DONOR_MATRIX;
                                    r.vel_in_a = r.vel_out_a = 0.;
                                }

                                r.gp_in_l = grad_pot_in_l();  r.gp_out_l = grad_pot_out_l();
                                r.gp_in_v = grad_pot_in_v();  r.gp_out_v = grad_pot_out_v();
                                r.gp_in_a = with_air_ ? double(grad_pot_in_a())  : 0.;
                                r.gp_out_a = with_air_ ? double(grad_pot_out_a()) : 0.;
                                // (stage-start holdup snapshots removed: they fed the
                                //  deleted drift diagnostic and were write-only since.
                                //  Drop the mv_*0 / ml_*0 fields from TripleRates too.)
                                stage_rates_.push_back(r);

                                // Per-PV drain rates for CheckDry voting: use the frozen
                                // velocity times the CURRENT-at-capture donor mass /pv.
                                // (Voting is a conservative sizing hint; capture-time
                                //  magnitude is the right scale — beats re-read holdup.)
                                PhaseDrainRates& d_mid = node_drain_rates_perpv_[mid_node];
                                d_mid.l += mass_l_removed_from_mid / pore_volume_mid();
                                d_mid.v += mass_v_removed_from_mid / pore_volume_mid();
                                PhaseDrainRates& d_in = node_drain_rates_perpv_[in_node];
                                d_in.l += mass_l_removed_from_in / pore_volume_in();
                                d_in.v += mass_v_removed_from_in / pore_volume_in();
                                PhaseDrainRates& d_out = node_drain_rates_perpv_[out_node];
                                d_out.l += mass_l_removed_from_out / pore_volume_out();
                                d_out.v += mass_v_removed_from_out / pore_volume_out();
                                if ( with_air_ )
                                {
                                    d_mid.a += mass_a_removed_from_mid / pore_volume_mid();
                                    d_in.a  += mass_a_removed_from_in  / pore_volume_in();
                                    d_out.a += mass_a_removed_from_out / pore_volume_out();
                                }
                            }
                            else
                            {
                                available_liquid_mid = ml_mid() * pore_volume_mid();
                                available_liquid_in = ml_in() * pore_volume_in();
                                available_liquid_out = ml_out() * pore_volume_out();

                                available_vapor_mid = mv_mid() * pore_volume_mid();
                                available_vapor_in = mv_in() * pore_volume_in();
                                available_vapor_out = mv_out() * pore_volume_out();

                                overdraining     = false;
                                overdrain_detail_ = "";   // phase+side(s) that tripped, for the log

                                // Drainage margin: overdraining fires when a node would drain to within
                                // drain_margin of EMPTY (not only past empty). Keeps fault nodes off
                                // exactly-zero, which otherwise seeds an mt<->phase desync in the
                                // equilibrator. Safe against dt-spiral: the drainage flux carries the
                                // donor concentration (net_Mflux ~ ml_donor), so the depletion FRACTION
                                // is independent of how full the node is — emptying proceeds geometrically
                                // over uncut steps; only a genuine over-the-margin drain (dt too big) cuts.
                                const double dm = 1. - sb_drain_margin;   // overdrain margin (namespace-scope knob)

                                // helper: record the first tripping phase/side without spamming
                                auto note_overdrain = [&](const char* phase, const char* side)
                                {
                                    overdraining = true;
                                    if (!overdrain_detail_.empty()) overdrain_detail_ += ",";
                                    overdrain_detail_ += std::string(phase) + "-" + side;
                                };

                                //-----------------------------------------------
                                //Check draining of mid node, phase by phase (margin-guarded)
                                if( definitelyGreaterThan(fabs(mass_l_removed_from_mid), dm*available_liquid_mid, numeric_limits<double>::epsilon()))
                                    note_overdrain("liquid","mid");
                                if( definitelyGreaterThan(fabs(mass_v_removed_from_mid), dm*available_vapor_mid, numeric_limits<double>::epsilon()))
                                    note_overdrain("vapor","mid");

                                //-----------------------------------------------
                                //Check draining of in node, phase by phase
                                if( definitelyGreaterThan(fabs(mass_l_removed_from_in), dm*available_liquid_in, numeric_limits<double>::epsilon()))
                                    note_overdrain("liquid","in");
                                if( definitelyGreaterThan(fabs(mass_v_removed_from_in), dm*available_vapor_in, numeric_limits<double>::epsilon()))
                                    note_overdrain("vapor","in");

                                //-----------------------------------------------
                                //Check draining of out node, phase by phase
                                if( definitelyGreaterThan(fabs(mass_l_removed_from_out), dm*available_liquid_out, numeric_limits<double>::epsilon()))
                                    note_overdrain("liquid","out");
                                if( definitelyGreaterThan(fabs(mass_v_removed_from_out), dm*available_vapor_out, numeric_limits<double>::epsilon()))
                                    note_overdrain("vapor","out");

                                // Air overdraining checks
                                if(with_air_)
                                {
                                    available_air_mid = ma_mid() * pore_volume_mid();
                                    available_air_in  = ma_in()  * pore_volume_in();
                                    available_air_out = ma_out() * pore_volume_out();

                                    if( definitelyGreaterThan(fabs(mass_a_removed_from_mid), dm*available_air_mid, numeric_limits<double>::epsilon()))
                                        note_overdrain("air","mid");
                                    if( definitelyGreaterThan(fabs(mass_a_removed_from_in), dm*available_air_in, numeric_limits<double>::epsilon()))
                                        note_overdrain("air","in");
                                    if( definitelyGreaterThan(fabs(mass_a_removed_from_out), dm*available_air_out, numeric_limits<double>::epsilon()))
                                        note_overdrain("air","out");
                                }

                                //==================================================================
                                //  OVERDRAIN CLAMP  (beat_mode_: per-beat budget clamp)
                                //==================================================================
                                //  Instead of cutting dt on overdrain, cap each donor node's
                                //  CUMULATIVE removal this leakage pass at (1-sb_drain_margin) of
                                //  its START-OF-PASS holdup, via a per-pass budget keyed by Node*
                                //  (drain_budgets_, cleared next to reported_nodes). Every visit
                                //  draws from the shared budget, so every node lands at exactly
                                //  margin*(start holdup) — never zero or negative — independent of
                                //  parent count (per-visit fraction clamping would land at
                                //  margin^P). Keying by Node* also shares one budget across the
                                //  in/out sides of a periphery doublet and across matrix nodes
                                //  common to adjacent triples.
                                //
                                //  SYMMETRY: a mid donating on BOTH legs is scaled symmetrically —
                                //  both legs are summed into one per-phase accumulator, so a single
                                //  draw yields ONE factor for both (proportional split; mid->in
                                //  cannot starve mid->out within a visit). Sequential-priority
                                //  asymmetry exists only across visits of a shared node and between
                                //  the in/out draws of a periphery doublet — endgame-only: it
                                //  redistributes which leg delivers the last of the budget, never
                                //  the total.
                                //
                                //  Legs are scaled per donor-phase by s = min(1, budget/requested),
                                //  budget decremented by the applied removal. E/X/Tr/Li legs are
                                //  net_Mflux * donor concentration ratios and MUST scale by the
                                //  same factor, or enthalpy/salt ships for mass that never moved.
                                //  Receipts are unconstrained, so one pass per visit is exact.
                                //
                                //  COST: mass moved is less than the pressure solve's frozen fluxes
                                //  implied; the mismatch lands in nQ next step, bounded by
                                //  margin*(start holdup) per node per step. Same lagged-closure
                                //  class as elsewhere; the in-step p_mid iteration is the planned
                                //  smooth replacement of this hard stop.
                                //
                                //  WARNING — STUFFING IS NOT GUARDED: the clamp bounds removals
                                //  only. A receiver can still be over-filled (up to the donors'
                                //  budgets), which is the over-fill -> lagged-nQ oscillation the
                                //  stuffing check used to watch. If it resurfaces, that is the
                                //  signal to prioritize the p_mid iteration, not to re-enable the
                                //  estimate-based guard.
                                //
                                //  CFL MACHINERY IS IRRELEVANT UNDER THE CLAMP: the drain/stuff
                                //  fractions below are gated off entirely. Drain cuts would guard
                                //  an outcome the budget already bounds, and the area-scaled
                                //  fractions describe pre-clamp DEMAND, not anything that will be
                                //  applied — logging them would only mislead.
                                //------------------------------------------------------------------
                                // beat mode -> clamp; one-shot legacy -> skip clamp,
                                // the CFL/overdrain detection below requests the cut.
                                if ( beat_mode_ )
                                    ApplyDrainBudgets();



                                //==================================================================
                                //  VOLUME-BASED STEP-SIZE CHECKS  (drainage + stuffing)
                                //==================================================================
                                // VOLUME-BASED CFL FRACTIONS (area-scaled)
                                //
                                // area_scale = tot_area/area lifts this cell-visit's contribution to
                                // the per-node total all parent cells will collectively apply.
                                //==================================================================

                                CFL_outflow_l_mid = CFL_outflow_v_mid = CFL_outflow_a_mid = 0.;
                                CFL_outflow_l_in  = CFL_outflow_v_in  = CFL_outflow_a_in  = 0.;
                                CFL_outflow_l_out = CFL_outflow_v_out = CFL_outflow_a_out = 0.;
                                CFL_inflow_l_mid  = CFL_inflow_v_mid  = CFL_inflow_a_mid  = 0.;
                                CFL_inflow_l_in   = CFL_inflow_v_in   = CFL_inflow_a_in   = 0.;
                                CFL_inflow_l_out  = CFL_inflow_v_out  = CFL_inflow_a_out  = 0.;
                                CFL_drain_max = CFL_stuff_max = 0.;

                                // Drain/stuff CFL fractions are only meaningful when overdrain
                                // CUTS the step. Under the budget clamp they would report
                                // area-scaled PRE-clamp demand — an outcome the budget already
                                // bounds — so the whole machinery is gated off: the maxima stay
                                // at their zero init, no cuts fire, no report triggers.
                                if ( !beat_mode_ )   // legacy one-shot: CFL cut machinery live
                                {
                                    // --- per-phase drainage fractions (volume drained / pore volume) ---
                                    CFL_outflow_l_mid = volume_l_removed_from_mid * area_scale / pore_volume_mid();
                                    CFL_outflow_v_mid = volume_v_removed_from_mid * area_scale / pore_volume_mid();
                                    CFL_outflow_l_in  = volume_l_removed_from_in  * area_scale / pore_volume_in();
                                    CFL_outflow_v_in  = volume_v_removed_from_in  * area_scale / pore_volume_in();
                                    CFL_outflow_l_out = volume_l_removed_from_out * area_scale / pore_volume_out();
                                    CFL_outflow_v_out = volume_v_removed_from_out * area_scale / pore_volume_out();
                                    if(with_air_)
                                    {
                                        CFL_outflow_a_mid = volume_a_removed_from_mid * area_scale / pore_volume_mid();
                                        CFL_outflow_a_in  = volume_a_removed_from_in  * area_scale / pore_volume_in();
                                        CFL_outflow_a_out = volume_a_removed_from_out * area_scale / pore_volume_out();
                                    }

                                    // --- per-phase stuffing fractions (volume added / pore volume) ---
                                    CFL_inflow_l_mid = volume_l_added_to_mid * area_scale / pore_volume_mid();
                                    CFL_inflow_v_mid = volume_v_added_to_mid * area_scale / pore_volume_mid();
                                    CFL_inflow_l_in  = volume_l_added_to_in  * area_scale / pore_volume_in();
                                    CFL_inflow_v_in  = volume_v_added_to_in  * area_scale / pore_volume_in();
                                    CFL_inflow_l_out = volume_l_added_to_out * area_scale / pore_volume_out();
                                    CFL_inflow_v_out = volume_v_added_to_out * area_scale / pore_volume_out();
                                    if(with_air_)
                                    {
                                        CFL_inflow_a_mid = volume_a_added_to_mid * area_scale / pore_volume_mid();
                                        CFL_inflow_a_in  = volume_a_added_to_in  * area_scale / pore_volume_in();
                                        CFL_inflow_a_out = volume_a_added_to_out * area_scale / pore_volume_out();
                                    }


                                    CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_l_mid);
                                    CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_v_mid);
                                    CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_l_in);
                                    CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_v_in);
                                    CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_l_out);
                                    CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_v_out);
                                    if(with_air_)
                                    {
                                        CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_a_mid);
                                        CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_a_in);
                                        CFL_stuff_max = std::max(CFL_stuff_max, CFL_inflow_a_out);
                                    }

                                    CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_l_mid);
                                    CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_v_mid);
                                    CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_l_in);
                                    CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_v_in);
                                    CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_l_out);
                                    CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_v_out);
                                    if(with_air_)
                                    {
                                        CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_a_mid);
                                        CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_a_in);
                                        CFL_drain_max = std::max(CFL_drain_max, CFL_outflow_a_out);
                                    }
                                }
                                //==================================================================

                                //Accumulate and store gravity source terms
                                grav_mass_source_in_tot()  += grav_mass_source_in_l  + grav_mass_source_in_v;
                                grav_mass_source_out_tot() += grav_mass_source_out_l + grav_mass_source_out_v;
                                grav_mass_source_mid_tot() -= grav_mass_source_in_l  + grav_mass_source_in_v + grav_mass_source_out_l + grav_mass_source_out_v;
                                if(with_air_)
                                {
                                    grav_mass_source_in_tot()  += grav_mass_source_in_a;
                                    grav_mass_source_out_tot() += grav_mass_source_out_a;
                                    grav_mass_source_mid_tot() -= grav_mass_source_in_a + grav_mass_source_out_a;
                                }

                                //Legacy mass source for visualization, not used in P equation anymore
                                in_node->Store(grav_mass_source_key, grav_mass_source_in_tot);
                                out_node->Store(grav_mass_source_key, grav_mass_source_out_tot);
                                mid_node->Store(grav_mass_source_key, grav_mass_source_mid_tot);

                                TotalNetFluxes_for_visualization();

                                UpdateTransportVariables();
                            } // end !capture_rates_

                        }
                    } //end nodes loop

                    // PATCH 16/09/2026 (D2.1): visualization only.
                    {
                        ScalarVariable hc(ANY, viz_hc_sum / static_cast<double>(n_nodes));
                        (*sb_c)->InterveningElement()->Store(viz_hc_key, hc);
                    }
                }
            }//end cells loop
        }
    }//end SBs loop

    cerr << "\nReturning cut_dt: ";
    if (cut_dt) cerr << "true\n";
    else cerr << "false\n";
    return cut_dt;
}

template<uint32_t dim>
void SplitRegionLeakage<dim>::UpdateTransportVariables()
{
    //Accumulate and store advection terms
    mass_exchange_term_in_l  = net_Mflux_in_l / pore_volume_in();
    mass_exchange_term_out_l = net_Mflux_out_l / pore_volume_out();
    mass_exchange_term_mid_l = -(net_Mflux_in_l + net_Mflux_out_l) / pore_volume_mid();

    mass_exchange_term_in_v  = net_Mflux_in_v / pore_volume_in();
    mass_exchange_term_out_v = net_Mflux_out_v / pore_volume_out();
    mass_exchange_term_mid_v = -(net_Mflux_in_v + net_Mflux_out_v) / pore_volume_mid();

    energy_exchange_term_in_l  = net_Eflux_in_l / pore_volume_in();
    energy_exchange_term_out_l = net_Eflux_out_l / pore_volume_out();
    energy_exchange_term_mid_l = -(net_Eflux_in_l + net_Eflux_out_l) / pore_volume_mid();

    energy_exchange_term_in_v  = net_Eflux_in_v / pore_volume_in();
    energy_exchange_term_out_v = net_Eflux_out_v / pore_volume_out();
    energy_exchange_term_mid_v = -(net_Eflux_in_v + net_Eflux_out_v) / pore_volume_mid();

    salt_mass_exchange_term_in_l  = net_Xflux_in_l / pore_volume_in();
    salt_mass_exchange_term_out_l = net_Xflux_out_l / pore_volume_out();
    salt_mass_exchange_term_mid_l = -(net_Xflux_in_l + net_Xflux_out_l) / pore_volume_mid();

    salt_mass_exchange_term_in_v  = net_Xflux_in_v / pore_volume_in();
    salt_mass_exchange_term_out_v = net_Xflux_out_v / pore_volume_out();
    salt_mass_exchange_term_mid_v = -(net_Xflux_in_v + net_Xflux_out_v) / pore_volume_mid();

    if (with_tracer_)
    {
        tracer_mass_exchange_term_in_l  = net_Trflux_in_l / pore_volume_in();
        tracer_mass_exchange_term_out_l = net_Trflux_out_l / pore_volume_out();
        tracer_mass_exchange_term_mid_l = -(net_Trflux_in_l + net_Trflux_out_l) / pore_volume_mid();

        tracer_mass_exchange_term_in_v  = net_Trflux_in_v / pore_volume_in();
        tracer_mass_exchange_term_out_v = net_Trflux_out_v / pore_volume_out();
        tracer_mass_exchange_term_mid_v = -(net_Trflux_in_v + net_Trflux_out_v) / pore_volume_mid();
    }

    if (with_lithium_)
    {
        lithium_mass_exchange_term_in_l  = net_Liflux_in_l / pore_volume_in();
        lithium_mass_exchange_term_out_l = net_Liflux_out_l / pore_volume_out();
        lithium_mass_exchange_term_mid_l = -(net_Liflux_in_l + net_Liflux_out_l) / pore_volume_mid();

        lithium_mass_exchange_term_in_v  = net_Liflux_in_v / pore_volume_in();
        lithium_mass_exchange_term_out_v = net_Liflux_out_v / pore_volume_out();
        lithium_mass_exchange_term_mid_v = -(net_Liflux_in_v + net_Liflux_out_v) / pore_volume_mid();
    }

    if (with_gold_)
    {
        so2_mass_exchange_term_in_l  = net_so2flux_in_l / pore_volume_in();
        so2_mass_exchange_term_out_l = net_so2flux_out_l / pore_volume_out();
        so2_mass_exchange_term_mid_l = -(net_so2flux_in_l + net_so2flux_out_l) / pore_volume_mid();

        so2_mass_exchange_term_in_v  = net_so2flux_in_v / pore_volume_in();
        so2_mass_exchange_term_out_v = net_so2flux_out_v / pore_volume_out();
        so2_mass_exchange_term_mid_v = -(net_so2flux_in_v + net_so2flux_out_v) / pore_volume_mid();

        auhs0_mass_exchange_term_in_l  = net_auhs0flux_in_l / pore_volume_in();
        auhs0_mass_exchange_term_out_l = net_auhs0flux_out_l / pore_volume_out();
        auhs0_mass_exchange_term_mid_l = -(net_auhs0flux_in_l + net_auhs0flux_out_l) / pore_volume_mid();

        auhs0_mass_exchange_term_in_v  = net_auhs0flux_in_v / pore_volume_in();
        auhs0_mass_exchange_term_out_v = net_auhs0flux_out_v / pore_volume_out();
        auhs0_mass_exchange_term_mid_v = -(net_auhs0flux_in_v + net_auhs0flux_out_v) / pore_volume_mid();
    }

    if(with_air_)
    {
        mass_exchange_term_in_a  = net_Mflux_in_a / pore_volume_in();
        mass_exchange_term_out_a = net_Mflux_out_a / pore_volume_out();
        mass_exchange_term_mid_a = -(net_Mflux_in_a + net_Mflux_out_a) / pore_volume_mid();

        energy_exchange_term_in_a  = net_Eflux_in_a / pore_volume_in();
        energy_exchange_term_out_a = net_Eflux_out_a / pore_volume_out();
        energy_exchange_term_mid_a = -(net_Eflux_in_a + net_Eflux_out_a) / pore_volume_mid();
    }

    double DRAIN_threshold = sb_drain_threshold;
    double STUFF_threshold = sb_stuff_threshold;
    {
        // Accumulate ALL active triggers into one reason (so a step that both
        // over-drains and over-stuffs reports both), rather than letting the last
        // condition silently overwrite. Format: "node=N | drain(..)=.. | stuff(..)=.. |
        // overdrain(phase-side,..)".
        std::string reasons;

        // Worst-contributor labels: drain and stuff are both reduced per-phase,
        // so the label is the worst phase-side (e.g. liquid-mid).
        // A labelled CFL contribution: a fraction and the phase-side it came from.
        struct PhaseSideCFL { double fraction; const char* label; };

        // Return the label of the largest fraction in a list of phase-side contributions.
        auto label_of_worst = [](std::initializer_list<PhaseSideCFL> contributions) -> std::string {
            const char* worst_label    = "?";
            double      worst_fraction = -1.;
            for (const auto& contribution : contributions)
                if (contribution.fraction > worst_fraction) {
                    worst_fraction = contribution.fraction;
                    worst_label    = contribution.label;
                }
            return worst_label;
        };

        // Drainage is always reduced per-phase -> worst contributor is a phase-side.
        auto worst_drain_label = [&]() -> std::string {
            return label_of_worst({
                                   {CFL_outflow_l_mid,"liquid-mid"},{CFL_outflow_v_mid,"vapor-mid"},{CFL_outflow_a_mid,"air-mid"},
                                   {CFL_outflow_l_in, "liquid-in"}, {CFL_outflow_v_in, "vapor-in"}, {CFL_outflow_a_in, "air-in"},
                                   {CFL_outflow_l_out,"liquid-out"},{CFL_outflow_v_out,"vapor-out"},{CFL_outflow_a_out,"air-out"} });
        };

        // Stuffing is checked per-phase -> worst contributor is a phase-side.
        auto worst_stuff_label = [&]() -> std::string {
            return label_of_worst({
                                   {CFL_inflow_l_mid,"liquid-mid"},{CFL_inflow_v_mid,"vapor-mid"},{CFL_inflow_a_mid,"air-mid"},
                                   {CFL_inflow_l_in, "liquid-in"}, {CFL_inflow_v_in, "vapor-in"}, {CFL_inflow_a_in, "air-in"},
                                   {CFL_inflow_l_out,"liquid-out"},{CFL_inflow_v_out,"vapor-out"},{CFL_inflow_a_out,"air-out"} });

        };

        if( CFL_drain_max > DRAIN_threshold )
        {
            cut_dt = true;
            reasons += " | drain(" + worst_drain_label() + ")=" + std::to_string(CFL_drain_max) + "/" + std::to_string(DRAIN_threshold);
            cerr<<endl<<" Cut dt!!!, drainage fraction: "<<CFL_drain_max;
        }
        if( CFL_stuff_max > STUFF_threshold )
        {
            cut_dt = true;
            reasons += " | stuff(" + worst_stuff_label() + ")=" + std::to_string(CFL_stuff_max) + "/" + std::to_string(STUFF_threshold);
            cerr<<endl<<" Cut dt!!!, stuffing fraction: "<<CFL_stuff_max;
        }
        if( overdraining )
        {
            cut_dt = true;
            reasons += " | overdrain(" + overdrain_detail_ + ")";
            cerr<<endl<<" Cut dt!!!, Overdraining!! ["<<overdrain_detail_<<"]";
        }
        if( cut_dt )
            sb_cut_reason_ = "node=" + std::to_string(mid_node->Idx()) + reasons;
    }

    if (reported_nodes.insert(mid_node).second)
        if (cut_dt or (CFL_drain_max > DRAIN_threshold) or (CFL_stuff_max > STUFF_threshold))
        {
            cerr << endl;
            cerr << endl << "================ SB LEAKAGE REPORT ================";
            // -- verdict first: why are we looking at this node? --------------------
            cerr << endl << "  mid node " << mid_node->Idx()
                 << (in_node == out_node ? "  [PERIPHERY]" : "")
                 << "   dt=" << dt;
            cerr << endl << "  CUT=" << (cut_dt ? "YES" : "no")
                 << "   drain=" << CFL_drain_max << "/" << DRAIN_threshold
                 << "   stuff=" << CFL_stuff_max << "/" << STUFF_threshold
                 << (overdraining ? "   OVERDRAINING" : "");

            // -- flow directions ----------------------------------------------------
            cerr << endl << "  -- flow (liquid potential) --";
            cerr << endl << "     in side:  "
                 << (grad_pot_in_l() < 0. ? "In->Mid" : grad_pot_in_l() > 0. ? "Mid->In" : "none")
                 << "     out side: "
                 << (grad_pot_out_l() < 0. ? "Out->Mid" : grad_pot_out_l() > 0. ? "Mid->Out" : "none");

            // -- geometry -----------------------------------------------------------
            cerr << endl << "  -- geometry  (in, mid, out) --";
            cerr << endl << "     parents=" << static_cast<double>(mid_node->Parents())
                 << "   dist in-mid=" << distance_in_mid << "  out-mid=" << distance_out_mid;
            cerr << endl << "     x: " << in_node->x() << ", " << mid_node->x() << ", " << out_node->x();
            cerr << endl << "     y: " << in_node->y() << ", " << mid_node->y() << ", " << out_node->y();
            cerr << endl << "     z: " << in_node->z() << ", " << mid_node->z() << ", " << out_node->z();
            cerr << endl << "     pore vol: " << pore_volume_in() << ", " << pore_volume_mid() << ", " << pore_volume_out();

            // -- pressure / driving state ------------------------------------------
            cerr << endl << "  -- pressure & drive  (in, mid, out) --";
            cerr << endl << "     p:      " << pressure_in()  << ", " << pressure_mid()  << ", " << pressure_out();
            cerr << endl << "     p_prev: " << ppressure_in() << ", " << ppressure_mid() << ", " << ppressure_out();
            cerr << endl << "     nfvs:   " << nfvs_in()      << ", " << nfvs_mid()      << ", " << nfvs_out();
            cerr << endl << "     perm:   " << permeability_i() << ", " << permeability_m() << ", " << permeability_o();
            cerr << endl << "     gradP in=" << gradP_in << " out=" << gradP_out
                 << "   grad_pot_l in=" << grad_pot_in_l() << " out=" << grad_pot_out_l();
            cerr << endl << "     K_l in=" << hydraulic_conductivity_in_l << " out=" << hydraulic_conductivity_out_l;
            cerr << endl << "     Sl: " << Sl_in() << ", " << Sl_mid() << ", " << Sl_out();

            // -- mass moved vs available  (the drainage safety check) --------------
            cerr << endl << "  -- mass removed vs available  (l | v) --";
            cerr << endl << "     mid: " << mass_l_removed_from_mid << " / " << available_liquid_mid
                 << "  |  " << mass_v_removed_from_mid << " / " << available_vapor_mid;
            cerr << endl << "     in:  " << mass_l_removed_from_in << " / " << available_liquid_in
                 << "  |  " << mass_v_removed_from_in << " / " << available_vapor_in;
            cerr << endl << "     out: " << mass_l_removed_from_out << " / " << available_liquid_out
                 << "  |  " << mass_v_removed_from_out << " / " << available_vapor_out;

            // // -- hydrostatic consistency check -------------------------------------
            // {
            //     const double exp_in  = rhol_mid() * 9.80665 * (in_node->y()  - mid_node->y());
            //     const double exp_out = rhol_mid() * 9.80665 * (out_node->y() - mid_node->y());
            //     cerr << endl << "  -- hydrostatic check (expected | actual | residual) --";
            //     cerr << endl << "     in:  " << exp_in  << " | " << (pressure_mid() - pressure_in())
            //          << " | " << ((pressure_mid() - pressure_in())  - exp_in);
            //     cerr << endl << "     out: " << exp_out << " | " << (pressure_mid() - pressure_out())
            //          << " | " << ((pressure_mid() - pressure_out()) - exp_out);
            // }

            // -- detailed state dump (verbose; uncomment as needed) ----------------
            cerr << endl << "  -- detail --";
            cerr << endl << "     mv:   " << mv_in()   << ", " << mv_mid()   << ", " << mv_out()
                 <<        "   rhov: " << rhov_in() << ", " << rhov_mid() << ", " << rhov_out();
            cerr << endl << "     ml:   " << ml_in()   << ", " << ml_mid()   << ", " << ml_out()
                 <<        "   rhol: " << rhol_in() << ", " << rhol_mid() << ", " << rhol_out();
            cerr << endl << "     net_Mflux_l in=" << net_Mflux_in_l << " out=" << net_Mflux_out_l
                 <<        "   net_Mflux_v in=" << net_Mflux_in_v << " out=" << net_Mflux_out_v;
            cerr << endl << "     mass_exch_l  (in,mid,out): " << mass_exchange_term_in_l
                 << ", " << mass_exchange_term_mid_l << ", " << mass_exchange_term_out_l;
            cerr << endl << "     energy_exch_l(in,mid,out): " << energy_exchange_term_in_l
                 << ", " << energy_exchange_term_mid_l << ", " << energy_exchange_term_out_l;
            cerr << endl << "==================================================" << endl;
        } // end SB Leakage Report


    if (!cut_dt)
    {
        // not used in pde-operator in next step anymore (computed inline there) kept here for visualization,
        in_node->Store(split_grad_pot_l_key, grad_pot_in_l);
        out_node->Store(split_grad_pot_l_key, grad_pot_out_l);
        in_node->Store(split_grad_pot_v_key, grad_pot_in_v);
        out_node->Store(split_grad_pot_v_key, grad_pot_out_v);

        if(with_air_)
        {
            in_node->Store(split_grad_pot_a_key, grad_pot_in_a);
            out_node->Store(split_grad_pot_a_key, grad_pot_out_a);
        }

        // We update transport variables

        ml_in()     += mass_exchange_term_in_l;
        ml_out()    += mass_exchange_term_out_l;
        ml_mid()    += mass_exchange_term_mid_l;

        mv_in()     += mass_exchange_term_in_v;
        mv_out()    += mass_exchange_term_out_v;
        mv_mid()    += mass_exchange_term_mid_v;

        hCl_in()    += energy_exchange_term_in_l;
        hCl_out()   += energy_exchange_term_out_l;
        hCl_mid()   += energy_exchange_term_mid_l;

        hCv_in()    += energy_exchange_term_in_v;
        hCv_out()   += energy_exchange_term_out_v;
        hCv_mid()   += energy_exchange_term_mid_v;

        xCl_in()    += salt_mass_exchange_term_in_l;
        xCl_out()   += salt_mass_exchange_term_out_l;
        xCl_mid()   += salt_mass_exchange_term_mid_l;

        xCv_in()    += salt_mass_exchange_term_in_v;
        xCv_out()   += salt_mass_exchange_term_out_v;
        xCv_mid()   += salt_mass_exchange_term_mid_v;

        if (with_tracer_)
        {
            TrCl_in()    += tracer_mass_exchange_term_in_l;
            TrCl_out()   += tracer_mass_exchange_term_out_l;
            TrCl_mid()   += tracer_mass_exchange_term_mid_l;

            TrCv_in()    += tracer_mass_exchange_term_in_v;
            TrCv_out()   += tracer_mass_exchange_term_out_v;
            TrCv_mid()   += tracer_mass_exchange_term_mid_v;
        }

        if (with_lithium_)
        {
            LiCl_in()    += lithium_mass_exchange_term_in_l;
            LiCl_out()   += lithium_mass_exchange_term_out_l;
            LiCl_mid()   += lithium_mass_exchange_term_mid_l;

            LiCv_in()    += lithium_mass_exchange_term_in_v;
            LiCv_out()   += lithium_mass_exchange_term_out_v;
            LiCv_mid()   += lithium_mass_exchange_term_mid_v;
        }

        if (with_gold_)
        {
            so2Cl_in()    += so2_mass_exchange_term_in_l;
            so2Cl_out()   += so2_mass_exchange_term_out_l;
            so2Cl_mid()   += so2_mass_exchange_term_mid_l;

            so2Cv_in()    += so2_mass_exchange_term_in_v;
            so2Cv_out()   += so2_mass_exchange_term_out_v;
            so2Cv_mid()   += so2_mass_exchange_term_mid_v;

            auhs0Cl_in()    += auhs0_mass_exchange_term_in_l;
            auhs0Cl_out()   += auhs0_mass_exchange_term_out_l;
            auhs0Cl_mid()   += auhs0_mass_exchange_term_mid_l;

            auhs0Cv_in()    += auhs0_mass_exchange_term_in_v;
            auhs0Cv_out()   += auhs0_mass_exchange_term_out_v;
            auhs0Cv_mid()   += auhs0_mass_exchange_term_mid_v;
        }

        if(with_air_)
        {
            ma_in()  += mass_exchange_term_in_a;
            ma_out() += mass_exchange_term_out_a;
            ma_mid() += mass_exchange_term_mid_a;

            hCa_in()  += energy_exchange_term_in_a;
            hCa_out() += energy_exchange_term_out_a;
            hCa_mid() += energy_exchange_term_mid_a;
        }

        // ── DUST-NEGATIVE SNAP (exact-landing residue, 2026-07) ─────────────
        // ═══════════════════════════════════════════════════════════════════
        //  DUST/CONSISTENCY RULE 1 of 3 — SBL beat landings (THE ONLY SBL RULE)
        //  (Rule 2: visitor mt-increment snap. Rule 3: visitor EMPTY-NODE
        //   override. Transport has NO rule — clean by construction, loud
        //   WriteResults abort kept as tripwire.)
        //
        //  The budget clamp scales removals multiplicatively, so zero-landings
        //  leave rounding residue: mass +-1e-13..1e-20, contents up to 1e-6
        //  (energy lives at O(1e6..1e10)). Applied HERE, on the local
        //  variables BEFORE they are stored (fix at the source, no re-reads).
        //  Logic is CONSISTENCY-GATED, not threshold-blanket: a near-empty
        //  phase is zeroed only if its state is IMPOSSIBLE (negative content,
        //  or content exceeding what the mass could carry) — a consistent
        //  tiny parcel (observed: m=2.2e-7 with h=2.8e6 J/kg, physical vapor)
        //  is REAL and passes untouched. A present phase gets only its
        //  negative rounding dust snapped; a large negative anywhere is a real
        //  bug and is WARNED, then zeroed (unphysical either way).
        // ═══════════════════════════════════════════════════════════════════
        {
            constexpr double M_EMPTY    = 1.e-6;  // per-PV: below = near-empty band
            constexpr double H_SPEC_MAX = 6.e6;   // J/kg: above any physical enthalpy
            auto enforce_phase =
                [&]( ScalarVariable& m, ScalarVariable* hC,
                    std::initializer_list<ScalarVariable*> fracs,
                    const char* phase, long node_idx )
            {
                if ( m() < M_EMPTY )
                {
                    // near-empty band: consistency gate
                    if ( m() < -M_EMPTY )
                        cerr << "\n[LeakLoop][phase WARNING] node=" << node_idx
                             << " " << phase << " mass NEGATIVE (m=" << m()
                             << ") — zeroed; investigate upstream";
                    const double carry = std::max( m(), 0. );
                    bool corrupt = ( m() < 0. );
                    if ( hC && ( (*hC)() < 0. || (*hC)() > H_SPEC_MAX * carry ) )
                        corrupt = true;
                    for ( ScalarVariable* c : fracs )
                        if ( (*c)() < 0. || (*c)() > carry ) corrupt = true;
                    if ( corrupt )
                    {
                        m() = 0.;
                        if ( hC ) (*hC)() = 0.;
                        for ( ScalarVariable* c : fracs ) (*c)() = 0.;
                    }
                    return;   // consistent tiny parcel: leave untouched
                }
                // present phase: snap only negative ROUNDING dust on contents
                // (tiny vs carrying capacity); larger negatives are real bugs
                // and stay visible for the loud downstream abort.
                if ( hC && (*hC)() < 0. && (*hC)() > -1.e-9 * H_SPEC_MAX * m() )
                    (*hC)() = 0.;
                for ( ScalarVariable* c : fracs )
                    if ( (*c)() < 0. && (*c)() > -1.e-9 * m() ) (*c)() = 0.;
            };

            enforce_phase( ml_in,  &hCl_in,  { &xCl_in,  &TrCl_in,  &LiCl_in  }, "liquid", in_node->Idx()  );
            enforce_phase( ml_out, &hCl_out, { &xCl_out, &TrCl_out, &LiCl_out }, "liquid", out_node->Idx() );
            enforce_phase( ml_mid, &hCl_mid, { &xCl_mid, &TrCl_mid, &LiCl_mid }, "liquid", mid_node->Idx() );
            enforce_phase( mv_in,  &hCv_in,  { &xCv_in  }, "vapor", in_node->Idx()  );
            enforce_phase( mv_out, &hCv_out, { &xCv_out }, "vapor", out_node->Idx() );
            enforce_phase( mv_mid, &hCv_mid, { &xCv_mid }, "vapor", mid_node->Idx() );
            if ( with_air_ )
            {
                enforce_phase( ma_in,  &hCa_in,  {}, "air", in_node->Idx()  );
                enforce_phase( ma_out, &hCa_out, {}, "air", out_node->Idx() );
                enforce_phase( ma_mid, &hCa_mid, {}, "air", mid_node->Idx() );
            }
        }

        if(
            definitelyLessThan( ml_in(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( ml_out(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( ml_mid(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( mv_in(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( mv_out(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( mv_mid(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( xCl_in(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( xCl_out(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( xCl_mid(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( xCv_in(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( xCv_out(), 0., numeric_limits<double>::epsilon())
            ||
            definitelyLessThan( xCv_mid(), 0., numeric_limits<double>::epsilon())
            ||
            (with_air_ && (
                 definitelyLessThan( ma_in(),  0., numeric_limits<double>::epsilon())
                 ||
                 definitelyLessThan( ma_out(), 0., numeric_limits<double>::epsilon())
                 ||
                 definitelyLessThan( ma_mid(), 0., numeric_limits<double>::epsilon())
                 )))
        {
            cerr<<endl<<"------------------";
            cerr<<endl<<"Mid node ID: "<<mid_node->Idx();
            cerr<<endl<<"Negative transport variable, salt or mass in split boundary!";
        }
    }


    if( in_node != out_node )
    {
        in_node->Store(ml_key, ml_in);
        in_node->Store(mv_key, mv_in);
        in_node->Store(hCl_key, hCl_in);
        in_node->Store(hCv_key, hCv_in);
        in_node->Store(xCl_key, xCl_in);
        in_node->Store(xCv_key, xCv_in);

        if (with_tracer_)
        {
            in_node->Store(TrCl_key, TrCl_in);
            in_node->Store(TrCv_key, TrCv_in);
        }

        if (with_lithium_)
        {
            in_node->Store(LiCl_key, LiCl_in);
            in_node->Store(LiCv_key, LiCv_in);
        }

        if (with_gold_)
        {
            in_node->Store(so2Cl_key, so2Cl_in);
            in_node->Store(so2Cv_key, so2Cv_in);
            in_node->Store(auhs0Cl_key, auhs0Cl_in);
            in_node->Store(auhs0Cv_key, auhs0Cv_in);
        }
        
        if(with_air_)
        {
            in_node->Store(ma_key, ma_in);
            in_node->Store(hCa_key, hCa_in);
        }
        
    }

    mid_node->Store(ml_key, ml_mid);
    mid_node->Store(mv_key, mv_mid);
    mid_node->Store(hCl_key, hCl_mid);
    mid_node->Store(hCv_key, hCv_mid);
    mid_node->Store(xCl_key, xCl_mid);
    mid_node->Store(xCv_key, xCv_mid);
    mid_node->Store(TrCl_key, TrCl_mid);
    if(with_lithium_) mid_node->Store(LiCl_key, LiCl_mid);
    if(with_air_) { mid_node->Store(ma_key, ma_mid); mid_node->Store(hCa_key, hCa_mid); }

    //if( in_node != out_node )
    {
        out_node->Store(ml_key, ml_out);
        out_node->Store(mv_key, mv_out);
        out_node->Store(hCl_key, hCl_out);
        out_node->Store(hCv_key, hCv_out);
        out_node->Store(xCl_key, xCl_out);
        out_node->Store(xCv_key, xCv_out);

        if (with_tracer_)
        {
            mid_node->Store(TrCl_key, TrCl_mid);
            mid_node->Store(TrCv_key, TrCv_mid);
        }

        if (with_lithium_)
        {
            mid_node->Store(LiCl_key, LiCl_mid);
            mid_node->Store(LiCv_key, LiCv_mid);
        }

        if (with_gold_)
        {
            mid_node->Store(so2Cl_key, so2Cl_mid);
            mid_node->Store(so2Cv_key, so2Cv_mid);
            mid_node->Store(auhs0Cl_key, auhs0Cl_mid);
            mid_node->Store(auhs0Cv_key, auhs0Cv_mid);
        }

        if(with_air_) { out_node->Store(ma_key, ma_out); out_node->Store(hCa_key, hCa_out); }
    }
}

template<uint32_t dim>
void SplitRegionLeakage<dim>::TotalNetFluxes_for_visualization() { //for visualization only
    //Accumulate and store
    qM_in()  += net_Mflux_in_l + net_Mflux_in_v;
    qM_out() += net_Mflux_out_l + net_Mflux_out_v;
    qM_mid() -= net_Mflux_in_l + net_Mflux_in_v + net_Mflux_out_l + net_Mflux_out_v;

    qE_in()  += net_Eflux_in_l + net_Eflux_in_v;
    qE_out() += net_Eflux_out_l + net_Eflux_out_v;
    qE_mid() -= net_Eflux_in_l + net_Eflux_in_v + net_Eflux_out_l + net_Eflux_out_v;

    qS_in()  += net_Xflux_in_l + net_Xflux_in_v;
    qS_out() += net_Xflux_out_l + net_Xflux_out_v;
    qS_mid() -= net_Xflux_in_l + net_Xflux_in_v + net_Xflux_out_l + net_Xflux_out_v;

    if(with_air_)
    {
        qM_in()  += net_Mflux_in_a;
        qM_out() += net_Mflux_out_a;
        qM_mid() -= net_Mflux_in_a + net_Mflux_out_a;

        qE_in()  += net_Eflux_in_a;
        qE_out() += net_Eflux_out_a;
        qE_mid() -= net_Eflux_in_a + net_Eflux_out_a;
    }

    if(with_gold_)
    {
        qSO2_in()  += net_so2flux_in_l + net_so2flux_in_v;
        qSO2_out() += net_so2flux_out_l + net_so2flux_out_v;
        qSO2_mid() -= net_so2flux_in_l + net_so2flux_in_v + net_so2flux_out_l + net_so2flux_out_v;

        qAuHS0_in()  += net_auhs0flux_in_l + net_auhs0flux_in_v;
        qAuHS0_out() += net_auhs0flux_out_l + net_auhs0flux_out_v;
        qAuHS0_mid() -= net_auhs0flux_in_l + net_auhs0flux_in_v + net_auhs0flux_out_l + net_auhs0flux_out_v;
    }

    in_node->Store(qM_key, qM_in);
    out_node->Store(qM_key, qM_out);
    mid_node->Store(qM_key, qM_mid);

    in_node->Store(qE_key, qE_in);
    out_node->Store(qE_key, qE_out);
    mid_node->Store(qE_key, qE_mid);

    in_node->Store(qS_key, qS_in);
    out_node->Store(qS_key, qS_out);
    mid_node->Store(qS_key, qS_mid);

    if(with_gold_)
    {
        in_node ->Store(qSO2_key, qSO2_in);
        out_node->Store(qSO2_key, qSO2_out);
        mid_node->Store(qSO2_key, qSO2_mid);

        in_node ->Store(qAuHS0_key, qAuHS0_in);
        out_node->Store(qAuHS0_key, qAuHS0_out);
        mid_node->Store(qAuHS0_key, qAuHS0_mid);
    }
}

template<uint32_t dim>
void SplitRegionLeakage<dim>::WithSpeciesTransport(Model<dim> &model,
                                                   bool with_air,
                                                   bool with_magmatic_fluids,
                                                   bool with_tracer,
                                                   bool with_gold,
                                                   bool with_lithium)
{
    if (with_air)
    {
        with_air_            = with_air;
        viscosity_a_key      = model.Database().StorageKey("viscosity air");
        sat_a_key            = model.Database().StorageKey("saturation air");
        density_a_key        = model.Database().StorageKey("density air");
        ma_key               = model.Database().StorageKey("fluid mass air");
        hCa_key              = model.Database().StorageKey("enthalpy content air");
        split_grad_pot_a_key = model.Database().StorageKey("split region potential gradient a");
        split_relperm_a_key  = model.Database().StorageKey("split region relperm a");
    }

    if (with_magmatic_fluids)
    {
        with_magmatic_fluids_ = with_magmatic_fluids;
    }

    if (with_tracer)
    {
        with_tracer_ = with_tracer;
        TrCl_key    = model.Database().StorageKey("tracer content liquid");
        TrCv_key    = model.Database().StorageKey("tracer content vapor");
    }

    if (with_gold)
    {
        with_gold_  = with_gold;
        so2Cl_key   = model.Database().StorageKey("so2 content liquid");
        so2Cv_key   = model.Database().StorageKey("so2 content vapor");
        qSO2_key    = model.Database().StorageKey("split region SO2 source");

        auhs0Cl_key = model.Database().StorageKey("auhs0 content liquid");
        auhs0Cv_key = model.Database().StorageKey("auhs0 content vapor");
        qAuHS0_key  = model.Database().StorageKey("split region AuHS0 source");
    }

    if (with_lithium)
    {
        with_lithium_= with_lithium;
        LiCl_key     = model.Database().StorageKey("li content liquid");
        LiCv_key     = model.Database().StorageKey("li content vapor");
    }
}

template class SplitRegionLeakage<1U>;
template class SplitRegionLeakage<2U>;
template class SplitRegionLeakage<3U>;

} // csmp
