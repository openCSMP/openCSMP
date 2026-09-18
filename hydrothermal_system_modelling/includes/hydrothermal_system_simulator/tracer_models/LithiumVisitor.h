// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LithiumVisitor_H
#define LithiumVisitor_H

#include "Visitor.h"
#include "Model.h"
#include "Node.h"
#include "compareFloats.h"

/** @file LithiumVisitor.h
 *  @author Nicolas Krattiger, Thierry Solms, Jonas Köpping
 *  @brief Simple and basic tracer visitor
 *  @date December 2023
 */

/* LithiumTracer based on the TracerVisitor template.
 * This visitor accounts for metal partitioning between
 * melt and fluid, melt and crystals, and liquid and vapor.
 *
 * Used in https://doi.org/10.3929/ethz-c-000790075.
 */

namespace csmp
{
  //  class Index;
  template<uint32_t> class PropertyDatabase;

  template <uint32_t dim>
  class LithiumVisitor : public Visitor<dim>
  {
    public:
      LithiumVisitor( Model<dim> &model);

      virtual ~LithiumVisitor();

      virtual void Visit( Region<dim> *r);
      virtual void Visit( Model<dim> *m);
      virtual void Visit( Element<dim> *e);
      virtual void Visit( Node<dim> *n);

      void CalculateMassFractionsAndMobilities(Model<dim> *model);
      void Partitioning(Model<dim> *model);
      void ApplyBoundaryFlowCorrections(Model<dim> *model);
      void InitializeLithiumVisitor(Model<dim> *model);
      void TranslateContentsToConcentrations(Model<dim> *model);
      double ReadPartitionCoefficient(double lv_partition_slope, double rho_liq,
                                      double rho_vap);

      void ReadInputArguments(double ext_partition_coefficient_lithium_fluid_melt,
                              double ext_partition_coefficient_lithium_xtel_melt,
                              double ext_partition_coefficient_lithium_vl_slope,
                              double ext_li_concentration_melt_initialize);

    private:

      LithiumVisitor();

      Model<dim> &model_ref_;

      double D_li_fm;
      double D_li_cm;
      double D_li_vl_slope;
      double li_concentration_melt_initialize;

      csmp::Index
      meltmass_key, crystalmass_key,
                    li_mass_liquid_key, li_mass_vapor_key, li_mass_melt_key, li_mass_crystals_key,
                    li_mass_fluid_key,
                    li_content_fluid_key, li_content_liquid_key, li_content_vapor_key, li_content_melt_key,
                    li_content_crystals_key,
                    massfrac_li_fluid_key, massfrac_li_fluid_rescaled_key, massfrac_li_liquid_key,
                    massfrac_li_vapor_key, mob_li_liquid_key, mob_li_vapor_key, mml_key, mmv_key,
                    temperature_key, fluid_pressure_key, lithostatic_pressure_key, ml_key, mv_key, mt_key,
                    bfm_key, pore_volume_key, bulk_volume_key, previous_pore_volume_key,
                    li_concentration_fluid_key, li_concentration_liquid_key, li_concentration_vapor_key,
                    li_concentration_melt_key, li_concentration_crystals_key,
                    fluid_state_key, rho_l_key, rho_v_key,
                    li_mass_out_total_key, bfm_mass_out_total_key,
                    D_li_vl_key,
                    li_mass_total_key;

      ScalarVariable
      meltmass, crystalmass,
                li_mass_liquid, li_mass_vapor, li_mass_melt, li_mass_crystals, li_mass_fluid,
                li_content_fluid, li_content_liquid, li_content_vapor, li_content_melt,
                li_content_crystals,
                massfrac_li_fluid, massfrac_li_fluid_rescaled, massfrac_li_liquid, massfrac_li_vapor,
                mob_li_liquid, mob_li_vapor, mml, mmv,
                temperature, fluid_pressure, lithostatic_pressure, ml, mv, mt, bfm, pore_volume,
                bulk_volume, previous_pore_volume,
                li_concentration_fluid, li_concentration_liquid, li_concentration_vapor,
                li_concentration_melt, li_concentration_crystals,
                fluid_state, rho_l, rho_v,
                li_mass_out_total, bfm_mass_out_total,
                D_li_vl,
                li_mass_total_output;

      double li_molar_mass = 6.941 * 1e-3; //kg mol-1
  };
}

#endif

