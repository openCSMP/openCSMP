// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TracerVisitor_H
#define TracerVisitor_H

#include "Visitor.h"
#include "Model.h"
#include "Node.h"
#include "compareFloats.h"

/** @file TracerVisitor.h
 *  @author Nicolas Krattiger, Jonas Köpping, Thierry Solms
 *  @brief Simple and basic tracer visitor
 *  @date December 2023
 */

/* Visitor template to write a basic tracer.
 * Can be modified to be applied in specific regions
 * defined by e.g., coordinates or property ranges.
 */

namespace csmp
{
  //  class Index;
  template<uint32_t> class PropertyDatabase;

  template <uint32_t dim>
  class TracerVisitor : public Visitor<dim>
  {
    public:
      TracerVisitor( Model<dim> &model);

      virtual ~TracerVisitor();

      virtual void Visit( Region<dim> *r);
      virtual void Visit( Model<dim> *m);
      virtual void Visit( Element<dim> *e);
      virtual void Visit( Node<dim> *n);

      void ReadVariables( Node<dim> *n );
      void StoreVariables( Node<dim> *n );
      void CalculateMobilities(Model<dim> *model);
      void CalculateMassFractions(Model<dim> *model);
      void SetTracerContent(Model<dim> *model, double timestep);
      void DistributeTracerBetweenPhases(Model<dim> *model);
      void ApplyBoundaryFlowCorrections(Model<dim> *model);
      void InitializeTracerVisitor(Model<dim> *model);
      void InitializeMagmaticConcentrations(Model<dim> *model, double timestep);
      void SetTimeIncrement( double time_increment, double model_time_);

      void ReadLiquidusTemperature(Model<dim> *model);
      void ReadSolidusTemperature(Model<dim> *model);

    private:

      TracerVisitor();

      Model<dim> &model_ref_;

      double dt, dt_max, model_time;

      double tracerscaling = 1.e-12;

      csmp::Index
      c_tr_fl_key, c_tr_l_key, c_tr_v_key, massfrac_tr_fl_key, massfrac_tr_fl_rescaled_key,
                   massfrac_tr_l_key, massfrac_tr_v_key, mob_tr_l_key, mob_tr_v_key, mml_key, mmv_key,
                   temperature_key, fluid_pressure_key, lithostatic_pressure_key, ml_key, mv_key, mt_key,
                   bfm_key, pore_volume_key, tr_mass_out_total_key,
                   coord_x_key, coord_y_key, fluid_enthalpy_key;

      ScalarVariable
      c_tr_fl, c_tr_l, c_tr_v, massfrac_tr_fl, massfrac_tr_fl_rescaled, massfrac_tr_l,
               massfrac_tr_v, mob_tr_l, mob_tr_v, mml, mmv,
               temperature, fluid_pressure, lithostatic_pressure, ml, mv, mt, bfm, pore_volume,
               tr_mass_out_total,
               coord_x, coord_y, fluid_enthalpy;
  };
}

#endif

