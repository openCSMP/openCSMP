// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FailureModeVisitor_h
#define FailureModeVisitor_h

#include "Visitor.h"
#include "Model.h"

#include "ErrorHandler.h"

namespace csmp
{

  template<uint32_t dim>
  class FailureModeVisitor : public Visitor<dim>
  {
    public:
      FailureModeVisitor( Model<dim> &sg, double cohesive_strength,
                          double differential_stress, bool near_critical );
      ~FailureModeVisitor();

      virtual void Visit(Node<dim> *n);
      virtual void Visit(Model<dim> *n);
      void SetTemperatureRelaxationTo(bool relax);
      void ChangeBrittleDuctileTransitionTemperature( double transition_start,
                                                      double transition_end );

    private:

      bool critical, temperature_relaxation;
      csmp::Index     failure_pressure_key, reference_pressure_key, lithostatic_pressure_key,
           temperature_key, thp_key, op_key, thop_key, fluid_pressure_key, mp_key;
      ScalarVariable  failure_pressure, reference_pressure, lithostatic_pressure, temperature,
                      thp, op, thop, fluid_pressure, mp;

      double C, diff_stress, normal_fault_offset;

      bool BD_changed;
      double   BD_low, BD_high, BD_diff;

  };

} // csmp

#endif
