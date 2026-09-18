// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FailureModeVisitor.h"

using namespace std;

namespace csmp
{

  template<uint32_t dim>
  FailureModeVisitor<dim>::FailureModeVisitor( Model<dim> &sg, double cohesive_strength,
                                               double differential_stress, bool near_critical )
    : C(cohesive_strength),
      diff_stress(differential_stress),
      normal_fault_offset(4.*(diff_stress - C) / 3.),
      critical(near_critical),
      BD_changed(false),
      BD_low(360.),
      BD_high(500.),
      BD_diff(BD_high - BD_low)

  {
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    failure_pressure_key        = sg.Database().StorageKey("failure pressure");
    lithostatic_pressure_key    = sg.Database().StorageKey("lithostatic pressure");
    reference_pressure_key      = sg.Database().StorageKey("reference pressure");
    temperature_key             = sg.Database().StorageKey("temperature");
    fluid_pressure_key          = sg.Database().StorageKey("fluid pressure");
    thp_key                     = sg.Database().StorageKey("threshold pressure");
    op_key                      = sg.Database().StorageKey("failure overpressure");
    thop_key                    = sg.Database().StorageKey("threshold overpressure");

    if (failure_pressure_key.type != SCALAR || failure_pressure_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "FailureModeVisitor::(constructor)",
                            "failure pressure", " must be a nodal scalar property." );

    if ( lithostatic_pressure_key.type != SCALAR || lithostatic_pressure_key.place != NODE )
      throw csmp::Exception(FATAL_ERROR, "FailureModeVisitor::(constructor)",
                            "lithostaic pressure", " must be a nodal scalar property." );

    if ( temperature_key.type != SCALAR || temperature_key.place != NODE )
      throw csmp::Exception(FATAL_ERROR, "FailureModeVisitor::(constructor)",
                            "temperature", " must be a nodal scalar property." );
  }

  template<uint32_t dim>
  FailureModeVisitor<dim>::~FailureModeVisitor()
  {}

  /** visit function for Region */
  template<uint32_t dim>
  void FailureModeVisitor<dim>::Visit(Model<dim> *n)
  {
    // no calcution for the region
  }

  template<uint32_t dim>
  void FailureModeVisitor<dim>::Visit(Node<dim> *n)
  {
    n->Read( lithostatic_pressure_key, lithostatic_pressure );
    n->Read( reference_pressure_key, reference_pressure );
    n->Read( temperature_key, temperature );
    n->Read( fluid_pressure_key, fluid_pressure );

    if ( !critical )
      {
        failure_pressure() = lithostatic_pressure() - normal_fault_offset;

        if ( failure_pressure() < reference_pressure())
          {
            failure_pressure() = reference_pressure();
          }
      }

    else
      {
        diff_stress = C + (1. - reference_pressure() / lithostatic_pressure()) * 3. / 4.*lithostatic_pressure();

        if ( diff_stress < 2.5 * C )
          {
            diff_stress = 0.5 * C + (1. - reference_pressure() / lithostatic_pressure()) * lithostatic_pressure();
          }

        diff_stress -= 0.5 * C; // 0.5*C = tensile strength

        if ( diff_stress < 0. )
          {
            diff_stress = 0.;
          }

        if ( temperature_relaxation )
          {
            if ( !BD_changed )
              {
                if ( temperature() > 360. && temperature() < 500. )
                  {
                    diff_stress *= pow((500. - temperature()) / 140., 2.);
                  }

                else if ( temperature() >= 500. )
                  {
                    diff_stress = 0.;
                  }
              }
            else
              {
                if ( temperature() > BD_low && temperature() < BD_high )
                  {
                    diff_stress *= pow((BD_high - temperature()) / BD_diff, 2.);
                  }

                else if ( temperature() >= BD_high )
                  {
                    diff_stress = 0.;
                  }
              }

          }

        failure_pressure() = (4.*C - lithostatic_pressure() + 4.*(lithostatic_pressure() - diff_stress)) / 3.;

        if ( diff_stress < 2.5 * C )
          {
            failure_pressure() = 0.5 * C + lithostatic_pressure() - diff_stress;
            //      if (failure_pressure() > 1.9e8)//BB what is this about?
            //        failure_pressure() = lithostatic_pressure();
          }
      }

    thp() = lithostatic_pressure() + 10.0e6; //BB, was 10e6
    op() = fluid_pressure() - failure_pressure(); //BB
    thop() = fluid_pressure() - thp(); //BB

    n->Store( thp_key, thp ); // threshold pressure
    n->Store( failure_pressure_key, failure_pressure );
    n->Store( op_key, op );//BB
    n->Store( thop_key, thop );//BB
  }

  template<uint32_t dim>
  void FailureModeVisitor<dim>::SetTemperatureRelaxationTo( bool relax )
  {
    temperature_relaxation = relax;
  }

  template<uint32_t dim>
  void FailureModeVisitor<dim>::ChangeBrittleDuctileTransitionTemperature( double transition_start, double transition_end )
  {
    BD_changed = true;
    BD_low = std::min(transition_start, transition_end);
    BD_high = std::max(transition_start, transition_end);
    BD_diff = BD_high - BD_low;
  }

  template class FailureModeVisitor<1U>;
  template class FailureModeVisitor<2U>;
  template class FailureModeVisitor<3U>;

} // csmp
