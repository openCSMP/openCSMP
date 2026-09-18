// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef AdditionalPorosityVisitor_h
#define AdditionalPorosityVisitor_h

#include "Visitor.h"
#include "Model.h"
#include "Brine.h"

#include "Exception.h"

namespace csmp
{

  template<int32_t dim>
  class AdditionalPorosityVisitor : public Visitor<dim>
  {
    public:
      AdditionalPorosityVisitor( Model<dim> &model );
      ~AdditionalPorosityVisitor();

      virtual void Visit(Node<dim> *n);
      virtual void Visit(Model<dim> *n);

    private:

      csmp::Index     p_key, mt_key, rho_key, lp_key, T_key, wt_key, phi_key, d_phi_key,
           dfp_key;
      ScalarVariable pressure, failure_pressure, temperature, salinity, porosity, mt, rho;
      ScalarVariable additional_porosity, dynamic_pressure;
      double current_density, target_density, factor, maximum_porosity;
      ;

      double tem, pre, sal;
      Brine      brine;

  };

} // csmp

#endif
