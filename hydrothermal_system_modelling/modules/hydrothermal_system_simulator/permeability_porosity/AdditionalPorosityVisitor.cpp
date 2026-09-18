// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "AdditionalPorosityVisitor.h"

using namespace std;

namespace csmp
{

  template<int32_t dim>
  AdditionalPorosityVisitor<dim>::AdditionalPorosityVisitor( Model<dim> &model )
    : brine(tem, pre, sal), maximum_porosity(0.0)
  {

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    p_key = model.Database().StorageKey("fluid pressure");
    mt_key = model.Database().StorageKey("fluid density");
    rho_key = model.Database().StorageKey("bulk fluid density");
    lp_key = model.Database().StorageKey("failure pressure");
    T_key = model.Database().StorageKey("temperature");
    wt_key = model.Database().StorageKey("salinity");
    phi_key = model.Database().StorageKey("nodal porosity");
    d_phi_key = model.Database().StorageKey("additional porosity");
    dfp_key = model.Database().StorageKey("dynamic fluid pressure");

    if (p_key.type != SCALAR || p_key.place != NODE)
      throw Exception(FATAL_ERROR, "AdditionalPorosityVisitor::(constructor)",
                      "fluid pressure", " must be a nodal scalar property." );

  }


  template<int32_t dim>
  AdditionalPorosityVisitor<dim>::~AdditionalPorosityVisitor()
  {}

  /** visit function for Region */
  template<int32_t dim>
  void AdditionalPorosityVisitor<dim>::Visit(Model<dim> *n)
  {
    // no calcution for the region
  }

  template<int32_t dim>
  void AdditionalPorosityVisitor<dim>::Visit(Node<dim> *n)
  {

    additional_porosity() = 0.0;

    n->Read(p_key, pressure);
    n->Read(lp_key, failure_pressure);

    if (pressure() >= failure_pressure())
      {
        dynamic_pressure() = failure_pressure();
        n->Read(T_key, temperature);
        n->Read(wt_key, salinity);
        n->Read(mt_key, mt);
        n->Read(rho_key, rho);
        n->Read(phi_key, porosity);

        tem = temperature();
        pre = pressure();
        sal = salinity() / 100.;
        current_density = brine.Density();

        pre = failure_pressure();
        target_density = brine.Density();

        factor = current_density / target_density - 1.0;
        additional_porosity() = porosity() * factor;

        if (additional_porosity() > maximum_porosity)
          maximum_porosity = additional_porosity();

      }
    else
      dynamic_pressure() = pressure();

    n->Store(d_phi_key, additional_porosity);
    n->Store(dfp_key, dynamic_pressure);

  }


  template class AdditionalPorosityVisitor<1U>;
  template class AdditionalPorosityVisitor<2U>;
  template class AdditionalPorosityVisitor<3U>;

} // csmp
