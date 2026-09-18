// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PorosityPermeabilityCouplingVisitor_H
#define PorosityPermeabilityCouplingVisitor_H

#include "Visitor.h"
#include "Model.h"
#include "Element.h"
#include "ErrorHandler.h"

#include "compareFloats.h"
#include <vector>

using namespace std;
using namespace csmp;

namespace csmp
{

  template<uint32_t dim>
  class PorosityPermeabilityCouplingVisitor: public Visitor<dim>
  {

    public:
      // with_pore_volume_fields (default true): if true, the constructor
      // looks up the five pore-volume fields used by
      // CalculatePoreVolumeChangeFactor(). If false, those lookups are
      // skipped and the "undefined property" error is raised only if that
      // method is later called. Pass false when the caller only uses
      // Visit(), so the model does not need to declare fields it will not
      // use.
      PorosityPermeabilityCouplingVisitor(Model<dim> &model,
                                          bool with_pore_volume_fields = true);

      ~PorosityPermeabilityCouplingVisitor();

      virtual void Visit(Region<dim> *r);
      virtual void Visit(Element<dim> *e);
      virtual void Visit(Model<dim> *m);

      void DepthDependent( bool d_dep );
      void CouplingOption( std::string coupling_option );
      void ChangeUpperLowerPorosityLimits( double lower_limit, double upper_limit );
      void CalculatePoreVolumeChangeFactor();

    private:

      PorosityPermeabilityCouplingVisitor();

      double CalculateDepthDependentPorosityPermeabilityCoupling_Costa (Element<dim> *e);
      double CalculateDepthDependentPorosityPermeabilityCoupling_VermaPruess( Element<dim> *e );

      Model<dim> &model;

      Index
      porosity_key,
      porosity_node_key,
      permeability_key,
      bulk_volume_key,
      pore_volume_before_key, pore_volume_key,
      volume_change_factor_key;

      bool
      depth_dependent,
      with_pore_volume_fields_;

      std::string
      coupling;

      double
      phi_D, phi_min, phi_max;

      ScalarVariable
      k, phi;

  };
} // end namespace csmp

#endif // PorosityPermeabilityCouplingVisitor_H
