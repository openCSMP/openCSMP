// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PARMIGIANI_VISITOR_H
#define PARMIGIANI_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"
#include "TensorVariable.h"
#include "VectorVariable.h"
#include "compareFloats.h"

#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"

namespace csmp
{

  template<uint32_t> class PropertyDatabase;
  template<uint32_t> class Model;

  /**
  @author BLC
  @date 2016
  edited NK 2024
  */

  template<size_t dim>
  class ParmigianiVisitor : public Visitor<dim>
  {

    public:
      ParmigianiVisitor(Model<dim> &model);

      virtual ~ParmigianiVisitor();
      //void CalculateVolatileSaturation(Model<dim>* model);
      virtual void Visit(Element<dim> *element);
      virtual void Visit(Model<dim> *element);
      virtual void Visit(Region<dim> *element);

      void GetIntrinsicPermeability();
      void GetRelativePermeability();
      void GetCriticalVolatileFraction();

    private:
      std::string                group_name;  // if vis is restricted to Region
      const PropertyDatabase<dim>  &prop_ref_;
      Model< dim>                  &model_ref_;
      csmp::Index

      volatile_fraction_key_,
      crystal_fraction_key_,
      crystallinity_key_,
      melt_fraction_key_,
      porous_flag_key_,
      permeability_key_,
      vertical_permeability_key_,
      horizontal_permeability_key_,
      porosity_key_,
      fracturing_ref_key_,
      relative_permeability_key_,
      intrinsic_permeability_key_,
      volatile_fraction_critical_key_,
      channels_permeability_key_,
      permeability_ID_key_;

      ScalarVariable             fracturing_ref;
      ScalarVariable             volatile_fraction;
      ScalarVariable             crystal_fraction;
      ScalarVariable             melt_fraction;
      ScalarVariable             crystallinity;
      ScalarVariable             porous_flag;
      ScalarVariable             porous_flag_before;
      ScalarVariable             permeability, vertical_permeability, horizontal_permeability;
      ScalarVariable             porosity;
      ScalarVariable             relative_permeability;
      ScalarVariable             intrinsic_permeability;
      ScalarVariable             volatile_fraction_critical;
      ScalarVariable             channels_permeability;
      ScalarVariable             permeability_ID;

  };

} // end namespace csmp

#endif
