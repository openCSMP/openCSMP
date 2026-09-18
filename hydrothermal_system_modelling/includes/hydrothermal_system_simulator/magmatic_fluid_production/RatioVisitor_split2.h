// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RatioVisitor_split2_h
#define RatioVisitor_split2_h

#include "Visitor.h"
#include "Model.h"

#include "ErrorHandler.h"

namespace csmp
{

  template<size_t dim>
  class RatioVisitor_split2 : public Visitor<dim>
  {
    public:
      RatioVisitor_split2( Model<dim> &sg, bool with_magmatic_fluid = true);
      ~RatioVisitor_split2();

      virtual void Visit(Node<dim> *n);
      virtual void Visit(Model<dim> *n);
      void SetOpenTopTo(bool open);
      void SetTimeIncrement( double time_increment );

    private:

      void ReadVariables( Node<dim> *n );
      void PerformMagmaticFluidAdvection();
      void StoreVariables( Node<dim> *n );

      bool open_boundaries;
      double dt;

      // mass
      csmp::Index        magmatic_fluid_mass_key, magmatic_fluid_mass_liquid_key,
           magmatic_fluid_mass_vapor_key, previous_magmatic_fluid_mass_liquid_key,
           previous_magmatic_fluid_mass_vapor_key;
      ScalarVariable     magmatic_fluid_mass, magmatic_fluid_mass_liquid,
                         magmatic_fluid_mass_vapor, previous_magmatic_fluid_mass_liquid,
                         previous_magmatic_fluid_mass_vapor;

      // salt
      csmp::Index         magmatic_mass_salt_key, magmatic_salt_content_liquid_key,
           magmatic_salt_content_vapor_key, previous_magmatic_salt_content_liquid_key,
           previous_magmatic_salt_content_vapor_key;
      ScalarVariable      magmatic_mass_salt, magmatic_salt_content_liquid,
                          magmatic_salt_content_vapor, previous_magmatic_salt_content_liquid,
                          previous_magmatic_salt_content_vapor;

      bool with_magmatic_fluid_;

  };

} // csmp

#endif
