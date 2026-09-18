// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RatioVisitor_split1_h
#define RatioVisitor_split1_h

#include "Visitor.h"
#include "Model.h"
#include "compareFloats.h"
#include "ErrorHandler.h"

namespace csmp
{

  template<size_t dim>
  class RatioVisitor_split1 : public Visitor<dim>
  {
    public:
      RatioVisitor_split1( Model<dim> &sg, bool with_magmatic_fluid = true);
      ~RatioVisitor_split1();

      virtual void Visit(Node<dim> *n);
      virtual void Visit(Model<dim> *n);
      void SetOpenTopTo(bool open);
      void SetTimeIncrement( double time_increment );

    private:

      void ReadVariables( Node<dim> *n );
      void CalculateMagmaticRatio();
      void CalculateMagmaticSaltRatio( );
      void CalculateMagmaticTransportVariables();
      void StoreVariables( Node<dim> *n );
      void AddBoundaryTerms( );
      void SubtractBoundaryTerms( );

      bool open_boundaries;
      double dt, new_mass;

      // mass
      csmp::Index         magmatic_ratio_key, magmatic_water_ratio_key, fluid_mass_total_key,
           fluid_mass_liquid_key, fluid_mass_vapor_key, liquid_mass_mobility_key,
           vapor_mass_mobility_key,
           magmatic_fluid_mass_key, magmatic_fluid_mass_liquid_key, magmatic_fluid_mass_vapor_key,
           magmatic_liquid_mass_mobility_key, magmatic_vapor_mass_mobility_key;
      ScalarVariable      fluid_mass_total, fluid_mass_liquid, fluid_mass_vapor,
                          liquid_mass_mobility, vapor_mass_mobility, magmatic_fluid_mass,
                          magmatic_fluid_mass_liquid, magmatic_fluid_mass_vapor, magmatic_liquid_mass_mobility,
                          magmatic_ratio, magmatic_water_ratio, magmatic_vapor_mass_mobility;

      // salt
      csmp::Index         previous_mass_salt_key, salt_content_liquid_key,
           salt_content_vapor_key, salt_content_halite_key, liquid_salt_mobility_key,
           vapor_salt_mobility_key,
           magmatic_mass_salt_key, magmatic_salt_content_liquid_key, magmatic_salt_content_vapor_key,
           magmatic_salt_content_halite_key, magmatic_liquid_salt_mobility_key,
           magmatic_vapor_salt_mobility_key,
           magmatic_salt_ratio_key;
      ScalarVariable      previous_mass_salt, salt_content_liquid, salt_content_vapor,
                          salt_content_halite, liquid_salt_mobility, vapor_salt_mobility,
                          magmatic_mass_salt, magmatic_salt_content_liquid, magmatic_salt_content_vapor,
                          magmatic_salt_content_halite, magmatic_liquid_salt_mobility, magmatic_vapor_salt_mobility,
                          magmatic_salt_ratio;

      // open top
      csmp::Index         boundary_flow_mass_key, boundary_flow_salt_key, pore_volume_key;
      ScalarVariable boundary_flow_mass, boundary_flow_salt, pore_volume;

      bool with_magmatic_fluid_;

  };

} // csmp

#endif
