// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "RatioVisitor_split1.h"

using namespace std;

namespace csmp
{

  template<size_t dim>
  RatioVisitor_split1<dim>::RatioVisitor_split1( Model<dim> &sg, bool with_magmatic_fluid )
    : open_boundaries(false),
      dt(0.),
      with_magmatic_fluid_(with_magmatic_fluid)
  {
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    // mass
    magmatic_ratio_key                  = sg.Database().StorageKey("magmatic ratio");
    magmatic_water_ratio_key            = sg.Database().StorageKey("magmatic water ratio");
    fluid_mass_total_key                = sg.Database().StorageKey("fluid density");
    fluid_mass_liquid_key               = sg.Database().StorageKey("fluid mass liquid");
    fluid_mass_vapor_key                = sg.Database().StorageKey("fluid mass vapor");
    liquid_mass_mobility_key            = sg.Database().StorageKey("liquid mass mobility");
    vapor_mass_mobility_key             = sg.Database().StorageKey("vapor mass mobility");
    magmatic_fluid_mass_key             = sg.Database().StorageKey("magmatic fluid mass");
    magmatic_fluid_mass_liquid_key      = sg.Database().StorageKey("magmatic fluid mass liquid");
    magmatic_fluid_mass_vapor_key       = sg.Database().StorageKey("magmatic fluid mass vapor");
    magmatic_liquid_mass_mobility_key   = sg.Database().StorageKey("magmatic liquid mass mobility");
    magmatic_vapor_mass_mobility_key    = sg.Database().StorageKey("magmatic vapor mass mobility");

    // salt
    magmatic_salt_ratio_key             = sg.Database().StorageKey("magmatic salt ratio");
    previous_mass_salt_key              = sg.Database().StorageKey("previous mass salt");
    salt_content_liquid_key             = sg.Database().StorageKey("salt content liquid");
    salt_content_vapor_key              = sg.Database().StorageKey("salt content vapor");
    salt_content_halite_key             = sg.Database().StorageKey("salt content halite");
    liquid_salt_mobility_key            = sg.Database().StorageKey("liquid salt mobility");
    vapor_salt_mobility_key             = sg.Database().StorageKey("vapor salt mobility");
    magmatic_mass_salt_key              = sg.Database().StorageKey("magmatic mass salt");
    magmatic_salt_content_liquid_key    = sg.Database().StorageKey("magmatic salt content liquid");
    magmatic_salt_content_vapor_key     = sg.Database().StorageKey("magmatic salt content vapor");
    magmatic_salt_content_halite_key    = sg.Database().StorageKey("magmatic salt content halite");
    magmatic_liquid_salt_mobility_key   = sg.Database().StorageKey("magmatic liquid salt mobility");
    magmatic_vapor_salt_mobility_key    = sg.Database().StorageKey("magmatic vapor salt mobility");

    // boundary variables
    boundary_flow_mass_key              = sg.Database().StorageKey("boundary flow mass");
    boundary_flow_salt_key              = sg.Database().StorageKey("boundary flow salt");

    pore_volume_key                     = sg.Database().StorageKey("pore volume");
  }


  template<size_t dim>
  RatioVisitor_split1<dim>::~RatioVisitor_split1()
  {}

  /** visit function for Region */
  template<size_t dim>
  void RatioVisitor_split1<dim>::Visit(Model<dim> *n)
  {
    // no calcution for the region
  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::Visit(Node<dim> *n)
  {
    ReadVariables( n );
    AddBoundaryTerms();

    if (with_magmatic_fluid_)
{
      CalculateMagmaticRatio();
      CalculateMagmaticSaltRatio( );
}

    SubtractBoundaryTerms(); // NK commented out on Nov 30 2023 // why?
    CalculateMagmaticTransportVariables();
    StoreVariables( n );

}


  template<size_t dim>
  void RatioVisitor_split1<dim>::ReadVariables(Node<dim> *n)
  {
    // mass advection variables
    n->Read(fluid_mass_total_key, fluid_mass_total ); //
    n->Read(fluid_mass_liquid_key, fluid_mass_liquid );
    n->Read(fluid_mass_vapor_key, fluid_mass_vapor );
    n->Read(liquid_mass_mobility_key, liquid_mass_mobility );
    n->Read(vapor_mass_mobility_key, vapor_mass_mobility );

    // salt advection variables
    n->Read(previous_mass_salt_key, previous_mass_salt );
    n->Read(salt_content_liquid_key, salt_content_liquid );
    n->Read(salt_content_vapor_key, salt_content_vapor );
    n->Read(salt_content_halite_key, salt_content_halite );
    n->Read(liquid_salt_mobility_key, liquid_salt_mobility );
    n->Read(vapor_salt_mobility_key, vapor_salt_mobility );

    if (with_magmatic_fluid_)
      {
        // ratios
        n->Read(magmatic_ratio_key, magmatic_ratio );
        n->Read(magmatic_water_ratio_key, magmatic_water_ratio);
        n->Read(magmatic_salt_ratio_key, magmatic_salt_ratio);

        // magmatic mass advection variables
        n->Read(magmatic_fluid_mass_key, magmatic_fluid_mass );

        // magmatic salt advection variables
        n->Read(magmatic_mass_salt_key, magmatic_mass_salt );
        n->Read(magmatic_salt_content_liquid_key, magmatic_salt_content_liquid );
        n->Read(magmatic_salt_content_vapor_key, magmatic_salt_content_vapor );
        n->Read(magmatic_salt_content_halite_key, magmatic_salt_content_halite );
      }




    // boundary flow
    n->Read(boundary_flow_mass_key, boundary_flow_mass );
    n->Read(boundary_flow_salt_key, boundary_flow_salt );

    n->Read(pore_volume_key, pore_volume );
  }


  template<size_t dim>
  void RatioVisitor_split1<dim>::CalculateMagmaticSaltRatio( )
  {

    extern ErrorHandler   skm_err;

    if (essentiallyEqual(magmatic_mass_salt(), 0., numeric_limits<double>::epsilon()))
      {
        magmatic_salt_ratio() = 0.0;
        magmatic_mass_salt() = 0.0;
      }
    else if (essentiallyEqual(magmatic_mass_salt(), previous_mass_salt(), numeric_limits<double>::epsilon()))
      {
        magmatic_salt_ratio() = 1.0;
        magmatic_mass_salt() = previous_mass_salt();
      }
    else if (definitelyGreaterThan( magmatic_mass_salt(), previous_mass_salt(), numeric_limits<double>::epsilon() ))
      {
        magmatic_mass_salt() = previous_mass_salt();
        //      throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticSaltRatio",
        //                     "magmatic mass salt is greater than total mass salt." );
        //      int stop;
        //      cout << "Node at " << n->x() << ", " << n->y() << endl;
        //      cout << "magmatic_mass_salt: " << magmatic_mass_salt << ", previous_mass_salt: " << previous_mass_salt << endl;
        //      cin >> stop;
      }
    else if (definitelyLessThan( magmatic_mass_salt(), 0.0, numeric_limits<double>::epsilon() ))
      {
        magmatic_mass_salt() = 0.0;
        //      throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticSaltRatio",
        //                     "magmatic mass salt is less than zero." );
      }
    else if (definitelyGreaterThan( previous_mass_salt(), 0.0, numeric_limits<double>::epsilon() ))
      {
        magmatic_salt_ratio() = magmatic_mass_salt() / previous_mass_salt();
      }
    else if (definitelyLessThan( previous_mass_salt(), 0.0, numeric_limits<double>::epsilon() ))
      {
        throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticSaltRatio",
                        "mass salt is less than zero." );
      }
    else
      {
        int temp;
        cerr << "\nmagmatic water ratio: magmatic_water_ratio " << magmatic_water_ratio() << endl;
        cerr << "\nmagmatic fluid mass: magmatic_fluid_mass " << magmatic_fluid_mass() << endl;
        cerr << "\nfluid density: fluid_mass_total " << fluid_mass_total() << endl;
        cerr << "\nmagmatic salt ratio: magmatic_salt_ratio " << magmatic_salt_ratio() << endl;
        cerr << "\nmagmatic salt content liquid: magmatic_salt_content_liquid " <<
             magmatic_salt_content_liquid() << endl;
        cerr << "\nmagmatic salt content vapor: magmatic_salt_content_vapor " <<
             magmatic_salt_content_vapor() << endl;
        cerr << "\nmagmatic mass salt: magmatic_mass_salt " << magmatic_mass_salt() << endl;
        cerr << "\nprevious mass salt: previous_mass_salt " << previous_mass_salt() << endl;
        cin >> temp;
        throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticSaltRatio",
                        "magmatic salt ratio could not be calculated." );
      }

  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::CalculateMagmaticRatio()
  {

    extern ErrorHandler   skm_err;

    if (essentiallyEqual(magmatic_fluid_mass(), 0., numeric_limits<double>::epsilon()))
      {
        magmatic_ratio() = 0.0;
        magmatic_water_ratio() = 0.0;
        magmatic_fluid_mass() = 0.0;
      }
    else if (essentiallyEqual(magmatic_fluid_mass(), fluid_mass_total(), numeric_limits<double>::epsilon()))
      {
        magmatic_ratio() = 1.0;
        magmatic_water_ratio() = 1.0;
        magmatic_fluid_mass() = fluid_mass_total();
      }
    else if (definitelyGreaterThan( magmatic_fluid_mass(), fluid_mass_total(), numeric_limits<double>::epsilon() ))
      {
        magmatic_fluid_mass() = fluid_mass_total();
      }
    else if (definitelyLessThan( magmatic_fluid_mass(), 0.0, numeric_limits<double>::epsilon() ))
      {
        magmatic_fluid_mass() = 0.0;
      }
    else if (definitelyGreaterThan( fluid_mass_total(), 0.0, numeric_limits<double>::epsilon() ))
      {
        // ratio of bulk fluid
        magmatic_ratio() = magmatic_fluid_mass() / fluid_mass_total();

        // ratio of H2O component
        if (definitelyGreaterThan(fluid_mass_total(), previous_mass_salt(), numeric_limits<double>::epsilon()))
          {
            magmatic_water_ratio() = (magmatic_fluid_mass() - magmatic_mass_salt()) / (fluid_mass_total() - previous_mass_salt());
            magmatic_water_ratio() = std::min(1.0, magmatic_water_ratio());
            magmatic_water_ratio() = std::max(0.0, magmatic_water_ratio());
          }
        else
          magmatic_water_ratio() = 0.0;

      }
    else if (essentiallyEqual(fluid_mass_total(), 0.0, numeric_limits<double>::epsilon()))
      {
        throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticRatio",
                        "fluid density is zero." );
      }
    else if (definitelyLessThan( fluid_mass_total(), 0.0, numeric_limits<double>::epsilon() ))
      {
        throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticRatio",
                        "fluid density is less than zero." );
      }
    else
      {
        int temp;
        cerr << "\nmagmatic water ratio: " << magmatic_water_ratio() << endl;
        cerr << "\nmagmatic fluid mass: " << magmatic_fluid_mass() << endl;
        cerr << "\nfluid density: " << fluid_mass_total() << endl;
        cerr << "\nmagmatic salt ratio: " << magmatic_salt_ratio() << endl;
        cerr << "\nmagmatic salt content liquid: " << magmatic_salt_content_liquid() << endl;
        cerr << "\nmagmatic salt content vapor: " << magmatic_salt_content_vapor() << endl;
        cerr << "\nmagmatic mass salt: " << magmatic_mass_salt() << endl;
        cerr << "\nprevious mass salt: " << previous_mass_salt() << endl;
        cin >> temp;
        throw Exception(FATAL_ERROR, "RatioVisitor_split1::CalculateMagmaticRatio",
                        "magmatic ratio could not be calculated." );
      }

  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::CalculateMagmaticTransportVariables()
  {
    if (with_magmatic_fluid_)
      {
        // salt
        magmatic_salt_content_liquid() = salt_content_liquid() * magmatic_salt_ratio();
        magmatic_salt_content_vapor() = salt_content_vapor() * magmatic_salt_ratio();
        magmatic_salt_content_halite() = salt_content_halite() * magmatic_salt_ratio();
        magmatic_liquid_salt_mobility() = liquid_salt_mobility() * magmatic_salt_ratio();
        magmatic_vapor_salt_mobility() = vapor_salt_mobility() * magmatic_salt_ratio();

        // mass
        magmatic_fluid_mass_liquid() = ( fluid_mass_liquid() - salt_content_liquid() ) * magmatic_water_ratio() + magmatic_salt_content_liquid();
        magmatic_fluid_mass_vapor() = ( fluid_mass_vapor() - salt_content_vapor() ) * magmatic_water_ratio() + magmatic_salt_content_vapor();

        if (fluid_mass_liquid() > 0.0)
          magmatic_liquid_mass_mobility() = liquid_mass_mobility() * magmatic_fluid_mass_liquid() / fluid_mass_liquid();
        else
          magmatic_liquid_mass_mobility() = 0.0;

        if (fluid_mass_vapor() > 0.0)
          magmatic_vapor_mass_mobility() = vapor_mass_mobility() * magmatic_fluid_mass_vapor() / fluid_mass_vapor();
        else
          magmatic_vapor_mass_mobility() = 0.0;
      }

  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::StoreVariables(Node<dim> *n)
  {

    if (with_magmatic_fluid_)
      {
        // ratios
        n->Store(magmatic_ratio_key, magmatic_ratio );
        n->Store(magmatic_water_ratio_key, magmatic_water_ratio);
        n->Store(magmatic_salt_ratio_key, magmatic_salt_ratio);

        // magmatic mass advection variables
        n->Store(magmatic_fluid_mass_key,   magmatic_fluid_mass );
        n->Store(magmatic_fluid_mass_liquid_key,  magmatic_fluid_mass_liquid );
        n->Store(magmatic_fluid_mass_vapor_key,  magmatic_fluid_mass_vapor );
        n->Store(magmatic_liquid_mass_mobility_key, magmatic_liquid_mass_mobility );
        n->Store(magmatic_vapor_mass_mobility_key, magmatic_vapor_mass_mobility );

        // magmatic salt advection variables
        n->Store(magmatic_mass_salt_key,   magmatic_mass_salt );
        n->Store(magmatic_salt_content_liquid_key,  magmatic_salt_content_liquid );
        n->Store(magmatic_salt_content_vapor_key,  magmatic_salt_content_vapor );
        n->Store(magmatic_salt_content_halite_key,  magmatic_salt_content_halite );
        n->Store(magmatic_liquid_salt_mobility_key,  magmatic_liquid_salt_mobility );
        n->Store(magmatic_vapor_salt_mobility_key,  magmatic_vapor_salt_mobility );
      }



  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::AddBoundaryTerms()
  {
    fluid_mass_total() += std::max(0.0, boundary_flow_mass());
    previous_mass_salt() += std::max(0.0, boundary_flow_salt());
  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::SubtractBoundaryTerms()
  {
    if (with_magmatic_fluid_)
      {
        magmatic_fluid_mass() -= std::max(0.0,
                                          (boundary_flow_mass() - boundary_flow_salt()) * magmatic_water_ratio() +
                                          boundary_flow_salt() * magmatic_salt_ratio());
        magmatic_mass_salt() -= std::max(0.0, boundary_flow_salt() * magmatic_salt_ratio());
      }

  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::SetOpenTopTo( bool open)
  {
    open_boundaries = open;
  }

  template<size_t dim>
  void RatioVisitor_split1<dim>::SetTimeIncrement( double time_increment )
  {
    dt = time_increment;
  }

  template class RatioVisitor_split1<1U>;
  template class RatioVisitor_split1<2U>;
  template class RatioVisitor_split1<3U>;

} // csmp
