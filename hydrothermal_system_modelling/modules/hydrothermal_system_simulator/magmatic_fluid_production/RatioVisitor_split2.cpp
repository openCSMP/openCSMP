// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "RatioVisitor_split2.h"

using namespace std;

namespace csmp
{

  template<size_t dim>
  RatioVisitor_split2<dim>::RatioVisitor_split2( Model<dim> &sg, bool with_magmatic_fluid )
    : open_boundaries(false),
      dt(0.),
      with_magmatic_fluid_(with_magmatic_fluid)
  {
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    // magmatic mass
    magmatic_fluid_mass_key                     = sg.Database().StorageKey("magmatic fluid mass");
    magmatic_fluid_mass_liquid_key              = sg.Database().StorageKey("magmatic fluid mass liquid");
    magmatic_fluid_mass_vapor_key               = sg.Database().StorageKey("magmatic fluid mass vapor");
    previous_magmatic_fluid_mass_liquid_key     = sg.Database().StorageKey("previous magmatic fluid mass liquid");
    previous_magmatic_fluid_mass_vapor_key      = sg.Database().StorageKey("previous magmatic fluid mass vapor");

    // magmatic salt
    magmatic_mass_salt_key                      = sg.Database().StorageKey("magmatic mass salt");
    magmatic_salt_content_liquid_key            = sg.Database().StorageKey("magmatic salt content liquid");
    magmatic_salt_content_vapor_key             = sg.Database().StorageKey("magmatic salt content vapor");
    previous_magmatic_salt_content_liquid_key   = sg.Database().StorageKey("previous magmatic salt content liquid");
    previous_magmatic_salt_content_vapor_key    = sg.Database().StorageKey("previous magmatic salt content vapor");

  }


  template<size_t dim>
  RatioVisitor_split2<dim>::~RatioVisitor_split2()
  {}

  /** visit function for Region */
  template<size_t dim>
  void RatioVisitor_split2<dim>::Visit(Model<dim> *n)
  {
    // no calcution for the region
  }

  template<size_t dim>
  void RatioVisitor_split2<dim>::Visit(Node<dim> *n)
  {
    if (n->Status( magmatic_fluid_mass_liquid_key ) != DIRICH || open_boundaries)
      {
        ReadVariables( n );
        PerformMagmaticFluidAdvection();
        StoreVariables( n );
      }
  }

  template<size_t dim>
  void RatioVisitor_split2<dim>::ReadVariables(Node<dim> *n)
  {
    if (with_magmatic_fluid_)
      {
        // magmatic mass advection variables
        n->Read(magmatic_fluid_mass_key, magmatic_fluid_mass );
        n->Read(magmatic_fluid_mass_liquid_key, magmatic_fluid_mass_liquid );
        n->Read(magmatic_fluid_mass_vapor_key, magmatic_fluid_mass_vapor );
        n->Read(previous_magmatic_fluid_mass_liquid_key, previous_magmatic_fluid_mass_liquid );
        n->Read(previous_magmatic_fluid_mass_vapor_key, previous_magmatic_fluid_mass_vapor );

        // magmatic salt advection variables
        n->Read(magmatic_mass_salt_key, magmatic_mass_salt );
        n->Read(magmatic_salt_content_liquid_key, magmatic_salt_content_liquid );
        n->Read(magmatic_salt_content_vapor_key, magmatic_salt_content_vapor );
        n->Read(previous_magmatic_salt_content_liquid_key, previous_magmatic_salt_content_liquid );
        n->Read(previous_magmatic_salt_content_vapor_key, previous_magmatic_salt_content_vapor );
      }

  }

  template<size_t dim>
  void RatioVisitor_split2<dim>::PerformMagmaticFluidAdvection()
  {

    if (with_magmatic_fluid_)
      {
        magmatic_fluid_mass()  += (magmatic_fluid_mass_liquid() - previous_magmatic_fluid_mass_liquid());
        magmatic_fluid_mass()  += (magmatic_fluid_mass_vapor() - previous_magmatic_fluid_mass_vapor());

        magmatic_mass_salt() += (magmatic_salt_content_liquid() - previous_magmatic_salt_content_liquid());
        magmatic_mass_salt() += (magmatic_salt_content_vapor() - previous_magmatic_salt_content_vapor());
      }

  }

  template<size_t dim>
  void RatioVisitor_split2<dim>::StoreVariables(Node<dim> *n)
  {
    if (with_magmatic_fluid_)
      {
        // magmatic mass advection variables
        n->Store(magmatic_fluid_mass_key,   magmatic_fluid_mass );
        // magmatic salt advection variables
        n->Store(magmatic_mass_salt_key,   magmatic_mass_salt );
      }

  }

  template<size_t dim>
  void RatioVisitor_split2<dim>::SetOpenTopTo( bool open)
  {
    open_boundaries = open;
  }

  template<size_t dim>
  void RatioVisitor_split2<dim>::SetTimeIncrement( double time_increment )
  {
    dt = time_increment;
  }



  template class RatioVisitor_split2<1U>;
  template class RatioVisitor_split2<2U>;
  template class RatioVisitor_split2<3U>;

} // csmp
