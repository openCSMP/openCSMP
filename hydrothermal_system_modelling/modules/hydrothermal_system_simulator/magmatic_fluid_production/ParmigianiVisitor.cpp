// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ParmigianiVisitor.h"

// Visitor meant to evaluate the permeability of a magma as a function of the crystal and fluid volume fractions
// Based on Parmigiani et al. (2017), Geochem. Geophys. Geosys. 18(8) 2887-2905

namespace csmp
{

  template<size_t dim>
  ParmigianiVisitor<dim>::ParmigianiVisitor(Model<dim> &model)
    : Visitor<dim>(MODEL, ELEMENT),

      model_ref_(model),
      prop_ref_(model.Database()),
      volatile_fraction_key_(model.Database().StorageKey("volatile fraction element")),
      volatile_fraction_critical_key_(model.Database().StorageKey("volatile fraction critical element")),
      crystal_fraction_key_(model.Database().StorageKey("crystal fraction element")),
      crystallinity_key_(model.Database().StorageKey("crystallinity element")),
      melt_fraction_key_(model.Database().StorageKey("melt fraction element")),
      porous_flag_key_(model.Database().StorageKey("porous flag element")),
      permeability_key_(model.Database().StorageKey("permeability")),
      vertical_permeability_key_(model.Database().StorageKey("vertical permeability")),
      horizontal_permeability_key_(model.Database().StorageKey("horizontal permeability")),
      porosity_key_(model.Database().StorageKey("porosity")),
      fracturing_ref_key_(model.Database().StorageKey("fracturing reference")),
      relative_permeability_key_(model.Database().StorageKey("relative permeability")),
      intrinsic_permeability_key_(model.Database().StorageKey("intrinsic permeability")),
      channels_permeability_key_(model.Database().StorageKey("channels permeability")),
      permeability_ID_key_( model.Database().StorageKey("permeability ID") )


  {} // end ParmigianiVisitor

  template<size_t dim>
  ParmigianiVisitor<dim>::~ParmigianiVisitor()

  {} // end ~ParmigianiVisitor

  template<size_t dim>
  void ParmigianiVisitor<dim>::Visit(Model<dim> *element)
  {}

  template<size_t dim>
  void ParmigianiVisitor<dim>::Visit(Region<dim> *element)
  {}

  template<size_t dim>
  void ParmigianiVisitor<dim>::Visit(Element<dim> *element)
  {
    element->Read(volatile_fraction_key_, volatile_fraction);
    element->Read(porosity_key_, porosity);
    element->Read(melt_fraction_key_, melt_fraction);
    element->Read(crystal_fraction_key_, crystal_fraction);
    element->Read(crystallinity_key_, crystallinity);
    element->Read(permeability_key_, permeability);
    element->Read(vertical_permeability_key_, vertical_permeability);
    element->Read(horizontal_permeability_key_, horizontal_permeability);
    element->Read(porous_flag_key_, porous_flag_before);

    intrinsic_permeability() = 1.0e-22;
    relative_permeability() = 0.0;
    channels_permeability() = 1.0e-22;

    porous_flag() = 0.0;
    fracturing_ref() = 0.0;

    permeability_ID() = -1;

    GetCriticalVolatileFraction();

    double scale_down_factor(1.0); // obsolete. TS
    double max_perm(1.0e-15); //Maximum allowed permeability


    //---------------------------------------------------------------------------------------------------------------------------
    //case "white" in Parmigiani et al. (2017), Fig. 9: "MVP is the carrier phase"
    if (volatile_fraction() >= 0.5)
      {
        /*        if (crystal_fraction()<0.4) //NK: fringe case within the field "MVP is the carrier phase" (Parmigiani et al. 2017, Fig. 9)
                //but why should this even be a special case?
                {
                    std::cerr<<std::endl<<"crystal_fraction()<0.4 && volatile_fraction()>0.5"<<std::endl; // "this is a fringe case" TS
                    std::cerr << " press ENTER to continue... ";
                    std::cin.get();

                    intrinsic_permeability()=1.0e-14;
                    relative_permeability()=2.0;

                    porous_flag()=1.0;
                    channels_permeability()=relative_permeability()*intrinsic_permeability()*scale_down_factor;

                    channels_permeability() = std::min(channels_permeability(),max_perm);
                    // to ensure that permeability does not exceed the maximum permeability that we allow for reasons of computational time

                    permeability()=channels_permeability();
                } */

      }

    //---------------------------------------------------------------------------------------------------------------------------
    //case "blue" in Parmigiani et al. (2017), Fig. 9: "bubbles & crystal suspension"
    else if (crystal_fraction() < 0.4 && volatile_fraction() < 0.5) // JK: is this crystallinity oder crystal fraction?
      {
        //for now there is no permeability rule that we apply for this case
      }

    //
    // else if (crystal_fraction()<0.4)
    // {
    //     // Don't do anything ... no channels have formed yet.
    //     // This is to ensure fracturing_ref() is not updated in "else".
    // }

    //---------------------------------------------------------------------------------------------------------------------------
    //case "black" in Parmigiani et al. (2017), Fig. 9: "MVP channels"
    else if (crystal_fraction() >= 0.4 && crystal_fraction() <= 0.7
             && volatile_fraction() >= volatile_fraction_critical()
             && volatile_fraction() < 0.5
             && !essentiallyEqual(crystallinity(), 1., std::numeric_limits<double>::epsilon()))
      //this is to ensure that the channels close again after full crystallisation, even though volatile fraction might be >0.3
      {
        GetIntrinsicPermeability();
        GetRelativePermeability();

        porous_flag() = 1.0;
        channels_permeability() = relative_permeability() * intrinsic_permeability() * scale_down_factor;

        channels_permeability() = std::min(channels_permeability(), max_perm);
        // to ensure that permeability does not exceed the maximum permeability that we allow for reasons of computational time

        permeability()              = channels_permeability();
        vertical_permeability()     = channels_permeability();
        horizontal_permeability()   = channels_permeability();

        // update permeability ID to indicate channels
        permeability_ID() = 6.;
        element->Store(permeability_ID_key_, permeability_ID);
      }

    //---------------------------------------------------------------------------------------------------------------------------
    //case "green" in Parmigiani et al. (2017), Fig. 9: "trapped bubbles; >potential for capillary fracturing"
    //including case where no melt is remaining (i.e. chamber is fully crystallized): This is to avoid channel opening in the "bug case"
    //where the initial volatile saturation is so high that after maximal exsolution, porosity is higher than 0.3. In this case, there
    //was a bug that crystal fraction was interpreted to be between 0.4 and 0.7
    //(because after full crystallisation: porosity + crystal fraction = 1)
    else
      {
        porous_flag() = 1.0;
        fracturing_ref() = 1.0; // Allow for fracturing

        relative_permeability() = 1.0;
        channels_permeability() = 1.0e-22;
      }

    //---------------------------------------------------------------------------------------------------------------------------

    // make sure to update fracturing reference to apply the permeability visitor once crystallinity is 1. or crystal fraction is >0.7
    if (crystallinity() >= 0.99 || crystal_fraction() >= 0.7 || (melt_fraction() <= 0.5 && volatile_fraction() >= 0.5))
      {
        fracturing_ref() = 2.0; // Allow for fracturing
      }


    // store variables
    element->Store(intrinsic_permeability_key_, intrinsic_permeability);
    element->Store(relative_permeability_key_, relative_permeability);
    element->Store(channels_permeability_key_, channels_permeability);
    element->Store(volatile_fraction_critical_key_, volatile_fraction_critical);
    element->Store(permeability_key_, permeability);
    element->Store(vertical_permeability_key_, vertical_permeability);
    element->Store(horizontal_permeability_key_, horizontal_permeability);
    element->Store(porous_flag_key_, porous_flag);
    element->Store(fracturing_ref_key_, fracturing_ref);


  } // end Visit

  template<size_t dim>
  void ParmigianiVisitor<dim>::GetIntrinsicPermeability()
  {
    intrinsic_permeability() = 1.0e-04 * (-0.0534 * pow(crystal_fraction(), 3) + 0.1083 * pow(crystal_fraction(), 2) - 0.0747 * crystal_fraction() + 0.0176);
    // Eq.2 Parmigiani et al. 2017 (10.1002/2017GC006912)
  }

  template<size_t dim>
  void ParmigianiVisitor<dim>::GetCriticalVolatileFraction()
  {
    volatile_fraction_critical() = 0.7495 * pow(crystal_fraction(), 3.) - 0.4268 * pow(crystal_fraction(), 2.) - 0.1626 * crystal_fraction() + 0.1478;
    // Eq.3.1 Degruyter et al. 2019 (NEWEST, updated from Parmigiani et al. 2017)
  }

  template<size_t dim>
  void ParmigianiVisitor<dim>::GetRelativePermeability()
  {
    relative_permeability() = -2.1778 * pow(crystal_fraction(), 4.) + 5.1511 * pow(crystal_fraction(), 3.) - 4.5199 * pow(crystal_fraction(),2.) + 1.7385 * crystal_fraction() - 0.2461;
    // Eq.3 Parmigiani et al 2017
  }


  template class ParmigianiVisitor<1U>;
  template class ParmigianiVisitor<2U>;
  template class ParmigianiVisitor<3U>;

} // end namespace csmp
