//
//  Fluid.h
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 19/01/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLUID_H
#define CSMP_FLUID_H

#include "CSMP_definitions.h"
#include "VariablePlacement.h"

namespace csmp {

//  unit conversions

    double64  molalNaClToMassFracNaClInAqueousPhase( double64 mSalt); // no CO2

    double64  massFracNaClToMolalNaClInAqueousPhase( double64 massFracSalt); // no CO2

    double64  massFracNaClToMolarFracNaClInAqueousPhase( double64 massFracSalt);// no CO2

    double64  molalNaClToMolarFracNaClInAqueousPhase( double64 mSalt);// no CO2

    double64  ppmNaClToMolalNaClInAqueousPhase( double64 ppmSalt );

    double64  molalNaClToPpmInAqueousPhase( double64 mSalt );

    double64  psiToPa( double64 pressureInPsi );

    double64  paToPsi( double64 pressureInPa );

    double64  paTobar( double64 pressureInPa );

    double64  barTopa( double64 pressureInbar );

    double64  degreeCToKelvin( double64 temperatureInC );

    double64  KelvinTodegreeC( double64 temperatureInK );

/**
    blueprint for any specific Fluid property class to be used in the generic transport scheme.
*/
template<size_t dim, template<size_t> class USER>
class Fluid {
  public:
    /// temperature oC at current initialisation point in Element
    template<class TARGET_PLACEMENT>
    double64 Temperature( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_T); }
  
    /// fluid pressure (Pa) at current initialisation point in Element
    template<class TARGET_PLACEMENT>
    double64 Pressure( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_pf); }
  
    /// salinity, mSalt (molality = moles/kg)
    template<class TARGET_PLACEMENT>
    double64 Salinity( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_msalt); }
  
    /// salinity, mass fraction
    template<class TARGET_PLACEMENT>
    double64 MassFractionNaCl( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_xSalt); }

    /// composition: mass fraction (0..1) of CO2 in the aqueous phase
    template<class TARGET_PLACEMENT>
    double64 XCO2_AqueousPhase( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_xCO2); }

    /// composition: mass fraction (0..1) of water in the aqueous phase
    template<class TARGET_PLACEMENT>
    double64 XH2O_AqueousPhase( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_xH2O); }

    /// composition: mass fraction (0..1) of CO2 in the carbonic phase
    template<class TARGET_PLACEMENT>
    double64 YCO2_CarbonicPhase( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_YCO2); }

    /// composition: mass fraction (0..1) of water in the carbonic phase
    template<class TARGET_PLACEMENT>
    double64 YH2O_CarbonicPhase( TARGET_PLACEMENT& p ) const { return p.Interpolate(User()->key_YH2O); }
  
    // with or without dissolved CO2
    template<class TARGET_PLACEMENT>
    double64 Viscosity( TARGET_PLACEMENT&, size_t phase=0 ) const;

    template<class TARGET_PLACEMENT>
    double64 Density( TARGET_PLACEMENT&, size_t phase=0 ) const;

    template<class TARGET_PLACEMENT>
    double64 DensityMixture( TARGET_PLACEMENT&, double64 salinity=0. ) const;
  
    template<class TARGET_PLACEMENT>
    double64 ViscosityRatio( TARGET_PLACEMENT&, double64 salinity=0. ) const;

  protected:
    Fluid() = delete;
    ~Fluid() = delete;
  
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    // member is EOS module
};

}

#endif /* CSMP_FLUID_H */
