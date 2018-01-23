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
    double64 Temperature() const { return User()->T; }
  
    /// fluid pressure (Pa) at current initialisation point in Element
    double64 Pressure() const { return User()->P; }
  
    /// salinity, mSalt (molality = moles/kg)
    double64 Salinity() const { return User()->salinity; }
  
    /// composition: mass fraction (0..1) of CO2 in the aqueous phase
    double64 XCO2_AqueousPhase() const { return User()->XCO2; }

    /// composition: mass fraction (0..1) of water in the aqueous phase
    double64 XH2O_AqueousPhase() const { return User()->XH2O; }

    /// composition: mass fraction (0..1) of CO2 in the carbonic phase
    double64 YCO2_CarbonicPhase() const { return User()->YCO2; }

    /// composition: mass fraction (0..1) of water in the carbonic phase
    double64 YH2O_CarbonicPhase() const { return User()->YH2O; }
  
    double64 Viscosity( double64 pf, double64 T, double64 salinity=0., size_t phase=0 ) const;
    double64 Density( double64 pf, double64 T, double64 salinity=0., size_t phase=0 ) const;
    double64 DensityMixture( double64 pf, double64 T, double64 sw, double64 salinity=0. ) const;

    // versions that account for dissolved CO2
    double64 Viscosity( double64 pf, double64 T, double64 salinity, double64 XCO2, double64 YH2O, size_t phase=0 ) const;
    double64 Density( double64 pf, double64 T, double64 salinity, double64 XCO2, double64 YH2O, size_t phase=0 ) const;
    double64 DensityMixture( double64 pf, double64 T, double64 sw, double64 salinity, double64 XCO2, double64 YH2O ) const;
  
    double64 ViscosityRatio( double64 pf, double64 T, double64 msalt=0. ) const;

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
