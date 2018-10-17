//
//  Fluid.h
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 19/01/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FLUID_H
#define CSMP_FLUID_H

#include "Node.h"

namespace csmp {

/**
    blueprint for any specific Fluid property class to be used in the generic transport scheme.
*/
template<size_t dim, template<size_t> class USER>
class Fluid {
  public:
  
 // TODO: these keys now are part of the new composition array
    /// salinity, mass fraction from the node
    double64 MassFractionNaCl( const Node<dim>* const n ) const; // fix { return n->Read(User()->key_xSalt); }

    /// composition: mass fraction (0..1) of CO2 in the aqueous phase
    double64 XCO2_AqueousPhase( const Node<dim>* const n ) const; // fix { return n->Read(User()->key_xCO2); }

    /// composition: mass fraction (0..1) of water in the aqueous phase
    double64 XH2O_AqueousPhase( const Node<dim>* const n) const; // fix { return n->Read(User()->key_xH2O); }

    /// composition: mass fraction (0..1) of CO2 in the carbonic phase
    double64 YCO2_CarbonicPhase( const Node<dim>* const n ) const; // fix { return n->Read(User()->key_YCO2); }

    /// composition: mass fraction (0..1) of water in the carbonic phase
    double64 YH2O_CarbonicPhase( const Node<dim>* const n ) const; // fix { return n->Read(User()->key_YH2O); }
  
    /// the node property fluid viscosity returned has been interpolated to the user-specified target placement (argument parameter)
    double64 Viscosity( const Node<dim>* const n, size_t phase ) const;
    double64 Viscosity( const Element<dim>* const n, size_t phase ) const;

    /// the returned node property fluid density (of phase) has been interpolated to the user-specified target placement (argument parameter)
    double64 Density( const Node<dim>* const n, size_t phase ) const;
    double64 Density( const Element<dim>* const n, size_t phase ) const;

    /// returns saturation-weighted density of the fluid mixture interpolated to the target placement
    double64 DensityMixture( const Node<dim>* const n, double64 salinity=0. ) const;
  
    /// returns the ratio of the phase viscosities at the target placement
    double64 ViscosityRatio( const Node<dim>* const n, double64 salinity=0. ) const;

  protected:
    Fluid() {}
    Fluid(bool PTX_dependent) {PTX_dependent_ = PTX_dependent;}
    
    bool PTX_dependent_ = false;
  
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    // member is EOS module
};

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

} // end csmp

#endif /* CSMP_FLUID_H */
