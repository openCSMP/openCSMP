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
#include "PVTX_Calculator_H2O_CO2_NaCl.h"

namespace csmp {

/**
    Blueprint for any specific Fluid property class to be used in the generic transport scheme.
    This generic one is not associated with any particular equation of state, but derived ones
    in this implementation of static polymorphism may well be.
 
    When properties are requested for the element, the methods of the class performs
    a smart interpolation of fluid properties dependent on their existance
    as indicated by the phase state.
    Else,the method just reads the properties from the model.
 
    @attention the fluid properties must be initialised elsewhere by an equation of state.
*/
template<size_t dim, template<size_t> class USER>
class Fluid {
  public:
    /// default constructor that tests the phase-state key
    Fluid();

    /// the node property fluid viscosity (Pa.s) returned has been interpolated to the user-specified target placement (argument parameter)
    double Viscosity( Node<dim>* const, size_t phase ) const;
    /// the node property fluid viscosity (Pa.s) interpolated to element barycentre
    double Viscosity( Element<dim>* const, size_t phase ) const;
    /// the node property fluid viscosity (Pa.s) at the node i of the element
    double Viscosity( Element<dim>* const, size_t node, size_t phase ) const;

    /// the returned node property fluid density (of phase) (kg/m3) interpolated to the user-specified target placement (argument parameter)
    double Density( Node<dim>* const, size_t phase ) const;
    double Density( Element<dim>* const, size_t phase ) const;
    double Density( Element<dim>* const, size_t node, size_t phase ) const;

    /// returns saturation-weighted density average (kg/m3) for the fluid mixture; properties are interpolated to the target placement
    double MixtureDensity( Node<dim>* const ) const;
    double MixtureDensity( Element<dim>* const ) const;
    double MixtureDensity( Element<dim>* const, size_t node ) const;

    /// returns the ratio of the phase viscosities at the target placement
    double ViscosityRatio( Node<dim>* const ) const;
    double ViscosityRatio( Element<dim>* const ) const;
    double ViscosityRatio( Element<dim>* const, size_t node ) const;

  protected:
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
};

//  unit conversions

inline double  molalNaClToMassFracNaClInAqueousPhase( double mSalt); // no CO2

inline double  massFracNaClToMolalNaClInAqueousPhase( double massFracSalt); // no CO2

inline double  massFracNaClToMolarFracNaClInAqueousPhase( double massFracSalt);// no CO2

inline double  molalNaClToMolarFracNaClInAqueousPhase( double mSalt);// no CO2

inline double  ppmNaClToMolalNaClInAqueousPhase( double ppmSalt );

inline double  molalNaClToPpmInAqueousPhase( double mSalt );

inline double  psiToPa( double pressureInPsi );

inline double  paToPsi( double pressureInPa );

inline double  paTobar( double pressureInPa );

inline double  barTopa( double pressureInbar );

inline double  degreeCToKelvin( double temperatureInC );

inline double  KelvinTodegreeC( double temperatureInK );

} // end csmp

#endif /* CSMP_FLUID_H */
