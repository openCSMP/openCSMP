#ifndef FLUIDPROPERTIES_H
#define FLUIDPROPERTIES_H

#include "States.h"
#include "CSMP_definitions.h"

/* Changelog

 13 July 2006, Thomas Driesner :
   - added dxdt, dxdp, smf
   - cp now has it's exact meaning (i.e., at constant p and x) only for single phases
   - a new member dhdt is an apparent heat capacity, i.e., the derivative at constant p and at constant phase state;
     for example, this would be the enthalpy change of liquid at twophase conditions (= non-constant composition) 
     per unit temperature
 17 August 2006, Thomas Driesner : added int "state", see definition in States.h
 2008-2011, Thomas Driesner : removed numerous properties (data members) that are not needed in standard csmp applications
 2012, Thomas Driesner : initial port to csmp++
 February 2014, Thomas Driesner : revising port to csmp++
*/


namespace csmp
{
  struct Fluidproperties 
  {
  public:
    Fluidproperties();
    ~Fluidproperties();

    Fluidproperties(const Fluidproperties&);
    Fluidproperties& operator=(const Fluidproperties&);
    double64 
      t,     // temperature [C]
      p,     // pressure    [Pa]
      x,     // mole fraction NaCl []
      wt,    // weight percent NaCl []
      smf,   // mass fraction NaCl, i.e. wt/100 []
      rho,   // density [kg/m^3]
      h,     // specific enthalpy [J/kg]
      cp,    // isobaric heat capacity [J/kg/C]
      beta,  // compressibility [Pa^-1]
      s,     // saturation []
      mf,    // mass fraction of the respective phase (i.e.: m/m_total) 
      mu;    // dynamic viscosity [Pa s]
    States state; // as defined in States.h
    void InitToZero();
    void InitToBogus();
  private:

  };

  inline void Fluidproperties::InitToZero()
  {
    t     = 0.0e0;
    p     = 0.0e0;
    x     = 0.0e0;
    wt    = 0.0e0;
    smf   = 0.0e0;
    rho   = 0.0e0;
    h     = 0.0e0;
    cp    = 0.0e0;
    beta  = 0.0e0;
    s     = 0.0e0;
    mf    = 0.0e0;
    mu    = 0.0e0;
    state = none;
    return;
  }

  inline void Fluidproperties::InitToBogus()
  {
    t     = 9.9e99;
    p     = 9.9e99;
    x     = 9.9e99;
    wt    = 9.9e99;
    smf   = 9.9e99;
    rho   = 9.9e99;
    h     = 9.9e99;
    cp    = 9.9e99;
    beta  = 9.9e99;
    s     = 9.9e99;
    mf    = 9.9e99;
    mu    = 9.9e99;
    state = none;
    return;
  }
  /**
     @struct Fluidproperties Fluidproperties.h "eos/h2o_nacl/Fluidproperties.h"

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact 
     thomas.driesner@erdw.ethz.ch

     @section motivation Motivation
     Fluidproperties is essentially a container for all kinds of fluid properties used in simulations.
   
     @section issues Known issues
     This version is legacy that is still used but the design should be improved and made more flexible. The current version contains way too many variables and cannot be constructed with fewer or more arguments. Also, converting to a struct might be better.

  */

}//csmp
#endif
// FLUIDPROPERTIES_H

