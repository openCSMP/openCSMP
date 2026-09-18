// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef WATER_H
#define WATER_H

#include "CSMP_definitions.h"

#include "steam4.h"

#include "CriticalPointH2O.h"

/* Changelog
   
   October 2014, Thomas Driesner: port to CSMP++, changed interface to comply with Pa as pressure unit, added doxygen documentation, adopted CSMP++ style guide
  
*/
namespace csmp
{
  /// Water properties as a function of T and P, provides interface to PROST implementation of IAPS84 equation of state
  class Water
  {

  public:
    Water( const double& externaltemperature,  // [C] 
           const double& externalpressure );   // [Pa]
    ~Water();

    double Pressure()        const ;
    double PressurePascal()  const ;
    double Temperature()     const ;
    double Density()               ;
    double MolarVolume()           ;
    double Enthalpy()              ;
    double FreeEnergy()            ;
    double HeatCapacity()          ;
    double Dvdppascal()            ;
    double Compressibilitypascal() ;
    double Dvdpbar()               ;
    double Dvdt()                  ;
    double Compressibilitybar()    ;
    double MinusDvdppascal()       ;
    double MinusDvdpbar()          ;

    double MolarVolumeAtPPlusDelP();
    double DensityAtPPlusDelP();

    double EnthalpyAtTPlusDelT()   ;
    double DensityAtTPlusDelT()    ;
    double LiquidSaturationDensityAtTPlusDelT() ;
    double VaporSaturationDensityAtTPlusDelT()  ;
    double SaturationPressureAtTPlusDelT()      ;

    double EnthalpyAtTMinusDelT()   ;
    double DensityAtTMinusDelT()    ;
    double LiquidSaturationDensityAtTMinusDelT() ;
    double VaporSaturationDensityAtTMinusDelT()  ;
    double SaturationPressureAtTMinusDelT()      ;

    double LiquidSaturationDensityForP()      ;
    double LiquidSaturationMolarVolumeForP()  ;
    double LiquidSaturationEnthalpyForP()     ;
    double LiquidSaturationHeatCapacityForP() ;
    double VaporSaturationDensityForP()       ;
    double VaporSaturationMolarVolumeForP()   ;
    double VaporSaturationEnthalpyForP()      ;
    double VaporSaturationHeatCapacityForP()  ;
    double LiquidSaturationDensityForT()      ;
    double LiquidSaturationMolarVolumeForT()  ;
    double LiquidSaturationEnthalpyForT()     ;
    double LiquidSaturationHeatCapacityForT() ;
    double VaporSaturationDensityForT()       ;
    double VaporSaturationMolarVolumeForT()   ;
    double VaporSaturationEnthalpyForT()      ;
    double VaporSaturationHeatCapacityForT()  ;
    double CriticalDensity()                  ;
    double CriticalMolarVolume()              ;
    double CriticalEnthalpy()                 ;
    double CriticalHeatCapacity()             ;
    double SaturationPressure()               ;
    double SaturationTemperature()            ;
    double SaturationTemperatureKelvin()      ;
    double Viscosity()                        ;
    double ViscosityFromTandP(const double& t, const double& p);
    double LiquidViscosityFromT(const double& t);
    
  private:
    
    Water();

    const double&    temperature;        ///< reference to temperature in flow code
    const double&    pressure;           ///< reference to pressure in flow code

    const double     kelvin;             ///< for C -> K conversion
    const double     mol_h2o;            ///< mole / mass conversion factor

    double           tcurrent;           ///< internal temperature variable [C]
    double           pcurrent;           ///< internal pressure variable [bar]
    double           d;                  ///< dummy variable, for density estimate when calling PROST functions
    double           dp;                 ///< dummy variable, for density tolerance when calling PROST functions
    double           dvdppascal;         ///< dv/dP in [(cm^3/kg)/pa]
    double           comprpascal;        ///< compressibility [Pa-1]
    double           delp;               ///< pressure delta [bar] for numerical evaluation of derivatives
    double           delt;               ///< temperature delta [C] for numerical evaluation of derivatives
    double           dvdt;               ///< dv/dT in [(cm^3/kg)/K]
    double           viscosity;          ///< dynamic viscosity [Pa s]
    double           ak[4];              ///< coefficients of viscosity formula
    double           bij[6][5];          ///< more coefficients of viscosity formula

    CriticalPointH2O   cp_h2o;

    Prop              *properties;         ///< PROST water properties, general
    Prop              *delPproperties;     ///< PROST water properties, at P+delp
    Prop              *delTproperties;     ///< PROST water properties, at T+delt
    Prop              *minTproperties;     ///< PROST water properties, at T-delt 
    Prop              *liqpropst;          ///< PROST water properties, liquid on saturation curve, for given T
    Prop              *vappropst;          ///< PROST water properties, vapor on saturation curve, for given T
    Prop              *liqpropstplusdelt;  ///< PROST water properties, liquid on saturation curve, at T+delt
    Prop              *vappropstplusdelt;  ///< PROST water properties, vapor on saturation curve, at T+delt
    Prop              *liqpropstminusdelt; ///< PROST water properties, liquid on saturation curve, at T-delt
    Prop              *vappropstminusdelt; ///< PROST water properties, vapor on saturation curve, at T-delt
    Prop              *liqpropsp;          ///< PROST water properties, liquid on saturation curve, for given P
    Prop              *vappropsp;          ///< PROST water properties, vapor on saturation curve, for given P

    void               CheckStatus();
    double           Viscosity(const double& T, const double& rho);

    double           T() const;
    double           P() const;

  };

  inline double Water::Pressure()                   const { return pcurrent; }
  inline double Water::PressurePascal()             const { return pcurrent*1.0e5; }
  inline double Water::Temperature()                const { return tcurrent; }

  inline double Water::T()                          const { return temperature+273.15e0; }
  inline double Water::P()                          const { return pressure*1.0e5; }

  inline double Water::Density()                          { CheckStatus(); return properties->d; }
  inline double Water::MolarVolume()                      { CheckStatus(); return mol_h2o/properties->d; }
  inline double Water::Enthalpy()                         { CheckStatus(); return properties->h; }
  inline double Water::FreeEnergy()                       { CheckStatus(); return properties->g; }
  inline double Water::HeatCapacity()                     { CheckStatus(); return properties->cp; }
  inline double Water::Dvdppascal()                       { CheckStatus(); return dvdppascal; }
  inline double Water::Compressibilitypascal()            { CheckStatus(); return comprpascal; }
  inline double Water::Dvdpbar()                          { CheckStatus(); return dvdppascal*1.0e5; }
  inline double Water::Compressibilitybar()               { CheckStatus(); return comprpascal*1.0e5; }
  inline double Water::MinusDvdppascal()                  { CheckStatus(); return -dvdppascal; }
  inline double Water::MinusDvdpbar()                     { CheckStatus(); return -dvdppascal*1.0e5; }
  inline double Water::Dvdt()                             { CheckStatus(); return dvdt; }

  inline double Water::LiquidSaturationDensityForP()      { CheckStatus(); return liqpropsp->d; }
  inline double Water::LiquidSaturationMolarVolumeForP()  { CheckStatus(); return mol_h2o/liqpropsp->d; }
  inline double Water::LiquidSaturationEnthalpyForP()     { CheckStatus(); return liqpropsp->h; }
  inline double Water::LiquidSaturationHeatCapacityForP() { CheckStatus(); return liqpropsp->cp; }

  inline double Water::VaporSaturationDensityForP()       { CheckStatus(); return vappropsp->d; }
  inline double Water::VaporSaturationMolarVolumeForP()   { CheckStatus(); return mol_h2o/vappropsp->d; }
  inline double Water::VaporSaturationEnthalpyForP()      { CheckStatus(); return vappropsp->h; }
  inline double Water::VaporSaturationHeatCapacityForP()  { CheckStatus(); return vappropsp->cp; }

  inline double Water::LiquidSaturationDensityForT()      { CheckStatus(); return liqpropst->d; }
  inline double Water::LiquidSaturationMolarVolumeForT()  { CheckStatus(); return mol_h2o/liqpropst->d; }
  inline double Water::LiquidSaturationEnthalpyForT()     { CheckStatus(); return liqpropst->h; }
  inline double Water::LiquidSaturationHeatCapacityForT() { CheckStatus(); return liqpropst->cp; }
  inline double Water::LiquidViscosityFromT(const double& t)
  {
    double t_backup=tcurrent;
    tcurrent = t;
    CheckStatus(); 
    double visc = Viscosity(t+273.15,liqpropst->d);
    tcurrent = t_backup;
    CheckStatus();
    return visc;
  }

  inline double Water::VaporSaturationDensityForT()         { CheckStatus(); return vappropst->d; }
  inline double Water::VaporSaturationMolarVolumeForT()     { CheckStatus(); return mol_h2o/vappropst->d; }
  inline double Water::VaporSaturationEnthalpyForT()        { CheckStatus(); return vappropst->h; }
  inline double Water::VaporSaturationHeatCapacityForT()    { CheckStatus(); return vappropst->cp; }

  inline double Water::SaturationPressure()                 { CheckStatus(); return liqpropst->p; }
  inline double Water::SaturationTemperature()              { CheckStatus(); return liqpropsp->T-kelvin; }
  inline double Water::SaturationTemperatureKelvin()        { CheckStatus(); return liqpropsp->T; }

  inline double Water::EnthalpyAtTPlusDelT()                { CheckStatus(); return delTproperties->h; }
  inline double Water::DensityAtTPlusDelT()                 { CheckStatus(); return delTproperties->d; }
  inline double Water::LiquidSaturationDensityAtTPlusDelT() { CheckStatus(); return liqpropstplusdelt->d; }
  inline double Water::VaporSaturationDensityAtTPlusDelT()  { CheckStatus(); return vappropstplusdelt->d; }
  inline double Water::SaturationPressureAtTPlusDelT()      { CheckStatus(); return liqpropstplusdelt->p; }

  inline double Water::EnthalpyAtTMinusDelT()               { CheckStatus(); return minTproperties->h; }
  inline double Water::DensityAtTMinusDelT()                { CheckStatus(); return minTproperties->d; }
  inline double Water::LiquidSaturationDensityAtTMinusDelT(){ CheckStatus(); return liqpropstminusdelt->d; }
  inline double Water::VaporSaturationDensityAtTMinusDelT() { CheckStatus(); return vappropstminusdelt->d; }
  inline double Water::SaturationPressureAtTMinusDelT()     { CheckStatus(); return liqpropstminusdelt->p; }

  inline double Water::DensityAtPPlusDelP()                 { CheckStatus(); return delPproperties->d; }
  inline double Water::MolarVolumeAtPPlusDelP()             { CheckStatus(); return mol_h2o/delPproperties->d; }

  inline double Water::Viscosity()                          { CheckStatus(); return viscosity; }
  inline double Water::ViscosityFromTandP(const double& t, 
                                            const double& p)
  {
    double t_backup=tcurrent;	
    double p_backup=pcurrent;
    tcurrent = t; pcurrent = p;
    CheckStatus(); 
    double visc = Viscosity(t+273.15,properties->d);
    tcurrent = t_backup; pcurrent = p_backup;
    CheckStatus();
    return visc;
  }


  /**
     author: Thomas Driesner, ETH Zuerich
     contact: thomas.driesner@erdw.ethz.ch
     latest modifications: 
     - 
     
     @section motivation Motivation
     Water provides convenient interface to the "PROST4" library (as seen from the "steam4.h", showing that PROST is rather non-straightforward to use for non-experts). 

     Water is legacy code coming from Thomas Driesner's developments for the H2O-NaCl system. It is REQUIRED to build accurate lookup tables for the H2O-NaCl system as it is being used by the Brine class. 


     @section licensing Licensing

     "Water" makes use of the "PROST4" library (as seen from the "steam4.h"). PROST is distributed under GPL (which version, please re-check) and, hence, PROST4 is distributed as stand-alone external support library.

     @section usage Usage:
     Construct an instance of "Water" with the temperature (in C) and pressure (in Pa) as constructor variables. "Water" has an internal mechanism to make sure that it always uses the current value of that temperature and pressure. Default construction has therefore been disabled. Public member names should be self-explanatory, I hope.
     
     @code
     double  temperature(100.0);           // define your temperature variable, in [C]
     double  pressure(10.0e6);             // define your pressure Variable, in [Pa]
     Water     water(temperature, pressure); // instantiate a Water object
     cout << water.Density() << endl;        // would now return the density [kg m-3] at 100C and 10MPa

     ... // continue code code

     temperature = 250.;
     pressure    = 50.e6;
     cout << water.Density() << endl;        // would now return the density [kg m-3] at 1250C and 50MPa; T and P were updated automatically
     @endcode

     @section requirements Requirements:
     The PROST4 library is needed, i.e., the compiled library and the header "steam4.h" are required.
   
     @test Tested manually, OK

     @section issues Known Issues
     - In the internal computations, pressure still in bars, a legacy from the SoWat development. This should, however, not affect API users or anybody not modifying the source code.
     @todo Implement float comparisons etc. from compareFloats.h
  */
}//csmp
#endif
