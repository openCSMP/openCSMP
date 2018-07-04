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
    Water( const double64& externaltemperature,  // [C] 
           const double64& externalpressure );   // [Pa]
    ~Water();

    double64 Pressure()        const ;
    double64 PressurePascal()  const ;
    double64 Temperature()     const ;
    double64 Density()               ;
    double64 MolarVolume()           ;
    double64 Enthalpy()              ;
    double64 FreeEnergy()            ;
    double64 HeatCapacity()          ;
    double64 Dvdppascal()            ;
    double64 Compressibilitypascal() ;
    double64 Dvdpbar()               ;
    double64 Dvdt()                  ;
    double64 Compressibilitybar()    ;
    double64 MinusDvdppascal()       ;
    double64 MinusDvdpbar()          ;

    double64 MolarVolumeAtPPlusDelP();
    double64 DensityAtPPlusDelP();

    double64 EnthalpyAtTPlusDelT()   ;
    double64 DensityAtTPlusDelT()    ;
    double64 LiquidSaturationDensityAtTPlusDelT() ;
    double64 VaporSaturationDensityAtTPlusDelT()  ;
    double64 SaturationPressureAtTPlusDelT()      ;

    double64 EnthalpyAtTMinusDelT()   ;
    double64 DensityAtTMinusDelT()    ;
    double64 LiquidSaturationDensityAtTMinusDelT() ;
    double64 VaporSaturationDensityAtTMinusDelT()  ;
    double64 SaturationPressureAtTMinusDelT()      ;

    double64 LiquidSaturationDensityForP()      ;
    double64 LiquidSaturationMolarVolumeForP()  ;
    double64 LiquidSaturationEnthalpyForP()     ;
    double64 LiquidSaturationHeatCapacityForP() ;
    double64 VaporSaturationDensityForP()       ;
    double64 VaporSaturationMolarVolumeForP()   ;
    double64 VaporSaturationEnthalpyForP()      ;
    double64 VaporSaturationHeatCapacityForP()  ;
    double64 LiquidSaturationDensityForT()      ;
    double64 LiquidSaturationMolarVolumeForT()  ;
    double64 LiquidSaturationEnthalpyForT()     ;
    double64 LiquidSaturationHeatCapacityForT() ;
    double64 VaporSaturationDensityForT()       ;
    double64 VaporSaturationMolarVolumeForT()   ;
    double64 VaporSaturationEnthalpyForT()      ;
    double64 VaporSaturationHeatCapacityForT()  ;
    double64 CriticalDensity()                  ;
    double64 CriticalMolarVolume()              ;
    double64 CriticalEnthalpy()                 ;
    double64 CriticalHeatCapacity()             ;
    double64 SaturationPressure()               ;
    double64 SaturationTemperature()            ;
    double64 SaturationTemperatureKelvin()      ;
    double64 Viscosity()                        ;
    double64 ViscosityFromTandP(const double64& t, const double64& p);
    double64 LiquidViscosityFromT(const double64& t);
    
  private:
    
    Water();

    const double64&    temperature;        ///< reference to temperature in flow code
    const double64&    pressure;           ///< reference to pressure in flow code

    const double64     kelvin;             ///< for C -> K conversion
    const double64     mol_h2o;            ///< mole / mass conversion factor

    double64           tcurrent;           ///< internal temperature variable [C]
    double64           pcurrent;           ///< internal pressure variable [bar]
    double64           d;                  ///< dummy variable, for density estimate when calling PROST functions
    double64           dp;                 ///< dummy variable, for density tolerance when calling PROST functions
    double64           dvdppascal;         ///< dv/dP in [(cm^3/kg)/pa]
    double64           comprpascal;        ///< compressibility [Pa-1]
    double64           delp;               ///< pressure delta [bar] for numerical evaluation of derivatives
    double64           delt;               ///< temperature delta [C] for numerical evaluation of derivatives
    double64           dvdt;               ///< dv/dT in [(cm^3/kg)/K]
    double64           viscosity;          ///< dynamic viscosity [Pa s]
    double64           ak[4];              ///< coefficients of viscosity formula
    double64           bij[6][5];          ///< more coefficients of viscosity formula

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
    double64           Viscosity(const double64& T, const double64& rho);

    double64           T() const;
    double64           P() const;

  };

  inline double64 Water::Pressure()                   const { return pcurrent; }
  inline double64 Water::PressurePascal()             const { return pcurrent*1.0e5; }
  inline double64 Water::Temperature()                const { return tcurrent; }

  inline double64 Water::T()                          const { return temperature+273.15e0; }
  inline double64 Water::P()                          const { return pressure*1.0e5; }

  inline double64 Water::Density()                          { CheckStatus(); return properties->d; }
  inline double64 Water::MolarVolume()                      { CheckStatus(); return mol_h2o/properties->d; }
  inline double64 Water::Enthalpy()                         { CheckStatus(); return properties->h; }
  inline double64 Water::FreeEnergy()                       { CheckStatus(); return properties->g; }
  inline double64 Water::HeatCapacity()                     { CheckStatus(); return properties->cp; }
  inline double64 Water::Dvdppascal()                       { CheckStatus(); return dvdppascal; }
  inline double64 Water::Compressibilitypascal()            { CheckStatus(); return comprpascal; }
  inline double64 Water::Dvdpbar()                          { CheckStatus(); return dvdppascal*1.0e5; }
  inline double64 Water::Compressibilitybar()               { CheckStatus(); return comprpascal*1.0e5; }
  inline double64 Water::MinusDvdppascal()                  { CheckStatus(); return -dvdppascal; }
  inline double64 Water::MinusDvdpbar()                     { CheckStatus(); return -dvdppascal*1.0e5; }
  inline double64 Water::Dvdt()                             { CheckStatus(); return dvdt; }

  inline double64 Water::LiquidSaturationDensityForP()      { CheckStatus(); return liqpropsp->d; }
  inline double64 Water::LiquidSaturationMolarVolumeForP()  { CheckStatus(); return mol_h2o/liqpropsp->d; }
  inline double64 Water::LiquidSaturationEnthalpyForP()     { CheckStatus(); return liqpropsp->h; }
  inline double64 Water::LiquidSaturationHeatCapacityForP() { CheckStatus(); return liqpropsp->cp; }

  inline double64 Water::VaporSaturationDensityForP()       { CheckStatus(); return vappropsp->d; }
  inline double64 Water::VaporSaturationMolarVolumeForP()   { CheckStatus(); return mol_h2o/vappropsp->d; }
  inline double64 Water::VaporSaturationEnthalpyForP()      { CheckStatus(); return vappropsp->h; }
  inline double64 Water::VaporSaturationHeatCapacityForP()  { CheckStatus(); return vappropsp->cp; }

  inline double64 Water::LiquidSaturationDensityForT()      { CheckStatus(); return liqpropst->d; }
  inline double64 Water::LiquidSaturationMolarVolumeForT()  { CheckStatus(); return mol_h2o/liqpropst->d; }
  inline double64 Water::LiquidSaturationEnthalpyForT()     { CheckStatus(); return liqpropst->h; }
  inline double64 Water::LiquidSaturationHeatCapacityForT() { CheckStatus(); return liqpropst->cp; }
  inline double64 Water::LiquidViscosityFromT(const double64& t)
  {
    double64 t_backup=tcurrent;
    tcurrent = t;
    CheckStatus(); 
    double64 visc = Viscosity(t+273.15,liqpropst->d);
    tcurrent = t_backup;
    CheckStatus();
    return visc;
  }

  inline double64 Water::VaporSaturationDensityForT()         { CheckStatus(); return vappropst->d; }
  inline double64 Water::VaporSaturationMolarVolumeForT()     { CheckStatus(); return mol_h2o/vappropst->d; }
  inline double64 Water::VaporSaturationEnthalpyForT()        { CheckStatus(); return vappropst->h; }
  inline double64 Water::VaporSaturationHeatCapacityForT()    { CheckStatus(); return vappropst->cp; }

  inline double64 Water::SaturationPressure()                 { CheckStatus(); return liqpropst->p; }
  inline double64 Water::SaturationTemperature()              { CheckStatus(); return liqpropsp->T-kelvin; }
  inline double64 Water::SaturationTemperatureKelvin()        { CheckStatus(); return liqpropsp->T; }

  inline double64 Water::EnthalpyAtTPlusDelT()                { CheckStatus(); return delTproperties->h; }
  inline double64 Water::DensityAtTPlusDelT()                 { CheckStatus(); return delTproperties->d; }
  inline double64 Water::LiquidSaturationDensityAtTPlusDelT() { CheckStatus(); return liqpropstplusdelt->d; }
  inline double64 Water::VaporSaturationDensityAtTPlusDelT()  { CheckStatus(); return vappropstplusdelt->d; }
  inline double64 Water::SaturationPressureAtTPlusDelT()      { CheckStatus(); return liqpropstplusdelt->p; }

  inline double64 Water::EnthalpyAtTMinusDelT()               { CheckStatus(); return minTproperties->h; }
  inline double64 Water::DensityAtTMinusDelT()                { CheckStatus(); return minTproperties->d; }
  inline double64 Water::LiquidSaturationDensityAtTMinusDelT(){ CheckStatus(); return liqpropstminusdelt->d; }
  inline double64 Water::VaporSaturationDensityAtTMinusDelT() { CheckStatus(); return vappropstminusdelt->d; }
  inline double64 Water::SaturationPressureAtTMinusDelT()     { CheckStatus(); return liqpropstminusdelt->p; }

  inline double64 Water::DensityAtPPlusDelP()                 { CheckStatus(); return delPproperties->d; }
  inline double64 Water::MolarVolumeAtPPlusDelP()             { CheckStatus(); return mol_h2o/delPproperties->d; }

  inline double64 Water::Viscosity()                          { CheckStatus(); return viscosity; }
  inline double64 Water::ViscosityFromTandP(const double64& t, 
                                            const double64& p)
  {
    double64 t_backup=tcurrent;	
    double64 p_backup=pcurrent;
    tcurrent = t; pcurrent = p;
    CheckStatus(); 
    double64 visc = Viscosity(t+273.15,properties->d);
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
     double64  temperature(100.0);           // define your temperature variable, in [C]
     double64  pressure(10.0e6);             // define your pressure Variable, in [Pa]
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
