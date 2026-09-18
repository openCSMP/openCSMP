// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef HALITELOOKUP_H
#define HALITELOOKUP_H

#include <cmath>

#include "CSMP_definitions.h"
#include "NaClMeltingCurveHaliteLookup.h"

namespace csmp
{
  /// Build and/or query lookup tables for thermodynamic properties of halite (crystalline NaCl)
  class HaliteLookup
  {
    
  public:
    HaliteLookup(const double& externaltemperature,
                 const double& externalpressure);
    ~HaliteLookup();

    double        MassFractionNaCl();// [mass fraction NaCl]
    double        Density();         // [kg m-3]
    double        Enthalpy();        // [J kg-3]
    double        HeatCapacity();    // [J kg-1 K-1]
    double        Compressibility(); // [Pa-1]
    double        Viscosity();       // [Pa-1], dummy for completeness
    double        ValueOf(const int& property_index);

  private:

    const double& temperature_;     ///< reference to temperature [C] in flow code
    const double& pressure_;        ///< reference to pressure [Pa] in flow code

    double        tcurrent_;        ///< temperature [C] for internal use
    double        pcurrent_;        ///< pressure [Pa], internal use
    double        tdummy_;          ///< another temperature variable [C] for internal use
    double        pdummy_;          ///< another pressure variable [Pa] for internal use
    double        t_res_;           ///< temperature interval [C] in lookup table
    double        p_res_;           ///< pressure interval [bar] in lookup table
    double        tnorm_;           ///< normalized temperature inside lookup cell
    double        pnorm_;           ///< normalized pressure inside lookup cell
    double        t_iA_;            ///< temperature at low-T, low-P corner of lookup cell
    double        t_iB_;            ///< temperature at high-T, low-P corner of lookup cell
    double        t_iC_;            ///< temperature at high-T, high-P corner of lookup cell
    double        t_iD_;            ///< temperature at low-T, high-P left corner of lookup cell
    double        p_iA_;            ///< pressure at low-T, low-P corner of lookup cell
    double        p_iB_;            ///< pressure at high-T, low-P corner of lookup cell
    double        p_iC_;            ///< pressure at high-T, high-P corner of lookup cell
    double        p_iD_;            ///< pressure at low-T, high-P left corner of lookup cell
    double        v_bottom_;        ///< value of interest at low-P endpoint of isothermal interpolation
    double        v_top_;           ///< value of interest at high-P endpoint of isothermal interpolation
    double        v_interpolated_;  ///< interpolated value of interest in lookup cell
    double        v_iA_;            ///< pressure at low-T, low-P corner of lookup cell
    double        v_iB_;            ///< pressure at high-T, low-P corner of lookup cell
    double        v_iC_;            ///< pressure at high-T, high-P corner of lookup cell
    double        v_iD_;            ///< pressure at low-T, high-P left corner of lookup cell
    double        v_before_;        ///< value of interest at low-T endpoint of isobaric interpolation
    double        v_behind_;        ///< value of interest at high-T endpoint of isobaric interpolation
    double        v_vlh_;           ///< value of interest on intersection with Vapor-Halite-Liquid coexistence surface
    double        v_vlh_behind_;    ///< value of interest on intersection with Vapor-Halite-Liquid coexistence surface if higher T
    double        v_vlh_before_;    ///< value of interest on intersection with Vapor-Halite-Liquid coexistence surface if lower T
    double        tvlh_;            ///< temperature at intersection with Vapor-Halite-Liquid coexistence surface

    long            it_;              ///< lookup grid line index along temperature axis
    long            ip_;              ///< lookup grid line index along pressure axis
    long            t_dim_;           ///< number of T-gridlines = maximum value of it_ + 1
    long            p_dim_;           ///< number of P-gridlines = maximum value of ip_ + 1
    long            iA_;              ///< lookup index of low-T, low-P corner of lookup cell
    long            iB_;              ///< lookup index of high-T, low-P corner of lookup cell  
    long            iC_;              ///< lookup index of high-T, high-P corner of lookup cell
    long            iD_;              ///< lookup index of low-T, high-P left corner of lookup cell
    long            i_dummy_;         ///< dummy index variable
    int             state_;           ///< phase state of system (always H for halite)
    int             state_iA_;        ///< phase state of system at low-T, low-P corner of lookup cell (always H for halite)
    int             state_iB_;        ///< phase state of system at high-T, low-P corner of lookup cell (always H for halite)
    int             state_iC_;        ///< phase state of system at high-T, high-P corner of lookup cell (always H for halite)
    int             state_iD_;        ///< phase state of system at low-T, high-P left corner of lookup cell (always H for halite)
    
    std::vector<double> storage_vector; ///< stores lookup data in sequence t-p-x-rho-h at each Lookup point
    // should the next one be <size_t> or <uint> rather than <int>?
    std::vector<int>      state_vector;   ///< stores phase state at each Lookup point

    double NormalInterpolation( const int& property_index );
    double NearVLHInterpolationLowT( const int& property_index );
    double NearVLHInterpolationHighT( const int& property_index );
    double NearNaClMeltInterpolation( const int& property_index );
   
    void     SetTemperatureAndPressure();
    void     GetIndex_iA(const int& property_index);
    void     GetTemperatureIndex(const double& t);
    void     GetPressureIndex(const double& p);

    NaClMeltingCurveHaliteLookup    naclmelt_h_lookup;
  };

  /**
     @class HaliteLookup HaliteLookup.h "eos/h2o_nacl/HaliteLookup.h"
                                                                                  
     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch

     @changes changes Latest Changes
     Jan 28, 2013: ported to CSMP++, changed all pressure units to Pa.
     Sep 10, 2014: adapted to CSMP++ coding style, updated doxygen comments

     @section motivation Motivation
     "HaliteLookup" stores various properties of the mineral halite (solid crystalline salt, NaCl) in lookup tables, as a function of temperature-pressure-composition [Celsius, Pa, mole fraction NaCl] according to the papers

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-composition space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.

     Driesner T. (2007): The system H2O-NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 to 1000oC, 1 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4902-4919.
                               
     @section usage Usage                                                                                             
     Construct an instance of "HaliteLookup" with temperature (in C) and pressure (in Pa) as they exist in the code that is supposed to use "HaliteLookup" as constructor variables. "HaliteLookup" has an internal mechanism to make sure that it always uses the current values of temperature and pressure . Public member names should be self-explanatory, I hope.

     Upon construction, HaliteLookup will check if the lookup tables do already exist. If not, they will be re-computed and written. The respective files are "HalitePropertiesLookupTable.bin" and "HaliteStateLookupTable.bin".

     @code                                                                                                            
     double t; // temperature [C] in user's application
     double p; // pressure [Pa] in user's application
     HaliteLookup haliteLookup(t,p);
     ...
     t = some_value;
     p = some_other_value;
     cout << haliteLookup.Density() << endl; // will return density at the new t and p conditions
     @endcode                                                                                                         
                                                                                                                      
     @section dependencies Dependencies                                                                               
     requires "ConvertConcentrationUnitsNaCl.h" to ensure consistent unit conversions   

     @section issues Known issues
     Composition() is in mole fraction NaCl, to convert to other units, the functions in "ConvertConcentrationUnitsNaCl.h" may be useful.

     To do: interpolation near the halite melting curve has not been coded optimally. Hence at extreme conditions (t > 800.7 C and salinities very close to 100%) tiny mismatches in fluid properties may appear and cause problems. 

     @section testing Testing
  */                                                                                                                  



}// namespace csmp
#endif






















