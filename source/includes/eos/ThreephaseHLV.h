// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef THREEPHASEHLV_H
#define THREEPHASEHLV_H

#include "CSMP_definitions.h"
#include "TriplePointNaCl.h"

/*   Changelog
     October 13, 2014, Thomas Driesner: porting to CSMP++, adjusting pressure to [Pa] units, doxygen documentation, style adjustments, basic manual testing against CSP version (OK)
*/

namespace csmp
{
  class ThreephaseHLV 
  {
  
  public:

    /// Pressure-temperature coordinates of threephase Halite+Liquid+Vapor coexistence
    ThreephaseHLV(const double& externaltemperature); // [C]
    ~ThreephaseHLV();

    double            Temperature() const;    // [C]
    double            Pressure();             // [Pa]
    double            DPressureDT();          // [Pa C-1]
    double            Tmax();                 // temperature for the pressure maximum near 590 C
    double            Pmax();                 // pressure at that maximum

  private:

    ThreephaseHLV();                            // disable default construction

    const double&     temperature_;           ///< const ref to temperature [C] in flow code

    const double      f0;                     ///< parameter f0  for eq. 10 in Driesner&Heinrich (2007)
    const double      f1;                     ///< parameter f1  for eq. 10 in Driesner&Heinrich (2007)
    const double      f2;                     ///< parameter f2  for eq. 10 in Driesner&Heinrich (2007)
    const double      f3;                     ///< parameter f3  for eq. 10 in Driesner&Heinrich (2007)
    const double      f4;                     ///< parameter f4  for eq. 10 in Driesner&Heinrich (2007)
    const double      f5;                     ///< parameter f5  for eq. 10 in Driesner&Heinrich (2007)
    const double      f6;                     ///< parameter f6  for eq. 10 in Driesner&Heinrich (2007)
    const double      f7;                     ///< parameter f7  for eq. 10 in Driesner&Heinrich (2007)
    const double      f8;                     ///< parameter f8  for eq. 10 in Driesner&Heinrich (2007)
    const double      f9;                     ///< parameter f9  for eq. 10 in Driesner&Heinrich (2007)
    const double      f10;                    ///< parameter f10 for eq. 10 in Driesner&Heinrich (2007)
 
    double            tcurrent_;              ///< temperature [C] for internal use
    double            myt_;                   ///< tcurrent_/nacl_triple.Temperature() both in [C]
    double            pressure_;              ///< pressure [bar], internal use
    double            dpressuredt_;           ///< pressure [bar] derivative w.r.t. temperature, internal use

    void                CheckStatus();
    double            Pressure(const double& myt);
    double            DPressureDT(const double& myt, const double& t);

    TriplePointNaCl     nacl_triple;
  };

  /**
     @class ThreephaeHLV ThreephaeHLV.h "eos/h2o_nacl/ThreephaeHLV.h"

     @author: Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation
     "ThreephaseHLV" computes the temperature-pressure-composition [Celsius, bar] coordinates of the vapor-liquid-halite threephase-coexistence surface of the system H2O-NaCl, according to the paper
     
     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-composition space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
     
     Other porperties such as density etc. are NOT provided.

     "ThreephaseHLV" is legacy code coming from Thomas Driesner's developments for the H2O-NaCl system. It is REQUIRED to build accurate lookup tables for the H2O-NaCl system. 
     
     @section usage Usage
     Construct an instance of "ThreephaseHLV" with the temperature (in C) as constructor variable. "ThreephaseHLV" has an internal mechanism to make sure that it always uses the current value of that temperature. Public member names should be self-explanatory, I hope.
     
     @section requirements Requirements:
     "ThreephaseHLV" needs "TriplePointNaCl"
     
     @section testing Testing:
     Manually tested against CSMP version
  */
  

}//csmp
#endif









