// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TRIPLEPOINTNACL_H
#define TRIPLEPOINTNACL_H

#include "CSMP_definitions.h"
/* Changelog

     Jan 29, 2013, Thomas Driesner: ported to CSMP++; disabled all functionality for volumetric and thermal properties as those are computed from functions that rather reside in Brine, Halite, etc.;
     September 9, 2014, Thomas Driesner: minor polishing, updated documentation and coding style to better match csmp++ conventions, converted pressure to output [Pa]
*/  

namespace csmp
{
  /// Temperature-Pressure-Composition coordinates of NaCl Triple Point  
  class TriplePointNaCl
  {

  public :
    TriplePointNaCl();
    ~TriplePointNaCl();
    
    double Temperature()       const;
    double Pressure()          const;

    double CompositionVapor()  const;
    double CompositionLiquid() const;
    double CompositionHalite() const;

  private :
	
    const double ttriple_nacl; 
    const double ptriple_nacl; 
    const double xtriple_nacl;

  };

  inline double TriplePointNaCl::Temperature()          const { return ttriple_nacl; }
  inline double TriplePointNaCl::Pressure()             const { return ptriple_nacl; }

  inline double TriplePointNaCl::CompositionVapor()     const { return xtriple_nacl; }
  inline double TriplePointNaCl::CompositionLiquid()    const { return xtriple_nacl; }
  inline double TriplePointNaCl::CompositionHalite()    const { return xtriple_nacl; }

  /**
     @class CriticalPointH2O CriticalPointH2O.h "eos/h2o_nacl/CriticalPointH2O.h"                                              
                                                                                                                      
     @author Thomas Driesner, ETH Zuerich                                                                             
     @section contact Contact                                                                                         
     thomas.driesner@erdw.ethz.ch                                                                                     
    
     @section motivation Motivation                                                                                   
     TriplePointNaCl contains the triple point parameters (temperature, pressure, and composition) of the pure NaCl phase halite, liquid, and vapor at its triple point according to the fits obtained in the paper

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-composition space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.

     @section dependencies Dependencies                                                                               

     @section testing Testing                                                                                         
  */
    

}//csmp
#endif
