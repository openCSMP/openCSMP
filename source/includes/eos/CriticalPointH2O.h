// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CRITICALPOINTH2O_H
#define CRITICALPOINTH2O_H

#include "CSMP_definitions.h"

namespace csmp
{
  class CriticalPointH2O
  {
    
  public :
    CriticalPointH2O();
    ~CriticalPointH2O();
    
    double Temperature()      const;
    double Pressure()         const;
    double MassFractionNaCl() const;
    double Density()          const;
    double MolarVolume()      const;
    double Enthalpy()         const;
  };
  


  /**
     @class CriticalPointH2O CriticalPointH2O.h "eos/h2o_nacl/CriticalPointH2O.h"                                              
                                                                                                                      
     @author Thomas Driesner, ETH Zuerich                                                                             
     @section contact Contact                                                                                         
     thomas.driesner@erdw.ethz.ch                                                                                     
    
     @changes changes Latest Changes                                                                                  
     Oct. 6, 2012: ported to CSMP++, changed Pressure() to Pa.               
    
     @section motivation Motivation                                                                                   
     CriticalPointH2O contains the critical parameters (temperature, pressure, composition, density, molar volume, and specific enthalpy) of pure H2O according the 1984 IAPS or Haar-Gallagher-Kell equation of state (check www.iapws.org for more information).

     @section dependencies Dependencies                                                                               
     requires "ConvertConcentrationUnitsNaCl.h"

     @section testing Testing                                                                                         
     testing was done in the period before publication in 2007, prior to writing this documentation, no details are a\
     vailable anymore                                                                                                      
  */                                                                                                                  

}
#endif
