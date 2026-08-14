#ifndef TRIPLEPOINTH2O_H
#define TRIPLEPOINTH2O_H

#include "CSMP_definitions.h"

/* Changelog

     October 15, 2014, Thomas Driesner: port to CSMP++, changed Pressure to output [Pa]
*/  

namespace csmp
{

  /// Temperature and pressure at the triple point of pure H2O
  class TriplePointH2O
  {
    
    public :

    TriplePointH2O();
    ~TriplePointH2O();
    
    double Temperature() const;
    double Pressure()    const;
    
  private :

  };

  /**
     @class CriticalPointH2O CriticalPointH2O.h "eos/h2o_h2o/CriticalPointH2O.h"                                              
                                                                                                                      
     @author Thomas Driesner, ETH Zuerich                                                                             
     @section contact Contact                                                                                         
     thomas.driesner@erdw.ethz.ch                                                                                     
    
     @section motivation Motivation                                                                                   
     TriplePointH2O contains the triple point parameters (temperature, pressure) of pure H2O.
  */
    

}//csmp
#endif
