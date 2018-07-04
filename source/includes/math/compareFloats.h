#ifndef COMPAREFLOATS_H
#define COMPAREFLOATS_H

#include <limits>
#include <cmath>

#include "CSMP_definitions.h"

namespace csmp
{

  // Functions that directly test within the numerical precision for double64 
  inline bool approximatelyEqual(double64 a, double64 b)
  {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double64>::epsilon());
  }
  
  inline bool essentiallyEqual(double64 a, double64 b)
  {
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double64>::epsilon());
  }
  
  inline bool definitelyGreaterThan(double64 a, double64 b)
  {
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double64>::epsilon());
  }
  
  inline bool definitelyLessThan(double64 a, double64 b)
  {
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double64>::epsilon());
  }


  // For user-defined tolerance
  inline bool approximatelyEqual(double64 a, double64 b, double64 epsilon)
  {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }
  
  inline bool essentiallyEqual(double64 a, double64 b, double64 epsilon )
  {
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }
  
  inline bool definitelyGreaterThan(double64 a, double64 b, double64 epsilon )
  {
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }
  
  inline bool definitelyLessThan(double64 a, double64 b, double64 epsilon )
  {
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }

  /**
     @file compareFloats.h "math/compareFloats.h"                                              
                                                                                                                      
     @author Thomas Driesner, ETH Zuerich                                                                             
     @section contact Contact                                                                                         
     thomas.driesner@erdw.ethz.ch                                                                                     
    
     @changes changes Latest Changes                                                                      
     Oct. 6, 2012: ported to CSMP++

     @section motivation Motivation                                                                                   
     A selection of inline functions for float comparisons that account for numerical precision. Function names should be self-explanatory. Functions can be traced back to Knuth's book and were found on: http://stackoverflow.com/questions/17333/most-effective-way-for-float-and-double-comparison

     @section testing Testing                                                                                         
     code was taken as is without further testing
  */                                                                                                                  
}
#endif
