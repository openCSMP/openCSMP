#ifndef CSMP_COMPARE_FLOATS_H
#define CSMP_COMPARE_FLOATS_H

#include <limits>
#include <cmath>

namespace csmp
{

  // Functions that directly test within the numerical precision for double 
  inline bool approximatelyEqual(double a, double b)
  {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double>::epsilon());
  }
  
  inline bool essentiallyEqual(double a, double b)
  {
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double>::epsilon());
  }
  
  inline bool definitelyGreaterThan(double a, double b)
  {
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double>::epsilon());
  }
  
  inline bool definitelyLessThan(double a, double b)
  {
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<double>::epsilon());
  }


  // For user-defined tolerance
  inline bool approximatelyEqual(double a, double b, double epsilon)
  {
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }
  
  // whichever is exceeded, the absolute or the relative tolerance
  inline bool approximatelyEqual2( double a, double b, double relTol = 1e-12, double absTol = 1e-14 )
  {
    return fabs(a - b) <= std::max( relTol * std::max( fabs(a), fabs(b)), absTol );
  }
  
  inline bool essentiallyEqual(double a, double b, double epsilon )
  {
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }
  
  inline bool definitelyGreaterThan(double a, double b, double epsilon )
  {
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon );
  }
  
  inline bool definitelyLessThan(double a, double b, double epsilon )
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
