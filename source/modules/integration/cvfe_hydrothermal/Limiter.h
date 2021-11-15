#ifndef LIMITER_H
#define LIMITER_H

#include "Interrelation.h"

/*   Changelog
     November 2014, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

 /// Limiting the value of a variable to a specified range

template<size_t dim>
class Limiter : public Interrelation<dim> {
    Operand<dim>&  op;
    double min, max;
    
  public:
    Limiter( const PropertyDatabase<dim>& p,
             const char* prop,              // variable name of property
			 double min_value,            // minimum value of specified range
			 double max_value );          // maximum value of specified range
                         
    ~Limiter() {};
    void Calculate();
};

  /**
     @class Limiter Limiter.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Designed to keep variable within the validity of equation of state.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).
	  Use with care, as the values are changed somewhat arbitrarily, which may lead to unwanted behavior.

     @code
	 The interralation reads in a variable and checks whether the value lies within the specified range.
	 If not, the value is corrected to the minimum or maximum value, respectively.
          
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues
	 Currently, the user does not receive any information whether the values are changed or not.
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
} // csmp

#endif

