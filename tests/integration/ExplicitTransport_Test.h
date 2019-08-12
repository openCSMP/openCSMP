//
//  ExplicitTransport_Test.h
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 2/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_EXPLICIT_TRANSPORT_TEST_H
#define CSMP_EXPLICIT_TRANSPORT_TEST_H

#include "Test.h"

namespace csmp {

template<size_t> class Model;

/**
     Testing programme
 
     1. Conservation on the initialised model
        - prescribed velocity
        - computed velocity
        - tracer transport in steady-state pressure field, 1st order
        - tracer transport 2nd order
        - transport for transient pressure case
        - transport in the presence of sources and sinks
*/
class ExplicitTransport_Test : public Test {
  public:
    ExplicitTransport_Test( const char* test_model, const char* test_variables );
    virtual ~ExplicitTransport_Test();
    
    virtual void run();
  
  private:
    const std::string test_model_, test_variable_file_;
};


} // end csmp

#endif /* defined(CSMP_EXPLICIT_TRANSPORT_TEST_H) */
