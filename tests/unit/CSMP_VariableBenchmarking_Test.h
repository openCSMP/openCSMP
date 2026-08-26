//
//  CSMP_VariableBenchmarking_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 1/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_VARIABLE_BENCHMARKING_TEST_H
#define CSMP_VARIABLE_BENCHMARKING_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "TensorVariable.h"

namespace csmp {

// comparing performance of different makeScalar() implementations

class VariableBenchmarking_Test : public Test {
  public:
    VariableBenchmarking_Test() = default;

    void run() override final;
 };


} // end csmp


#endif /* CSMP_VARIABLE_BENCHMARKING_TEST_H */
