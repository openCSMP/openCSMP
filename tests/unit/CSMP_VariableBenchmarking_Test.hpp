//
//  CSMP_VariableBenchmarking_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 1/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_VariableBenchmarking_Test_hpp
#define CSMP_VariableBenchmarking_Test_hpp

#include "Test.h"
#include "TensorVariable.h"

namespace csmp {

// comparing performance of different makeScalar() implementations

class VariableBenchmarking_Test : public Test {
    virtual void run();
 };


} // end csmp


#endif /* CSMP_VariableBenchmarking_Test_hpp */
