//
//  TRIANGLE_Interface_Test.h
//  Open CSMP++
//
//  Verifies the input generated reading file suite from TRIANGLE mesher
//  (J. Shewchuk, UC Berkeley)
//
//  Created by Stephan Matthai on 14/11/2024.
//

#ifndef CSMP_TRIANGLE_INTERFACE_TEST_H
#define CSMP_TRIANGLE_INTERFACE_TEST_H

#include "Test.h"

namespace csmp {

// TODO: flex out to include extra checks
class TRIANGLE_Interface_Test : public Test {
  public:
    virtual void run();
    
    void Test_BOX_BOUNDARY_NodeFlags();

   const bool verbose_ = false;
};


} // csmp

#endif /* CSMP_TRIANGLE_INTERFACE_TEST_H */
