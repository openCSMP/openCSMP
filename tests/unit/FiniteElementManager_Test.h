//
//  FiniteElementManager_Test.hpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 29/6/2024.
//

#ifndef CSMP_FINITE_ELEMENT_MANAGER_TEST_H
#define CSMP_FINITE_ELEMENT_MANAGER_TEST_H

#include "Test.h"
#include "FiniteElementManager.h"

namespace csmp {

/**
    Verifies functionality of FiniteElement manager and then tests its performance with a pressure diffusion model based on prism_test
*/
class FiniteElementManager_Test : public Test {
  public:
    virtual void run();
    
    /// builds different managers and demonstrates that they contain the correct element types
    void TestBasicFunctionality();
    void TestBasicFunctionality_FiniteElementManager1();
    
    /// uses a transient pressure diffusion calculation with tetra-triangle model to measure the speed of the accumulation
    void TestAccumulationSpeed();
    
  private:
    bool verbose_ = true;
};


} // end csmp

#endif /* CSMP_FINITE_ELEMENT_MANAGER_TEST_H */
