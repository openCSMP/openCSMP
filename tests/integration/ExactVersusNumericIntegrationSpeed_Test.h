//
//  ExactVersusNumericIntegrationSpeed_Test.hpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 7/12/2024.
//

#ifndef EXACT_VERSUS_NUMERIC_INTEGRATION_SPEED_TEST_H
#define EXACT_VERSUS_NUMERIC_INTEGRATION_SPEED_TEST_H

#include "Test.h"

namespace csmp {

class ExactVersusNumericIntegrationSpeed_Test : public Test {
  public:
    void run() override;
    
    ///  LocalCoordinateFiniteElement vs GlobalCoordinateFiniteElement pressure computation
    void CompareSpeed();
  
  private:
    const bool verbose_ = true;
};

} // end csmp

#endif /* EXACT_VERSUS_NUMERIC_INTEGRATION_SPEED_TEST_H */
