//
//  IsoparametricLinearPyramid_Test.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 15/10/2024.
//

#ifndef ISOPARAMETRIC_LINEAR_PYRAMID_TEST_H
#define ISOPARAMETRIC_LINEAR_PYRAMID_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Element;

class IsoparametricLinearPyramid_Test : public Test {
  public:
    virtual void run() override;
    
    Element<3> CreateLocalCoordinatePyramidElement( int n_integration_points=1 );
    Element<3> CreateDistortedPyramidElement(  int n_integration_points=1 );
    
    bool InterpolationFunctionAtIntegrationPointsTest();
    bool SumOfInterpolationFunctionTest(); // at integration points
    bool InterpolationFunctionDerivative_Test(); // at the nodes
    bool Volume_Test();

 private:
   const bool verbose_ = true;
};

} // end csmp

#endif /* ISOPARAMETRIC_LINEAR_PYRAMID_TEST_H */
