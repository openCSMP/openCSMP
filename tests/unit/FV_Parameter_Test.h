#ifndef CSMP_FV_PARAMETER_TEST_H
#define CSMP_FV_PARAMETER_TEST_H

#include "Test.h"
#include "CSMP_number_types.h"

namespace csmp 
{

/// test does not generate any console output
class FV_Parameter_Test : public Test {
  public:
    FV_Parameter_Test();
    ~FV_Parameter_Test();
  
    void run(); // runs all the tests for the class (register other methods)
    
    void FV_ParameterCtor();
    void FV_ParameterCopyCtor();
    void FV_ParameterEqual();
    void FV_ParameterInitialize();
    void FV_ParameterResize();
    void FV_ParameterSectorVolume();
    void FV_ParameterFacetArea();
    void FV_ParameterFacetNormal();
    void FV_ParameterFacetNormalVelocity();
    void FV_ParameterFacetNormalProjection();

  private:
    double64 fTolerance;
  
}; // end class

} // end csmp

#endif
