#ifndef TWO_PHASE_MODEL_WITH_HYSTERESIS_TEST_H
#define TWO_PHASE_MODEL_WITH_HYSTERESIS_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

class TwoPhaseModelwithHysteresis_Test : public Test {
  public:
    TwoPhaseModelwithHysteresis_Test();
    ~TwoPhaseModelwithHysteresis_Test();
  
    virtual void run();
    
  private:
  
    enum TestCases : std::int8_t { CapillaryPressure, DerivativeOfCapillaryPressure, krw, krn };
  
    std::string parseTestCases( TestCases cases)
    {
      if ( cases == CapillaryPressure ) return "CapillaryPressure";
      if ( cases == DerivativeOfCapillaryPressure ) return "DerivativeOfCapillaryPressure";
      if ( cases == krw ) return "krw";
      if ( cases == krn ) return "krn";
      return "NONE" ;
    }
  
    void runOverSaturationRange( TestCases  );
    std::vector <std::pair< double64,double64 > > Extract_data( Model<1U>&, TestCases  ) ;
    
};

} // csmp

#endif // TWO_PHASE_MODEL_WITH_HYSTERESIS_TEST_H
