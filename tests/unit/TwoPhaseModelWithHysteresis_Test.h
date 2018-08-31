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

    /// do drainage and imbibition cycles, move along scanning curves
  
    bool Out_of_Bond(double64 , double64 , double64 );
  
    void UpdatePseduoResidualSaturations(Model<1U>& , double64 , double64 );
  
    void runOverSaturationRange(Model<1U>& );

    void CreateArtifitialDrainageProcess(Model<1U>&) ;
  
    void CreateArtifitialImbibitionProcess(Model<1U>&);
  
    std::vector <std::pair< double64,double64 > > Process_Test(Model<1U>&) ;
  
    std::vector <std::pair< double64,double64 > > Imbibition_Test(Model<1U>&, const std::array<std::array<double64,2>,2>& , const std::array<std::array<double64,2>,2>& ) ;
  
    std::vector <std::pair< double64,double64 > > Drainage_Test(Model<1U>&,  const std::array<std::array<double64,2>,2>& , const std::array<std::array<double64,2>,2>& ) ;
  
    std::vector <std::pair< double64,double64 > > TransitionProcess_Test(Model<1U>&) ;
  
    std::vector <std::pair< double64,double64 > > DrainageToImbibition_Test(Model<1U>&, const std::array<std::array<double64,2>,2>& , const std::array<std::array<double64,2>,2>& ) ;
  
    std::vector <std::pair< double64,double64 > > ImbibitionToDrainage_Test(Model<1U>&, const std::array<std::array<double64,2>,2>&, const std::array<std::array<double64,2>,2>& ) ;

    
};

} // csmp

#endif // TWO_PHASE_MODEL_WITH_HYSTERESIS_TEST_H
