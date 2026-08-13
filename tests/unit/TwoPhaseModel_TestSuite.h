#ifndef TWO_PHASE_MODEL_TESTS_H
#define TWO_PHASE_MODEL_TESTS_H

#include "Test.h"
#include "TestSuite.h"

#include "Model1D.h"

namespace csmp {

/// PL Nov 2010
class TwoPhaseModel_TestSuite : public Test
{
public:
  TwoPhaseModel_TestSuite( TestSuite& suite ) : suite_(suite), free_(true) {}
  virtual ~TwoPhaseModel_TestSuite();
  
  virtual void run();

private:
  void AssignSaturationValues( Model<1U>* );

  TestSuite&  suite_;
  bool        free_;

  Model<1U>*  rock_model_          = nullptr;
  Model<1U>*  fracture_rock_model_ = nullptr;

};

} // csmp

#endif // TWOPHASEMODEL_TESTS_H
