#ifndef TWOPHASEMODEL_TEST_H
#define TWOPHASEMODEL_TEST_H

#include "Test.h"

#include "VTU_Interface.h"

#include "TwoPhaseModel.h"

#include "BrooksCorey.h"
#include "FourarLenormand.h"
#include "VanGenuchten.h"
#include "TwoPhaseFileBased.h"
#include "Experimental2PhaseModel.h"
#include "ExperimentalRT.h"
#include "FractureMatrixUpscaled.h"
#include "BrooksCoreyWithHysteresis.h"
#include "LinearTwoPhaseModel.h"


namespace csmp{

template<size_t> class Model;

/// PL Nov 2010
class TwoPhaseModel_Test : public Test
{
public:

  explicit TwoPhaseModel_Test( Model<1U>*,
                               TwoPhaseModel<1U>*,
                               const char* testName,
                               const char* krn,
                               const char* krw,
                               const char* pc,
                               size_t upToElementNumber = 0 );
  ~TwoPhaseModel_Test();
  virtual void run();

private:

   Model<1U>* model_;

   TwoPhaseModel<1U>* twoPhaseModel_;

   size_t upToElementNumber_;

   Index krn_key_,
         krw_key_,
         pc_key_;

   std::list<string> vtuOutputProps_;
};

} // csmp

#endif // TWOPHASEMODEL_TEST_H
