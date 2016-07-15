#ifndef PARALLEL_TESTCASE_H
#define PARALLEL_TESTCASE_H

#include "Test.h"
#include "Region.h"

namespace csmp {

  template<size_t> class Model;
  template<size_t> class TwoPhaseModel;

  class Parallel_TestCase : public Test
    {
    public:
      virtual void run();

    private:
      void UpdateFlowProps( Model<3>& model, TwoPhaseModel<3>& saturationFunctions );
      void ComputeTotalVelocity( csmp::Model<3U>& model ) const;
    };

  } // csmp

#endif