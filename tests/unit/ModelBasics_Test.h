/*
 *  ModelBasics_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *
 */

#ifndef CSMP_MODEL_BASICS_TEST_H
#define CSMP_MODEL_BASICS_TEST_H

#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

/**
       Estabilishes that basic functionality of the Model works:
        - Model can be built from VSet
        - Model with Boundaries can be built from ModelTopology and VSet, written to disk and read again.
*/
class ModelBasics_Test: public Test {
  public:
      virtual void run();
      
  private:
     // from VSet and variables file only
     bool TestModelConstructionFromVSet();
  
     bool TestWriteModelToDiskAndReadBack();
     bool TestWriteModelToDiskAndReadBackWithInterfaces();
  
  private:
    const static bool verbose_ = true; // TURN VERBOSE ON HERE

};

} // end csmp

#endif /* CSMP_MODEL_BASICS_TEST_H */
