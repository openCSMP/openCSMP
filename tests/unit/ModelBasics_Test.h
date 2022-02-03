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

template<size_t> class Model;

class ModelBasics_Test: public Test {
  public:
      explicit ModelBasics_Test( bool verbose=false );
      ~ModelBasics_Test();
      
      virtual void run();
  
  private:
    bool verbose_;

};

} // end csmp

#endif /* CSMP_MODEL_BASICS_TEST_H */
