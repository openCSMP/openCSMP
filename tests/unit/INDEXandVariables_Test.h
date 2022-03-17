#ifndef CSMP_INDEX_AND_VARIABLES_TEST_H
#define CSMP_INDEX_AND_VARIABLES_TEST_H

#include "Test.h"

namespace csmp  {

template<uint32_t> class Model;

/**
same as Variables_Test, but instead of using runtime csmp::Index variables, the typed compile time versions are used 
and the corresponding Read(), Store() and Status() operations.

@todo test new methods that allow the selective reading of ArrayVariable members.

@author SKM
@date 21/6/2020

*/
class INDEXandVariables_Test : public Test {
  public:
    explicit INDEXandVariables_Test( const char* );
    virtual void run();

  private:
    void runModel( Model<3>& );
    const char* prefix_;
};

  } // csmp

#endif /* CSMP_INDEX_AND_VARIABLES_TEST_H */
