#ifndef VARIABLES_TEST_H
#define VARIABLES_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp  {

template<uint32_t> class Model;

/**

  Tests Read(), Store() and Status() operations for all variable placements and types with the FracBox model which contains Boundary objects.
  The operations are performed twice: 1) on the live model and 2) on a model recovered form disk.
  
  @todo Add tests for variables placed on SplitBoundary objects.
  
@author Phillip Lang
@date 9/24/2012

*/
class Variables_Test : public Test {
  public:
    /// takes the name of a test model as argument; this model should contain Boundary and SplitBoundary objects
    explicit Variables_Test( const char* );
    virtual void run();

  private:
    void runModel( Model<3U>& );
    const char* prefix_ = nullptr;
};

  } // csmp

#endif /* VARIABLES_TEST_H */
