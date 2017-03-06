#ifndef VARIABLES_TEST_H
#define VARIABLES_TEST_H

#include "Test.h"

namespace csmp  {

template<size_t> class Model;

class Variables_Test : public Test {
  public:
    explicit Variables_Test( const char* );
    virtual void run();

  private:
    void runModel( Model<3>& );
    const char* prefix_;
};

  } // csmp

#endif /* VARIABLES_TEST_H */
