#ifndef VARIABLES_TEST_CASE_H
#define VARIABLES_TEST_CASE_H

#include "Test.h"

namespace csmp
  {
    template<size_t> class Model;

  class Variables_TestCase : public Test
    {
    const char* prefix_;
    public:
      explicit Variables_TestCase( const char* );
      virtual void run();
    private:
      virtual void runModel( Model<3>& );
    };

  } // csmp

#endif