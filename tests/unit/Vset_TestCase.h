#ifndef VSET_TESTCASE_H
#define VSET_TESTCASE_H

#include <iostream>
#include "Model.h"
#include "Test.h"

namespace csmp {

template<size_t> class Model;

class Vset_TestCase : public Test
  {
    public:
    explicit Vset_TestCase( const char* prefix);
        ~Vset_TestCase();
        virtual void run();
  };

} // csmp

#endif // VSET_TESTCASE_H
