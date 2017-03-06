#ifndef VSET_TESTCASE_H
#define VSET_TESTCASE_H

#include "Model.h"
#include "Test.h"

namespace csmp {

template<size_t> class Model;

class Vset_TestCase : public Test
  {
    public:
      explicit Vset_TestCase( const char* prefix,
                              bool verbose );
      ~Vset_TestCase();
      virtual void run();
    
    private:
      const std::string model_file_;
      const bool verbose_;
  };

} // csmp

#endif // VSET_TESTCASE_H
