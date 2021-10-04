#ifndef VSET_TESTCASE_H
#define VSET_TESTCASE_H

#include "Model.h"
#include "Test.h"

namespace csmp {

template<size_t> class Model;

/**
       Tests the use of VSet in model construction 
       storage to disk,
       and reconstruction from binary file.
*/
class VSet_TestCase : public Test
  {
    public:
      explicit VSet_TestCase( const char* prefix,
                              bool verbose );
      ~VSet_TestCase();
      virtual void run();
    
    private:
      void TestModelConstructionAndSaving2D();
      void TestModelConstructionAndSaving3D();
      
      // uses specific model with line elements from vset_makers
      bool Test_EstablishElementConnectivity2D();
      
    private:
      const std::string model_file_;
      const bool verbose_;
  };

} // csmp

#endif // VSET_TESTCASE_H
