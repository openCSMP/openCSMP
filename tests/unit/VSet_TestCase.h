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
      explicit VSet_TestCase( bool verbose=false );
      ~VSet_TestCase();
      
      virtual void run();
    
    private:
    
      /// using VSet from VSet makers
      bool Test_ModelConstructionAndSaving2D();
      void Test_ANSYS_ModelConstructionAndSaving2D( const std::string& input_file_set="HorFracs2D" );
      void Test_ANSYS_ModelConstructionAndSaving3D( const std::string& input_file_set );
      
      /// uses  VSet maker-made model with line elements from vset_makers
      bool Test_EstablishElementConnectivity2D();
      
    private:
      const bool verbose_;
  };

} // csmp

#endif // VSET_TESTCASE_H
