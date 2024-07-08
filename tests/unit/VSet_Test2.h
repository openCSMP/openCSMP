#ifndef CSMP_VSET_TEST_2_H
#define CSMP_VSET_TEST_2_H

#include "Test.h"

namespace csmp {

template<uint32_t> class VSet;
template<uint32_t> class Model;

/**  Second test focussing on the VSet's extended functionality (higher level)
       Tests the use of VSet in model construction
       storage to disk,
       and reconstruction from binary file.
*/
class VSet_Test2 : public Test
  {
    public:
      explicit VSet_Test2( bool verbose=false );
      ~VSet_Test2();
      
      virtual void run();
    
    private:
    
      /// using VSet from VSet makers
      bool Test_ModelConstructionAndSaving2D();
      void Test_ModelConstructionAndSaving3D();
      
      /// uses  VSet maker-made model with line elements from vset_makers
      bool Test_EstablishElementConnectivity2D();
      bool Test_EstablishElementConnectivity3D();
      
      // helpers
      void BoundaryFlagsToVTK( VSet<3>& vset );
      
    private:
      const bool verbose_;
  };

} // csmp

#endif // CSMP_VSET_TEST_2_H
