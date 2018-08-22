#ifndef REGION_TEST_H
#define REGION_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp
{
 // class Index;
  class Region_Test: public Test
  {
  public:
      Region_Test( bool verbose=false );
      ~Region_Test();
      void run();
    
  private:
     bool TestBoundaryFaceFunctionality( const std::string& );
     bool TestBoundaryFaceFunctionality( /* "VSET_MAKER" */ );

    const bool verbose_;
};
  
/// checks whether the neighbor information matches the boundary face info for region "Model"
bool consistencyCheckNeighborVersusPerimeterFaces( const Model<3U>& );
  
} // end csmp

#endif // REGION_TEST_H
