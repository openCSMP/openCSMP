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
     bool TestBoundaryFaceFunctionality( /* "VSET_MAKER" */ );
     bool TestBoundaryFaceFunctionality( const std::string& );
    
     /// test covers functionality of the base class ModelSubDomain
     bool TestRegionFileInputOutput( Model<3U>& model, const char* region );

    const bool verbose_;
};
  
/// checks whether the neighbor information matches the boundary face info for region "Model"
bool consistencyCheckNeighborVersusPerimeterFaces( const Model<3U>& );
  
} // end csmp

#endif // REGION_TEST_H
