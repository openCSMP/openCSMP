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
     void run();
    
  private:
     bool TestBoundaryFaceFunctionality( /* "VSET_MAKER" */ );
     bool TestBoundaryFaceFunctionality( const std::string& );
    
     /// test covers functionality of the base class ModelSubDomain
     bool TestRegionFileInputOutput( Model<3U>& model, const char* region );

    const static bool verbose_ = false;
};
  
/// checks whether the neighbor information matches the boundary face info for region "Model"
template<uint32_t dim>
bool consistencyCheckNeighborVersusPerimeterFaces( const Model<dim>& );
  
} // end csmp

#endif // REGION_TEST_H
