#ifndef BOX_TEST_H
#define BOX_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "Box.h"

namespace csmp {

template<uint32_t> class Element;

class Box_Test : public Test
{
  public:
    explicit Box_Test( const char* name="cube_flag" )
      : model_name_(name) {}

    virtual void run();
    
    void TestParsingOfFlags();
  
    bool TestBoundaryFlagAssigment2D();
    
    bool TestWhetherSideBoundaryFlagsArePresent();
    
    bool TestWhetherAllBoxFlagsArePresent();
    
    bool TestBoundaryFlagging();
  
    bool TestBoundaryVersusBOX_BOUNDARY_Flagging();

    bool TestBoundaryFlagRecreation();
    
    bool TestWhetherBoundaryFlagsArePreservedInBinaryFile();
  
    // new format of binary IO
    bool TestWhetherBoundaryFlagsArePreservedInBinaryFile1();

    bool TestConsistencyOfBoxFlaggingWithBoundaryIdentification();
  
    /// testing the normals of the volumetric elements and the boundary faces
    void TestWhetherElementNormalsAreOutwardPointing();

  private:
    std::map<std::string,std::vector<Point<3>>> FormPointCloudsFromBOX_BOUNDARY_Flags( Model<3>& );
  
    std::string model_name_;
    static const bool verbose_ = true;
};

bool writePrismWithNormalsToVTK( const csmp::Element<3>* eptr, const std::string& filename );
bool writePyramidWithNormalsToVTK( const csmp::Element<3>* eptr, const std::string& filename );

} // csmp

#endif // BOX_TEST_H
