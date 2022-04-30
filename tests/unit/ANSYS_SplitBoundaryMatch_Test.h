#ifndef ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
#define ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H

#include "Test.h"

namespace csmp {

/**
    Tests high level utility function that matches node-matched elements
    across SplitBoundaries created in ANSYS.
*/
class ANSYS_SplitBoundaryMatch_Test : public Test {
  public:
    virtual void run();

  private:
    const static bool verbose_ = true;
    
    bool TestForContiguousModel( /* "Split_Edges" */ );
    bool TestForDiscontiguousModel();
};

} // csmp

#endif // ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
