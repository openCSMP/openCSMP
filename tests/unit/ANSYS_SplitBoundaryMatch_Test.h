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
    ANSYS_SplitBoundaryMatch_Test( bool verbose=false ) : verbose_(verbose) {}

    virtual void run();

  private:
    bool verbose_;
};

} // csmp

#endif // ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
