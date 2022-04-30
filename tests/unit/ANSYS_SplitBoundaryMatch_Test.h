#ifndef ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
#define ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H

#include "Test.h"
#include "CSMP_definitions.h"

namespace csmp {

class Index;
template<uint32_t> class Node;

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
    
    /// for testing how the nodes are organised in the node vector
    void PrintSplitBoundaryNodeVariableVector( const csmp::Index& var_key,
                                               std::vector<Node<3U>*>::const_iterator begin,
                                               std::vector<Node<3U>*>::const_iterator end,
                                               INTERFACE_SIDE );
};

} // csmp

#endif // ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
