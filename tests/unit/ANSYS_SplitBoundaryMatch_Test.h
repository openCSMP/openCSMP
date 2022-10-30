#ifndef ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
#define ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H

#include "Test.h"
#include "CSMP_definitions.h"

namespace csmp {

class Index;
template<uint32_t> class Node;
template<uint32_t> class SplitBoundary;

/**
    Tests high level utility function that matches node-matched elements
    across SplitBoundaries created in ANSYS.
    
    @todo insert _test() macros into the two test functions deciding on quantities to analyse
    @todo function that matches node coordinates across interfaces fails in some cases although the interfaces themselves were matched before; fix
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
                                               
    void TestThatManifoldNodesAreCollocated( const SplitBoundary<3U>& );
};

} // csmp

#endif // ANSYS_SPLIT_BOUNDARY_MATCH_TEST_H
