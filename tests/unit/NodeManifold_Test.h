#ifndef NODE_MANIFOLD_TEST_H
#define NODE_MANIFOLD_TEST_H

#include "Test.h"

namespace csmp {

class NodeManifold_Test : public Test {
  public:
    virtual void run();
    
    static const bool verbose_ = false;
    
  private:
    void Test_Classify2D();
    void Test_Classify3D();
};

} // csmp

#endif // NODE_MANIFOLD_TEST_H
