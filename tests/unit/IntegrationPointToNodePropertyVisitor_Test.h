#ifndef CSMP_INTEGRATION_POINT_TO_NODE_PROPERTY_VISITOR_TEST_H
#define CSMP_INTEGRATION_POINT_TO_NODE_PROPERTY_VISITOR_TEST_H

#include "Test.h"

namespace csmp{

  template<uint32_t> class Model;

class IntegrationPointToNodePropertyVisitor_Test : public Test
{
public:
  
  // uses model prism test
  virtual void run();

private:
  
  void Test_NodeToIntegrationPointInterpolation( Model<3U>& model );
  void Test_IntegrationPointToNodePropertyVisitor( Model<3U>& model );
  void Test_CoordinateInterpolationToBaryCenter( Model<3U>& model );

};


} // csmp

#endif // CSMP_INTEGRATION_POINT_TO_NODE_PROPERTY_VISITOR_TEST_H
