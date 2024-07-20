#ifndef CONSTRAINTPOINTTONODEPROPERTYVISITOR_TEST_H
#define CONSTRAINTPOINTTONODEPROPERTYVISITOR_TEST_H

#include "Test.h"

namespace csmp{

  template<uint32_t> class Model;

class IntegrationPointToNodePropertyVisitor_Test : public Test
{
public:
  explicit IntegrationPointToNodePropertyVisitor_Test( Model<3U>& model );
  virtual void run();

private:
  Model<3U>& model_;

};


} // csmp

#endif // CONSTRAINTPOINTTONODEPROPERTYVISITOR_TEST_H
