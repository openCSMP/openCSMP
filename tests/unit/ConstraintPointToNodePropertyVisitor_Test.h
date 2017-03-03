#ifndef CONSTRAINTPOINTTONODEPROPERTYVISITOR_TEST_H
#define CONSTRAINTPOINTTONODEPROPERTYVISITOR_TEST_H

#include "Test.h"

namespace csmp{

  template<size_t> class Model;

class ConstraintPointToNodePropertyVisitor_Test : public Test
{
public:
  explicit ConstraintPointToNodePropertyVisitor_Test( Model<3U>& model );
  virtual void run();

private:
  Model<3U>& model_;

};


} // csmp

#endif // CONSTRAINTPOINTTONODEPROPERTYVISITOR_TEST_H
