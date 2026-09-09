#ifndef CSMP_COPY_REPLACE_VISITOR_TEST_H
#define CSMP_COPY_REPLACE_VISITOR_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

  template<uint32_t> class Model;

/// PL Nov 2010
class CopyReplaceVisitor_Test : public Test
{
public:
  CopyReplaceVisitor_Test() {}

  virtual void run();

private:

  template<class V>
  void TestNodes( Model<2>& model, const V& value, const char* propertyName );

  template<class V>
  void TestElements( Model<2>& model, const V& value, const char* propertyName );

  template<class V>
  void TestElementIntegrationPoints( Model<2>& model, const V& value, const char* propertyName );
};


} // csmp

#endif // CSMP_COPY_REPLACE_VISITOR_TEST_H
