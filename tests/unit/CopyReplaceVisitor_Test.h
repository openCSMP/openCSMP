#ifndef COPYREPLACEVISITOR_TEST_H
#define COPYREPLACEVISITOR_TEST_H

#include "Test.h"


namespace csmp{

  template<size_t> class Model;

/// PL Nov 2010
class CopyReplaceVisitor_Test : public Test
{
public:
  explicit CopyReplaceVisitor_Test( Model<3U>* model) : model_( model ) {}
  ~CopyReplaceVisitor_Test();
  virtual void run();

private:

  template<class V>
  void testNodes( Model<3>* model, const V& value, const char* propertyName );

  template<class V>
  void testElements( Model<3>* model, const V& value, const char* propertyName );

  template<class V>
  void testElementIntegrationPoints( Model<3>* model, const V& value, const char* propertyName );

  Model<3U>* model_;
};


} // csmp

#endif // COPYREPLACEVISITOR_TEST_H
