#ifndef VISITOR_TESTS_H
#define VISITOR_TESTS_H

#include "Test.h"
#include "TestSuite.h"

namespace csmp{

  class ANSYS_Model3D;

/// PL Nov 2010
class Visitor_TestSuite : public Test
{
public:

  explicit Visitor_TestSuite( TestSuite& suite ) : suite_( suite ), free_( true ){}
  ~Visitor_TestSuite();

  virtual void run();

private:
  TestSuite& suite_;
  ANSYS_Model3D* model_;
  bool free_;
};


} // csmp

#endif // VISITOR_TESTS_H
