#ifndef GENERIC_SINGLETON_TEST_H
#define GENERIC_SINGLETON_TEST_H

#include "Test.h"
#include "GenericSingleton.h"

namespace csmp{

class GenericSingleton_Test : public Test
{
public:
  virtual void run();
};


 

class TestSingleton : public GenericSingleton<TestSingleton>
{
  friend class GenericSingleton<TestSingleton>;
 
  private:
    TestSingleton() : i(3) {}

  public:
    int i;
};

} // csmp

#endif // GENERIC_SINGLETON_TEST_H
