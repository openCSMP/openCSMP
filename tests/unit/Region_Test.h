#ifndef REGION_TEST_H
#define REGION_TEST_H

#include "Test.h"



namespace csmp
{
 // class Index;
  class Region_Test: public Test
  {
  public:
      Region_Test( bool verbose=false );
      ~Region_Test();
      void run();
  private:
      const bool verbose_;
  };
}
#endif // REGION_TEST_H
