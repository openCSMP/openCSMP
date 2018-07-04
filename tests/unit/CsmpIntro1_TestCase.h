#ifndef CSMP_INTRO_1_TEST_CASE_H
#define CSMP_INTRO_1_TEST_CASE_H

#include "Test.h"

namespace csmp {

class  CsmpIntro1_TestCase : public Test {
  public:
    CsmpIntro1_TestCase( bool verbose ) : verbose_(verbose) {}
  
    virtual void run();
  
  private:
    const bool verbose_;
};

} // csmp


#endif // CSMP_INTRO_1_TEST_CASE_H
