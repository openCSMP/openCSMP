#ifndef EXPLICIT_ADVECTION_2D_VVCASE_H
#define EXPLICIT_ADVECTION_2D_VVCASE_H

#include "Test.h"

namespace csmp {

  /** Explicit Advection test!: More description to come.

*/
  class ExplicitAdvection2D_VVCase : public Test {
  public:
    explicit ExplicitAdvection2D_VVCase();
    explicit ExplicitAdvection2D_VVCase(const char* prefix);
    ~ExplicitAdvection2D_VVCase();
    void run(); // runs all the tests for the class

  };
}

#endif // EXPLICIT_ADVECTION_2D_VVCASE_H
