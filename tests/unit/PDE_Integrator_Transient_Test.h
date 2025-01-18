#ifndef CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H
#define CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

class PDE_Integrator_Transient_Test : public Test {
  public:
   void run();
   void outVector(const std::vector<double>& vector, std::string file);
   /*void oldIntegrate(PDE_Integrator_UoM_Mock<2U, Region>& pde, Region<2U>& domain);
   void newIntegrate(PDE_Integrator_UoM_Mock<2U, Region>& pde, Region<2U>& domain);*/

  };

} // end namespace csmp

#endif /* CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H */
