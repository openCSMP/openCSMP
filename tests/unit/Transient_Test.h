#ifndef TRANSIENT_TEST_H
#define TRANSIENT_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "VSet.h"

#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "PDE_Integrator.h"
#include "PDE_Integrator_UoM_Mock.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "CSMP_highLevelUtilities.h"

#include "LinearSolver.h"

#include "PropertyHandle.h"


namespace csmp {

  class Transient_Test : public Test {
  public:
   void run();
   void outVector(const std::vector<double>& vector, std::string file);
   /*void oldIntegrate(PDE_Integrator_UoM_Mock<2U, Region>& pde, Region<2U>& domain);
   void newIntegrate(PDE_Integrator_UoM_Mock<2U, Region>& pde, Region<2U>& domain);*/

  };

} // end namespace csmp

#endif
