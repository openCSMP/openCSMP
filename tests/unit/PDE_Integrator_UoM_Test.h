#ifndef PDE_INTEGRATOR_UOM_TEST_H
#define PDE_INTEGRATOR_UOM_TEST_H

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

  class PDE_Integrator_UoM_Test : public Test {
  public:
    explicit PDE_Integrator_UoM_Test(std::string);
    ~PDE_Integrator_UoM_Test();
    void run();

  private:
    Model<2U>* model                                        = nullptr;
    PDE_Integrator_UoM_Mock<2U, Region>* pde_validate       = nullptr;         // for validating
    PDE_Integrator_UoM_Mock<2U, Region>* pde_test           = nullptr;         // for testing

    /*
    =================> LHS:
    NumIntegral_BT_D_B_dV
    NumIntegral_BT_D_op_dV
    NumIntegral_dNT_dN_dV
    NumIntegral_dNT_mixed_op_dN_dV
    NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV
    NumIntegral_dNT_op_dN_dV                                      1
    NumIntegral_dNT_op_dN_dV_NT_v_dN_dV
    NumIntegral_dNT_op_dN_NT_op_dop_dN_dV
    NumIntegral_DNT_op_DN_NT_v_DN_dV
    NumIntegral_NT_dNi_dV
    NumIntegral_NT_dNi_dV_sc
    NumIntegral_NT_lhsop_N_dV
    NumIntegral_PT_lhsop_P_dV        

    ==================> RHS:
    MathOperatorRHS
    NumIntegral_DNi_rhsop_dV
    NumIntegral_dNT_op_dV                                         1
    NumIntegral_DNT_rhsop_DN_dV
    NumIntegral_DNT_v_dV
    NumIntegral_NT_mixed_op_dNi_dV
    NumIntegral_NT_op_dNi_dV
    NumIntegral_NT_op_N_dS
    NumIntegral_NT_op_N_dV                                        1
    NumIntegral_NT_op1_op2_dNi_dV
    NumIntegral_op_NT_dN_orthogonal_dV
    NumIntegral_op_NT_N_dV
    NumIntegral_op_PT_P_dV
    NumIntegral_PT_op_dS
    NumIntegral_PT_op_dV
    NumIntegral_PT_op_P_dV
    PointSource_rhsop                                             1
    */


    // pressure equation
    NumIntegral_dNT_op_dN_dV<2U, Element<2U> >* pressureLHS   = nullptr;
    NumIntegral_NT_op_N_dV<2U, Element<2U> >* sourceVolume    = nullptr;
    PointSource_rhsop<2U, Element<2U>>*       sourcePoint     = nullptr;
    NumIntegral_dNT_op_dV<2U, Element<2U>>*   gravityTerm     = nullptr;
    VelocityAndVolumeFlux<2U, Element<2U>> *  fluid_velocity  = nullptr;

    // thermal equation
    NumIntegral_dNT_op_dN_dV<2U, Element<2U> >* temperatureLHS = nullptr;
    NumIntegral_NT_op_N_dV<2U, Element<2U> >* sourceHeat       = nullptr;

  private:
    void Reset();
    void TestSingleVariable();
    void TestTwoScalarVariables();
    void TestOutputSingleVariable();
  };

} // end namespace csmp

#endif
