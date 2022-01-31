#ifndef CSMP_PDE_INTEGRATOR_TEST_H
#define CSMP_PDE_INTEGRATOR_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "VSet.h"

#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "CSMP_highLevelUtilities.h"

#include "LinearSolver.h"

#include "PropertyHandle.h"

namespace csmp {

class PDE_Integrator_Test : public Test {
public:
    /// builds a test model that gets assigned to sg_
    PDE_Integrator_Test();
    ~PDE_Integrator_Test();
    
    // TODO: test is not developed yet, nothing gets tested for
    virtual void run();
    
    void SolveMatrixEquationWithSAMG();

private:

  Model<2U>* sg_                       = nullptr;
  PDE_Integrator<2U,Region>* alg_      = nullptr;
  NumIntegral_dNT_op_dN_dV<2U>* stiff_ = nullptr;
  NumIntegral_NT_op_N_dV<2U>* source_  = nullptr;
};

} // end namespace csmp

#endif /* CSMP_PDE_INTEGRATOR_TEST_H */
