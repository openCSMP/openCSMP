#ifndef ALGORITHM_TEST_H
#define ALGORITHM_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "VSet.h"

#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "PDE_Integrator.h"
#include "Integral_dNT_op_dN_dV.h"
#include "Integral_NT_op_N_dV.h"
#include "CSMP_highLevelUtilities.h"

#include "LinearSolver.h"

#include "PropertyHandle.h"
namespace csmp {

class PDE_Integrator_Test : public Test {
public:
    PDE_Integrator_Test();
    ~PDE_Integrator_Test();
    void run();
    //void exchangeSolverTest();
    void sameSolverTest();

private:

  Model<2U>* sg_;
  PDE_Integrator<2U,Region>* alg_;
  Integral_dNT_op_dN_dV<2U,Element<2U> >* stiff_;
  Integral_NT_op_N_dV<2U,Element<2U> >* source_;
};

} // end namespace csmp

#endif
