#ifndef PDE_INTEGRATOR_CRM_TEST_H
#define PDE_INTEGRATOR_CRM_TEST_H

// FE algorithm
#include "PDE_Integrator_CRM.h"

#include "CSMP_definitions.h"
#include "Test.h"

// CSMP model
#include "ANSYS_Model2D.h"
#include "Model.h"
#include "VSet.h"

// PDE operators building the FE algorithm
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_dNi_dV_sc.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "PointSource_rhsop.h"

#include "LinearSolver.h"

// interfaces
#include "VTU_Interface.h"

// utility functions
#include "CSMP_highLevelUtilities.h"

// visitors
#include "StabilizationParameterVisitor.h"
namespace csmp {

class PDE_Integrator_CRM_Test : public Test {
public:
  PDE_Integrator_CRM_Test();
  ~PDE_Integrator_CRM_Test();
  void run();

private:
  Model<2U>* sg_;
  PDE_Integrator<2U, Region>* alg_;
  NumIntegral_dNT_op_dN_dV<2U, Element<2U> >* stiff_;
  NumIntegral_NT_op_N_dV<2U, Element<2U> >* source_;

  void scaleRegion( Model<2U>& sg, double64 scale_factor );
  void constructVelocityVector( Model<2U>& mdl );
  void assignFluxToPointSource( Model<2U>& mdl, const char* flux );
};

} // end namespace csmp

#endif
