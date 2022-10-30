#ifndef SNEDDON_CRACK_COUPLED_3D_VVCASE_H
#define SNEDDON_CRACK_COUPLED_3D_VVCASE_H

#include "Test.h"
#include "GlobalVerbose.h"

#include "CSMP_highLevelUtilities.h"
#include "ANSYS_Model3D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
//#include "NumIntegral_BT_C_B_dV.h"
//#include "StressAndStrainOutput.h"
#include "ExtractTensorVariableComponent.h"
#include "NumIntegral_PT_op_dS.h"
#include "ModelSubDomain.h"
#include "Fracture.h"
#include <stdio.h>
#include <functional>
#include <cmath>

//Interrellations for conductivity operator
#include "InterFaceFractureVisitor.h"

//Numerical Integrals Lubricaiton Equation
#include "NumIntegral_dNT_mixed_op_dN_dV.h"           // to interpolate the aperture
#include "NumIntegral_dNT_op_dN_dV.h"                 // to interpolate the aperture cubed  (try with quadratic base functions...)
#include "NumIntegral_NT_lhsop_N_dV.h"                // Mass Matrix LHS
#include "NumIntegral_SetRHS_to_Zero.h"               // Zero right hand side
#include "NumIntegral_NT_op_N_dV.h"                   // Lumped Mass Matrix RHS or Source term with AccumulateLater()
#include "NumIntegral_PT_op_dV.h"
#include "PT_op.h"                                    //b force num int

//Coupled Case
#include "NumIntegral_BT_D_B_dV.h"                    //stiffness
#include "NumIntegral_PT_op_dS.h"                     //rhs traction terms (split)
//#include "NumIntegral_PT_n_N_dS.h"                    //lhs pressure traction term (splitboundary)
#include "NumIntegral_dNT_op_dN_dV.h"                 //conductance (fluid pressure)
#include "PointSource_rhsop.h"                        //fluid rhs point source (on all nodes of element)
//#include "NumIntegral_NT_n_P_dS.h"                    // Aperture change LHS
//#include "NumIntegral_NT_n_rhs_P_dS.h"                // Aperture change RHS
//#include "NumIntegral_PT_n_N_dS.h"

//#include "RecoveryBasedOnDisplacement.h"

namespace csmp
  {

  class SneddonCrackCoupled3D_VVCase: public Test
    {
    public:
      SneddonCrackCoupled3D_VVCase(const char* prefix);

      void set_min_avg_error(double epsilon){ min_avg_error = epsilon; };

      virtual void run();

    private:

     double min_avg_error = 8.0;


  };




  } // csmp

#endif
