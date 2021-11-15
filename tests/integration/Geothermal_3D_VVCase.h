#ifndef GEOTHERMAL_3D_VVCASE_H
#define GEOTHERMAL_3D_VVCASE_H

#include "Test.h"

#include "Model.h"
#include "VTU_Interface.h"
#include "ANSYS_Model3D.h"
#include "CSMP_highLevelUtilities.h"
#include "LinearSolver.h"
#include "LUdcmp_Solver.h"
#include "CSMP_definitions.h"
#include "GlobalVerbose.h"

// finite volumes 
#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"

// finite elements
#include "PDE_Integrator.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "PointSource_rhsop.h"
#include "VelocityAndVolumeFlux.h"

// visitors
#include "ComputeGravityTermVisitor.h"
#include "ThermalVisitor.h"
#include "SourceVisitor.h"
#include "ConductivityVisitor.h"

// utilities and monitoring
#include "InputDataManager.h"
#include "ComputationalSettings.h"

using namespace std;
using namespace csmp;


namespace csmp
  {

  /** Geothermal 3D Test Case
  =================================
  Mesh:       Linear Triangles
  Model:       m 3D
  Test:       Temperature
  BC:         constant pressure, constant temperature
  Criterion:
  =================================
  */
  class Geothermal_3D_VVCase : public Test
    {
      public:
        Geothermal_3D_VVCase(const char* prefix);
        virtual void run();
      private:
        void outputToVTU( Model<3U>& model,std::string model_name, const list<string>& props, size_t timestep );
        void ComputeMassConductivity (Model<3U>& model);
        void ComputeMassGravityTerm (Model<3U>& model);
    };

  } // csmp

#endif // GEOTHERMAL_3D_VVCASE_H
