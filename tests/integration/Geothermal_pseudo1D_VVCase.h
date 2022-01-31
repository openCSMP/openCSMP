#ifndef GEOTHERMAL_PSEUDO1D_VVCASE_H
#define GEOTHERMAL_PSEUDO1D_VVCASE_H

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
//#include "NumIntegral_dNT_op_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "PointSource_rhsop.h"
#include "VelocityAndVolumeFlux.h"

// visitors
//#include "ComputeGravityTermVisitor.h"
#include "ThermalVisitor.h"
#include "SourceVisitor.h"
#include "ConductivityVisitor.h"
#include "PropertyAtPointVisitor.h"

// utilities and monitoring
#include "InputDataManager.h"
#include "ComputationalSettings.h"

using namespace std;
using namespace csmp;


namespace csmp
  {

  /** Geothermal pseudo 1D Test Case
  =================================
  Mesh:       Linear Triangles
  Model:      2000x1000 m 2D
  Test:       Temperature
  BC:         constant pressure, constant temperature
  Criterion:
  =================================
  */
  class Geothermal_pseudo1D_VVCase : public Test
    {
      public:
        Geothermal_pseudo1D_VVCase(const char* prefix);
        virtual void run();
      private:
        void outputToVTU( Model<3U>& model,string model_name, const list<string>& props, size_t timestep );
        void ComputeMassConductivity (Model<3U>& model);
		bool Compare (Model<3U>& model, string file);
    };

  } // csmp

#endif // GEOTHERMAL_PSEUDO1D_VVCASE_H
