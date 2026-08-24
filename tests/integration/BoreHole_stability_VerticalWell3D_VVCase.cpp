#include "ANSYS_Model3D.h"
#include "Model.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "NumIntegral_PT_op_dV.h"
#include "StressesAndStrains.h"
#include "ExtractTensorVariableComponent.h"
#include "BoreHole_stability_VerticalWell3D_VVCase.h"
#include "Test.h"
#include "VSet.h"
#include "VSetConverter.h"
#include "ANSYS_Interface.h"
#include "MohrCoulombFailure_Visitor.h"
#include "VTK_Interface.h"
#include"ScalarVariable.h"
#include "PropertyHandle.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "SAMG_Exception.h"
#else
#include "LinearSolver.h"
#endif

using namespace std;

namespace csmp {

BoreHole_stability_VerticalWell3D_VVCase::BoreHole_stability_VerticalWell3D_VVCase(const char* prefix)
  {
    this->setName("BoreHole_stability_VerticalWell3D_VVCase");
    prefix_=prefix;
  }



void BoreHole_stability_VerticalWell3D_VVCase::run()
  {
     // Rw = 13.33
  const double SHmax( 90.0e+06 );
  const double SHmin( 51.5e+06);
  const double Pmud( 31.5e+06 );
  const double SV( 88.2e+06 );
  const double Pp( 31.5e+06);
  const ScalarVariable zeroScalar( PLAIN, 0. );
  const double displacement( 0. );
  const VectorVariable<3U> zeroVector( PLAIN, 0. );
  //const VectorVariable<3U> zeroVectorDirichlet( DIRICH, 0. );
  const VectorVariable<3U> DisplacementVectorBOTTOM(PLAIN, DIRICH ,PLAIN,0., displacement, 0.); // it can not move in the Y direction
  const VectorVariable<3U> Rho_g( PLAIN, PLAIN ,PLAIN, 0. , - 2.4e+4 , 0.);

  // establishing model & output facility
  string input_file_name(prefix_);

  //  string modelName( "BoreHoleVertical_stability3D" );
  ANSYS_Model3D model( input_file_name.data(),this->getName().c_str(), "BoreHole_stability_VerticalWell3D_VVCase-variables.txt" );

  VTU_Interface<3U> vtu( model );
  vtu.OmitZeroInFileName(true);

  // output properties & initial output
  list<string> outputProps;

  //    outputProps.push_back( "mean stress" );
  model.InputPropertyValue( "mean stress", zeroScalar );
  model.InputPropertyValue( "displacement", zeroVector );
  model.InputPropertyValue( "force", zeroVector );
  model.InputPropertyValue( "Neumann stress", zeroVector );
  model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN, Pp) );
  vtu.OutputDataToVTU( ((string)(input_file_name.data())+"UnLoaded").data(), outputProps, "Model", static_cast<int>(0) );

  // model configuration (rock is ZONE0), hole volume is WELL
  model.InputPropertyValue( "Poisson's ratio", makeScalar( PLAIN, 0.29 ) ); //carbonate
  model.InputPropertyValue( "Young's modulus", makeScalar( PLAIN, 24.0e+9 ) ); //carbonate
  model.InputPropertyValue("gravity force", Rho_g );
  Boundary<3U>& rightBoundary( model.Boundary( "RIGHT" ) );
  Boundary<3U>& leftBoundary( model.Boundary( "LEFT" ) );
  Boundary<3U>& topBoundary( model.Boundary( "TOP" ) );
  Boundary<3U>& bottomBoundary( model.Boundary( "BOTTOM" ) );
  Boundary<3U>& frontBoundary( model.Boundary( "FRONT" ) );
  Boundary<3U>& backBoundary( model.Boundary( "BACK" ) );
  
  // BOUNDARY_WELL is the corresponding surface
  Region<3U>& wellbore( model.Region("WELL") );

  const VectorVariable<3U> SressRight( DIRICH, DIRICH,DIRICH, -SHmax, 0.,0. );
  const VectorVariable<3U> SressLeft( DIRICH, DIRICH,DIRICH, SHmax, 0.,0. );
  const VectorVariable<3U> SressTop( DIRICH, DIRICH,DIRICH, 0.,-SV,0. );
  const VectorVariable<3U> Sressfront( DIRICH, DIRICH,DIRICH,  0.,0.,-SHmin );
  const VectorVariable<3U> Sressback( DIRICH, DIRICH,DIRICH, 0.,0.,SHmin);
  const ScalarVariable Pressure( DIRICH,Pmud);

  rightBoundary.InputPropertyValue( "Neumann stress", SressRight );
  leftBoundary.InputPropertyValue( "Neumann stress", SressLeft );
  topBoundary.InputPropertyValue( "Neumann stress", SressTop );
  frontBoundary.InputPropertyValue( "Neumann stress", Sressfront );
  backBoundary.InputPropertyValue( "Neumann stress", Sressback );

  wellbore.InputPropertyValue( "fluid pressure", Pressure );

  bottomBoundary.InputPropertyValue("displacement",DisplacementVectorBOTTOM);// BOTTOM Boundary can not move in the Y direction

  // input properties & initial output
  list<string> inputProps = {"Young's modulus", "Poisson's ratio", "displacement",
                             "gravity force", "Neumann stress", "mean stress", "force", "fluid pressure", "displacement"};


  printRangeOfVariable( model, "Neumann stress");
  // will not print Neumann stress as it is a FACE variable
  vtu.OutputDataToVTU( "BoreHole_stability_VerticalWell3D_VVCase_input", inputProps, "ZONE_0", 0 );


  // setting up & solving linear elasticity fea problem
  #ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings settings;
  settings.Set_napproach(2);
  SAMG_Solver solver(&settings);
  PDE_Integrator<3U,Element> deformation( solver  );
  #else
  CSMP_DEFAULT_LINEAR_SOLVER solver;
  PDE_Integrator<3U,Element> deformation( solver );
  #endif

  PT_op<3U> bforces( model.Database(), "force", "displacement" );
  NumIntegral_BT_D_B_dV<3U>  stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
  NumIntegral_PT_op_dV<3U>   bodyforce(  model.Database(), "gravity force", "displacement");
  NumIntegral_PT_op_dV<3U>   AppliedStress( model.Database(), "Neumann stress", "displacement");
  NumIntegral_BT_op_dV<3U>   WellBorePressure( model.Database(),"fluid pressure", "displacement");

  deformation.Add( &stiffness );
  deformation.Add( &bforces );
  deformation.Add( &bodyforce );
  deformation.Add( &AppliedStress );
  deformation.Add( &WellBorePressure );

  constexpr bool output_principal_vectors{true}, extrapolate_results_to_nodes{false};
  StressesAndStrains<3U>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", output_principal_vectors, extrapolate_results_to_nodes );
  deformation.AddPostProcess( &postpro );

  model.Apply( deformation );

  // Failure mohr coulomb
  MohrCoulombFailure_Visitor<3U> failureTest( model, ELEMENT ); //friction angle 31 degree -->0.6 friction coefficient

  model.InputPropertyValue( "cohesion", makeScalar( PLAIN, 3000000. ) ); // 14 MPa --->UCS=50 MPa (unconfined compression strenght)
  model.InputPropertyValue( "friction angle", makeScalar( PLAIN, 31. ) ); //friction angle 31 degree -->0.6 friction coefficient

  model.Accept(failureTest);

  // apply resulting displacement
  model.MoveNodeCoordinatesBy("displacement");

  // collect all output variables into a single set and write one VTU file
  const set<string> output_variables {
      "strain1",
      "strain2",
      "strain3",
      "sigma1",
      "sigma2",
      "sigma3",
      "mean stress",
      "displacement",
      "failure",
      "stress",
      "strain",
      "failure01",
      "stress node",
      "strain node"
  };

  vtu.OutputDataToVTU( "BoreholeStability", output_variables, "Model", 0 );}

} // csmp
