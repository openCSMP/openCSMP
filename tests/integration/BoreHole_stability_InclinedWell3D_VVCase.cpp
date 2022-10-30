#include "ANSYS_Model3D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "SAMG_Solver.h"
#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include"NumIntegral_BT_op_dV.h"
#include"NumIntegral_PT_op_dV.h"
#include "StressesAndStrains.h"
#include "ExtractTensorVariableComponent.h"
#include "BoreHole_stability_InclinedWell3D_VVCase.h"
#include "Test.h"
#include "VSet.h"
#include "VSetConverter.h"
#include "ANSYS_Interface.h"
#include "MohrCoulombFailure_Visitor.h"
//#include "RegionInterface.h"
//#include "Point.h"
#include "VTK_Interface.h"
#include"ScalarVariable.h"
#include "PropertyHandle.h"
#include "NodeManifold.h"


using namespace std;


namespace csmp
{

BoreHole_stability_InclinedWell3D_VVCase::BoreHole_stability_InclinedWell3D_VVCase(const char* prefix)

{
  this->setName("BoreHole_stability_InclinedWell3D_VVCase");
  prefix_=prefix;
}

void BoreHole_stability_InclinedWell3D_VVCase::run()
{

// some constants

//overbalanced drilling

// Normal fault regime

// SHmax=57000000 Pa
// SHmin=38000000 Pa
// Sv=67000000 Pa
// Pmud= 32000000 Pa
// Pp= 28000000 Pa hydrostatic pressure for a depth of 2800 m

// Strike-slip fault regime

// SHmax=85000000 Pa
// SHmin=45000000 Pa
// Sv=67000000 Pa
// Pmud= 32000000 Pa --> 9.1 ppg
// Pp= 28000000 Pa

// Reverse fault regime

// SHmax=105000000 Pa
// SHmin=85000000 Pa
// Sv=67000000 Pa
// Pmud= 32000000 Pa
// Pp= 28000000 Pa

// Rw=0.5 m

const double SHmax( 57000000. );
const double SHmin( 38000000.);
const double Pmud( 32000000. );
const double SV( 67000000. );
const double Pp( 28000000. );

const ScalarVariable zeroScalar( PLAIN, 0. );
const double displacement( 0. );
const VectorVariable<3U> zeroVector( PLAIN, 0. );
//    const VectorVariable<3U> zeroVectorDirichlet( DIRICH, 0. );
const VectorVariable<3U> DisplacementVectorBOTTOM(PLAIN, DIRICH ,PLAIN,0., displacement, 0.); // it can not move in the Y direction
const VectorVariable<3U> Rho_g( PLAIN, PLAIN ,PLAIN, 0. , - 2.4e+4 , 0.);


// establishing model & output facility
 string input_file_name(prefix_);

 //  string modelName( "BoreHoleInclined_stability3D" );

  ANSYS_Model3D model( input_file_name.data(),this->getName().c_str(),(this->getName()+".txt").c_str() );

  VTU_Interface<3U> vtu( model );
  vtu.OmitZeroInFileName(true);

  VTK_Interface<3U> vtk;


  // output properties & initial output
  list<string> outputProps;
  //    outputProps.push_back( "mean stress" );
   model.InputPropertyValue( "mean stress", zeroScalar );
   model.InputPropertyValue( "displacement", zeroVector );
   model.InputPropertyValue( "force", zeroVector );
   model.InputPropertyValue( "Neumann stress", zeroVector );
   model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN, Pp) );


   vtu.OutputDataToVTU( ((string)(input_file_name.data())+"UnLoaded").data(), outputProps, "Model", static_cast<int>(0) );

    // model configuration

    model.InputPropertyValue( "Poisson's ratio", makeScalar( PLAIN, 0.29 ) ); //carbonate
    model.InputPropertyValue( "Young's modulus", makeScalar( PLAIN, 24.0e+9 ) ); //carbonate
    model.InputPropertyValue("gravity force", Rho_g );

    Boundary<3U>& rightBoundary( model.Boundary( "RIGHT" ) );
    Boundary<3U>& leftBoundary( model.Boundary( "LEFT" ) );
    Boundary<3U>& topBoundary( model.Boundary( "TOP" ) );
    Boundary<3U>& bottomBoundary( model.Boundary( "BOTTOM" ) );
    Boundary<3U>& frontBoundary( model.Boundary( "FRONT" ) );
    Boundary<3U>& backBoundary( model.Boundary( "BACK" ) );
    Region<3U>& wellbore( model.Region("BOREHOLE") );

    Region<3U>& TOPwellbore( model.Region("TOPBOREHOLE") );
    Region<3U>& BOTTOMwellbore( model.Region("BOTTOMBOREHOLE") );

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
//  bottomBoundary.InputPropertyValue( "displacement", zeroVectorDirichlet ); // BOTTOM Boundary can not move in the 3 directions direction
    wellbore.InputPropertyValue( "fluid pressure", Pressure );
    bottomBoundary.InputPropertyValue("displacement",DisplacementVectorBOTTOM);// BOTTOM Boundary can not move in the Y direction

    TOPwellbore.InputPropertyValue( "fluid pressure", Pressure );
    BOTTOMwellbore.InputPropertyValue( "fluid pressure", Pressure );

    vtu.OutputDataToVTU( "borehole_fluid pressure", "fluid pressure",    "Model", 0 );

    vtu.OutputDataToVTU( "borehole_Neumann stress", "Neumann stress",    "Model", 0 );

    // setting up & solving linear elasticity fea problem

    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_napproach(2);
    SAMG_Solver solver(&settings);
    PDE_Integrator<3U,Element> deformation( solver );
    #else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    PDE_Integrator<3U,Element> deformation( solver );
    #endif

    PT_op<3U> bforces( model.Database(), "force", "displacement" );

    NumIntegral_BT_D_B_dV<3U> stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
    NumIntegral_PT_op_dV<3U>  bodyforce(  model.Database(), "gravity force", "displacement");
    NumIntegral_PT_op_dV<3U>  AppliedStress( model.Database(), "Neumann stress", "displacement");
    NumIntegral_BT_op_dV<3U>  WellBorePressure( model.Database(),"fluid pressure", "displacement");


    deformation.Add( &stiffness );
    deformation.Add( &bforces );
    deformation.Add( &bodyforce );
    deformation.Add( &AppliedStress );
    deformation.Add( &WellBorePressure );

    StressesAndStrains<3U>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", true,true);

    deformation.AddPostProcess( &postpro );

    model.Apply( deformation );

//    Failure mohr coulomb

    printRangeOfVariable(model,"strain1");

    cout << model.Database().StorageKey("stress") << endl;


     MohrCoulombFailure_Visitor<3U> failureTest(model,true /* positive compressive stress convention */, false/* verbose mode */);

     model.InputPropertyValue( "cohesion", makeScalar( PLAIN, 1400000. ) ); // 14 MPa --->UCS=50 MPa (unconfined compression strenght)
     model.InputPropertyValue( "friction angle", makeScalar( PLAIN, 31. ) ); //friction angle 31 degree -->0.6 friction coefficient

     model.Accept(failureTest);


// apply resulting displacement
      model.MoveNodeCoordinatesBy("displacement");


//    vtk.OutputDataToVTK(model,"Model_strain1","strain1",0);
//    vtk.OutputDataToVTK(model,"Model_strain2","strain2",0);
//    vtk.OutputDataToVTK(model,"Model_strain3","strain3",0);
//    vtk.OutputDataToVTK(model,"Model_sigma1","sigma1",0);
//    vtk.OutputDataToVTK(model,"Model_sigma2","sigma2",0);
//    vtk.OutputDataToVTK(model,"Model_sigma3","sigma3",0);
    vtk.OutputDataToVTK(model,"Model_meanstress","mean stress",0);
//    vtk.OutputDataToVTK(model,"Model_displacement","displacement",0);
//    vtk.OutputDataToVTK(model,"Model_failure","failure",0);
    vtk.OutputDataToVTK(model,"Model_stress","stress",0);
//    vtk.OutputDataToVTK(model,"Model_strain","strain",0);
//    vtk.OutputDataToVTK(model,"Model_stress_node","stress node",0);
//    vtk.OutputDataToVTK(model,"Model_strain_node","strain node",0);
    vtk.OutputDataToVTK(model,"Model_failure01","failure01",0);
//    vtk.OutputDataToVTK(model,"Model_principal stress","principal stress",0);



      return;
     }

} // csmp
