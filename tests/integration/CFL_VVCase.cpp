#include "CFL_VVCase.h"
#include "ANSYS_Model3D.h"
#include "Model.h"
#include "LinearSolver.h"
#include "PropertyHandle.h"
#include "PDE_Integrator.h"
#include "ScalarVariable.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "PL_Utilities.h"
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "CSMP_highLevelUtilities.h"
#include "StencilProcessor.h"
#include "BrooksCorey.h"
#include "VTU_Interface.h"
#include "Timer.hpp"

#include <fstream>


using namespace std;

namespace csmp
{


CFL_TestCase::CFL_TestCase()
  {
  }

CFL_TestCase::~CFL_TestCase()
  {
  }

struct CFL_TestData
{
  double time, CFLmultiplier, mobilityRatio, phi, velX, frontSat, frontX, lenghtInFlowDir, frontSat_fw_Derivative, testCFL, anisotropicCFL, difference;
  string testName, testVersion, inflowBoundary;
  bool passed;
};

void CFL_TestCase::run()
  {

    enum{DIM=3};
    //const
    const  double satOil(1.), tolCFL(.1);

    //test container
    vector<CFL_TestData> tests;


    // test 1.U
    CFL_TestData data11;
    data11.time = 1.; data11.CFLmultiplier = 1.; data11.mobilityRatio = 3.; data11.velX = 0.01; data11.phi = 0.25;
    data11.frontSat = 0.4; data11.frontX = 2.0; data11.frontSat_fw_Derivative = 3.918396790; data11.lenghtInFlowDir = 0.1;
    data11.testName = "SmallRectangle"; data11.testVersion = "U"; data11.inflowBoundary = "RIGHT";
    //Test( data11, satOil, tolCFL );

    // test 1.V
    CFL_TestData data12;
    data12.time = 1.; data12.CFLmultiplier = 1.; data12.mobilityRatio = 3.; data12.velX = 0.01; data12.phi = 0.25;
    data12.frontSat = 0.4; data12.frontX = 2.0; data12.frontSat_fw_Derivative = 3.918396790; data12.lenghtInFlowDir = 0.1;
    data12.testName = "SmallRectangleQuad"; data12.testVersion = "V"; data12.inflowBoundary = "RIGHT";
    Test( data12, satOil, tolCFL );


  }//run

void CFL_TestCase::Test( CFL_TestData& data, double satOil, double tolCFL )
{
  ANSYS_Model3D model( data.testName.data(), "CSMP-transport-variables.txt" );
  printModelDimensions( model );
  SetModel( model, data.mobilityRatio, data.velX, data.phi );
  SetFront( model, data.frontSat, data.frontX );
  OutputVTU( model, (data.testName + data.testVersion).c_str(), (long)(data.time) );
  TwoPhaseImplicitNodeCenteredFVTransport<3U,StencilProcessor> TPINCVT( "Model", model );
  BrooksCorey<3> brooksCorey( model.Database(), "permeability", "viscosity oil", "viscosity water",
                              "density oil", "density water",
                              "brooks corey parameter", "entry pressure",
                              "saturation water","residual saturation oil", "residual saturation water");
  TPINCVT.CFL_Multiplier( data.CFLmultiplier );
  data.anisotropicCFL = TPINCVT.AnisotropicCourantIncrement( brooksCorey, 1.0E+05 );


  data.testCFL =  data.lenghtInFlowDir * data.phi / ( data.velX * data.frontSat_fw_Derivative );
  data.difference = data.testCFL - data.anisotropicCFL;
  data.passed = ( abs(data.difference) <= tolCFL * data.testCFL );

  cout << "\n\nCFL TestCase test results:"
       <<   "\n========================\n";
  cout << "Test CFL: " << data.testCFL << "\nAnisotropicCourantIncrement: " << data.anisotropicCFL << "\n";

  _equal( data.anisotropicCFL, data.testCFL, data.testCFL*tolCFL );

}//Test


void CFL_TestCase::OutputVTU( Model<3>& model,
               const char* fileName,
               long time )
{
  VTU_Interface<3> vtu( model );
  model.ExtrapolateElementToNodeProperty( "velocity", "nodal velocity" );
  list<string> props;
  props.push_back( "saturation water" );
  props.push_back( "nodal velocity" );
  vtu.OutputDataToVTU( fileName, "saturation water","Model", static_cast<size_t>(time) );
}//OutputVTU


void CFL_TestCase::SetModel( Model<3>& model,
                             double mobilityRatio,
                             double velX,
                             double phi )
{
    const ScalarVariable viscosityWater( PLAIN, 1.0E-03 );
    const ScalarVariable viscosityOil( PLAIN, viscosityWater()*mobilityRatio );
    const ScalarVariable permeability( PLAIN, 1.0E-12 );
    const ScalarVariable densityWater( PLAIN, 1.0E+03 );
    const ScalarVariable densityOil( PLAIN, 1.0E+03 );
    const ScalarVariable saturationWater( PLAIN, 0. );
    const ScalarVariable saturationOil( PLAIN, 1.0 );
    const ScalarVariable porosity( PLAIN, phi );
    const ScalarVariable conductivity( PLAIN, permeability() / viscosityWater() );
    const ScalarVariable lambda( PLAIN, 3.0 );
    const ScalarVariable residualWater( PLAIN, 0.0 );
    const ScalarVariable residualOil( PLAIN, 2.0E-01 );
    const ScalarVariable entryPressure( PLAIN, 1.0E+03 );
    const VectorVariable<3> velocity( PLAIN, PLAIN, PLAIN, velX, 0., 0. );

    model.InputPropertyValue( "fluid volume source", makeScalar( PLAIN, 0. ) );
    model.InputPropertyValue( "viscosity water",     viscosityWater );
    model.InputPropertyValue( "viscosity oil",       viscosityOil );
    model.InputPropertyValue( "density water",       densityWater );
    model.InputPropertyValue( "density oil",         densityOil );
    model.InputPropertyValue( "saturation water",    saturationWater );
    model.InputPropertyValue( "saturation oil",      saturationOil );
    model.InputPropertyValue( "residual saturation water",    residualWater );
    model.InputPropertyValue( "residual saturation oil",      residualOil );
    model.InputPropertyValue( "porosity",            porosity );
    model.InputPropertyValue( "permeability",        permeability );
    model.InputPropertyValue( "entry pressure",      entryPressure );
    model.InputPropertyValue( "brooks corey parameter",        lambda );
    model.InputPropertyValue( "conductivity",        conductivity );
    model.InputPropertyValue( "velocity",            velocity );
    model.InputPropertyValue( "nodal fluid volume source", makeScalar(PLAIN, 0.) );
} //SetModel


void CFL_TestCase::SetFront( Model<3>& model,
                               double satWater,
                               double xFront
                           )
  {
    Point<3> min, max;
    minMaxXYZ( min, max, model );

    Index satOilKey( model.Database().StorageKey("saturation oil" ) );
    Index satWaterKey( model.Database().StorageKey("saturation water" ) );
    Point<3> satFront( min );
    satFront[0] += xFront;

    for( vector<Node<3>*>::const_iterator node(model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
    {
      if( (*node)->x() < satFront[0] ){
        const ScalarVariable satOilScalar( PLAIN, (1 - satWater) );
        (*node)->Store( satOilKey, satOilScalar ) ;
        const ScalarVariable satWaterScalar( PLAIN, satWater );
        (*node)->Store( satWaterKey, satWaterScalar );
      }
    }
  } //SetFront


} //csmp
