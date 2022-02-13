#include "TPINCFVT3D_VVCase.h"
#include "ANSYS_Model3D.h"
#include "LinearSolver.h"
#include "PropertyHandle.h"
#include "PDE_Integrator.h"
#include "ScalarVariable.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "PL_Utilities.h"
#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "StencilProcessor.h"
#include "BrooksCorey.h"
#include "FourarLenormand.h"
#include "LinearTwoPhaseModel.h"
#include "VTU_Interface.h"
#include "Timer.hpp"

#include <fstream>


using namespace std;

namespace csmp
{


TPINCFVT3D_VVCase::TPINCFVT3D_VVCase()
  {
  }
TPINCFVT3D_VVCase::TPINCFVT3D_VVCase( const char* prefix )
{
  name_="TPINCFVT3D_VVCase";
  prefix_=prefix;
}
TPINCFVT3D_VVCase::~TPINCFVT3D_VVCase()
  {
  }

struct TPINCFVT3D_TestData
{
  TPINCFVT3D_TestData() : time(1.), CFLmultiplier(1.), mobilityRatio(1.), porosity(1.), velX(1.), analyticFront(1.), difference(1.), numericFront(1.), timeElapsed(1.),
    timeMultiplier(1.) {}
  enum SATFCTS{BC, FL, LIN};
  SATFCTS relpermModel;
  double time, CFLmultiplier, mobilityRatio, porosity, velX, analyticFront, difference, numericFront,
         timeElapsed, timeMultiplier;
  string testName, testVersion, inflowBoundary;
  bool passed;
};

void TPINCFVT3D_VVCase::run()
  {
    string input_file_name(prefix_);

    enum{DIM=3};
    //const
    const  double tolSat(0.05), satOil(1.), tolFront(.1);
    //for each test
    double time(0.), CFLmultiplier(1.), mobilityRatio(0.), velX(0.), analyticFront(0.);

    //test container
    vector<TPINCFVT3D_TestData> tests;


    /// TestCase for Brooks-Corey two-phase flow
    /// running a 1.0 x 1.0 x 1.0 m cube model at different transport velocities

    // test 1.A // failing CFL multiplier
    TPINCFVT3D_TestData data11;
    data11.time = 100.; data11.CFLmultiplier = 1.; data11.mobilityRatio = 1.; data11.velX = 0.001; data11.analyticFront = 0.626;
    data11.testName = "UnitCube"; data11.testVersion = "A"; data11.porosity = 2.5E-01; data11.relpermModel = TPINCFVT3D_TestData::BC; data11.inflowBoundary = "RIGHT";
    //Test( data11, tolSat, satOil, tolFront );
    //tests.push_back(data11);

    // test 1.B // failing CFL multiplier
    TPINCFVT3D_TestData data12;
    data12.time = 10.; data12.CFLmultiplier = 1.; data12.mobilityRatio = 1.; data12.velX = 0.01; data12.analyticFront = 0.626;
    data12.testName = "UnitCube"; data12.testVersion = "B"; data12.porosity = 2.5E-01; data12.relpermModel = TPINCFVT3D_TestData::BC; data12.inflowBoundary = "RIGHT";
    //Test( data12, tolSat, satOil, tolFront );
    //tests.push_back(data12);

    // test 1.C
    TPINCFVT3D_TestData data13;
    data13.time = 100.; data13.CFLmultiplier = 0.1; data13.mobilityRatio = 1.; data13.velX = 0.001; data13.analyticFront = 0.626;
    data13.testName = "UnitCube"; data13.testVersion = "C"; data13.porosity = 2.5E-01; data13.relpermModel = TPINCFVT3D_TestData::BC; data13.inflowBoundary = "RIGHT";
    //Test( data13, tolSat, satOil, tolFront );
    //tests.push_back(data13);

    // test 1.D
    TPINCFVT3D_TestData data14;
    data14.time = 10.; data14.CFLmultiplier = 0.1; data14.mobilityRatio = 1.; data14.velX = 0.01; data14.analyticFront = 0.626;
    data14.testName = "UnitCube"; data14.testVersion = "D"; data14.porosity = 2.5E-01; data14.relpermModel = TPINCFVT3D_TestData::BC; data14.inflowBoundary = "RIGHT";
    //Test( data14, tolSat, satOil, tolFront );
    //tests.push_back(data14);

    /// TestCase for Brooks-Corey two-phase flow
    /// running a 40.0 x 10.0 m rectangular slab model with different viscosity ratios and at different transport velocities

    // test 2.A // failing CFL multiplier
    TPINCFVT3D_TestData data21;
    data21.time = 5000.; data21.CFLmultiplier = 1.; data21.mobilityRatio = 1.; data21.velX = 0.001; data21.analyticFront = 31.2;
    data21.testName = "Square40x10m"; data21.testVersion = "A"; data21.porosity = 2.5E-01; data21.relpermModel = TPINCFVT3D_TestData::BC; data21.inflowBoundary = "LEFT";
    //Test( data21, tolSat, satOil, tolFront );
    //tests.push_back(data21);

    // test 2.B // failing CFL multiplier
    TPINCFVT3D_TestData data22;
    data22.time = 500.; data22.CFLmultiplier = 1.; data22.mobilityRatio = 1.; data22.velX = 0.01; data22.analyticFront = 31.2;
    data22.testName = "Square40x10m"; data22.testVersion = "B"; data22.porosity = 2.5E-01; data22.relpermModel = TPINCFVT3D_TestData::BC; data22.inflowBoundary = "LEFT";
    //Test( data22, tolSat, satOil, tolFront );
    //tests.push_back(data22);

    // test 2.C
    TPINCFVT3D_TestData data23;
    data23.time = 5000.; data23.CFLmultiplier = 0.1; data23.mobilityRatio = 1.; data23.velX = 0.001; data23.analyticFront = 31.2;
    data23.testName = "Square40x10m"; data23.testVersion = "C"; data23.porosity = 2.5E-01; data23.relpermModel = TPINCFVT3D_TestData::BC; data23.inflowBoundary = "LEFT";
    //Test( data23, tolSat, satOil, tolFront );
    //tests.push_back(data23);

    // test 2.D
    TPINCFVT3D_TestData data24;
    data24.time = 500.; data24.CFLmultiplier = 0.1; data24.mobilityRatio = 1.; data24.velX = 0.01; data24.analyticFront = 31.2;
    data24.testName = "Square40x10m"; data24.testVersion = "D"; data24.porosity = 2.5E-01; data24.relpermModel = TPINCFVT3D_TestData::BC; data24.inflowBoundary = "LEFT";
    //Test( data24, tolSat, satOil, tolFront );
    //tests.push_back(data24);

    // test 2.E // failing CFL multiplier
    TPINCFVT3D_TestData data25;
    data25.time = 5000.; data25.CFLmultiplier = 1.; data25.mobilityRatio = 3.; data25.velX = 0.001; data25.analyticFront = 35.0;
    data25.testName = "Square40x10m"; data25.testVersion = "E"; data25.porosity = 2.5E-01; data25.relpermModel = TPINCFVT3D_TestData::BC; data25.inflowBoundary = "LEFT";
    //Test( data25, tolSat, satOil, tolFront );
    //tests.push_back(data25);

    // test 2.F // failing CFL multiplier
    TPINCFVT3D_TestData data26;
    data26.time = 500.; data26.CFLmultiplier = 1.; data26.mobilityRatio = 3.; data26.velX = 0.01; data26.analyticFront = 35.0;
    data26.testName = "Square40x10m"; data26.testVersion = "F"; data26.porosity = 2.5E-01; data26.relpermModel = TPINCFVT3D_TestData::BC; data26.inflowBoundary = "LEFT";
    //Test( data26, tolSat, satOil, tolFront );
    //tests.push_back(data26);

    // test 2.G
    TPINCFVT3D_TestData data27;
    data27.time = 5000.; data27.CFLmultiplier = 0.1; data27.mobilityRatio = 3.; data27.velX = 0.001; data27.analyticFront = 35.0;
    data27.testName = "Square40x10m"; data27.testVersion = "G"; data27.porosity = 2.5E-01; data27.relpermModel = TPINCFVT3D_TestData::BC; data27.inflowBoundary = "LEFT";
    //Test( data27, tolSat, satOil, tolFront );
    //tests.push_back(data27);

    // test 2.H
    TPINCFVT3D_TestData data28;
    data28.time = 500.; data28.CFLmultiplier = 0.1; data28.mobilityRatio = 3.; data28.velX = 0.01; data28.analyticFront = 35.0;
    data28.testName = "Square40x10m"; data28.testVersion = "H"; data28.porosity = 2.5E-01; data28.relpermModel = TPINCFVT3D_TestData::BC; data28.inflowBoundary = "LEFT";
    //Test( data28, tolSat, satOil, tolFront );
    //tests.push_back(data28);

    // test 2.I // failing CFL multiplier
    TPINCFVT3D_TestData data29;
    data29.time = 5000.; data29.CFLmultiplier = 1.; data29.mobilityRatio = 5.; data29.velX = 0.001; data29.analyticFront = 37.316;
    data29.testName = "Square40x10m"; data29.testVersion = "I"; data29.porosity = 2.5E-01; data29.relpermModel = TPINCFVT3D_TestData::BC; data29.inflowBoundary = "LEFT";
    //Test( data29, tolSat, satOil, tolFront );
    //tests.push_back(data29);

    // test 2.J // failing CFL multiplier
    TPINCFVT3D_TestData data30;
    data30.time = 500.; data30.CFLmultiplier = 1.; data30.mobilityRatio = 5.; data30.velX = 0.01; data30.analyticFront = 37.316;
    data30.testName = "Square40x10m"; data30.testVersion = "J"; data30.porosity = 2.5E-01; data30.relpermModel = TPINCFVT3D_TestData::BC; data30.inflowBoundary = "LEFT";
    //Test( data30, tolSat, satOil, tolFront );
    //tests.push_back(data30);

    // test 2.K
    TPINCFVT3D_TestData data31;
    data31.time = 5000.; data31.CFLmultiplier = 0.1; data31.mobilityRatio = 5.; data31.velX = 0.001; data31.analyticFront = 37.316;
    data31.testName = "Square40x10m"; data31.testVersion = "K"; data31.porosity = 2.5E-01; data31.relpermModel = TPINCFVT3D_TestData::BC; data31.inflowBoundary = "LEFT";
    //Test( data31, tolSat, satOil, tolFront );
    //tests.push_back(data31);

    // test 2.M
    TPINCFVT3D_TestData data32;
    data32.time = 500.; data32.CFLmultiplier = 0.1; data32.mobilityRatio = 5.; data32.velX = 0.01; data32.analyticFront = 37.316;
    data32.testName = "Square40x10m"; data32.testVersion = "M"; data32.porosity = 2.5E-01; data32.relpermModel = TPINCFVT3D_TestData::BC; data32.inflowBoundary = "LEFT";
    //Test( data32, tolSat, satOil, tolFront );
    //tests.push_back(data32);

    /// TestCase for Brooks-Corey two-phase flow
    /// running a 40.0 x 10.0 m rectangular slab model of gradually changing cell sizes at different transport velocities

    // test 3.N // failing CFL multiplier
    TPINCFVT3D_TestData data33;
    data33.time = 5000.; data33.CFLmultiplier = 1.; data33.mobilityRatio = 3.; data33.velX = 0.001; data33.analyticFront = 35.;
    data33.testName = "Square40x10GradientGrid"; data33.testVersion = "N"; data33.porosity = 2.5E-01; data33.relpermModel = TPINCFVT3D_TestData::BC; data33.inflowBoundary = "LEFT";
    //Test( data33, tolSat, satOil, tolFront );
    //tests.push_back(data33);

    // test 3.O // failing CFL multiplier
    TPINCFVT3D_TestData data34;
    data34.time = 500.; data34.CFLmultiplier = 1.; data34.mobilityRatio = 3.; data34.velX = 0.01; data34.analyticFront = 35.;
    data34.testName = "Square40x10GradientGrid"; data34.testVersion = "O"; data34.porosity = 2.5E-01; data34.relpermModel = TPINCFVT3D_TestData::BC; data34.inflowBoundary = "LEFT";
    //Test( data34, tolSat, satOil, tolFront );
    //tests.push_back(data34);

    // test 3.P
    TPINCFVT3D_TestData data35;
    data35.time = 5000.; data35.CFLmultiplier = 0.1; data35.mobilityRatio = 3.; data35.velX = 0.001; data35.analyticFront = 35.;
    data35.testName = "Square40x10GradientGrid"; data35.testVersion = "P"; data35.porosity = 2.5E-01; data35.relpermModel = TPINCFVT3D_TestData::BC; data35.inflowBoundary = "LEFT";
    //Test( data35, tolSat, satOil, tolFront );
    //tests.push_back(data35);

    // test 3.R
    TPINCFVT3D_TestData data36;
    data36.time = 500.; data36.CFLmultiplier = 0.1; data36.mobilityRatio = 3.; data36.velX = 0.01; data36.analyticFront = 35.;
    data36.testName = "Square40x10GradientGrid"; data36.testVersion = "R"; data36.porosity = 2.5E-01; data36.relpermModel = TPINCFVT3D_TestData::BC; data36.inflowBoundary = "LEFT";
    //Test( data36, tolSat, satOil, tolFront );
    ////tests.push_back(data36);

    // test 3.S
    TPINCFVT3D_TestData data37;
    data37.time = 5000.; data37.CFLmultiplier = 0.05; data37.mobilityRatio = 3.; data37.velX = 0.001; data37.analyticFront = 35.;
    data37.testName = "Square40x10GradientGrid"; data37.testVersion = "S"; data37.porosity = 2.5E-01; data37.relpermModel = TPINCFVT3D_TestData::BC; data37.inflowBoundary = "LEFT";
    //Test( data37, tolSat, satOil, tolFront );
    //tests.push_back(data37);

    // test 3.T
    TPINCFVT3D_TestData data38;
    data38.time = 500.; data38.CFLmultiplier = 0.05; data38.mobilityRatio = 3.; data38.velX = 0.01; data38.analyticFront = 35.;
    data38.testName = "Square40x10GradientGrid"; data38.testVersion = "T"; data38.porosity = 2.5E-01; data38.relpermModel = TPINCFVT3D_TestData::BC; data38.inflowBoundary = "LEFT";
    //Test( data38, tolSat, satOil, tolFront );
    //tests.push_back(data38);

    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab model of quadrilateral elements

    // test 4.U
    TPINCFVT3D_TestData data41;
    data41.time = 100.; data41.CFLmultiplier = 0.1 ; data41.mobilityRatio = 3.; data41.velX = 0.01; data41.analyticFront = 7.;
    data41.testName = "SmallRectangleQuad"; data41.testVersion = "U"; data41.porosity = 2.5E-01; data41.relpermModel = TPINCFVT3D_TestData::BC; data41.inflowBoundary = "RIGHT";
    //Test( data41, tolSat, satOil, tolFront );
    //tests.push_back(data41);

    // test 4.V // failing CFL multiplier
    TPINCFVT3D_TestData data42;
    data42.time = 100.; data42.CFLmultiplier = 1.; data42.mobilityRatio = 3.; data42.velX = 0.01; data42.analyticFront = 7.;
    data42.testName = "SmallRectangleQuad"; data42.testVersion = "V"; data42.porosity = 2.5E-01; data42.relpermModel = TPINCFVT3D_TestData::BC; data42.inflowBoundary = "RIGHT";
    //Test( data42, tolSat, satOil, tolFront );
    //tests.push_back(data42);


    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 10.0 x 1.0 x 0.25 m rectangular slab model

    // test 5.W
    TPINCFVT3D_TestData data51;
    data51.time = 50.; data51.CFLmultiplier = 0.1; data51.mobilityRatio = 2.; data51.velX = 0.01; data51.analyticFront = 6.;
    data51.testName = "SmallRectangle"; data51.testVersion = "W"; data51.porosity = 2.5E-01; data51.relpermModel = TPINCFVT3D_TestData::LIN; data51.inflowBoundary = "RIGHT";
    //Test( data51, tolSat, satOil, tolFront );
    //tests.push_back(data51);

    // test 5.X
    TPINCFVT3D_TestData data52;
    data52.time = 50.; data52.CFLmultiplier = 1.0; data52.mobilityRatio = 2.; data52.velX = 0.01; data52.analyticFront = 6.;
    data52.testName = "SmallRectangle"; data52.testVersion = "X"; data52.porosity = 2.5E-01; data52.relpermModel = TPINCFVT3D_TestData::LIN; data52.inflowBoundary = "RIGHT";
    //Test( data52, tolSat, satOil, tolFront );
    //tests.push_back(data52);

    // test 5.Y // failing CFL multiplier
    TPINCFVT3D_TestData data53;
    data53.time = 50.; data53.CFLmultiplier = 10.; data53.mobilityRatio = 2.; data53.velX = 0.01; data53.analyticFront = 6.;
    data53.testName = "SmallRectangle"; data53.testVersion = "Y"; data53.porosity = 2.5E-01; data53.relpermModel = TPINCFVT3D_TestData::LIN; data53.inflowBoundary = "RIGHT";
    //Test( data53, tolSat, satOil, tolFront );
    //tests.push_back(data53);


/////////////////////////////


    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab (fracture) model accounting for realistic fracture aperture conditions

    // test 6.1.1.BC-A
    TPINCFVT3D_TestData data611;
    data611.time = 50.; data611.CFLmultiplier = 0.01; data611.mobilityRatio = 3.; data611.velX = 0.001; data611.analyticFront = 35.;
    data611.testName = "Square40x10m"; data611.testVersion = "BC-A-FRAC-CFL0.01_"; data611.porosity = 2.5E-03; data611.relpermModel = TPINCFVT3D_TestData::BC; data611.inflowBoundary = "LEFT";
    //Test( data611, tolSat, satOil, tolFront );
    //tests.push_back(data611);

    // test 6.1.2.BC-B
    TPINCFVT3D_TestData data612;
    data612.time = 50.; data612.CFLmultiplier = 0.05; data612.mobilityRatio = 3.; data612.velX = 0.001; data612.analyticFront = 35.;
    data612.testName = "Square40x10m"; data612.testVersion = "BC-B-FRAC-CFL0.05_"; data612.porosity = 2.5E-03; data612.relpermModel = TPINCFVT3D_TestData::BC; data612.inflowBoundary = "LEFT";
    //Test( data612, tolSat, satOil, tolFront );
    //tests.push_back(data612);

    // test 6.1.3.BC-C
    TPINCFVT3D_TestData data613;
    data613.time = 50.; data613.CFLmultiplier = 0.1; data613.mobilityRatio = 3.; data613.velX = 0.001; data613.analyticFront = 35.;
    data613.testName = "Square40x10m"; data613.testVersion = "BC-C-FRAC-CFL0.1_"; data613.porosity = 2.5E-03; data613.relpermModel = TPINCFVT3D_TestData::BC; data613.inflowBoundary = "LEFT";
    Test( data613, tolSat, satOil, tolFront );
    tests.push_back(data613);

    // test 6.1.4.BC-D // failing CFL multiplier
    TPINCFVT3D_TestData data614;
    data614.time = 50.; data614.CFLmultiplier = 0.33; data614.mobilityRatio = 3.; data614.velX = 0.001; data614.analyticFront = 35.;
    data614.testName = "Square40x10m"; data614.testVersion = "BC-D-FRAC-CFL0.33_"; data614.porosity = 2.5E-03; data614.relpermModel = TPINCFVT3D_TestData::BC; data614.inflowBoundary = "LEFT";
    //Test( data614, tolSat, satOil, tolFront );
    //tests.push_back(data614);

    // test 6.1.5.BC-E // failing CFL multiplier
    TPINCFVT3D_TestData data615;
    data615.time = 50.; data615.CFLmultiplier = 0.66; data615.mobilityRatio = 3.; data615.velX = 0.001; data615.analyticFront = 35.;
    data615.testName = "Square40x10m"; data615.testVersion = "BC-E-FRAC-CFL0.66_"; data615.porosity = 2.5E-03; data615.relpermModel = TPINCFVT3D_TestData::BC; data615.inflowBoundary = "LEFT";
    //Test( data615, tolSat, satOil, tolFront );
    //tests.push_back(data615);

    // test 6.1.6.BC-F // failing CFL multiplier
    TPINCFVT3D_TestData data616;
    data616.time = 50.; data616.CFLmultiplier = 1.0; data616.mobilityRatio = 3.; data616.velX = 0.001; data616.analyticFront = 35.;
    data616.testName = "Square40x10m"; data616.testVersion = "BC-F-FRAC-CFL1.0_"; data616.porosity = 2.5E-03; data616.relpermModel = TPINCFVT3D_TestData::BC; data616.inflowBoundary = "LEFT";
    //Test( data616, tolSat, satOil, tolFront );
    //tests.push_back(data616);


    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab (fracture) model accounting for realistic fracture aperture conditions
    /// Discontinous transport using time multiplier 0.1

    // test 6.1.7.BC-A
    TPINCFVT3D_TestData data617;
    data617.time = 50.; data617.timeMultiplier = 0.1; data617.CFLmultiplier = 0.01; data617.mobilityRatio = 3.; data617.velX = 0.001; data617.analyticFront = 35.;
    data617.testName = "Square40x10m"; data617.testVersion = "BC-A-FRAC-CFL0.01_TMP0.1_"; data617.porosity = 2.5E-03; data617.relpermModel = TPINCFVT3D_TestData::BC; data617.inflowBoundary = "LEFT";
    //Test( data617, tolSat, satOil, tolFront );
    //tests.push_back(data617);

    // test 6.1.8.BC-B
    TPINCFVT3D_TestData data618;
    data618.time = 50.; data618.timeMultiplier = 0.1; data618.CFLmultiplier = 0.05; data618.mobilityRatio = 3.; data618.velX = 0.001; data618.analyticFront = 35.;
    data618.testName = "Square40x10m"; data618.testVersion = "BC-B-FRAC-CFL0.05_TMP0.1_"; data618.porosity = 2.5E-03; data618.relpermModel = TPINCFVT3D_TestData::BC; data618.inflowBoundary = "LEFT";
    //Test( data618, tolSat, satOil, tolFront );
    //tests.push_back(data618);

    // test 6.1.9.BC-C
    TPINCFVT3D_TestData data619;
    data619.time = 50.; data619.timeMultiplier = 0.1; data619.CFLmultiplier = 0.1; data619.mobilityRatio = 3.; data619.velX = 0.001; data619.analyticFront = 35.;
    data619.testName = "Square40x10m"; data619.testVersion = "BC-C-FRAC-CFL0.1_TMP0.1_"; data619.porosity = 2.5E-03; data619.relpermModel = TPINCFVT3D_TestData::BC; data619.inflowBoundary = "LEFT";
    Test( data619, tolSat, satOil, tolFront );
    tests.push_back(data619);

    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab (fracture) model accounting for realistic fracture aperture conditions
    /// Discontinous transport using time multiplier 0.01

    // test 6.2.1.BC-A
    TPINCFVT3D_TestData data621;
    data621.time = 50.; data621.timeMultiplier = 0.01; data621.CFLmultiplier = 0.01; data621.mobilityRatio = 3.; data621.velX = 0.001; data621.analyticFront = 35.;
    data621.testName = "Square40x10m"; data621.testVersion = "BC-A-FRAC-CFL0.01_TMP0.01_"; data621.porosity = 2.5E-03; data621.relpermModel = TPINCFVT3D_TestData::BC; data621.inflowBoundary = "LEFT";
    //Test( data621, tolSat, satOil, tolFront );
    //tests.push_back(data621);

    // test 6.2.2.BC-B
    TPINCFVT3D_TestData data622;
    data622.time = 50.; data622.timeMultiplier = 0.01; data622.CFLmultiplier = 0.05; data622.mobilityRatio = 3.; data622.velX = 0.001; data622.analyticFront = 35.;
    data622.testName = "Square40x10m"; data622.testVersion = "BC-B-FRAC-CFL0.05_TMP0.01_"; data622.porosity = 2.5E-03; data622.relpermModel = TPINCFVT3D_TestData::BC; data622.inflowBoundary = "LEFT";
    //Test( data622, tolSat, satOil, tolFront );
    //tests.push_back(data622);

    // test 6.2.3.BC-C
    TPINCFVT3D_TestData data623;
    data623.time = 50.; data623.timeMultiplier = 0.01; data623.CFLmultiplier = 0.1; data623.mobilityRatio = 3.; data623.velX = 0.001; data623.analyticFront = 35.;
    data623.testName = "Square40x10m"; data623.testVersion = "BC-C-FRAC-CFL0.1_TMP0.01_"; data623.porosity = 2.5E-03; data623.relpermModel = TPINCFVT3D_TestData::BC; data623.inflowBoundary = "LEFT";
    Test( data623, tolSat, satOil, tolFront );
    tests.push_back(data623);


    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab model

    // test 6.3.1.BC-A
    TPINCFVT3D_TestData data631;
    data631.time = 500.; data631.CFLmultiplier = 0.01; data631.mobilityRatio = 3.; data631.velX = 0.01; data631.analyticFront = 35.;
    data631.testName = "Square40x10m"; data631.testVersion = "BC-A-CFL0.01_"; data631.porosity = 2.5E-01; data631.relpermModel = TPINCFVT3D_TestData::BC; data631.inflowBoundary = "LEFT";
    //Test( data631, tolSat, satOil, tolFront );
    //tests.push_back(data631);

    // test 6.3.2.BC-B
    TPINCFVT3D_TestData data632;
    data632.time = 500.; data632.CFLmultiplier = 0.05; data632.mobilityRatio = 3.; data632.velX = 0.01; data632.analyticFront = 35.;
    data632.testName = "Square40x10m"; data632.testVersion = "BC-B-CFL0.05_"; data632.porosity = 2.5E-01; data632.relpermModel = TPINCFVT3D_TestData::BC; data632.inflowBoundary = "LEFT";
    //Test( data632, tolSat, satOil, tolFront );
    //tests.push_back(data632);

    // test 6.3.3.BC-C
    TPINCFVT3D_TestData data633;
    data633.time = 500.; data633.CFLmultiplier = 0.1; data633.mobilityRatio = 3.; data633.velX = 0.01; data633.analyticFront = 35.;
    data633.testName = "Square40x10m"; data633.testVersion = "BC-C-CFL0.1_"; data633.porosity = 2.5E-01; data633.relpermModel = TPINCFVT3D_TestData::BC; data633.inflowBoundary = "LEFT";
    Test( data633, tolSat, satOil, tolFront );
    tests.push_back(data633);

    // test 6.3.4.BC-D // failing CFL multiplier
    TPINCFVT3D_TestData data634;
    data634.time = 500.; data634.CFLmultiplier = 0.33; data634.mobilityRatio = 3.; data634.velX = 0.01; data634.analyticFront = 35.;
    data634.testName = "Square40x10m"; data634.testVersion = "BC-D-CFL0.33_"; data634.porosity = 2.5E-01; data634.relpermModel = TPINCFVT3D_TestData::BC; data634.inflowBoundary = "LEFT";
    //Test( data634, tolSat, satOil, tolFront );
    //tests.push_back(data634);

    // test 6.3.5.BC-E // failing CFL multiplier
    TPINCFVT3D_TestData data635;
    data635.time = 500.; data635.CFLmultiplier = 0.66; data635.mobilityRatio = 3.; data635.velX = 0.01; data635.analyticFront = 35.;
    data635.testName = "Square40x10m"; data635.testVersion = "BC-E-CFL0.66_"; data635.porosity = 2.5E-01; data635.relpermModel = TPINCFVT3D_TestData::BC; data635.inflowBoundary = "LEFT";
    //Test( data635, tolSat, satOil, tolFront );
    //tests.push_back(data635);

    // test 6.3.6.BC-F // failing CFL multiplier
    TPINCFVT3D_TestData data636;
    data636.time = 500.; data636.CFLmultiplier = 1.0; data636.mobilityRatio = 3.; data636.velX = 0.01; data636.analyticFront = 35.;
    data636.testName = "Square40x10m"; data636.testVersion = "BC-F-CFL1.0_"; data636.porosity = 2.5E-01; data636.relpermModel = TPINCFVT3D_TestData::BC; data636.inflowBoundary = "LEFT";
    //Test( data636, tolSat, satOil, tolFront );
    //tests.push_back(data636);

    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab model
    /// Discontinous transport using time multiplier 0.1

    // test 6.3.7.BC-A
    TPINCFVT3D_TestData data637;
    data637.time = 500.; data637.timeMultiplier = 0.1; data637.CFLmultiplier = 0.01; data637.mobilityRatio = 3.; data637.velX = 0.01; data637.analyticFront = 35.;
    data637.testName = "Square40x10m"; data637.testVersion = "BC-A-CFL0.01_TMP0.1_"; data637.porosity = 2.5E-01; data637.relpermModel = TPINCFVT3D_TestData::BC; data637.inflowBoundary = "LEFT";
    //Test( data637, tolSat, satOil, tolFront );
    //tests.push_back(data637);

    // test 6.3.8.BC-B
    TPINCFVT3D_TestData data638;
    data638.time = 500.; data638.timeMultiplier = 0.1; data638.CFLmultiplier = 0.05; data638.mobilityRatio = 3.; data638.velX = 0.01; data638.analyticFront = 35.;
    data638.testName = "Square40x10m"; data638.testVersion = "BC-B-CFL0.05_TMP0.1_"; data638.porosity = 2.5E-01; data638.relpermModel = TPINCFVT3D_TestData::BC; data638.inflowBoundary = "LEFT";
    //Test( data638, tolSat, satOil, tolFront );
    //tests.push_back(data638);

    // test 6.3.9.BC-C
    TPINCFVT3D_TestData data639;
    data639.time = 500.; data639.timeMultiplier = 0.1; data639.CFLmultiplier = 0.1; data639.mobilityRatio = 3.; data639.velX = 0.01; data639.analyticFront = 35.;
    data639.testName = "Square40x10m"; data639.testVersion = "BC-C-CFL0.1_TMP0.1_"; data639.porosity = 2.5E-01; data639.relpermModel = TPINCFVT3D_TestData::BC; data639.inflowBoundary = "LEFT";
    Test( data639, tolSat, satOil, tolFront );
    tests.push_back(data639);

    /// TestCase for Brooks-Corey two-phase flow
    /// running a 10.0 x 1.0 m rectangular slab model
    /// Discontinous transport using time multiplier 0.01

    // test 6.4.1.BC-A
    TPINCFVT3D_TestData data641;
    data641.time = 500.; data641.timeMultiplier = 0.01; data641.CFLmultiplier = 0.01; data641.mobilityRatio = 3.; data641.velX = 0.01; data641.analyticFront = 35.;
    data641.testName = "Square40x10m"; data641.testVersion = "BC-A-CFL0.01_TMP0.01_"; data641.porosity = 2.5E-01; data641.relpermModel = TPINCFVT3D_TestData::BC; data641.inflowBoundary = "LEFT";
    //Test( data641, tolSat, satOil, tolFront );
    //tests.push_back(data641);

    // test 6.4.2.BC-B
    TPINCFVT3D_TestData data642;
    data642.time = 500.; data642.timeMultiplier = 0.01; data642.CFLmultiplier = 0.05; data642.mobilityRatio = 3.; data642.velX = 0.01; data642.analyticFront = 35.;
    data642.testName = "Square40x10m"; data642.testVersion = "BC-B-CFL0.05_TMP0.01_"; data642.porosity = 2.5E-01; data642.relpermModel = TPINCFVT3D_TestData::BC; data642.inflowBoundary = "LEFT";
    //Test( data642, tolSat, satOil, tolFront );
    //tests.push_back(data642);

    // test 6.4.3.BC-C
    TPINCFVT3D_TestData data643;
    data643.time = 500.; data643.timeMultiplier = 0.01; data643.CFLmultiplier = 0.1; data643.mobilityRatio = 3.; data643.velX = 0.01; data643.analyticFront = 35.;
    data643.testName = "Square40x10m"; data643.testVersion = "BC-C-CFL0.1_TMP0.01_"; data643.porosity = 2.5E-01; data643.relpermModel = TPINCFVT3D_TestData::BC; data643.inflowBoundary = "LEFT";
    Test( data643, tolSat, satOil, tolFront );
    tests.push_back(data643);


    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 40.0 x 10.0 m rectangular slab (fracture) model accounting for realistic fracture aperture conditions

    // test 7.1.1.LIN-A
    TPINCFVT3D_TestData data711;
    data711.time = 20.; data711.CFLmultiplier = 0.1; data711.mobilityRatio = 2.; data711.velX = 0.001; data711.analyticFront = 25.0;
    data711.testName = "Square40x10m"; data711.testVersion = "LIN-A-FRAC-CFL0.1_"; data711.porosity = 2.5E-03; data711.relpermModel = TPINCFVT3D_TestData::LIN; data711.inflowBoundary = "LEFT";
    //Test( data711, tolSat, satOil, tolFront );
    //tests.push_back(data711);

    // test 7.1.2.LIN-B
    TPINCFVT3D_TestData data712;
    data712.time = 20.; data712.CFLmultiplier = 0.5; data712.mobilityRatio = 2.; data712.velX = 0.001; data712.analyticFront = 25.0;
    data712.testName = "Square40x10m"; data712.testVersion = "LIN-B-FRAC-CFL0.5_"; data712.porosity = 2.5E-03; data712.relpermModel = TPINCFVT3D_TestData::LIN; data712.inflowBoundary = "LEFT";
    //Test( data712, tolSat, satOil, tolFront );
    //tests.push_back(data712);

    // test 7.1.3.LIN-C
    TPINCFVT3D_TestData data713;
    data713.time = 20.; data713.CFLmultiplier = 1.0; data713.mobilityRatio = 2.; data713.velX = 0.001; data713.analyticFront = 25.0;
    data713.testName = "Square40x10m"; data713.testVersion = "LIN-C-FRAC-CFL1.0_"; data713.porosity = 2.5E-03; data713.relpermModel = TPINCFVT3D_TestData::LIN; data713.inflowBoundary = "LEFT";
    Test( data713, tolSat, satOil, tolFront );
    tests.push_back(data713);

    // test 7.1.4.LIN-D // failing CFL multiplier
    TPINCFVT3D_TestData data714;
    data714.time = 20.; data714.CFLmultiplier = 2.5; data714.mobilityRatio = 2.; data714.velX = 0.001; data714.analyticFront = 25.0;
    data714.testName = "Square40x10m"; data714.testVersion = "LIN-D-FRAC-CFL2.5_"; data714.porosity = 2.5E-03; data714.relpermModel = TPINCFVT3D_TestData::LIN; data714.inflowBoundary = "LEFT";
    //Test( data714, tolSat, satOil, tolFront );
    //tests.push_back(data714);

    // test 7.1.5.LIN-E // failing CFL multiplier
    TPINCFVT3D_TestData data715;
    data715.time = 20.; data715.CFLmultiplier = 5.; data715.mobilityRatio = 2.; data715.velX = 0.001; data715.analyticFront = 25.0;
    data715.testName = "Square40x10m"; data715.testVersion = "LIN-E-FRAC-CFL5.0_"; data715.porosity = 2.5E-03; data715.relpermModel = TPINCFVT3D_TestData::LIN; data715.inflowBoundary = "LEFT";
    //Test( data715, tolSat, satOil, tolFront );
    //tests.push_back(data715);


    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 40.0 x 10.0 m rectangular slab (fracture) model accounting for realistic fracture aperture conditions
    /// Discontinous transport using time multiplier 0.1

    // test 7.1.6.LIN-A
    TPINCFVT3D_TestData data716;
    data716.time = 20.; data716.timeMultiplier = 0.1; data716.CFLmultiplier = 0.1; data716.mobilityRatio = 2.; data716.velX = 0.001; data716.analyticFront = 25.0;
    data716.testName = "Square40x10m"; data716.testVersion = "LIN-A-FRAC-CFL0.1_TMP0.1_"; data716.porosity = 2.5E-03; data716.relpermModel = TPINCFVT3D_TestData::LIN; data716.inflowBoundary = "LEFT";
    //Test( data716, tolSat, satOil, tolFront );
    //tests.push_back(data716);

    // test 7.1.7.LIN-B
    TPINCFVT3D_TestData data717;
    data717.time = 20.; data717.timeMultiplier = 0.1; data717.CFLmultiplier = 0.5; data717.mobilityRatio = 2.; data717.velX = 0.001; data717.analyticFront = 25.0;
    data717.testName = "Square40x10m"; data717.testVersion = "LIN-B-FRAC-CFL0.5_TMP0.1_"; data717.porosity = 2.5E-03; data717.relpermModel = TPINCFVT3D_TestData::LIN; data717.inflowBoundary = "LEFT";
    //Test( data717, tolSat, satOil, tolFront );
    //tests.push_back(data717);

    // test 7.1.8.LIN-C
    TPINCFVT3D_TestData data718;
    data718.time = 20.; data718.timeMultiplier = 0.1; data718.CFLmultiplier = 1.0; data718.mobilityRatio = 2.; data718.velX = 0.001; data718.analyticFront = 25.0;
    data718.testName = "Square40x10m"; data718.testVersion = "LIN-C-FRAC-CFL0.1_TMP0.1_"; data718.porosity = 2.5E-03; data718.relpermModel = TPINCFVT3D_TestData::LIN; data718.inflowBoundary = "LEFT";
    Test( data718, tolSat, satOil, tolFront );
    tests.push_back(data718);


    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 40.0 x 10.0 m rectangular slab (fracture) model accounting for realistic fracture aperture conditions
    /// Discontinous transport using time multiplier 0.01

    // test 7.1.9.LIN-A
    TPINCFVT3D_TestData data719;
    data719.time = 20.; data719.timeMultiplier = 0.01; data719.CFLmultiplier = 0.1; data719.mobilityRatio = 2.; data719.velX = 0.001; data719.analyticFront = 25.0;
    data719.testName = "Square40x10m"; data719.testVersion = "LIN-A-FRAC-CFL0.1_TMP0.01_"; data719.porosity = 2.5E-03; data719.relpermModel = TPINCFVT3D_TestData::LIN; data719.inflowBoundary = "LEFT";
    //Test( data719, tolSat, satOil, tolFront );
    //tests.push_back(data719);

    // test 7.2.1.LIN-B
    TPINCFVT3D_TestData data721;
    data721.time = 20.; data721.timeMultiplier = 0.01; data721.CFLmultiplier = 0.5; data721.mobilityRatio = 2.; data721.velX = 0.001; data721.analyticFront = 25.0;
    data721.testName = "Square40x10m"; data721.testVersion = "LIN-B-FRAC-CFL0.5_TMP0.01_"; data721.porosity = 2.5E-03; data721.relpermModel = TPINCFVT3D_TestData::LIN; data721.inflowBoundary = "LEFT";
    //Test( data721, tolSat, satOil, tolFront );
    //tests.push_back(data721);

    // test 7.2.2.LIN-C // failing time multiplier
    TPINCFVT3D_TestData data722;
    data722.time = 20.; data722.timeMultiplier = 0.01; data722.CFLmultiplier = 1.0; data722.mobilityRatio = 2.; data722.velX = 0.001; data722.analyticFront = 25.0;
    data722.testName = "Square40x10m"; data722.testVersion = "LIN-C-FRAC-CFL1.0_TMP0.01_"; data722.porosity = 2.5E-03; data722.relpermModel = TPINCFVT3D_TestData::LIN; data722.inflowBoundary = "LEFT";
    //Test( data722, tolSat, satOil, tolFront );
    //tests.push_back(data722);


    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 40.0 x 10.0 m rectangular slab model

    // test 7.3.1.LIN-A
    TPINCFVT3D_TestData data731;
    data731.time = 200.; data731.CFLmultiplier = 0.1; data731.mobilityRatio = 2.; data731.velX = 0.01; data731.analyticFront = 25.0;
    data731.testName = "Square40x10m"; data731.testVersion = "LIN-A-CFL0.1_"; data731.porosity = 2.5E-01; data731.relpermModel = TPINCFVT3D_TestData::LIN; data731.inflowBoundary = "LEFT";
    //Test( data731, tolSat, satOil, tolFront );
    //tests.push_back(data731);


    // test 7.3.2LIN-B
    TPINCFVT3D_TestData data732;
    data732.time = 200.; data732.CFLmultiplier = 0.5; data732.mobilityRatio = 2.; data732.velX = 0.01; data732.analyticFront = 25.0;
    data732.testName = "Square40x10m"; data732.testVersion = "LIN-B-CFL0.5_"; data732.porosity = 2.5E-01; data732.relpermModel = TPINCFVT3D_TestData::LIN; data732.inflowBoundary = "LEFT";
    //Test( data732, tolSat, satOil, tolFront );
    //tests.push_back(data732);

    // test 7.3.3.LIN-C
    TPINCFVT3D_TestData data733;
    data733.time = 200.; data733.CFLmultiplier = 1.0; data733.mobilityRatio = 2.; data733.velX = 0.01; data733.analyticFront = 25.0;
    data733.testName = "Square40x10m"; data733.testVersion = "LIN-C-CFL1.0_"; data733.porosity = 2.5E-01; data733.relpermModel = TPINCFVT3D_TestData::LIN; data733.inflowBoundary = "LEFT";
    Test( data733, tolSat, satOil, tolFront );
    tests.push_back(data733);

    // test 7.3.4.LIN-D // failing CFL multiplier
    TPINCFVT3D_TestData data734;
    data734.time = 200.; data734.CFLmultiplier = 2.5; data734.mobilityRatio = 2.; data734.velX = 0.01; data734.analyticFront = 25.0;
    data734.testName = "Square40x10m"; data734.testVersion = "LIN-D-CFL2.5_"; data734.porosity = 2.5E-01; data734.relpermModel = TPINCFVT3D_TestData::LIN; data734.inflowBoundary = "LEFT";
    //Test( data734, tolSat, satOil, tolFront );
    //tests.push_back(data734);

    // test 7.3.5.LIN-E // failing CFL multiplier
    TPINCFVT3D_TestData data735;
    data735.time = 200.; data735.CFLmultiplier = 5.; data735.mobilityRatio = 2.; data735.velX = 0.01; data735.analyticFront = 25.0;
    data735.testName = "Square40x10m"; data735.testVersion = "LIN-E-CFL5.0_"; data735.porosity = 2.5E-01; data735.relpermModel = TPINCFVT3D_TestData::LIN; data735.inflowBoundary = "LEFT";
    //Test( data735, tolSat, satOil, tolFront );
    //tests.push_back(data735);

    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 40.0 x 10.0 m rectangular slab model
    /// Discontinous transport using time multiplier 0.1

    // test 7.3.6.LIN-A
    TPINCFVT3D_TestData data736;
    data736.time = 200.; data736.timeMultiplier = 0.1; data736.CFLmultiplier = 0.1; data736.mobilityRatio = 2.; data736.velX = 0.01; data736.analyticFront = 25.0;
    data736.testName = "Square40x10m"; data736.testVersion = "LIN-A-CFL0.1_TMP0.1_"; data736.porosity = 2.5E-01; data736.relpermModel = TPINCFVT3D_TestData::LIN; data736.inflowBoundary = "LEFT";
    //Test( data736, tolSat, satOil, tolFront );
    //tests.push_back(data736);

    // test 7.3.7.LIN-B
    TPINCFVT3D_TestData data737;
    data737.time = 200.; data737.timeMultiplier = 0.1; data737.CFLmultiplier = 0.5; data737.mobilityRatio = 2.; data737.velX = 0.01; data737.analyticFront = 25.0;
    data737.testName = "Square40x10m"; data737.testVersion = "LIN-B-CFL0.5_TMP0.1_"; data737.porosity = 2.5E-01; data737.relpermModel = TPINCFVT3D_TestData::LIN; data737.inflowBoundary = "LEFT";
    //Test( data737, tolSat, satOil, tolFront );
    //tests.push_back(data737);

    // test 7.3.8.LIN-C
    TPINCFVT3D_TestData data738;
    data738.time = 200.; data738.timeMultiplier = 0.1; data738.CFLmultiplier = 1.0; data738.mobilityRatio = 2.; data738.velX = 0.01; data738.analyticFront = 25.0;
    data738.testName = "Square40x10m"; data738.testVersion = "LIN-C-CFL1.0_TMP0.1_"; data738.porosity = 2.5E-01; data738.relpermModel = TPINCFVT3D_TestData::LIN; data738.inflowBoundary = "LEFT";
    Test( data738, tolSat, satOil, tolFront );
    tests.push_back(data738);

    /// TestCase for Linear Two-Phase model (Warning: viscosities are hard-coded in the two-phase contructor)
    /// running a 40.0 x 10.0 m rectangular slab model
    /// Discontinous transport using time multiplier 0.01

    // test 7.3.9.LIN-A
    TPINCFVT3D_TestData data739;
    data739.time = 200.; data739.timeMultiplier = 0.01; data739.CFLmultiplier = 0.1; data739.mobilityRatio = 2.; data739.velX = 0.01; data739.analyticFront = 25.0;
    data739.testName = "Square40x10m"; data739.testVersion = "LIN-A-CFL0.1_TMP0.01_"; data739.porosity = 2.5E-01; data739.relpermModel = TPINCFVT3D_TestData::LIN; data739.inflowBoundary = "LEFT";
    //Test( data739, tolSat, satOil, tolFront );
    //tests.push_back(data739);

    // test 7.4.1.LIN-B
    TPINCFVT3D_TestData data741;
    data741.time = 200.; data741.timeMultiplier = 0.01; data741.CFLmultiplier = 0.5; data741.mobilityRatio = 2.; data741.velX = 0.01; data741.analyticFront = 25.0;
    data741.testName = "Square40x10m"; data741.testVersion = "LIN-B-CFL0.5_TMP0.01_"; data741.porosity = 2.5E-01; data741.relpermModel = TPINCFVT3D_TestData::LIN; data741.inflowBoundary = "LEFT";
    //Test( data741, tolSat, satOil, tolFront );
    //tests.push_back(data741);

    // test 7.4.2.LIN-C // failing time multiplier
    TPINCFVT3D_TestData data742;
    data742.time = 200.; data742.timeMultiplier = 0.01; data742.CFLmultiplier = 1.0; data742.mobilityRatio = 2.; data742.velX = 0.01; data742.analyticFront = 25.0;
    data742.testName = "Square40x10m"; data742.testVersion = "LIN-C-CFL1.0_TMP0.01_"; data742.porosity = 2.5E-01; data742.relpermModel = TPINCFVT3D_TestData::LIN; data742.inflowBoundary = "LEFT";
    //Test( data742, tolSat, satOil, tolFront );
    //tests.push_back(data742);


    /// Report
    ofstream outfile( "TPINCFVT3DtestResults.txt");
    cout << "\n\nTPINCFVT 3D TestCase test results:"
         <<   "\n=================================\n";
    for( size_t i(0); i < tests.size(); ++i )
      cout << "Test:\t"<< tests.at(i).testName << " " << tests.at(i).testVersion << "\nM:\t" << tests.at(i).mobilityRatio << "\nTime:\t" << tests.at(i).time << "\nTime Multiplier:\t" << tests.at(i).timeMultiplier << "\nCFLmultiplier:\t" << tests.at(i).CFLmultiplier << "\nPorosity:\t" << tests.at(i).porosity
           << "\nAnalytic:\t"<< tests.at(i).analyticFront << "\nNumeric:\t" << tests.at(i).numericFront
           << "\nVelX:\t" << tests.at(i).velX << "\nDiff [m]:\t" << tests.at(i).difference << "\nDiff [%]:\t" << tests.at(i).difference / (tests.at(i).analyticFront/100.)
           << "\nTime elapsed:\t" << tests.at(i).timeElapsed << "\nPassed:\t" << tests.at(i).passed << endl << endl;
    for( size_t i(0); i < tests.size(); ++i )
      outfile << "Test:\t"<< tests.at(i).testName << " " << tests.at(i).testVersion << "\nM:\t" << tests.at(i).mobilityRatio << "\nTime:\t" << tests.at(i).time << "\nPorosity:\t"<< tests.at(i).porosity << "\nCFLmultiplier:\t" << tests.at(i).CFLmultiplier
           << "\nAnalytic:\t"<< tests.at(i).analyticFront << "\nNumeric:\t" << tests.at(i).numericFront
           << "\nVelX:\t" << tests.at(i).velX << "\nDiff [m]:\t" << tests.at(i).difference << "\nDiff [%]:\t" << tests.at(i).difference / (tests.at(i).analyticFront/100.)
           << "\nTime elapsed:\t" << tests.at(i).timeElapsed << "\nPassed:\t" << tests.at(i).passed << endl << endl;
    outfile.close();


    /// Spreadsheet
    ofstream spreadsheet( "TPINCFVT3DtestResultsSpreadSheet.txt");
    cout << "\n\nTPINCFVT 3D test results:"
         <<   "\n========================\n";

    spreadsheet << "Test";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).testName << " " << tests.at(i).testVersion;
    spreadsheet << "\nM";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).mobilityRatio;
    spreadsheet << "\nTime";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).time;
    spreadsheet << "\nTime Multiplier";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).timeMultiplier;
    spreadsheet << "\nCFLmultiplier";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).CFLmultiplier;
    spreadsheet << "\nPorosity";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).porosity;
    spreadsheet << "\nAnalytic";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t"<< tests.at(i).analyticFront;
    spreadsheet << "\nNumeric";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).numericFront;
    spreadsheet << "\nVelX";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).velX;
    spreadsheet << "\nDiff [m]";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).difference;
    spreadsheet << "\nDiff [%]";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).difference / (tests.at(i).analyticFront/100.);
    spreadsheet << "\nTime elapsed";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).timeElapsed;
    spreadsheet << "\nPassed";
    for( size_t i(0); i < tests.size(); ++i )
      spreadsheet << "\t" << tests.at(i).passed;
    spreadsheet << endl << endl;
    spreadsheet.close();
  }//run

void TPINCFVT3D_VVCase::Test( TPINCFVT3D_TestData& data, double tolSat, double satOil, double tolFront )
{
  ANSYS_Model3D model( data.testName.data(), "CSMP-transport-variables.txt" );
  printModelDimensions( model );
  SetModel( model, data.mobilityRatio, data.velX, data.porosity );

  model.Boundary( data.inflowBoundary.data() ).InputPropertyValue( "saturation oil", makeScalar( DIRICH, 0.2 ) );
  model.Boundary( data.inflowBoundary.data() ).InputPropertyValue( "saturation water", makeScalar( DIRICH, 0.8 ) );
  TwoPhaseImplicitNodeCenteredFVTransport<3U,StencilProcessor> TPINCVT( "Model", model );

  TwoPhaseModel<3>* relPerm(NULL);

  switch( data.relpermModel )
  {
    case TPINCFVT3D_TestData::BC: relPerm = new BrooksCorey<3>( model.Database(), "permeability", "viscosity oil", "viscosity water",
                                                               "density oil", "density water",
                                                               "brooks corey parameter", "entry pressure",
                                                               "saturation water","residual saturation oil", "residual saturation water");
    break;

    case TPINCFVT3D_TestData::FL: relPerm = new FourarLenormand<3> ( model.Database(), "fracture aperture", "permeability", "viscosity oil", "viscosity water",
                                                                    "density oil", "density water",
                                                                    "saturation water","residual saturation oil", "residual saturation water" );
    break;

    case TPINCFVT3D_TestData::LIN: relPerm = new LinearTwoPhaseModel<3> ( model.Database(), "permeability", 2.0E-03, 1.0E-03,
                                                                        1.0E+3, 1.0E+3, "entry pressure",
                                                                        "saturation water","residual saturation oil", "residual saturation water");
    break;


    default: throw csmp::Exception( ERROR, "TPINCFVT3D_VVCase::Test", "Internal Error finding active relperm model");
  }
  assert( relPerm );

  TPINCVT.CFL_Multiplier( data.CFLmultiplier );
  Timer runTime;
  runTime.Start();

  /// Checking numerical dispersivity

  double currentTime(0.);
  double simulationTime( data.time );

  if( data.timeMultiplier < 1.)
    while( simulationTime >= currentTime )
    {
     TPINCVT.TransportPhase( *relPerm, data.time * data.timeMultiplier );
     currentTime += data.time * data.timeMultiplier;
     OutputVTU( model, (data.testName + data.testVersion ).c_str(), (long)(currentTime) );
    }
  else
  TPINCVT.TransportPhase( *relPerm, data.time );


  data.timeElapsed = runTime.Stop();
  data.numericFront = TestFront( model, tolSat, satOil, tolFront, data.analyticFront );
  data.difference = data.analyticFront-data.numericFront;
  data.passed = ( abs(data.difference) <= tolFront*data.analyticFront );
  OutputVTU( model, (data.testName + data.testVersion ).c_str(), (long)(data.time) );
  delete relPerm;
}//Test


void TPINCFVT3D_VVCase::OutputVTU( Model<3>& model,
               const char* fileName,
               long time )
{
  VTU_Interface<3> vtu( model );
  model.ExtrapolateElementToNodeProperty( "velocity", "nodal velocity" );
  list<string> props;
  props.push_back( "saturation water" );
  props.push_back( "nodal velocity" );
  vtu.OutputDataToVTU( fileName, "saturation water","Model", static_cast<uint32_t>(time) );
}//OutputVTU


void TPINCFVT3D_VVCase::SetModel( Model<3>& model,
               double mobilityRatio,
               double velX,
               double poro )
{
    const ScalarVariable viscosityWater( PLAIN, 1.0E-03 );
    const ScalarVariable viscosityOil( PLAIN, viscosityWater()*mobilityRatio );
    const ScalarVariable permeability( PLAIN, 1.0E-12 );
    const ScalarVariable densityWater( PLAIN, 1.0E+03 );
    const ScalarVariable densityOil( PLAIN, 1.0E+03 );
    const ScalarVariable saturationWater( PLAIN, 0. );
    const ScalarVariable saturationOil( PLAIN, 1.0 );
    const ScalarVariable porosity( PLAIN, poro );
    const ScalarVariable conductivity( PLAIN, permeability() / viscosityWater() );
    const ScalarVariable lambda( PLAIN, 3.0 );
    const ScalarVariable residualWater( PLAIN, 0.0 );
    const ScalarVariable residualOil( PLAIN, 2.0E-01 );
    const ScalarVariable entryPressure( PLAIN, 1.0E+03 );
    const ScalarVariable fractureAperture( PLAIN, porosity() );
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
    model.InputPropertyValue( "fracture aperture",   fractureAperture );
    model.InputPropertyValue( "nodal fluid volume source", makeScalar(PLAIN, 0.) );

} //SetModel


double TPINCFVT3D_VVCase::TestFront( Model<3>& model,
                                      double tolSat,
                                      double satOil,
                                      double tolFront,
                                      double analyticFront
                                     )
  {
    Point<3> min, max;
    minMaxXYZ( min, max, model );

    Index satKey( model.Database().StorageKey("saturation oil" ) );
    Point<3> firstUntouchedNode( min );

    for( vector<Node<3>*>::const_iterator node(model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
      if( (*node)->Read(satKey) < (satOil-tolSat) && (*node)->x() > firstUntouchedNode[0] )
        firstUntouchedNode = (*node)->Coordinate();

    double numericFront( firstUntouchedNode[0]-min[0] );
    double difference( (analyticFront-numericFront)/(analyticFront/100) );
    cout << "\n\nFront advancement: " << numericFront
         << "\nDifference to analytic [%]: "<< difference << endl << endl;

    _equal( numericFront, analyticFront, analyticFront*tolFront );

    return numericFront;
  }


} //csmp
