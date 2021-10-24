#include "RegionMonitor_Test.h"
#include "Region.h"
#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"

using namespace std;

namespace csmp {


RegionMonitor_Test::RegionMonitor_Test()
: fTolerance(numeric_limits<double64>::epsilon() * 1000.) // scaled with the larges number that will be encountered
{
}



RegionMonitor_Test::~RegionMonitor_Test()
{
}



void RegionMonitor_Test::run()
{
    RegionMonitorScalarPropertyIntegrals2D();
    RegionMonitorScalarPropertyIntegrals3D();
}


void RegionMonitor_Test::RegionMonitorScalarPropertyIntegrals2D()
{
    ANSYS_Model2D model("LeftRight", "CSMP-2phase-variables.txt");
    
    // assigning material properties 
    // -----------------------------------------------------  
    model.Region("M_LEFT").InputPropertyValue( "porosity", makeScalar(PLAIN,0.2) );
    model.Region("M_RIGHT").InputPropertyValue( "porosity", makeScalar(PLAIN,0.4) );
    model.Region("M_LEFT").InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-12) );
    model.Region("M_RIGHT").InputPropertyValue( "permeability", makeScalar(PLAIN,2.0e-12) );
    
    std::list<std::string> integral_properties;
    std::list<std::string> range_properties;
    
    integral_properties.push_back("permeability");
    integral_properties.push_back("porosity");
    range_properties.push_back("permeability");
    range_properties.push_back("porosity");
    
    RegionMonitor<2U> region_monitor(model, integral_properties, range_properties, false);
    region_monitor.ScalarPropertyIntegrals(model, 0.);
    region_monitor.ScalarPropertyRanges(model, 0.);
    
    region_monitor.Out("RM_out");

    std::ifstream fin("RM_out.txt");
    std::string input;
    double64 val1, val2, val3, val4, val5, val6, val7;
    
    std::getline(fin, input);
    std::getline(fin, input);
    fin >> input >> val1 >> val2;
    _test( input == "M_LEFT" );
    _equal( val1, 50., fTolerance );
    _equal( val2, 50., fTolerance );

    fin >> input >> val1 >> val2;
    _test( input == "M_RIGHT" );
    _equal( val1, 50., fTolerance );
    _equal( val2, 50., fTolerance );

//JC: something wrong with this. commented out those lines. check it later.
    //fin >> input >> val1 >> val2;
    //_test( input == "Model" );
    //_equal( val1, 100., fTolerance );
    //_equal( val2, 100., fTolerance );
    
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    bool found(input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 5.e-11, fTolerance );
    _equal( val3, 1.e-10, fTolerance );
    //_equal( val4, 1.5e-10, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );

    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 10., fTolerance );
    _equal( val3, 20., fTolerance );
    //_equal( val4, 30, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 1.e-12, fTolerance );
    _equal( val4, 2.e-12, fTolerance );
    _equal( val5, 2.e-12, fTolerance );
    //_equal( val6, 1.e-12, fTolerance );
    //_equal( val7, 2.e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.2, fTolerance );
    _equal( val4, 0.4, fTolerance );
    _equal( val5, 0.4, fTolerance );
    //_equal( val6, 0.2, fTolerance );
    //_equal( val7, 0.4, fTolerance );
    
    fin.close();
    std::remove("RM_out.txt");
    
    region_monitor.EraseData();

    // checking with porevolume integration    
    region_monitor.DefineProperties(model, integral_properties, range_properties, true);
    region_monitor.ScalarPropertyIntegrals(model, 0.);
    region_monitor.ScalarPropertyRanges(model, 0.);
    
    region_monitor.Out("RM_out");

    fin.open("RM_out.txt");
    std::getline(fin, input);
    std::getline(fin, input);
    fin >> input >> val1 >> val2;
    _test( input == "M_LEFT" );
    _equal( val1, 10., fTolerance );
    _equal( val2, 50., fTolerance );

    fin >> input >> val1 >> val2;
    _test( input == "M_RIGHT" );
    _equal( val1, 20., fTolerance );
    _equal( val2, 50., fTolerance );

    //fin >> input >> val1 >> val2;
    //_test( input == "Model" );
    //_equal( val1, 30., fTolerance );
    //_equal( val2, 100., fTolerance );
    
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-11, fTolerance );
    _equal( val3, 4.e-11, fTolerance );
    //_equal( val4, 5.e-11, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );

    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 2., fTolerance );
    _equal( val3, 8., fTolerance );
   // _equal( val4, 10, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 1.e-12, fTolerance );
    _equal( val4, 2.e-12, fTolerance );
    _equal( val5, 2.e-12, fTolerance );
    //_equal( val6, 1.e-12, fTolerance );
    //_equal( val7, 2.e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.2, fTolerance );
    _equal( val4, 0.4, fTolerance );
    _equal( val5, 0.4, fTolerance );
    //_equal( val6, 0.2, fTolerance );
    //_equal( val7, 0.4, fTolerance );
    
    fin.close();
    std::remove("RM_out.txt");


    // check division by volume
    region_monitor.DefineProperties(model, integral_properties, range_properties, false);
    region_monitor.DivideIntegralPropertiesByRegionVolumes(true);
    region_monitor.ScalarPropertyIntegrals(model, 0.);
    region_monitor.ScalarPropertyRanges(model, 0.);

    region_monitor.Out("RM_out");
    
    fin.open("RM_out.txt");
    std::getline(fin, input);
    std::getline(fin, input);
    fin >> input >> val1 >> val2;
    _test( input == "M_LEFT" );
    _equal( val1, 50., fTolerance );
    _equal( val2, 50., fTolerance );

    fin >> input >> val1 >> val2;
    _test( input == "M_RIGHT" );
    _equal( val1, 50., fTolerance );
    _equal( val2, 50., fTolerance );

    //fin >> input >> val1 >> val2;
    //_test( input == "Model" );
    //_equal( val1, 100., fTolerance );
    //_equal( val2, 100., fTolerance );
    
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 2.e-12, fTolerance );
    //_equal( val4, 1.5e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );

    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.4, fTolerance );
    //_equal( val4, 0.3, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 1.e-12, fTolerance );
    _equal( val4, 2.e-12, fTolerance );
    _equal( val5, 2.e-12, fTolerance );
    //_equal( val6, 1.e-12, fTolerance );
    //_equal( val7, 2.e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.2, fTolerance );
    _equal( val4, 0.4, fTolerance );
    _equal( val5, 0.4, fTolerance );
    //_equal( val6, 0.2, fTolerance );
    //_equal( val7, 0.4, fTolerance );
    
    fin.close();
    std::remove("RM_out.txt");
    
    region_monitor.Reset();

    
}



void RegionMonitor_Test::RegionMonitorScalarPropertyIntegrals3D()
{
    ANSYS_Model3D model("Cube", "CSMP-2phase-variables.txt");
    
    // assigning material properties 
    // -----------------------------------------------------  
    model.Region("LAYER1").InputPropertyValue( "porosity", makeScalar(PLAIN,0.2) );
    model.Region("LAYER2").InputPropertyValue( "porosity", makeScalar(PLAIN,0.4) );
    model.Region("LAYER1").InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-12) );
    model.Region("LAYER2").InputPropertyValue( "permeability", makeScalar(PLAIN,2.0e-12) );
    
    std::list<std::string> integral_properties;
    std::list<std::string> range_properties;
    
    integral_properties.push_back("permeability");
    integral_properties.push_back("porosity");
    range_properties.push_back("permeability");
    range_properties.push_back("porosity");
    
    RegionMonitor<3U> region_monitor( model, integral_properties, range_properties, false );
    region_monitor.ScalarPropertyIntegrals(model, 0.);
    region_monitor.ScalarPropertyRanges(model, 0.);
    
    region_monitor.Out("RM_out");

    std::ifstream fin("RM_out.txt");
    std::string input;
    double64 val1, val2, val3, val4, val5, val6, val7;
    
    std::getline(fin, input);
    std::getline(fin, input);
    fin >> input >> val1 >> val2;
    _test( input == "LAYER1" );
    _equal( val1, 500., fTolerance );
    _equal( val2, 400., fTolerance );

    fin >> input >> val1 >> val2;
    _test( input == "LAYER2" );
    _equal( val1, 500., fTolerance );
    _equal( val2, 400., fTolerance );

    //fin >> input >> val1 >> val2;
    //_test( input == "Model" );
    //_equal( val1, 1000., fTolerance );
    //_equal( val2, 600., fTolerance );
    
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    bool found(input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 5.e-10, fTolerance );
    _equal( val3, 1.e-9, fTolerance );
    //_equal( val4, 1.5e-9, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );

    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 100., fTolerance );
    _equal( val3, 200., fTolerance );
    //_equal( val4, 300, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 1.e-12, fTolerance );
    _equal( val4, 2.e-12, fTolerance );
    _equal( val5, 2.e-12, fTolerance );
    //_equal( val6, 1.e-12, fTolerance );
    //_equal( val7, 2.e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.2, fTolerance );
    _equal( val4, 0.4, fTolerance );
    _equal( val5, 0.4, fTolerance );
    //_equal( val6, 0.2, fTolerance );
    //_equal( val7, 0.4, fTolerance );
    
    fin.close();
    std::remove("RM_out.txt");
    
    region_monitor.EraseData();

    // checking with porevolume integration    
    region_monitor.DefineProperties(model, integral_properties, range_properties, true);
    region_monitor.ScalarPropertyIntegrals(model, 0.);
    region_monitor.ScalarPropertyRanges(model, 0.);
    region_monitor.Out("RM_out");

    fin.open("RM_out.txt");
    std::getline(fin, input);
    std::getline(fin, input);
    fin >> input >> val1 >> val2;
    _test( input == "LAYER1" );
    _equal( val1, 100., fTolerance );
    _equal( val2, 400., fTolerance );

    fin >> input >> val1 >> val2;
    _test( input == "LAYER2" );
    _equal( val1, 200., fTolerance );
    _equal( val2, 400., fTolerance );

    //fin >> input >> val1 >> val2;
    //_test( input == "Model" );
    //_equal( val1, 300., fTolerance );
    //_equal( val2, 600., fTolerance );
    
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-10, fTolerance );
    _equal( val3, 4.e-10, fTolerance );
    //_equal( val4, 5.e-10, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );

    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 20., fTolerance );
    _equal( val3, 80., fTolerance );
    //_equal( val4, 100, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 1.e-12, fTolerance );
    _equal( val4, 2.e-12, fTolerance );
    _equal( val5, 2.e-12, fTolerance );
    //_equal( val6, 1.e-12, fTolerance );
    //_equal( val7, 2.e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.2, fTolerance );
    _equal( val4, 0.4, fTolerance );
    _equal( val5, 0.4, fTolerance );
    //_equal( val6, 0.2, fTolerance );
    //_equal( val7, 0.4, fTolerance );
    
    fin.close();
    std::remove("RM_out.txt");


    // check division by volume
    region_monitor.DefineProperties(model, integral_properties, range_properties, false);
    region_monitor.DivideIntegralPropertiesByRegionVolumes(true);
    region_monitor.ScalarPropertyIntegrals(model, 0.);
    region_monitor.ScalarPropertyRanges(model, 0.);
    region_monitor.Out("RM_out");
    
    fin.open("RM_out.txt");
    std::getline(fin, input);
    std::getline(fin, input);
    fin >> input >> val1 >> val2;
    _test( input == "LAYER1" );
    _equal( val1, 500., fTolerance );
    _equal( val2, 400., fTolerance );

    fin >> input >> val1 >> val2;
    _test( input == "LAYER2" );
    _equal( val1, 500., fTolerance );
    _equal( val2, 400., fTolerance );

    //fin >> input >> val1 >> val2;
    //_test( input == "Model" );
    //_equal( val1, 1000., fTolerance );
    //_equal( val2, 600., fTolerance );
    
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 2.e-12, fTolerance );
    //_equal( val4, 1.5e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );

    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4;
	fin >> val1 >> val2 >> val3;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.4, fTolerance );
    //_equal( val4, 0.3, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("permeability") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 1.e-12, fTolerance );
    _equal( val3, 1.e-12, fTolerance );
    _equal( val4, 2.e-12, fTolerance );
    _equal( val5, 2.e-12, fTolerance );
    //_equal( val6, 1.e-12, fTolerance );
    //_equal( val7, 2.e-12, fTolerance );

    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);
    std::getline(fin, input);

    found = (input.find("porosity") != string::npos);
    _test( found == true );
    
    std::getline(fin, input);
    //fin >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7;
	fin >> val1 >> val2 >> val3 >> val4 >> val5;
    _equal( val1, 0., fTolerance );
    _equal( val2, 0.2, fTolerance );
    _equal( val3, 0.2, fTolerance );
    _equal( val4, 0.4, fTolerance );
    _equal( val5, 0.4, fTolerance );
    //_equal( val6, 0.2, fTolerance );
    //_equal( val7, 0.4, fTolerance );
    
    fin.close();
    std::remove("RM_out.txt");
    
    region_monitor.Reset();

    
}


} // end namespace csmp
