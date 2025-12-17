//
//  ExactVersusNumericIntegrationSpeed_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 7/12/2024.
//

#include "ExactVersusNumericIntegrationSpeed_Test.h"
//#include "ANSYS_Model2D.h"
//#include "ANSYS_Model3D.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "InputDataManager.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"

#include "vsetMakers.h"
#include "compareFloats.h"

#include "VSet.h"
#include "ModelTopology.h"

// exact integration
#include "Integral_dNT_op_dN_dV.h"
#include "Integral_NT_op_N_dV.h"

// numeric integration
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"

using namespace std;

namespace csmp {

void ExactVersusNumericIntegrationSpeed_Test::run()
 {
    // best result Dec. 8, 24: 2729 (exact) vs 2878 (numeric) microseconds (serial run M1MAX).
    CompareSpeed();
 }


void ExactVersusNumericIntegrationSpeed_Test::CompareSpeed()
 {
    // result comparison
    double max_value_exact_integration, max_value_numeric_integration;
    
  // 1. Building model of elements with local coordinate system
  // ----------------------------------------------------------
  {
    enum{DIM=3};
    VSet<DIM> vset;
    ModelTopology topo = create_FracBox( vset );

    // 1. local coordinate system computation
    cout <<"Building Model with global coordinate elements..."<<endl;
    const bool use_regions_file{ true };
    Model<DIM> model( topo, vset, "CSMP-1phase-variables.txt", use_regions_file );
    cout <<"Finished reading mesh..."<<endl;
    printModelDimensions( model, true );
    double domain_volume = model.Region("Model").Volume();
    cout <<"\nThe model has a volume of: "<< domain_volume <<" m^3."<< endl;

    // left-right pressure gradient, some fluid volue source
    //InputDataManager<DIM>  model_configuration;
    //model_configuration.ConfigureFromFile( model, this->getName().c_str(),false, true, true, true, false );
    model.InputPropertyValue( "conductivity", makeScalar(ANY,1.0e-12) );
    model.Region("FRACTURE").InputPropertyValue( "conductivity", makeScalar(ANY,1.0e-10) );
    model.InputPropertyValue( "fluid volume source", makeScalar(ANY,0.) );
    model.Region("FRACTURE").InputPropertyValue( "fluid volume source", makeScalar(ANY,1.0e-5) );
    model.InputPropertyValue( "fluid pressure", makeScalar(ANY,0.) );
    // not box-shaped
    model.Boundary("BOUNDARY1").InputPropertyValue( "fluid pressure", makeScalar(DIRICH,1.0e6) ); // LEFT
    model.Boundary("BOUNDARY2").InputPropertyValue( "fluid pressure", makeScalar(DIRICH,1.0e5) ); // RIGHT
    
    printRangeOfVariable(model,"fluid pressure");
    printRangeOfVariable(model,"conductivity");
    printRangeOfVariable(model,"fluid volume source");

    cout <<"\nRunning Test Case Simulation - "<< model.Name()<<endl;
    // setting up & solving linear pressure diffusion
    CSMP_DEFAULT_LINEAR_SOLVER    solver;
    PDE_Integrator<DIM,Element>   pressure_diffusion(solver);
    NumIntegral_dNT_op_dN_dV<DIM> conductance( model.Database(), "conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<DIM>   fluid_src( model.Database(),  "fluid volume source", "fluid pressure" );
    //Assemble the integrator
    pressure_diffusion.Add( &conductance );
    pressure_diffusion.Add( &fluid_src );
    // computation
    auto start = std::chrono::high_resolution_clock::now();
    model.Apply( pressure_diffusion );
    auto stop  = std::chrono::high_resolution_clock::now();
    max_value_numeric_integration = printRangeOfVariable(model,"fluid pressure");
    std::cout << "\n"<<"ExactVersusNumericIntegrationSpeed_Test::CompareSpeed::Computation using LocalCoordinateFiniteElement: "
              << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << " microseconds."<< endl;
    // check
    VTU_Interface<DIM> vtu( model );
    vtu.OmitZeroInFileName(true);
    list<string> outputProps;
    outputProps.push_back( "fluid pressure" );
    vtu.OutputDataToVTU( "FracBox_isoparametric_result", outputProps, "Model", static_cast<int>(0) );
  }


  // 2. Building model of elements with global coordinate system
  // -----------------------------------------------------------
  {
    enum{DIM=3};
    VSet<DIM> vset;
    const bool isoparametric_elmts{ false };
    ModelTopology topo = create_FracBox( vset /* isoparametric_elmts */ );

    // 1. local coordinate system computation
    cout <<"Building Model with global coordinate elements..."<<endl;
    const bool use_regions_file{ true };
    Model<DIM> model( topo, vset, "CSMP-1phase-variables.txt", use_regions_file );
    cout <<"Finished reading mesh..."<<endl;
    printModelDimensions( model, true );
    double domain_volume = model.Region("Model").Volume();
    cout <<"\nThe model has a volume of: "<< domain_volume <<" m^3."<< endl;

    // left-right pressure gradient, some fluid volue source
    //InputDataManager<DIM>  model_configuration;
    //model_configuration.ConfigureFromFile( model, this->getName().c_str(),false, true, true, true, false );
    model.InputPropertyValue( "conductivity", makeScalar(ANY,1.0e-12) );
    model.Region("FRACTURE").InputPropertyValue( "conductivity", makeScalar(ANY,1.0e-10) );
    model.InputPropertyValue( "fluid volume source", makeScalar(ANY,0.) );
    model.Region("FRACTURE").InputPropertyValue( "fluid volume source", makeScalar(ANY,1.0e-5) );
    model.InputPropertyValue( "fluid pressure", makeScalar(ANY,0.) );
    // not box-shaped
    // model.InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH,1.0e6) );
    // model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.0e5) );
    model.Boundary("BOUNDARY1").InputPropertyValue( "fluid pressure", makeScalar(DIRICH,1.0e6) );
    model.Boundary("BOUNDARY2").InputPropertyValue( "fluid pressure", makeScalar(DIRICH,1.0e5) );

    printRangeOfVariable(model,"conductivity");
    printRangeOfVariable(model,"fluid volume source");
    printRangeOfVariable(model,"fluid pressure");

    cout <<"\nRunning Test Case Simulation - "<< model.Name()<<endl;
    // setting up & solving linear pressure diffusion
    CSMP_DEFAULT_LINEAR_SOLVER    solver;
    PDE_Integrator<DIM,Element>   pressure_diffusion(solver);
    NumIntegral_dNT_op_dN_dV<DIM> conductance( model.Database(), "conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<DIM>   fluid_src( model.Database(),  "fluid volume source", "fluid pressure" );
    //Assemble the integrator
    pressure_diffusion.Add( &conductance );
    pressure_diffusion.Add( &fluid_src );
    // computation
    auto start = std::chrono::high_resolution_clock::now();
    model.Apply( pressure_diffusion );
    auto stop  = std::chrono::high_resolution_clock::now();
    max_value_exact_integration = printRangeOfVariable(model,"fluid pressure");
    std::cout << "\n"<<"ExactVersusNumericIntegrationSpeed_Test::CompareSpeed::Computation using GlobalCoordinateFiniteElement: "
              << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << " microseconds."<< endl;
    // check
    const double tolerance_relax_factor{ 1. }; // no scaling
    _test( approximatelyEqual(max_value_numeric_integration,max_value_exact_integration,tolerance_relax_factor) );
    VTU_Interface<DIM> vtu( model );
    vtu.OmitZeroInFileName(true);
    list<string> outputProps;
    outputProps.push_back( "fluid pressure" );
    vtu.OutputDataToVTU( "FracBox_exact_result", outputProps, "Model", static_cast<int>(0) );
  }

} // end CompareSpeed

} // end
