#include "PDE_Integrator_Computation_Test.h"

#include "vsetMakers.h"
#include "VSetConverter.h"
#include "geometricCalculations.h"
#include "CSMP_physical_constants.h"

// the CSMP model
#include "Model.h"
#include "ANSYS_Model3D.h"
#include "ModelTime.h"
#include "ModelTopology.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "FiniteElement.h"
#include "meshManagementUtilities.h"

// the FE algorithm
#include "PDE_Integrator.h"

// PDE operators building the FE algorithm
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "VelocityAndVolumeFlux.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

// output interfaces
#include "VTU_Interface.h"

// utility functions
//#include "CSMP_highLevelUtilities.h"
#include "ConstantFactor.h"

using namespace std;

namespace csmp {

// must delete the model
PDE_Integrator_Computation_Test::~PDE_Integrator_Computation_Test() { delete model_ptr_; }



void PDE_Integrator_Computation_Test::run()
{
    // -----------------------------------------------------------------------
    // 1. Mesh and model setup
    // -----------------------------------------------------------------------
    double& model_time(ModelTime::Instance().modelTime);
    model_time = 0.0;
    
    /// model with poly- element types including prisms and pyramids
    model_ptr_ = BuildAnsysModel3D("prism_test");
    //model_ptr_ = BuildFracBoxModel3D();
    //model_ptr_ = BuildAnsysModel3D("Shuaiba");

    VTU_Interface<3> vtu_output( *model_ptr_ );
    list<string> input_props = {"fluid pressure","conductivity"};
    vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_input" ).c_str(), input_props, string("Model"), 0 );

    /// angle criteria tests for simplex and poly-element meshes
    auto skewed_elmts = VerifyMeshQuality();
    //                  ^^^^^^^^^^^^^^^^^^^^
    if ( !skewed_elmts.empty() && verbose_ ) {
       model_ptr_->FormRegionFrom( "skewed-elements", skewed_elmts.begin(), skewed_elmts.end() );
       vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_input" ).c_str(),
                                   input_props, string("skewed-elements"), 0 );
    }
    /// testing for which skewed elements a singular value decomposition of dNT_dN reveals true degeneracy
    if ( !skewed_elmts.empty() ) {
        vector<Element<3U>*> degenerate_elmts = testCellSkewing<3,Element>( skewed_elmts, 1.0e10, true );

        // --- form diagnostic region and write to VTK ---
        if ( !degenerate_elmts.empty() ) {
            std::vector<size_t> degenerate_ids;
            for ( const Element<3U>* e : degenerate_elmts )
                degenerate_ids.push_back( e->Idx() );

            model_ptr_->FormRegionFrom( "degenerate-elements", degenerate_elmts.begin(), degenerate_elmts.end() );
            vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_input" ).c_str(), input_props, string("degenerate-elements"), 0 );
        }
        _test( degenerate_elmts.empty() == true );
    }
  
    // -----------------------------------------------------------------------
    // 2. Material properties and boundary conditions
    // -----------------------------------------------------------------------
    AssignProperties_PressureDiffusionAndFlow();
    AssignEssentialConditions_PressureDiffusionAndFlow();
    // make fracture highly permeable
//    model_ptr_->Region("FRACTURE").InputPropertyValue("permeability", makeScalar(ANY,1.0e-9) );
    model_ptr_->Region("FRAC_VOLUMES").InputPropertyValue("permeability", makeScalar(ANY,1.0e-9) );
    model_ptr_->Region("FRAC_VOLUMES").InputPropertyValue("porosity", makeScalar(ANY,1.0) );

    // -----------------------------------------------------------------------
    // 3. Compute hydraulic conductivity, K
    // -----------------------------------------------------------------------
    constexpr double dynamic_viscosity{ 0.001 };
    ConstantFactor<3,divides> conductivity( model_ptr_->Database(),
                                           "conductivity", "permeability", dynamic_viscosity );
    model_ptr_->Apply(conductivity);
    printRangeOfVariable( *model_ptr_, "conductivity" );


    // -----------------------------------------------------------------------
    // 4. Single PDE_Integrator — running test
    // -----------------------------------------------------------------------
    list<string> output_props = {"fluid pressure","velocity","volume flux"};

    // --- Steady state solution test (scalar equation) ---
    //     --------------------------------------------
    //     NumIntegral_dNT_lhsop_dN_dV<3>
    //     NumIntegral_NT_rhsop_N_dV<3>
    //     VelocityAndVolumeFlux<3>
    TestSteadyState_PressureDiffusionAndFlow();
    if ( verbose_ ) {
         printRangeOfVariable( *model_ptr_, "fluid pressure" );
         vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_output" ).c_str(), output_props, string("Model"), 0 );
         vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_output" ).c_str(), output_props, string("FRAC_VOLUMES"), 0 );
      }


    // TODO: test that the solution is divergence free for source=0
    
    // TODO: test performance for quadratic elements and ones with VECTOR and TENSOR properties

    // --- Transient solution test ---
    // --------------------------------
    //  extra operators
    //     NumIntegral_NT_lhsop_N_dV<3>
    //     NumIntegral_NT_rhsop_N_dV<3>
    double time_interval_tested = TestTransient_PressureDiffusionAndFlow();

    if ( verbose_ ) {
         vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_output" ).c_str(), output_props, string("Model"),       static_cast<long>(time_interval_tested) );
         vtu_output.OutputDataToVTU( string( string( model_ptr_->Name() ) + "-test_output" ).c_str(), output_props, string("FRAC_VOLUMES"), static_cast<long>(time_interval_tested) );
      }

    cout << "\nPDE_Integrator_Computation_Test: all tests run.\n";

} // end run



 
 // from vsetMakers
Model<3U>* PDE_Integrator_Computation_Test::BuildFracBoxModel3D()
 {
    // creating a valid CSMP model
    string  model_name("FracBox3D");
    {
      VSet<3> vset;
      
      ModelTopology subdivision_in_regions = create_FracBox( vset );
      
      // refining the model
      VSet<3>       refined_lin_mesh;
      ModelTopology refined_linear_topo;
      VSet<3>       quadratic_mesh;
      ModelTopology quadratic_topo;
//      transformMeshIntoRefinedLinearAndQuadraticMeshes( const VSet<3>& linear_mesh, const ModelTopology& linear_topo,
//                                                      VSet<3>& refined_lin_mesh, ModelTopology& refined_linear_topo,
//                                                       VSet<3>& quadratic_mesh, ModelTopology& quadratic_topo )
//      transformMeshIntoRefinedLinearAndQuadraticMeshes( vset, subdivision_in_regions,
//                                                        refined_lin_mesh,  refined_linear_topo,
//                                                        quadratic_mesh,  quadratic_topo );
      const bool use_regions_file{ true };
      Model<3> model( subdivision_in_regions, vset, "CSMP-1phase-variables.txt", use_regions_file );
//      Model<3> model( refined_linear_topo, refined_lin_mesh, "CSMP-1phase-variables.txt", use_regions_file );
      model.Name( (model_name + "_PDE_Integrator_Computation_Test").c_str() );
      printModelDimensions(model);

      // renaming the boundaries
      model.RenameBoundary( "BOUNDARY1", "LEFT" );
      model.RenameBoundary( "BOUNDARY2", "RIGHT" );
      model.RenameBoundary( "BOUNDARY3", "FRONT" );
      model.RenameBoundary( "BOUNDARY4", "BACK" );
      model.RenameBoundary( "BOUNDARY5", "BOTTOM" );
      model.RenameBoundary( "BOUNDARY6", "TOP" );
      
      model.BoundariesOut();
      
      model.OutputToBinaryFile( model.Name() );
    }
    model_ptr_ = new Model<3U>( model_name + "_PDE_Integrator_Computation_Test" );
    return model_ptr_;

 } // end BuildFracBoxModel3D






Model<3U>* PDE_Integrator_Computation_Test::BuildAnsysModel3D( const std::string& model_name )
 {
    ANSYS_Model3D model( model_name.c_str(), model_name.c_str(), "CSMP-1phase-variables.txt", true );
    printModelDimensions(model);

    model.OutputToBinaryFile( (model_name + "_PDE_Integrator_Computation_Test").c_str() );
    
    model_ptr_ = new Model<3U>( model_name + "_PDE_Integrator_Computation_Test" );
    return model_ptr_;

 } // end BuildModel3D





void PDE_Integrator_Computation_Test::AssignProperties_PressureDiffusionAndFlow()
 {
    // setting up necessary variables
    model_ptr_->InputPropertyValue( "permeability", makeScalar( PLAIN, 1.0e-14 ) );
    model_ptr_->InputPropertyValue( "porosity", makeScalar( PLAIN, 0.25 ) );
    model_ptr_->InputPropertyValue( "fluid volume source", makeScalar( PLAIN, 0.0 ) );
    model_ptr_->InputPropertyValue( "velocity", makeVector( ANY, ANY, ANY, 0., 0., 0. ) );
    model_ptr_->InputPropertyValue( "pore velocity", makeVector( ANY, ANY, ANY, 0., 0., 0. )  );
    model_ptr_->InputPropertyValue( "volume flux", makeScalar( PLAIN, 0. ) );
    model_ptr_->InputPropertyValue( "storativity", makeScalar( PLAIN, 1.0e-8 ) );
 }
 
 
 
 
 
 // left-right pressure gradient
 void PDE_Integrator_Computation_Test::AssignEssentialConditions_PressureDiffusionAndFlow()
  {
    constexpr double patm{ 100325. };
    constexpr double hydrostatic_grad{ 1000. * csmp::ACC_GRAVITY };
    
    const double max_length = printModelDimensions( *model_ptr_ );
    // setting background fluid pressure value as it is used in the elimination process
    model_ptr_->InputPropertyValue( "fluid pressure", makeScalar(ANY,0.0) );
    // assuming that x-length is equal to max length
    model_ptr_->Boundary("LEFT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * max_length) );
    model_ptr_->Boundary("RIGHT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
  }




/**
  Tests for small-angle = degenerate finite elements
*/
vector<Element<3>*> PDE_Integrator_Computation_Test::VerifyMeshQuality()
 {
    // are there any issues?
    vector<Element<3>*> skewed_elements;
    skewed_elements.reserve(1000);
    for ( const auto& cell : model_ptr_->Region("Model").CellVector() ) {
          if ( !isValidElement( cell ) ) skewed_elements.push_back(cell);
          // Test fails correctly for skewed elements:
          // _test( isValidElement( cell ) == true );
      }

    // reporting issues to User
    if ( verbose_ && skewed_elements.size() > 0 )
      for ( const auto& cell : model_ptr_->Region("Model").CellVector() )
        if ( !isValidElement( cell ) ) {
             cout <<"\n"<<"Cell "<< cell->Idx() <<": "<< parseFiniteElementType(cell->FE_Type()) <<" is overly skewed.";
             for ( auto nit=cell->NodesBegin(); nit!=cell->NodesEnd(); ++nit )
               cout <<" "<< (*nit)->Coordinate();
          }
          
    return skewed_elements;
 }



  
void PDE_Integrator_Computation_Test::TestSteadyState_PressureDiffusionAndFlow()
{
  cout <<"\n"<<"PDE_Integrator_Computation_Test::TestSteadyState_PressureDiffusionAndFlow: testing..."<< endl;

  // -----------------------------------------------------------------------
  // 4. PDE operators
  // -----------------------------------------------------------------------
  NumIntegral_dNT_lhsop_dN_dV<3>  stiffness_matrix( model_ptr_->Database(), "conductivity", "fluid pressure", "fluid pressure");

  NumIntegral_NT_rhsop_N_dV<3> source_term( model_ptr_->Database(), "fluid volume source", "fluid pressure");

  VelocityAndVolumeFlux<3>  velo( *model_ptr_, "conductivity", "porosity", "fluid pressure", false );

  // -----------------------------------------------------------------------
  // 5. Single PDE_Integrator — setup and computation
  // -----------------------------------------------------------------------
  PDE_Integrator<3,Element> pde_integrator;
#ifdef CSMP_WITH_SAMG_SOLVER
  SAMG_Settings settings;
  SAMG_Solver   samg_solver( &settings );
  pde_integrator.SetSolver( samg_solver );
#else
  CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
  pde_integrator.SetSolver( linear_solver );
#endif

  pde_integrator.Add(&stiffness_matrix);
  pde_integrator.Add(&source_term);
  pde_integrator.AddPostProcess(&velo);

  // -------------------------------------------------------------------
  // Robust performance measurement
  // -------------------------------------------------------------------
  {
      // 1. Warm-up pass — populates caches, triggers frequency scaling,
      //    runs any lazy initialisations inside Apply.
      //    Not measured.
      model_ptr_->Apply( pde_integrator );

      // turning excessive diagnostic and output off, speeding up SAMG so that accumulation plays a bigger role in test
#ifdef CSMP_WITH_SAMG_SOLVER
      settings.Set_iout1( 0 );
      settings.Set_iout2( 0 );
      settings.Set_idmp( -1 );
#endif
      // 2. Repeated measurements (30 is minimum to achieve some reproducibility)
      constexpr int RUNS = 30;
      std::vector<double> times_ms;
      times_ms.reserve( RUNS );

      for ( int run = 0; run < RUNS; ++run )
      {
          auto t0 = std::chrono::high_resolution_clock::now();
          model_ptr_->Apply( pde_integrator );
          // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          auto t1 = std::chrono::high_resolution_clock::now();

          times_ms.push_back(
              std::chrono::duration<double, std::milli>( t1 - t0 ).count() );
      }

      // 3. Statistics
      const double mean = std::accumulate( times_ms.begin(),
                                           times_ms.end(), 0.0 ) / RUNS;

      const double sq_sum = std::inner_product( times_ms.begin(),
                                                 times_ms.end(),
                                                 times_ms.begin(), 0.0 );
      const double stddev = std::sqrt( sq_sum / RUNS - mean * mean );

      const double best  = *std::min_element( times_ms.begin(), times_ms.end() );
      const double worst = *std::max_element( times_ms.begin(), times_ms.end() );

      std::cout << "\nTestSteadyState_PressureDiffusionAndFlow performance ("
                << RUNS << " runs, Release build):\n"
                << "\t best   : " << best   << " ms\n"
                << "\t worst  : " << worst  << " ms\n"
                << "\t mean   : " << mean   << " ms\n"
                << "\t std dev: " << stddev << " ms\n";
  }
  // -----------------------------------------------------------------------
  // 12. VTK output (non-critical, just verify no crash)
  // -----------------------------------------------------------------------

  cout << "\nPDE_Integrator_Transient_Test: all tests passed.\n";

} // end TestSteadyState_PressureDiffusionAndFlow




/*
    Averages the time spent on each step of a pressure diffusion run with 30 steps and reports the result.
    A final snapshot of the pressure distribution and the other features with flux down to ~zero once
    equilibrium is reached is also output
    
    returns final time reached
*/
double PDE_Integrator_Computation_Test::TestTransient_PressureDiffusionAndFlow()
{
    // -----------------------------------------------------------------------
    // 1. PDE operators for initial steady-state pressure computation
    // -----------------------------------------------------------------------
    NumIntegral_dNT_dN_dV<3>    stiffness_matrix( model_ptr_->Database(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<3>   source_term(      model_ptr_->Database(), "fluid volume source", "fluid pressure" );
    VelocityAndVolumeFlux<3>    velo( *model_ptr_, "conductivity", "porosity", "fluid pressure", false );

    // 1.1 Setting up PDE_Integrator and computation
    // -----------------------------------------------------------------------
    {
        PDE_Integrator<3,Element> pde_integrator;

    #ifdef CSMP_WITH_SAMG_SOLVER
        SAMG_Settings settings;
        SAMG_Solver   samg_solver( &settings );
        pde_integrator.SetSolver( samg_solver );
    #else
        CSMP_DEFAULT_LINEAR_SOLVER linear_solver;
        pde_integrator.SetSolver( linear_solver );
    #endif

        pde_integrator.Add( &stiffness_matrix );
        pde_integrator.Add( &source_term      );
        pde_integrator.AddPostProcess( &velo  );

        model_ptr_->Apply( pde_integrator );
    }
    

    // -----------------------------------------------------------------------
    // 2. PDE operators for transient computation
    // -----------------------------------------------------------------------
    // ([C] + dt[K]){p}t+dt = [C]{p}t + dt{Q}t+dt
    NumIntegral_NT_lhsop_N_dV<3> mass_matrix_lhs( model_ptr_->Database(), "storativity", "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<3> mass_matrix_rhs( model_ptr_->Database(), "storativity", "fluid pressure" );

    mass_matrix_lhs.MultiplyWithTimeIncrement( true );
    mass_matrix_rhs.MultiplyWithTimeIncrement( true );
    mass_matrix_lhs.LumpedFormulation( true );
    mass_matrix_rhs.LumpedFormulation( true );
    source_term.LumpedFormulation( true );
    source_term.AddAccumulateLater();

    PDE_Integrator<3,Element> pde_integrator;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    SAMG_Solver   samg_solver( &settings );
    pde_integrator.SetSolver( samg_solver );
#else
    CSMP_DEFAULT_LINEAR_SOLVER linear_solver;
    pde_integrator.SetSolver( linear_solver );
#endif

    pde_integrator.Add( &stiffness_matrix );
    pde_integrator.Add( &source_term      );
    pde_integrator.Add( &mass_matrix_lhs );
    pde_integrator.Add( &mass_matrix_rhs );
    pde_integrator.AddPostProcess( &velo  );

#ifdef CSMP_WITH_SAMG_SOLVER
    settings.Set_iout1(   0       );
    settings.Set_iout2(   0       );
    settings.Set_idmp(   -1       );
//    settings.Set_iswit(   4       );  // reuse solver setup from last call
    settings.Set_itypu(   0       );  // use solution from last step as initial guess
//    settings.Set_rel_eps( 1.0e-14 );  // relative tolerance
#endif

    // -----------------------------------------------------------------------
    // 8. Timestepping with performance measurement
    // -----------------------------------------------------------------------
    constexpr int    N_STEPS{ 30 };
    constexpr double DT_INITIAL{ 0.001 };
    constexpr double DT_GROWTH{  1.1 };

    std::vector<double> times_ms;
    times_ms.reserve( N_STEPS );

    double dt(       DT_INITIAL );
    double duration( 0.0        );

    // drop pressure on the left boundary to atmospheric
    model_ptr_->Boundary("LEFT").InputPropertyValue(
        "fluid pressure", makeScalar( DIRICH, 100325. ) );


    for ( int tstep{1}; tstep <= N_STEPS; ++tstep )
      {
        pde_integrator.TimeIncrement( 1.0 / dt );

        const auto t0 = std::chrono::high_resolution_clock::now();
        model_ptr_->Apply( pde_integrator );
        const auto t1 = std::chrono::high_resolution_clock::now();

        times_ms.push_back(
            std::chrono::duration<double, std::milli>( t1 - t0 ).count() );

        duration += dt;
        dt       *= DT_GROWTH;
      }

    // -----------------------------------------------------------------------
    // 9. Statistics over all timesteps
    // -----------------------------------------------------------------------
    const double mean = std::accumulate( times_ms.begin(),
                                         times_ms.end(), 0.0 )
                        / static_cast<double>( N_STEPS );

    const double sq_sum = std::inner_product( times_ms.begin(),
                                              times_ms.end(),
                                              times_ms.begin(), 0.0 );
    const double stddev = std::sqrt( sq_sum / static_cast<double>( N_STEPS )
                                     - mean * mean );

    const double best  = *std::min_element( times_ms.begin(), times_ms.end() );
    const double worst = *std::max_element( times_ms.begin(), times_ms.end() );

    std::cout << "\nTestTransient_PressureDiffusionAndFlow performance ("
              << N_STEPS << " timesteps, Release build):\n"
              << "\t best   : " << best   << " ms\n"
              << "\t worst  : " << worst  << " ms\n"
              << "\t mean   : " << mean   << " ms\n"
              << "\t std dev: " << stddev << " ms\n"
              << "\t total  : " << std::accumulate( times_ms.begin(),
                                                    times_ms.end(), 0.0 )
              << " ms\n"
              << "\t duration simulated: " << duration << " s\n";

    std::cout << "\nTestTransient_PressureDiffusionAndFlow: all tests passed.\n";
    
    return duration;

} // end TestTransient_PressureDiffusionAndFlow


/*
{
    // Use huge time increment -> mass term negligible
    constexpr double dt_large = 1.0e+20;
    attorney.TimeIncrement(1.0 / dt_large);
    attorney.EstablishMatrixSetup(region);
    attorney.Accumulate(region);
    attorney.AssignInitialConditions(region);
    attorney.LateAccumulate(region);
    attorney.AssignEssentialConditions(region);

    // 7. Solve linear algebraic system of equations
    attorney.Solve();

    // 8. Write results from the solution vector back to Model
    attorney.OutputResults( region );
                           
    // 9. Calculation of result-dependent properties
//    attorney.PostProcess( region );

    // Check solution is linear between boundary values
    const Index pKey = model.Database().StorageKey("fluid pressure");
    const double p_left  = 3.0e+07;
    const double p_right = 1.0e+07;
    const double L       = 10.0;   // model length

    for ( auto n=region.NodesBegin(); n!=region.NodesEnd(); ++n )
    {
        if ( (*n)->Status(pKey) == DIRICH ) continue;
        const double x        = (*n)->x();
        const double p_exact  = p_left + (p_right - p_left) * x / L;
        const double p_solved = (*n)->Read(pKey);
        _equal( p_solved, p_exact, 1.0e+02 );   // 100 Pa tolerance
    }
}
*/


} // csmp
