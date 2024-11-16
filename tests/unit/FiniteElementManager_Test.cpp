//
//  FiniteElementManager_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 29/6/2024.
//

#include <chrono>
#include <random>
#include "FiniteElementManager_Test.h"
#include "FiniteElementManager.h"
#include "FiniteElementManager_stack_version.h"
#include "vsetMakers.h"
#include "VSet.h"
#include "ModelTopology.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "ConstantFactor.h"
#include "VTK_Interface.h"

#include "PDE_Integrator.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"


using namespace std;

namespace csmp {

void FiniteElementManager_Test::run()
 {
    TestBasicFunctionality();
    TestBasicFunctionality_FiniteElementManager1();
//    TestAccumulationSpeed();
 
 } // end run



/**
   Tests methods of public interface of FiniteElementManager,
   demonstrating that it generates the correct element types
   and that they are accessible
*/
void FiniteElementManager_Test::TestBasicFunctionality()
 {
    // TODO: test InterpolationOrder( uint32_t interpolation_order ); // needs refactoring because it does not reset element types
        
    // 2D MODEL: LINEAR ISOPARAMETRIC ELEMENTS
    uint32_t dimensions{2u}, interpolation_order{1};
    bool isoparametric{ true };
    {
        // calls InitializeElements()
        FiniteElementManager mgr( dimensions, interpolation_order, isoparametric );
        _test( dimensions == mgr.Dimensions() );
        _test( interpolation_order == mgr.InterpolationOrder() );
        _test( isoparametric == mgr.UsesElementsWithLocalCoordinateSystem() );
        // testing that the right elements are assigned and can be retrieved
        _test( mgr.E(ISOPARAMETRIC_LINEAR_BAR)->ElementType() == ISOPARAMETRIC_LINEAR_BAR );
        _test( mgr.LinearBarElement()->ElementType() == ISOPARAMETRIC_LINEAR_BAR );
        _test( mgr.LinearTriangleElement()->ElementType() == ISOPARAMETRIC_LINEAR_TRIANGLE );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_LINEAR_TRIANGLE ) == 3u );
        _test( mgr.ContainsElementType(ISOPARAMETRIC_LINEAR_TRIANGLE) );
        list<CSMP_FEM_TYPE> etypes;
        mgr.CurrentElementTypes(etypes);
        _test( etypes.size() == 3 ); // line, tria, quad
        if ( verbose_ ) mgr.Out();
    }

    // 3D MODEL: QUADRATIC ISOPARAMETRIC ELEMENTS
    dimensions = 3u;
    interpolation_order = 2;
    isoparametric = true;
    {
        // calls InitializeElements()
        FiniteElementManager mgr( dimensions, interpolation_order, isoparametric );
        _test( dimensions == mgr.Dimensions() );
        _test( interpolation_order == mgr.InterpolationOrder() );
        _test( isoparametric == mgr.UsesElementsWithLocalCoordinateSystem() );
        // testing that the right elements are assigned and can be retrieved
        _test( mgr.E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON)->ElementType() == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON );
        const int8_t prism_type = ISOPARAMETRIC_QUADRATIC_PRISM18;
        _test( mgr.E(prism_type)->ElementType() == ISOPARAMETRIC_QUADRATIC_PRISM18 );
        _test( mgr.LinearBarElement() == nullptr );
        _test( mgr.LinearTriangleElement() == nullptr );
        _test( mgr.LinearTetrahedronElement() == nullptr );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) == 27u );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_QUADRATIC_PRISM18 ) == 18u );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) == 14u ); // one node in centre of base plane
        _test( mgr.ContainsElementType(ISOPARAMETRIC_QUADRATIC_TRIANGLE) );
        list<CSMP_FEM_TYPE> etypes;
        mgr.CurrentElementTypes(etypes);
        _test( etypes.size() == 7 ); // line, tria, quad, tet, hex, prism, pyra
        if ( verbose_ ) mgr.Out();
    }

    // 3D MODEL: STANDARD FINITE ELEMENTS
    dimensions = 3u;
    interpolation_order = 1;
    isoparametric = false;
    {
        // calls InitializeElements()
        FiniteElementManager mgr( dimensions, interpolation_order, isoparametric );
        _test( dimensions == mgr.Dimensions() );
        _test( interpolation_order == mgr.InterpolationOrder() );
        _test( isoparametric == mgr.UsesElementsWithLocalCoordinateSystem() );
        // testing that the right elements are assigned and can be retrieved
        _test( mgr.E(LINEAR_TETRAHEDRON)->ElementType() == LINEAR_TETRAHEDRON );
        _test( mgr.LinearBarElement()->ElementType() == LINEAR_BAR );
        _test( mgr.LinearTriangleElement()->ElementType() == LINEAR_TRIANGLE3D );
        _test( mgr.LinearTetrahedronElement()->ElementType() == LINEAR_TETRAHEDRON );
        _test( mgr.NodesOfElementType( LINEAR_CUBOID ) == 8u );
        _test( mgr.NodesOfElementType( LINEAR_TETRAHEDRON ) == 4u );
        list<CSMP_FEM_TYPE> etypes;
        mgr.CurrentElementTypes(etypes);
        _test( etypes.size() == 5 ); // line, tria, quad, tet, hex
        if ( verbose_ ) mgr.Out();
        // testing copy construction and assignment
        FiniteElementManager mgr2(mgr);
        FiniteElementManager mgr3;
        mgr3 = mgr;
        list<CSMP_FEM_TYPE> etypes2, etypes3;
        mgr2.CurrentElementTypes(etypes2);
        mgr3.CurrentElementTypes(etypes3);
        _test( etypes2 == etypes3 );
    }
     
 } // end TestBasicFunctionality(FiniteElementManager)



/**
     Builds templatiised version of finite element manager (FiniteElementManager1) and demonstrates that it contains the correct element types
 */
void FiniteElementManager_Test::TestBasicFunctionality_FiniteElementManager1()
 {
    // 2D MODEL: LINEAR ISOPARAMETRIC ELEMENTS
    {
        const uint32_t dimensions{2u}, interpolation_order{1};
        const bool isoparametric{ true };
        FiniteElementManager1<dimensions,interpolation_order,isoparametric> mgr;
        // testing that the right elements are assigned and can be retrieved
        _test( mgr.E(ISOPARAMETRIC_LINEAR_BAR)->ElementType() == ISOPARAMETRIC_LINEAR_BAR );
        _test( mgr.LinearBarElement()->ElementType() == ISOPARAMETRIC_LINEAR_BAR );
        _test( mgr.LinearTriangleElement()->ElementType() == ISOPARAMETRIC_LINEAR_TRIANGLE );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_LINEAR_TRIANGLE ) == 3u );
        _test( mgr.ContainsElementType(ISOPARAMETRIC_LINEAR_TRIANGLE) );
        list<CSMP_FEM_TYPE> etypes;
        mgr.CurrentElementTypes(etypes);
        _test( etypes.size() == 3 ); // line, tria, quad
        if ( verbose_ ) mgr.Out();
    }

    // 3D MODEL: QUADRATIC ISOPARAMETRIC ELEMENTS
    {
        const uint32_t dimensions = 3u;
        const uint32_t interpolation_order = 2;
        const bool isoparametric = true;
        // calls InitializeElements()
        FiniteElementManager mgr( dimensions, interpolation_order, isoparametric );
        _test( dimensions == mgr.Dimensions() );
        _test( interpolation_order == mgr.InterpolationOrder() );
        _test( isoparametric == mgr.UsesElementsWithLocalCoordinateSystem() );
        // testing that the right elements are assigned and can be retrieved
        _test( mgr.E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON)->ElementType() == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON );
        const int8_t prism_type = ISOPARAMETRIC_QUADRATIC_PRISM18;
        _test( mgr.E(prism_type)->ElementType() == ISOPARAMETRIC_QUADRATIC_PRISM18 );
        _test( mgr.LinearBarElement() == nullptr );
        _test( mgr.LinearTriangleElement() == nullptr );
        _test( mgr.LinearTetrahedronElement() == nullptr );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) == 27u );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_QUADRATIC_PRISM18 ) == 18u );
        _test( mgr.NodesOfElementType( ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) == 14u ); // one node in centre of base plane
        _test( mgr.ContainsElementType(ISOPARAMETRIC_QUADRATIC_TRIANGLE) );
        list<CSMP_FEM_TYPE> etypes;
        mgr.CurrentElementTypes(etypes);
        _test( etypes.size() == 7 ); // line, tria, quad, tet, hex, prism, pyra
        if ( verbose_ ) mgr.Out();
    }

    // 3D MODEL: STANDARD FINITE ELEMENTS
    {
        const uint32_t dimensions = 3u;
        const uint32_t interpolation_order = 1;
        const bool isoparametric = false;
        // calls InitializeElements()
        FiniteElementManager mgr( dimensions, interpolation_order, isoparametric );
        _test( dimensions == mgr.Dimensions() );
        _test( interpolation_order == mgr.InterpolationOrder() );
        _test( isoparametric == mgr.UsesElementsWithLocalCoordinateSystem() );
        // testing that the right elements are assigned and can be retrieved
        _test( mgr.E(LINEAR_TETRAHEDRON)->ElementType() == LINEAR_TETRAHEDRON );
        _test( mgr.LinearBarElement()->ElementType() == LINEAR_BAR );
        _test( mgr.LinearTriangleElement()->ElementType() == LINEAR_TRIANGLE3D );
        _test( mgr.LinearTetrahedronElement()->ElementType() == LINEAR_TETRAHEDRON );
        _test( mgr.NodesOfElementType( LINEAR_CUBOID ) == 8u );
        _test( mgr.NodesOfElementType( LINEAR_TETRAHEDRON ) == 4u );
        list<CSMP_FEM_TYPE> etypes;
        mgr.CurrentElementTypes(etypes);
        _test( etypes.size() == 5 ); // line, tria, quad, tet, hex
        if ( verbose_ ) mgr.Out();
        // testing copy construction and assignment
        FiniteElementManager mgr2(mgr);
        FiniteElementManager mgr3;
        mgr3 = mgr;
        list<CSMP_FEM_TYPE> etypes2, etypes3;
        mgr2.CurrentElementTypes(etypes2);
        mgr3.CurrentElementTypes(etypes3);
        _test( etypes2 == etypes3 );
    }
     
 } // end TestBasicFunctionality_FiniteElementManager1






    /// uses a transient pressure diffusion calculation with tetra-triangle model to measure the speed of the accumulation
 void FiniteElementManager_Test::TestAccumulationSpeed()
   {
      // 0. model construction
      // ---------------------
      VSet<3U>      vset;
      ModelTopology topology = create_FracBox( vset );
      // creating model
      const bool do_not_use_regions_file{true};
      Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", do_not_use_regions_file );
      model.Name("FracBox");

      // 1. Input material properties and initial conditions
      // ---------------------------------------------------
      model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.0) );
      model.InputPropertyValue( "permeability",        makeScalar(PLAIN,1.0e-14) );
      model.InputPropertyValue( "porosity",            makeScalar(PLAIN,0.2) );
      model.InputPropertyValue( "storativity",         makeScalar(PLAIN,2.0e-9) );
      model.InputPropertyValue( "fluid pressure",      makeScalar(PLAIN,1.0e+7) );
      
      Region<3U> fracdomain = model.Region("FRACTURE");
      fracdomain.InputPropertyValue( "permeability", makeScalar(PLAIN,1.0e-8 * 1.0e-3) );

      printRangeOfVariable( model, "permeability" );
      printRangeOfVariable( model, "porosity" );
      printRangeOfVariable( model, "storativity" );

      const double fluid_viscosity(1.0e-03);
      ConstantFactor<3U,divides>  conductivity( model.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model.Apply( conductivity );
      printRangeOfVariable( model, "conductivity" );


      // 2. Assign atmospheric pressure Dirichlet boundary conditions for fluid pressure
      // -------------------------------------------------------------------------------
      // these also serve as initial condition
      for ( auto bit=model.BoundariesBegin(); bit!=model.BoundariesEnd(); ++bit )
        (*bit).second.InputPropertyValue( "fluid pressure", makeScalar(DIRICH,1.0e+5) );
      printRangeOfVariable( model, "fluid pressure" );

      // 3. Variables for transient loop
      // --------------------------------
      const double day(86400.);  // 1 year in seconds
      double       model_time(0.), maxtime(10. * day), time_increment(0.01 * day);
      size_t       timestep(1);

      // 4. Build a transient fluid pressure algorithm using Backward-Euler time-stepping
      // ----------------------------------------------------------------------------------
      // ([C] + dt[K]){p}t+dt = [C]{p}t + dt {Q}t+dt
      PDE_Integrator<3U,Element>  transient_pressure;
    #ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Solver  samg_solver;
      transient_pressure.SetSolver( samg_solver );
    #else
      CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
      transient_pressure.SetSolver( linear_solver );
    #endif

      NumIntegral_dNT_op_dN_dV<3U> conductance( model.Database(), "conductivity",  "fluid pressure", "fluid pressure" );
                                                conductance.MultiplyWithTimeIncrement(true);

      NumIntegral_NT_lhsop_N_dV<3U> capacitance_lhs( model.Database(), "storativity",  "fluid pressure", "fluid pressure" );
                                                 capacitance_lhs.LumpedFormulation(true);

      NumIntegral_NT_op_N_dV<3U> capacitance_rhs( model.Database(), "storativity",  "fluid pressure" );
                                              capacitance_rhs.LumpedFormulation(true);

      NumIntegral_NT_op_N_dV<3U> source( model.Database(), "fluid volume source",  "fluid pressure" );
                                              source.MultiplyWithTimeIncrement(true);
                                              source.AddAccumulateLater();
                                              source.LumpedFormulation(true);

      transient_pressure.Add( &conductance );
      transient_pressure.Add( &capacitance_lhs );
      transient_pressure.Add( &capacitance_rhs );
      transient_pressure.Add( &source );
      transient_pressure.TimeIncrement( time_increment );


      // 11. Transient loop: Compute fluid pressure during each time-step and output the results for each time step
      // ----------------------------------------------------------------------------------------------------------
      cout <<"\n"<<"FiniteElementManager::TestAccumulationSpeed: running pressure diffusion simulation..."<< endl;
 		  auto t0 = chrono::high_resolution_clock::now();
      while ( model_time <= maxtime )
        {
          cout << "\n\nmain: COMPUTING TIMESTEP " << timestep << endl;

          // transient pressure
          model.Apply( transient_pressure );

          // output variables screen
          if ( verbose_ ) printRangeOfVariable( model, "fluid pressure" );

          // Preparing next Time Step
          time_increment *= 1.2;
          model_time     += time_increment;
          timestep++;

          cout << "\nmain: ELAPSED TIME " << model_time / day << " days " << endl;
        }
    
 		  auto t1 = chrono::high_resolution_clock::now();
	    cerr <<"\n\tCPU clock ticks used for pressure diffusion: " << chrono::duration_cast<chrono::nanoseconds>(t1-t0).count();
      cerr  << " nanoseconds." << endl;

      cout <<"\nmain: That's it..."<< endl;




  } // end TestAccumulationSpeed


} // end csmp

