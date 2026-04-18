/**
 * @file StaggeredGridStokesSolver_Example
 * @brief Former part of CSMP++ pore-scale modelling application (2-stage Stokes flow).
 *
 * @details
 * Created by Stephan Matthai on 6/12/2016.  
 * © 2016 Stephan Matthai. All rights reserved.
 *
 * This solver was successfully tested with the following input models:
 *
 * - bcc_0366_5x5x5_raw_quadra_linearized
 * - bcc_0366_5x5x5_raw_quadra
 * - bcc_0366_10x10x10_coarse_quadra_linearized
 * - bcc_0366_10x10x10_coarse_quadra
 * - bcc_0366_10x10x10_raw_quadra_linearized
 * - bcc_0366_10x10x10_raw_quadra
 * - bcc_0366_25x25x25_coarse_quadra_linearized
 * - bcc_0366_25x25x25_coarse_quadra
 * - bcc_0366_25x25x25_raw_quadra_linearized
 * - bcc_0366_25x25x25_raw_quadra
 * - bcc_0366_50x50x50_coarse_quadra_linearized
 * - bcc_0366_50x50x50_coarse_quadra
 * - bcc_0366_50x50x50_raw_quadra_linearized
 * - bcc_0366_50x50x50_raw_quadra
 *
 * - bcc_0366_quadratic_fine_linarized
 * - bcc_0366_quadratic_fine
 * - bcc_0366_basecase_quadra_linearized
 * - bcc_0366_basecase_quadra
 * - bcc_0366_coarse_quadra_linearized_new (== basecase)
 * - bcc_0366_coarse_quadra_new (== basecase)
 * - bcc_0366_coarse2_quadra_linearized
 * - bcc_0366_coarse2_quadra
 * - bcc_0366_2elements_channel_quadra_linearized
 * - bcc_0366_2elements_channel_quadra
 *
 * @note 5 June 2017 — Apoorv Jyoti models:
 * - 500Voxel_Berea_Quadratic_linear
 * - 500Voxel_Berea_Quadratic
 * - MtGambier_500Voxel_Quadratic_Linear
 * - MtGambier_500Voxel_Quadratic
 *
 * @note Berea 100 — Caroline Milliotte:
 * - berea_100x100x100_MP_in_porespace_quadra_linearized
 * - berea_100x100x100_MP_in_porespace_quadra
 *
 * @note Berea subvolume 100 — Chloe Burney:
 * - berea_sub100_quadra_linearized
 * - berea_sub100_quadra
 * - berea_sub100_20_quadra_linearized (seed elements: 20)
 * - berea_sub100_20_quadra
 *
 * @note ANLEC 100:
 * - anlec100_quadra_linearized
 * - anlec100_quadra
 *
 * @note Tubes (~60MB):
 * - tubes_all_quadra_linearized
 * - tubes_all_quadra
 */
#include "StaggeredGridStokesSolver_Example.h"

#include "ModelTopology.h"
#include "VSet.h"
#include "ANSYS_Model3D.h"
#include "ANSYS_Interface.h"
#include "VTU_Interface.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "IsoparametricQuadraticTetrahedron.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

#include "PDE_Integrator.h"
#include "IsoparametricQuadraticTetrahedron.h"
// analytic integration
#include "LHS_Integral_dNT_dN_dV.h"
#include "Integral_NT_op_N_dV.h"
#include "Integral_SetRHS_to_Zero.h"
// numeric integration
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_DNi_rhsop_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "StatisticalAnalyzer.h"

using namespace std;
using namespace csmp;

namespace csmp {

void StaggeredGridStokesSolver_Example::Specifications() {
  SetTitle("Staggered-grid 2-step Stokes equation solver" );
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "Stephan Matthai, 15 Dec. 2016." );
  AddDescription( "Sequential pressure (piecewise linear) - velocity (piecewise quadratic) FEM solver for incompressible Stokes lubrication equation." );
  AddDescription( "Note: uses node-matched staggered tetrahedral element meshes as input; each quadratic tetrahedron consists of 10 linear tetrahedra." );
  AddRequirement( "1. input:  mesh-files from ANSYS (*.asc and *.dat), created using CSP-output interface;" );
  AddRequirement( "2 sets of files containing linear and quadratic meshes, respectively." );
  AddRequirement( "3. Input models: must be box-shaped, consisting of contiguous mesh domains. \
                      Side boundaries must be labeled LEFT, RIGHT etc., see CSMP User's guide." );
  AddRequirement( "Models must contain a region that denotes the pore space called 'PORES' \
                   and the intersection curves of the grain boundaries with the box boundaries." );
  AddRequirement( "This line-element region must be called 'PORE_EDGE'." );
  AddRequirement( "Input model in example data that fulfills these needs is 'tubes'." );
}


static size_t countDirichletNodes( const Model<3U>& model, const char* var ) {
   const csmp::Index var_key = model.Database().StorageKey(var);
   assert( var_key.place == NODE );
   size_t counter{0lu};
   const Region<3U>& model_domain = model.Region("Model");
   for ( const auto& nd : model_domain.NodeVector() )
     if ( nd->Status(var_key) == DIRICH )
       counter++;
   return counter;
}

void StaggeredGridStokesSolver_Example::Run()
 {
    auto& csmp_error = ErrorHandler::Instance();
    
    const int  dim(3U);
    const bool verbose(true);

// EXAMPLE MODELS
// Fails as CNR1 cannot be identified
//    const string input_model_linear("tubes_all_quadra_linearized");
//    const string input_model_quadratic("tubes_all_quadra");
  
// Fails as CNR1 cannot be identified
    const string input_model_linear("bcc_0366_coarse_quadra_linearized_new");
    const string input_model_quadratic("bcc_0366_coarse_quadra_new");

// Only has BACK and FRONT (Z-axis flow from back to front)
//    string input_model_linear("spiral_quadra_linearized");
//    const string input_model_quadratic("spiral_quadra");

    // VTU-file output variables of to document initial configuration
    list<string> outputProps_before;
    outputProps_before.push_back("fluid pressure");
    outputProps_before.push_back("nodal velocity x");
    outputProps_before.push_back("nodal velocity y");
    outputProps_before.push_back("nodal velocity z");
    // boundary flags (enum VARIABLE_FLAG) of the nodes output to integer values
    if ( verbose ) {
         outputProps_before.push_back("node-flag velocity x");
         outputProps_before.push_back("node-flag velocity y");
         outputProps_before.push_back("node-flag velocity z");
      }
    
    // VTU-file output variables for results
    list<string> outputProps_after;
    outputProps_after.push_back("fluid pressure");
    outputProps_after.push_back("nodal velocity");
    outputProps_after.push_back("nodal velocity x");
    outputProps_after.push_back("nodal velocity y");
    outputProps_after.push_back("nodal velocity z");


    const bool        binary_file(true), isoparametric(true); // analytically integrated LFEM
    vector<Point<3> > linmodel_node_coords; // node coordinates in VSet order to re-establish original node numbering if necessary
    Point<3>          xyz_min, xyz_max;
    const string      pore_space("PORES");
    const string      pore_edge("PORE_EDGE");
    const double      viscosity(1.0e-3);
    const double      patm(100325.);
    const double      hydrostatic_grad(1.e+6 * 9.80665); // hydrostatic pressure gradient to drive flow

    // SAMG solver and integrator for linear algebraic system
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    SAMG_Solver   solver(&settings);
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
    // pressure map of linear model
    map<Point<3U>,array<double,3>>  pressure_linmodel;


    // ===========================================================================================================================
    // A. First linear FEM model (LFEM) to solve the pressure equation
    // ===========================================================================================================================
    {
      // including the reduction of regions according to the regions file (vars file must contain node-number for reconstruction
      ANSYS_Model3D   linear_model( isoparametric, input_model_linear.c_str(), "staggered_grid_stokes_variables_p.txt" );
      linear_model.MinMaxCoordinates( xyz_min, xyz_max ); // of the bounding box
      printModelDimensions( linear_model );
      linear_model.RestoreOriginalNodeNumbering(true);
      linmodel_node_coords.assign( linear_model.VerticesBegin(), linear_model.VerticesEnd() );
      
      linear_model.RegionsOut();
      linear_model.BoundariesOut();
     
      // FEM assembly for pressure equation: Laplacian = 0
      PDE_Integrator<dim,Element>  pressure_solver(solver);
      /* exactly integrated elements
      LHS_Integral_dNT_dN_dV<dim>  laplacian( linear_model.Database(), "fluid pressure", "fluid pressure" );
      Integral_SetRHS_to_Zero<dim> zero( linear_model.Database(), "fluid pressure" );
      */
      /* isoparametric elements */
      NumIntegral_dNT_dN_dV<dim>       laplacian( linear_model.Database(), "fluid pressure", "fluid pressure" );
      NumIntegral_SetRHS_to_Zero<dim>  zero( linear_model.Database(), "fluid pressure" );
      pressure_solver.Add( &laplacian );
      pressure_solver.Add( &zero );
     
      // flow domain (subregion of the model that must be called PORES)
      Region<dim>& pore_domain = linear_model.Region(pore_space);

      VTU_Interface<3> vtu_linear( linear_model );

      // 1. LEFT-RIGHT fluid-pressure computation: pressure gradient in X direction
      // --------------------------------------------------------------------------
      double length_x = xyz_max[0] - xyz_min[0];
      // boundary conditions
      linear_model.InputPropertyValue("fluid pressure", makeScalar(ANY,0.)); // to visualise model with boundary conditions only
      linear_model.Boundary("LEFT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * length_x) );
      linear_model.Boundary("RIGHT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
      vtu_linear.OutputDataToVTU( string("fluid-pressure-x_pboundaries"), string("fluid pressure"), string("Model"), 1 /* 1 = solution stage of 2-stage solver */ );
      cout <<"\n\t"<<"px: Nodes with Dirichlet constraints on pressure: "<< countDirichletNodes( linear_model, "fluid pressure") << endl;
      
      // computation of 'fluid pressure' on 'pore_domain'
      pressure_solver.IntegrateOver( pore_domain );
      // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pressure_solver.Reset();
      pore_domain.ChangePropertyStatus( "fluid pressure", ANY, PERIMETER );
      
      // output
      printRangeOfVariable( linear_model, pore_space.c_str(), "fluid pressure" );
      vtu_linear.OutputDataToVTU( string("fluid-pressure-x"), string("fluid pressure"), string("PORES"), 1 );
      
      // output to map of arrays 'pressure_linmodel'
      ///// Pressure mapping of linear model
      const auto nodesEnd_linmodel(linear_model.Mesh().NodesEnd());
       const csmp::Index p_key (linear_model.Database().StorageKey("fluid pressure"));
      for ( auto nit=linear_model.Mesh().NodesBegin(); nit!=nodesEnd_linmodel; ++nit  ) {
        pressure_linmodel.insert( make_pair( nit->Coordinate(), array<double,3>{ nit->Read(p_key),0.,0.} ) );
      }
      assert( linear_model.Mesh().Nodes() == pressure_linmodel.size() );
      cerr  <<"\nnumber of nodes in pressure map of linear model, X direction: "<< pressure_linmodel.size() <<endl;
      cout<< LinearModelElementVolume(linear_model, pore_space.c_str())<<endl;
     
     
      // 2. BOTTOM-TOP fluid-pressure computation: pressure gradient in Y direction
      // --------------------------------------------------------------------------
      double length_y = xyz_max[1] - xyz_min[1];
      // boundary conditions
      linear_model.InputPropertyValue("fluid pressure", makeScalar(ANY,0.));
      linear_model.Boundary("BOTTOM").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * length_y) );
      linear_model.Boundary("TOP").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
      vtu_linear.OutputDataToVTU( string("fluid-pressure-y_pboundaries"), string("fluid pressure"), string("Model"), 1 );
      cout <<"\n\t"<<"py: Nodes with Dirichlet constraints on pressure: "<< countDirichletNodes( linear_model, "fluid pressure") << endl;
      // PDE integrator was reset because boundary conditions & DOF changed, so it needs to be rebuild here
      pressure_solver.Add( &laplacian );
      pressure_solver.Add( &zero );
      // computation
      pressure_solver.IntegrateOver( pore_domain );
      // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pressure_solver.Reset();
      pore_domain.ChangePropertyStatus( "fluid pressure", ANY, PERIMETER );
      // output
      printRangeOfVariable( linear_model, pore_space.c_str(), "fluid pressure" );
      vtu_linear.OutputDataToVTU( string("fluid-pressure-y"), string("fluid pressure"), string("PORES"), 1 );
      // output to FEM_Data
      ///// Pressure mapping of linear model
      for ( auto nit=linear_model.Mesh().NodesBegin(); nit!=nodesEnd_linmodel; ++nit  ) {
          auto it = pressure_linmodel.find( nit->Coordinate() );
          if ( it != pressure_linmodel.end() ) (*it).second[1] = nit->Read(p_key);
      }
      assert( linear_model.Mesh().Nodes() == pressure_linmodel.size() );
      cerr  <<"\nnumber of nodes in pressure map of linear model, Y direction: "<< pressure_linmodel.size();
     
     
      // 3. BACK-FRONT fluid-pressure computation: pressure gradient in Z direction
      // --------------------------------------------------------------------------
      double length_z = xyz_max[2] - xyz_min[2];
      // boundary conditions
      linear_model.InputPropertyValue("fluid pressure", makeScalar(ANY,0.));
      linear_model.Boundary("BACK").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * length_z) );
      linear_model.Boundary("FRONT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
      vtu_linear.OutputDataToVTU( string("fluid-pressure-z_pboundaries"), string("fluid pressure"), string("Model"), 1 );
      cout <<"\n\t"<<"pz: Nodes with Dirichlet constraints on pressure: "<< countDirichletNodes( linear_model, "fluid pressure") << endl;
      // PDE integrator was reset because boundary conditions & DOF changed, so it needs to be rebuild here
      pressure_solver.Add( &laplacian );
      pressure_solver.Add( &zero );
      // computation
      pressure_solver.IntegrateOver( pore_domain );
      // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      pressure_solver.Reset();
      pore_domain.ChangePropertyStatus( "fluid pressure", ANY, PERIMETER );
      // output
      printRangeOfVariable( linear_model, pore_space.c_str(), "fluid pressure" );
      vtu_linear.OutputDataToVTU( string("fluid-pressure-z"), string("fluid pressure"), string("PORES"), 1 );
      // output to FEM_Data
      ///// Pressure mapping of linear model
      for ( auto nit=linear_model.Mesh().NodesBegin(); nit!=nodesEnd_linmodel; ++nit  ) {
          auto it = pressure_linmodel.find( nit->Coordinate() );
          if ( it != pressure_linmodel.end() ) (*it).second[2] = nit->Read(p_key);
      }
      assert( linear_model.Mesh().Nodes() == pressure_linmodel.size() );
      cerr  <<"\nnumber of nodes in pressure map of linear model, Z direction: "<< pressure_linmodel.size()<<endl;
     
      cout<<"Linear model, Number of interior nodes= "<<pore_domain.InteriorNodes()<<endl;
      // output of the model to CSMP native format (for future runs)
      string binary_output_file_set_name(linear_model.Name());
      binary_output_file_set_name += "_csmp_LFEM";
      
      // TODO: work with the binary after doing the conversion only once
      linear_model.OutputToBinaryFile( binary_output_file_set_name.c_str() );
           
    } // end of LFEM computations
   
   
   
    // ===========================================================================================================================
    // B. Second model (quadratic FEM) to solve the Stokes-flow equation
    // ===========================================================================================================================
    // output from the Ansys mesh generator via the CSMP interface
    ANSYS_Model3D  quadratic_model( input_model_quadratic.c_str(),
                                   "staggered_grid_stokes_variables_u.txt",
                                    binary_file );

    recreateBoxBoundaryFlags( quadratic_model );
   
    Region<dim>&  pore_domain2(quadratic_model.Region(pore_space));
    Region<dim>&  grain_edges(quadratic_model.Region(pore_edge));
    assert( quadratic_model.Mesh().Nodes() == linmodel_node_coords.size() );
    
    VTU_Interface<3> vtu_quadratic( quadratic_model );

    // 0. Placing the element integration points inside of the quadratic tetrahedra
    // ----------------------------------------------------------------------------
    // (by default 4-point rule, Zienkiewicz and Taylor, Vol I, page 223, the gauss-quadrature points are already away from the corners)
    //dynamic_cast<IsoparametricQuadraticTetrahedron*>( quadratic_model.FE_Manager().E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON) )->QuadratureRules_11Points();

    // 1. QFEM algorithm: mu laplacian u = grad p
    // ------------------------------------------  
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings2;
    // settings usually used to compute solution for quadratic FEM
    settings.Set_napproach(2); // interpolation separate for each unknown
    settings2.Set_ndefault(26);
    SAMG_Solver   solver2(&settings2);
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver2;
#endif
    PDE_Integrator<dim,Element>  velocity_solver( solver2 );

/*  
    // settings for variable-by-variable smoothing
    // DOES NOT CONVERGENCE
    settings.Set_napproach(3); // interpolation separate for each unknown
    settings.Set_nxtyp(1);     // ILU relaxation
    settings.Set_internal(0);  // primary matrix is user defined
    settings.Set_nprim(1);
    // ncyc
    settings.Set_igam(1);
    settings.Set_ncgrad(3);
    settings.Set_nkdim(9);
    settings.Set_ncycle(50);
*/

    NumIntegral_dNT_op_dN_dV<dim> viscosity_matr_x( quadratic_model.Database(), "viscosity", "nodal velocity x", "nodal velocity x" ),
                                  viscosity_matr_y( quadratic_model.Database(), "viscosity", "nodal velocity y", "nodal velocity y" ),
                                  viscosity_matr_z( quadratic_model.Database(), "viscosity", "nodal velocity z", "nodal velocity z" );
                                  
    NumIntegral_DNi_rhsop_dV<dim> gradient_x( quadratic_model.Database(), "fluid pressure", "nodal velocity x"),
                                  gradient_y( quadratic_model.Database(), "fluid pressure", "nodal velocity y"),
                                  gradient_z( quadratic_model.Database(), "fluid pressure", "nodal velocity z");

    gradient_x.SpatialDerivative(X_DIRECTION);
    gradient_y.SpatialDerivative(Y_DIRECTION);
    gradient_z.SpatialDerivative(Z_DIRECTION);
    // required to get flow from high to low pressure
    gradient_x.MultiplyBy( -1. );
    gradient_y.MultiplyBy( -1. );
    gradient_z.MultiplyBy( -1. );
    // PDE operators are added within x,y,z loop

    // to get velocity histograms for each flow direction
    StatisticalAnalyzer<3U>  modelAnalyzer( quadratic_model );
    HistogramBins            velocityRanges;
    // specify via ascii input file, the x-axis range of the histogram columns, i.e. the bin size
    modelAnalyzer.DefineBins( "pore_velocity_bins.bins", velocityRanges );
    // histogram output
    map<string,pair<HistogramBins,size_t> > velo_histograms;

    // fluid viscosity
    quadratic_model.InputPropertyValue( "viscosity", makeScalar(ANY,viscosity) );
    // model porosity
    const double  porosity = pore_domain2.Volume() / quadratic_model.Region("Model").Volume();
    cout <<"\n\n\nmodel: '"<<  quadratic_model.Name() <<"', porosity(computed): "<< std::fixed << std::setprecision(6) << porosity <<"\n";
    cout << std::scientific << std::setprecision(5);
    // permeability as measured in the 3 coordinate directions
    vector<double> k_equiv(3,0.);


    // 2. Computation of flow velocities from the XYZ axes-aligned pressure gradients
    // ------------------------------------------------------------------------------
    const size_t X_DIRECTION{0u}, Y_DIRECTION{1u}, Z_DIRECTION{2u};
    BOX_BOUNDARY inlet(LEFT), outlet(RIGHT);
    const auto nodesEnd_quadrmodel(quadratic_model.Mesh().NodesEnd());
    const Index nodeKey_quadrmod( quadratic_model.Database().StorageKey( "fluid pressure") );
    const auto nodesEnd_pressure_linmodel = pressure_linmodel.end();
    assert( quadratic_model.Mesh().Nodes() == pressure_linmodel.size() );

    for ( size_t flow_direction = 0; flow_direction < 3; flow_direction++ )
      {
        // 2.1 assigning pre-calculated pressure values from linear to quadratic model
        // ---------------------------------------------------------------------------
        // by lambda function
        auto assign_pressure_component = [&](size_t component, const std::string& suffix) {
            for ( auto nit = quadratic_model.Mesh().NodesBegin(); nit != nodesEnd_quadrmodel; ++nit ) {
                auto lin_it = pressure_linmodel.find( (*nit).Coordinate() );
                if ( lin_it != nodesEnd_pressure_linmodel ) {
                    (*nit).Store( nodeKey_quadrmod, makeScalar((*nit).Status(nodeKey_quadrmod), (*lin_it).second[component] ) );
                } else {
                    (*nit).Coordinate().Out();
                    csmp_error.Note( ERROR, "quadratic_model:", "node from linear model not found in quadratric model; do you want to continue?" );
                    continue;
                }
            }
            vtu_quadratic.OutputDataToVTU( std::string("fluid-pressure-quadratic-") + suffix, "fluid pressure", "PORES", 2 /* solution stage */ );
        };

        // dispatch lambda based on flow direction
        if      (flow_direction == X_DIRECTION) assign_pressure_component(0, "x");
        else if (flow_direction == Y_DIRECTION) assign_pressure_component(1, "y");
        else                                    assign_pressure_component(2, "z");

       // output assigned values for checking
       cout <<"\n\n"<<"StaggeredGridStokesSolver_Example::Run: analysing permeability / flow velocity distribution in "<< parse(static_cast<SPATIAL_DERIVATIVE>(flow_direction));
       cout <<" direction for the given boundary pressure differential...\n";
       printRangeOfVariable( quadratic_model, pore_space.c_str(), "fluid pressure" );
       // output for testing
       if ( verbose ) vtu_quadratic.OutputDataToVTU( string("quadratic_model-interpolated-fluid-pressure"),
                                                     string("fluid pressure"),
                                                     string(pore_space), flow_direction+1 );

       // 2.2 no-slip boundary conditions except for inlet and outlet
       // ---------------------------------------------------------------------------
       AssignNoSlipBoundaryConditions( quadratic_model, pore_domain2, grain_edges, flow_direction );

       // testing of boundary conditions
       if ( verbose ) {
            flagToNumber( quadratic_model, "nodal velocity x", "node-flag velocity x" );
            flagToNumber( quadratic_model, "nodal velocity y", "node-flag velocity y" );
            flagToNumber( quadratic_model, "nodal velocity z", "node-flag velocity z" );
            vtu_quadratic.OutputDataToVTU( "input", outputProps_before, pore_space.c_str(), static_cast<long>(flow_direction+1) );
         }
       // adding integrals to PDE_Integrator
       velocity_solver.Add( &viscosity_matr_x );
       velocity_solver.Add( &viscosity_matr_y );
       velocity_solver.Add( &viscosity_matr_z );
       velocity_solver.Add( &gradient_x );
       velocity_solver.Add( &gradient_y );
       velocity_solver.Add( &gradient_z );
       // computation
       velocity_solver.IntegrateOver( pore_domain2 );
       // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       velocity_solver.Reset();
       ConstructVelocityVector( quadratic_model );
       printRangeOfVariable( quadratic_model, pore_space.c_str(), "nodal velocity" );

       // 2.3 analysis of velocity distribution and permeability in flow direction
       // ---------------------------------------------------------------------------
       if      ( flow_direction == Y_DIRECTION ) {  inlet = BOTTOM; outlet = TOP; }
       else if ( flow_direction == Z_DIRECTION ) { inlet = BACK; outlet = FRONT; }
       double flux_through_model;
       k_equiv[flow_direction] = EquivalentPermeabilityQuadraticFEM( quadratic_model, pore_space, inlet, outlet, flux_through_model );
       //                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       cout <<"\n\n"<<"StaggeredGridStokesSolver_Example::Run: k_equivalent (old-style computation) ("<< flow_direction <<"): "<< k_equiv[flow_direction] <<" m2.\n";

       modelAnalyzer.RegionPropertyHistogramsIntegrationPoint("nodal velocity", velocityRanges, velo_histograms, pore_space );
       // writes Maple statlist format with weights; results must be copied and pasted into Maple workbook
       //const bool log10_of_bin_values(true);
       const bool log10_of_bin_values(false);
       modelAnalyzer.OutputRegionPropertyHistogramsMaple( (string("nodal velocity") + to_string(flow_direction+1)).c_str(),
                                                          velocityRanges, velo_histograms, log10_of_bin_values );
 
       vtu_quadratic.OutputDataToVTU( string("velocity"), outputProps_after, string("Model"), static_cast<long>(flow_direction+1) );
        
      } // end flow_direction loop


      cout<<"Quadratic model, Number of interior nodes= "<<pore_domain2.InteriorNodes()<<endl;
      PrintAnalysisResults( quadratic_model, porosity, k_equiv );

      // output of quadratic model to CSMP native format (for future runs)
      string binary_output_file_set_name2(quadratic_model.Name());
      binary_output_file_set_name2 += "_csmp_QFEM";
      quadratic_model.OutputToBinaryFile( binary_output_file_set_name2.c_str() );
   
 } // end Run



// debugging code: OK
/*
{
double gmin, gmax;
set<VARIABLE_FLAG> bflags;
const csmp::Index var_key = linear_model.Database().StorageKey("fluid pressure");
{
  Boundary<3U> boundary = linear_model.Boundary("BOTTOM");
  boundary.MinMaxOf( "fluid pressure", gmin, gmax );
  cout <<"\nrun: pf at BOTTOM: "<< gmin <<" - "<< gmax << endl;
  for( const auto& bit : boundary.NodeVector() ) bflags.insert( bit->Status( var_key ) );
  cout <<"\nrun: variable flags at bottom: ";
  for ( const auto& cit : bflags ) cout <<"  "<< parseStatus(cit);
  cout << endl;
  bflags.clear();
}
{
  Boundary<3U> boundary = linear_model.Boundary("TOP");
  boundary.MinMaxOf( "fluid pressure", gmin, gmax );
  cout <<"\nrun: pf at TOP: "<< gmin <<" - "<< gmax << endl;
  for( const auto& bit : boundary.NodeVector() ) bflags.insert( bit->Status( var_key ) );
  cout <<"\nrun: variable flags at top: ";
  for ( const auto& cit : bflags ) cout <<"  "<< parseStatus(cit);
  cout << endl;
  bflags.clear();
}
}
*/



/**
    Composes the 3D vector variable velocity from its components all across the
    entire domain "Model"
    
    @attention velocity components x, y, z must have been initialised everywhere
    before this function is called.
*/
void StaggeredGridStokesSolver_Example::ConstructVelocityVector( Model<3>& mdl )
{
	VectorVariable<3> v;
	const csmp::Index vx_key(mdl.Database().StorageKey("nodal velocity x")),
		                vy_key(mdl.Database().StorageKey("nodal velocity y")),
		                vz_key(mdl.Database().StorageKey("nodal velocity z")),
		                v_key(mdl.Database().StorageKey("nodal velocity"));

	Region<3>& mref = mdl.Region("Model");

	for ( auto& nit : mref.NodeVector() )
    {
       v(0) = (*nit).Read(vx_key);
       v(1) = (*nit).Read(vy_key);
       v(2) = (*nit).Read(vz_key);
       (*nit).Store(v_key, v);
    }
}



/**
     AssignNoSlipBoundaryConditions:  Assigns
     
     1. No slip conditions to grain edges
 
     2. Box-boundary parallel slip conditions to model side boundaries and edges
 
     @param flow_direction X=0, Y=1, Z=2
 
     @attention Method does not correct for ill-configured meshes where all the nodes
     of boundary triangles touch grain edges.
     
     @todo SKM distinguish properly between value and flag changes
     
     tested:  22-11-2016 SKM
*/
void StaggeredGridStokesSolver_Example::AssignNoSlipBoundaryConditions( Model<3U>& model,
                                                                         Region<3U>& stokes_flow_domain,
                                                                         Region<3U>& grain_edges,
                                                                         size_t flow_direction )
 {
    // flow direction 0=x
    string inlet("LEFT"), outlet("RIGHT");
    if ( flow_direction == 1 ) { // Y-direction
         inlet="BOTTOM"; outlet="TOP";
      }
    else if ( flow_direction == 2 ) { // Z-direction with increasing coordinates
         inlet="BACK"; outlet="FRONT";
      }
    assert( flow_direction <=2 );
   
    // --------------------------------------------------------
    // 1. Zapping pre-existing conditions
    // --------------------------------------------------------
    // by default all velocity scalars are flagged ANY
    const double  slip(0.);
    stokes_flow_domain.InputPropertyValue("nodal velocity x", makeScalar(DIRICH,slip), PERIMETER );
    stokes_flow_domain.InputPropertyValue("nodal velocity y", makeScalar(DIRICH,slip), PERIMETER );
    stokes_flow_domain.InputPropertyValue("nodal velocity z", makeScalar(DIRICH,slip), PERIMETER );

    // opening the inflow and outflow boundaries
    // -----------------------------------------
    ChangeBoxBoundaryFlags( model, parseBoundary(inlet), "nodal velocity x", ANY );
    ChangeBoxBoundaryFlags( model, parseBoundary(inlet), "nodal velocity y", ANY );
    ChangeBoxBoundaryFlags( model, parseBoundary(inlet), "nodal velocity z", ANY );
    ChangeBoxBoundaryFlags( model, parseBoundary(outlet), "nodal velocity x", ANY );
    ChangeBoxBoundaryFlags( model, parseBoundary(outlet), "nodal velocity y", ANY );
    ChangeBoxBoundaryFlags( model, parseBoundary(outlet), "nodal velocity z", ANY );
    // tested: OK
    // nodesToFile( model, inlet.c_str(), PERIMETER );
    // nodesToFile( model, inlet.c_str(), INTERIOR );

    // ---------------------------------------------------------
    // 2. Assigning flow-direction dependent boundary conditions
    // ---------------------------------------------------------
    const csmp::Index  nvx_key(model.Database().StorageKey("nodal velocity x"));
    const csmp::Index  nvy_key(model.Database().StorageKey("nodal velocity y"));
    const csmp::Index  nvz_key(model.Database().StorageKey("nodal velocity z"));

    // LEFT->RIGHT flow, LEFT = inflow boundary = inlet
    // --------------------------------------------------------
    // setting up the side boundaries TOP, BOTTOM, FRONT BACK to allow boundary-parallel flow

    if ( inlet == "LEFT" ) {
         // no-flow perpendicular to these boundaries
         // -----------------------------------------
         // top (only y-Dirich constraint is retained)
         model.Boundary("TOP").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("TOP").ChangePropertyStatus( "nodal velocity z", ANY );
         // bottom
         model.Boundary("BOTTOM").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("BOTTOM").ChangePropertyStatus( "nodal velocity z", ANY );
         // front
         model.Boundary("FRONT").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("FRONT").ChangePropertyStatus( "nodal velocity y", ANY );
         // back
         model.Boundary("BACK").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("BACK").ChangePropertyStatus( "nodal velocity y", ANY );
         // re-imposing no-slip constraints on side edges
         for ( auto it=stokes_flow_domain.PerimeterNodesBegin(); it!=stokes_flow_domain.NodesEnd(); ++it ) {
              const BOX_BOUNDARY edge( (*it)->AtBoundary() );
              if ( edge == EDGE1 || edge == EDGE3 || edge == EDGE9 || edge == EDGE11 || (edge <= CNR1 &&edge >= CNR8) ) {
                   (*it)->Status( nvy_key, DIRICH );
                   (*it)->Status( nvz_key, DIRICH );
               }
           }
      }

    else if ( inlet == "BOTTOM" ) {
         model.Boundary("LEFT").ChangePropertyStatus( "nodal velocity y", ANY );
         model.Boundary("LEFT").ChangePropertyStatus( "nodal velocity z", ANY );
         model.Boundary("RIGHT").ChangePropertyStatus( "nodal velocity y", ANY );
         model.Boundary("RIGHT").ChangePropertyStatus( "nodal velocity z", ANY );
         model.Boundary("FRONT").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("FRONT").ChangePropertyStatus( "nodal velocity y", ANY );
         model.Boundary("BACK").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("BACK").ChangePropertyStatus( "nodal velocity y", ANY );
         for ( auto it=stokes_flow_domain.PerimeterNodesBegin(); it!=stokes_flow_domain.NodesEnd(); ++it ) {
              const BOX_BOUNDARY edge( (*it)->AtBoundary() );
              if ( edge == EDGE2 || edge == EDGE4 || edge == EDGE10 || edge == EDGE12 || (edge <= CNR1 &&edge >= CNR8) ) {
                   (*it)->Status( nvx_key, DIRICH );
                   (*it)->Status( nvz_key, DIRICH );
                }
           }
      }

    else if ( inlet == "BACK" ) {
         model.Boundary("TOP").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("TOP").ChangePropertyStatus( "nodal velocity z", ANY );
         model.Boundary("BOTTOM").ChangePropertyStatus( "nodal velocity x", ANY );
         model.Boundary("BOTTOM").ChangePropertyStatus( "nodal velocity z", ANY );
         model.Boundary("LEFT").ChangePropertyStatus( "nodal velocity y", ANY );
         model.Boundary("LEFT").ChangePropertyStatus( "nodal velocity z", ANY );
         model.Boundary("RIGHT").ChangePropertyStatus( "nodal velocity y", ANY );
         model.Boundary("RIGHT").ChangePropertyStatus( "nodal velocity z", ANY );
         for ( auto it=stokes_flow_domain.PerimeterNodesBegin(); it!=stokes_flow_domain.NodesEnd(); ++it ) {
              const BOX_BOUNDARY edge( (*it)->AtBoundary() );
              if ( edge == EDGE5 || edge == EDGE6 || edge == EDGE7 || edge == EDGE8 || (edge <= CNR1 &&edge >= CNR8) ) {
                   (*it)->Status( nvx_key, DIRICH );
                   (*it)->Status( nvy_key, DIRICH );
                }
           }
      }

    // ------------------------------------------------------------
    // 3. No-slip boundary conditions at grain edges
    // ------------------------------------------------------------
    grain_edges.ChangePropertyStatus( "nodal velocity x", DIRICH );
    grain_edges.ChangePropertyStatus( "nodal velocity y", DIRICH );
    grain_edges.ChangePropertyStatus( "nodal velocity z", DIRICH );

/* TODO: deal with isolated points touching the boundary which belong to mineral surfaces
    throw csmp::Exception( ERROR, "assignNeumannDirichletConditions", "other cases not handled yet.");

    if ( !stokes_domain_points.empty() ) {
        const csmp::Index vx_key(model.Database().StorageKey("nodal velocity x"));
        const csmp::Index vy_key(model.Database().StorageKey("nodal velocity y"));
        const csmp::Index vz_key(model.Database().StorageKey("nodal velocity z"));
        for ( auto nit=stokes_domain_points.begin(); nit!=stokes_domain_points.end(); ++nit ) {
             (*nit)->Store( vx_key, makeScalar(DIRICH, 0.) );
             (*nit)->Store( vy_key, makeScalar(DIRICH, 0.) );
             (*nit)->Store( vz_key, makeScalar(DIRICH, 0.) );
          }
      }
*/
   
 } // end assignNoSlipBoundaryConditions





/**
    Assuming a box-shaped model, function calculates the equivalent permeability parallel to 
    its side boundaries.
*/
double StaggeredGridStokesSolver_Example::EquivalentPermeabilityQuadraticFEM( Model<3U>& model,
                                                                              const string& region,
                                                                              BOX_BOUNDARY inflow_boundary,
                                                                              BOX_BOUNDARY outflow_boundary,
                                                                              double& flux_through_model )
 {
    const uint32_t DIM(3U);

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    const Region<DIM>&  flow_domain(model.Region(region));
   
    // 1. calculating the cross-sectional area of the model
    // ---------------------------------------------------
    // (2 options: left-right and top bottom)
    Point<DIM>  xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max );
    const double length_x = xyz_max[0] - xyz_min[0];
    const double length_y = xyz_max[1] - xyz_min[1];
    const double length_z = xyz_max[2] - xyz_min[2];
    // assumption is that flow boundaries are opposite to one another
    double length{ numeric_limits<double>::quiet_NaN() }, area{ numeric_limits<double>::quiet_NaN() };
    if ( inflow_boundary == LEFT ) {
         length = length_x;
         area   = length_y * length_z;
      }
    else if ( inflow_boundary == BOTTOM ) {
         length = length_y;
         area   = length_x * length_z;
      }
    else if ( inflow_boundary == BACK ) {
         length = length_z;
         area   = length_x * length_y;
      }
    double xsect_area(length * length);
    cout<<"\n\t"<<" length (m): "<< length <<", cross-sectional area (m2) of inflow boundary."<< endl;
   
    // 2. calculating the far-field fluid pressure gradient
    // -----------------------------------------------------
    // O.K.
    const csmp::Index pf_key = model.Database().StorageKey("fluid pressure");
    double inflow_p(0.), outflow_p(0.);
    size_t in_nodes(0U), out_nodes(0U);
    // for the entire flow domain
    for ( auto nit=flow_domain.PerimeterNodesBegin(); nit!=flow_domain.NodesEnd(); ++nit )
      {
         if ( (*nit)->AtBoundary() == inflow_boundary ) {
              inflow_p += (*nit)->Read( pf_key );
              in_nodes++;
           }
         else if ( (*nit)->AtBoundary() == outflow_boundary ) {
              outflow_p += (*nit)->Read( pf_key );
              out_nodes++;
           }
      }
    inflow_p  /= static_cast<double>(in_nodes);
    outflow_p /= static_cast<double>(out_nodes);
   
    double farfield_pf_gradient = (inflow_p - outflow_p) / length;
    cout<<"\n\t"<<"inflow pressure: "<<inflow_p<<", outflow pressure: "<<outflow_p<<endl;

	  // 3a. fluid throughput through model as determined by velocity projections on model boundaries
    // -------------------------------------------------------------------------------------------
    const csmp::Index vt_key = model.Database().StorageKey("nodal velocity");
    double inlet_area, outlet_area;
    const double  inflow  = fabs(boxBoundaryFluxFEM_Quadratic<3>( model.Boundary(parseBoundary(inflow_boundary)), vt_key, inlet_area ));
    const double  outflow = fabs(boxBoundaryFluxFEM_Quadratic<3>( model.Boundary(parseBoundary(outflow_boundary)), vt_key, outlet_area ));

    // computing the mean of the in- and outflow
    flux_through_model = (fabs(inflow) + fabs(outflow)) / 2.;
    // making sure that the difference between in- and outflux from the model is less than 1%
    const double flux_discrepancy_tolerance(0.01);
    if ( fabs(fabs(inflow)-fabs(outflow)) / flux_through_model > flux_discrepancy_tolerance ) {
        cerr <<"\n\t"<<"averaged influx versus averaged outflux: "<< inflow <<" vs. "<< outflow <<" m3/s.";
        cout<<"relative flux error = "<<fabs(fabs(inflow)-fabs(outflow)) / flux_through_model <<endl;
        if ( fabs(fabs(inflow)-fabs(outflow)) / flux_through_model > 1.0e-2 ) {
             csmp_error.Note( ERROR, "equivalentPermeabilityQuadraticFEM:", "flux analysis is inaccurate; check integrity of boundary mesh." );
             // using the maximum flux estimate in this case
             flux_through_model = std::max( fabs(inflow), fabs(outflow) );
          }
      }
   
	  cout <<"\nequivalentPermeabilityQuadraticFEM: average flux through inlet and outlet boundaries: "<< flux_through_model;
    cout <<" m3/s, flux per m2: "<< 2* flux_through_model / (inlet_area+outlet_area) <<" m/s."<< endl;
   
    // recovering fluid viscosity
    double  mumin, fluid_viscosity;
    flow_domain.MinMaxOf( "viscosity", mumin, fluid_viscosity );
    // asserting that there is only a single viscosity value
    assert( fabs(fluid_viscosity - mumin) <= numeric_limits<double>::epsilon() );
   
    cout<<"\n\t"<<"permeability 1 (over inlet/outlet area) = "<<(2*flux_through_model * fluid_viscosity) / ((inlet_area+outlet_area) * farfield_pf_gradient)<<" m2"<<endl;
   
   // 3b. fluid throughput through model, averaged over pore space
   // -------------------------------------------------------------------------------------------
   double velo_integral(0.);
   double avg_velocity = volumeAveragedVelocity( model, region, velo_integral );
   
   // directly from Darcy's law where v_avg = Darcy's velocity
   cout<<"averaged velocity over pore space, m/s= "<<avg_velocity<<endl;
   cout <<"tube model: permeability 2 (from volume-averaged velocity) = "<< (avg_velocity * fluid_viscosity) / farfield_pf_gradient <<" m2"<<endl;
   cout <<"not for tubes, permeability 3( from volume-averaged velocity/unit cell volume) = "<< (velo_integral * fluid_viscosity) / (farfield_pf_gradient* length* length *length ) <<" m2"<<endl;

   // TODO: analyse velocity at the element integration points using their weights (rather than at the element barycentre)
   
    // equivalent permeability, old version
	  return (flux_through_model * fluid_viscosity) / (xsect_area * farfield_pf_gradient);

 } // end equivalentPermeabilityQuadraticFEM (using box-model boundaries)




/**
    For the computation of the seepage forces the no-slip walls of the soil or rock skeleton are identified and converted into Boundary objects
    compute seepage forces and fluid induced tractions acting on pore walls / solid skeleton
    
    @attention since the computation of  the pressure force requires consideration of Face normals it is performed on the linear mesh
    
    @attention since velocities are accurate only and stored at the integration points of the quadratic faces, (flow velocity related) tractions are computed on the quadratic mesh
*/
void StaggeredGridStokesSolver_Example::SeepageForces( Model<3U>& model )
 {
    // auto& csmp_error( ErrorHandler::Instance() );
    auto interpolation_order = model.Mesh().OrderOfShapeFunctions();
    if ( interpolation_order.size() > 1 /* p-refinement */ or  (*interpolation_order.begin()) != 1 )
      throw csmp::Exception( ERROR, "StaggeredGridStokesSolver_Example::SeepageForces", "function linear finite element interpolaton (LFEM)");
    
    // 1. create a csmp::Boundary between "PORES" and "GRAINS" for the integration of surface tractions
    // use the appropriate method from csmp::BoundaryInterface
    //   bname  success (y/n)
    pair<string,bool>    boundary_name1 = model.CreateBoundaryBetween( "PORES", "GRAINS" );
    assert( boundary_name1.second );
    const Boundary<3U>&  surf_integration_domain = model.Boundary( boundary_name1.first );

    // 1b. alternative: create a csmp::Boundary around "PARTICLE" for the integration of surface tractions
    // use the appropriate method from csmp::BoundaryInterface
    //   bname  success (y/n)
    bool created_successfully = model.CreateBoundaryAround( "PARTICLE" );
    assert( created_successfully );
    string boundary_name2 = "BOUNDARY_PARTICLE_HULL"; // BOUNDARY_region_HULL
    const Boundary<3U>&  particle_surface = model.Boundary( boundary_name2 );

    // 2. integration of pressure over boundary surfaces yields force vector acting on region hull
    //   (the unit normals of the elements have to be considered in this operation)
    const Index& p_key = model. Database().StorageKey("fluid pressure");
    Point<3U>    pforce_on_object(0.,0.,0.);
    
    // summing the face forces over the object of interest
    for ( const auto& fit : surf_integration_domain.CellVector() ) {
          // interpolating pressure to element barycentre
          // const double barycenter_pressure = fit->PropertyValueAtBaryCenter(p_key);
          // integrating face-normal pressure to obtain force acting on face
          auto pforce_face  = fit->UnitNormal() * fit->Area() * fit->PropertyValueAtBaryCenter(p_key);
          pforce_on_object += pforce_face;
      }
      
    cout <<"\n"<<"StaggeredGridStokesSolver_Example::SeepageForces: "<< pforce_on_object <<" force (N) acting on object '"<< boundary_name1.first <<"'"<< endl;
    
 } //end SeepageForces



/**
    For the computation of the seepage forces the no-slip walls of the soil or rock skeleton are identified and converted into Boundary objects
    compute seepage forces and fluid induced tractions acting on pore walls / solid skeleton
    
    @attention since the computation of  the pressure force requires consideration of Face normals it is performed on the linear mesh
    
    @attention since velocities are accurate only and stored at the integration points of the quadratic faces, (flow velocity related) tractions are computed on the quadratic mesh
*/
void StaggeredGridStokesSolver_Example::SeepageTractions( Model<3U>& model )
 {
    // 1. create a csmp::Boundary between "PORES" and "GRAINS" for the integration of surface tractions
    // use the appropriate method from csmp::BoundaryInterface
    //   bname  success (y/n)
    pair<string,bool>    boundary_name = model.CreateBoundaryBetween( "PORES", "GRAINS" );
    assert( boundary_name.second );
    const Boundary<3U>&  surf_integration_domain = model.Boundary( boundary_name.first );
    
    // 2. integration of pressure over boundary surfaces yields force vector acting on region hull
    //   (the unit normals of the elements have to be considered in this operation)
    for ( const auto& fit : surf_integration_domain.CellVector() ) {
          // interpolating pressure to the element integration points
      }
    
    
    // 3. integration of the shear stress acting on the object due to the passage of a viscous fluid
    
 } //end SeepageFTractions




/**
     Reports phi, k and averages to file
*/
void StaggeredGridStokesSolver_Example::PrintAnalysisResults( const Model<3>& model, double porosity, const vector<double>&  k_equiv )
 {
    assert( k_equiv.size() == 3U );
    ofstream  ofs( string(model.Name()) + "-simulation-results.txt" );
  
    ofs << "model: '"<<  model.Name() <<"', ";
    ofs << model.Mesh().Nodes() <<" nodes, "<< model.Mesh().Elements() <<" elements\n";
    ofs <<"porosity(computed): "<< porosity <<"\nequivalent permeability (kx,ky,kz):";
    double average_k(0.);
    for ( size_t i=0U; i<3U; ++i ) {
         ofs <<" "<< k_equiv[i];
         average_k += k_equiv[i];
      }
    average_k /= static_cast<double>(3U);
    ofs <<"\naverage k: "<< average_k;
 }



/**
    Sets the flags of scalar node variables on specified boundary / subregion thereof
*/
void StaggeredGridStokesSolver_Example::ChangeBoundaryFlags( Model<3U>& model,
                                                             const char* boundary,
                                                             SUBDOMAIN_PART part,
                                                             const char* node_variable,
                                                             VARIABLE_FLAG flag )
 {
    const csmp::Index  key(model.Database().StorageKey(node_variable));
    Boundary<3U>&      domain(model.Boundary(boundary));
    assert( key.place == NODE );
    assert( key.type  == SCALAR );
 
      if ( part == PERIMETER ) {
          for ( auto nit=domain.PerimeterNodesBegin(); nit!=domain.NodesEnd(); ++nit )
            {
               cerr << (*nit)->Idx() <<":"<< parseStatus((*nit)->Status(key)) <<" ";
               (*nit)->Status( key, flag );
               //(*nit)->Store( key, makeScalar(flag,(*nit)->Read(key)) );
            }
       }
     else if ( part == INTERIOR ) {
          for ( auto nit=domain.NodesBegin(); nit!=domain.PerimeterNodesBegin(); ++nit )
            {
               cerr << (*nit)->Idx() <<":"<< parseStatus((*nit)->Status(key)) <<" ";
               (*nit)->Status( key, flag );
               //(*nit)->Store( key, makeScalar(flag,(*nit)->Read(key)) );
            }
       }
      else {
          for ( auto nit=domain.NodesBegin(); nit!=domain.NodesEnd(); ++nit )
            {
               cerr << (*nit)->Idx() <<":"<< parseStatus((*nit)->Status(key)) <<" ";
               (*nit)->Status( key, flag );
               //(*nit)->Store( key, makeScalar(flag,(*nit)->Read(key)) );
            }
       }

 } // end changeBoundaryFlags



/// only for those nodes for which have the correct AtBoundary() flag, the VARIABLE_FLAG changes are applied
void StaggeredGridStokesSolver_Example::ChangeBoxBoundaryFlags( Model<3U>& model,
                                                                BOX_BOUNDARY boundary,
                                                                const char* node_variable,
                                                                VARIABLE_FLAG flag )
 {
    const csmp::Index  key(model.Database().StorageKey(node_variable));
    Region<3U>&        domain(model.Region("Model"));
    assert( key.place == NODE );
    assert( key.type  == SCALAR );
 
    for ( auto nit=domain.NodesBegin(); nit!=domain.NodesEnd(); ++nit )
      if ( (*nit)->AtBoundary() == boundary ) {
           //cerr << (*nit)->Idx() <<":"<< parseStatus((*nit)->Status(key)) <<" ";
           (*nit)->Status( key, flag );
           //(*nit)->Store( key, makeScalar(flag,(*nit)->Read(key)) );
        }

 } // end changeBoxBoundaryFlags




  double  StaggeredGridStokesSolver_Example::LinearModelElementVolume( const Model<3U>& model, const std::string& flow_domain )
  {
    double tot_volume_element(0.);
    const Region<3U>& model_domain(model.Region(flow_domain));
    const csmp::Index p_key(model.Database().StorageKey("fluid pressure"));
    
    for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it ) {
   
        for ( auto it2=(*it)->NodesBegin(); it2!=(*it)->NodesEnd(); ++it2 ){
          // Volume check
          if ((*it)->Volume() < 0. ) {
            (*it2)->Coordinate().Out();
            cout<<"\t"<< (*it2)->Read(p_key) << endl;
          }
          // BC check
          
        }
      tot_volume_element +=(*it)->Volume();


    }
    cout<<"linearModelElementVolume: tot_volume_element= "<< tot_volume_element << " volume of flow_domain = "<<model_domain.Volume()<<endl;

    return ((tot_volume_element / model_domain.Volume()) - 1);
  }



// =====================================================================================================================
//
// NON-MEMBER - AUXILIARY FUNCTIONS
//
// ======================================================================================================================

/**
    Computes flux through the boundary of the model using quadratic finite-element integration.
    
    1. we start with velocities on the nodes of the quadratic elements
    
    2. we integrate these using the quadrature rules of the boundary (lower-dimensional) faces
    
    A single unit normal is used for the flux calculation. 
    It is taken from the Box class.
*/
template<size_t dim>
double boxBoundaryFluxFEM_Quadratic( const Boundary<dim>& boundary, const csmp::Index& velo_key, double& boundary_area )
 {
   // todo to find face area and maximum velocity magnitude at boundary
    vector<double>     nrml;
    Box().UnitNormalTo( parseBoundary(boundary.Name()), dim, nrml );

    VectorVariable<dim>  velo;
   double             boundary_flux(0.);
   // verification
   double maxflux(0.);
  
    for ( auto it=boundary.CellsBegin(); it!=boundary.CellsEnd(); ++it )
      {
         double face_flux(0.), face_area(0.);
         for ( uint32_t i=0U; i<(*it)->IntegrationPoints(); ++i )
           {
              (*it)->PropertyValueAtIntegrationPoint( velo_key, i, velo );
              // project velocity onto boundary normal
              double face_normal_flux(0.);
              for ( uint32_t j=0U; j<dim; ++j )
                face_normal_flux += nrml[j] * velo[j];
               // integrate flux over face
               face_flux += face_normal_flux * (*it)->det_JINV_AtIntegrationPoint(i) * (*it)->WeightAtIntegrationPoint(i);
               face_area += (*it)->det_JINV_AtIntegrationPoint(i) * (*it)->WeightAtIntegrationPoint(i);
              // max flux
              maxflux = (maxflux < fabs(face_normal_flux)) ? fabs(face_normal_flux) : maxflux;
           }
        boundary_flux += face_flux;
        boundary_area += face_area;
//cout<<"(*it)->Area() = "<< (*it)->Area()<<" face_area = "<<face_area<<endl;
      }
     cout<<"boxBoundaryFluxFEM_Quadratic: area of "<<boundary.Name()<<" = "<<boundary_area<< " max flux = "<<maxflux<<endl;

    return boundary_flux;

 } // end boxBoundaryFluxFEM_Quadratic

template double boxBoundaryFluxFEM_Quadratic<3>( const Boundary<3U>&, const csmp::Index&, double& );












/**
    Computation using velocity magnitude in the averaging.
*/
double  volumeAveragedVelocity( const Model<3>& model, const std::string& flow_domain, double& velo_integral )
 {
   double volume_integral(0.);
//double tot_velocity_element(0.), tot_volume_element(0.);
     const Region<3U>& model_domain(model.Region(flow_domain));
     const csmp::Index v_key(model.Database().StorageKey("nodal velocity"));
     VectorVariable<3> velocity;

//   double patm(100325.);
   
   const csmp::Index p_key(model.Database().StorageKey("fluid pressure"));

     for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it ) {
         for ( uint32_t i=0u; i<(*it)->IntegrationPoints(); ++i )
           {
   //          if ((*it)->PropertyValueAtIntegrationPoint( p_key, i) >= patm)
   //          {
              (*it)->PropertyValueAtIntegrationPoint( v_key, i, velocity );
              double velo_magnitude = velocity.Length();
              double det_J((*it)->det_JINV_AtIntegrationPoint(i));
              velo_integral   += velo_magnitude * det_J * (*it)->WeightAtIntegrationPoint(i);
              // test: error
              volume_integral += det_J * (*it)->WeightAtIntegrationPoint(i);
  //            }
           }
       // diff way to calc
//       tot_volume_element +=(*it)->Volume();
 //      (*it)->Read(v_key,velocity);
 //      tot_velocity_element += velocity.Length() * (*it)->Volume();
       }
   cout<<"volumeAveragedVelocity: velo_integral= "<< velo_integral << " volume = "<<volume_integral<<endl;
// alternative way to calc
//cout
//   <<"tot_velocity_element= "<< tot_velocity_element
//   << " tot_volume_element = "<<tot_volume_element
//   <<" average velocity over whole pore volume= "<< tot_velocity_element/tot_volume_element
//   <<endl;
   
/*   for ( auto nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
   {
     (*nit)->Read(v_key,velocity);
     (*nit)->Read(p_key,pressure);
     
   }
 */
     return velo_integral / volume_integral;
 }





// *****************************************************************************************************************

// CODE FOR TESTING

// *****************************************************************************************************************


// FOR TESTING
void StaggeredGridStokesSolver_Example::ComputePressureLEFT_RIGHT( Model<3>& model )
 {
    Boundary<3>&  inlet(model.Boundary("LEFT"));
    Boundary<3>&  outlet(model.Boundary("RIGHT"));
    const double patm(100325.);
    inlet.InputPropertyValue( "fluid pressure", makeScalar(DIRICH,2. * patm) );
    outlet.InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
   
    // 1. Set-up to compute fluid pressure
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_ncgtyp(5);    // strong connectivity
    settings.Set_nxtyp(0);     // Gauss-Seidl relaxation
    settings.Set_ndefault(40); // standard coarsening strategy
    SAMG_Solver   solver(&settings);
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
    PDE_Integrator<3,Element>       pressure_solver(solver);
    NumIntegral_dNT_dN_dV<3>        laplacian( model.Database(), "fluid pressure", "fluid pressure" );
    NumIntegral_SetRHS_to_Zero<3>   zero( model.Database(), "fluid pressure" );
    pressure_solver.Add( &laplacian );
    pressure_solver.Add( &zero );
   
    Region<3>& pore_domain(model.Region("PORES"));
    pressure_solver.IntegrateOver( pore_domain );
    printRangeOfVariable( model, "PORES", "fluid pressure" );
 }




/// TESTING helper for function below
void nodeCoordinatesToText( const char* file, const set<Node<3U>*>& points )
 {
    ofstream  ofs(file);
    //ofs <<"Point data representing node locations (x,y,z); "<< file << endl;
    ofs <<"x,y,z,node-id\n";
    for ( set<Node<3U>*>::const_iterator pt=points.begin(); pt!=points.end(); ++pt )
      ofs << (*pt)->x() <<","<< (*pt)->y() <<","<< (*pt)->z() <<","<< (*pt)->Idx() << endl;
   
    ofs.close();
 }



/// @test OK - no issues found
void flagReadWriteTest( const Model<3U>& model, const char* variable )
 {
    const Region<3U>& domain(model.Region("Model"));
    const csmp::Index  key(model.Database().StorageKey(variable));
    assert( key.place == NODE );
 
    cerr <<"\n\nflagReadWriteTest: node variable status: 1) original, 2) var written 3) flag written and retrieved:";
    for ( auto nit=domain.PerimeterNodesBegin(); nit!=domain.NodesEnd(); ++nit )
      {
         VARIABLE_FLAG flag1 = (*nit)->Status( key );
         (*nit)->Store( key, makeScalar(flag1,0.) );
         VARIABLE_FLAG flag2 = (*nit)->Status( key );
         (*nit)->Status( key, flag2 );
         VARIABLE_FLAG flag3 = (*nit)->Status( key );
         cerr <<"\n\t"<< (*nit)->Idx() <<": "<< parseStatus(flag1) <<" "<< parseStatus(flag2) <<" "<< parseStatus(flag3) <<" ";
      }
      
} // end flagReadWriteTest





void nodesToFile( const Model<3U>& model, const char* region, SUBDOMAIN_PART part )
 {
     const Region<3U>& domain(model.Region(region));
     set<Node<3U>*>    points;
   
     if ( part == PERIMETER ) {
          for ( auto nit=domain.PerimeterNodesBegin(); nit!=domain.NodesEnd(); ++nit )
            points.insert( (*nit) );
       }
     else if ( part == INTERIOR ) {
          for ( auto nit=domain.NodesBegin(); nit!=domain.PerimeterNodesBegin(); ++nit )
            points.insert( (*nit) );
       }
      else {
          for ( auto nit=domain.NodesBegin(); nit!=domain.NodesEnd(); ++nit )
            points.insert( (*nit) );
       }
   
      string file_name(model.Name());
      file_name += "_" + parseSubdomainPart(part);
      file_name += ".txt";
      nodeCoordinatesToText( file_name.c_str(), points );
   
 } // end nodesToFile







/*  CODE FRAGMENTS FOR CHECKING

    // Extract pressure values at the nodes from the first mesh
    // --------------------------------------------------------
    if ( linear_model.RestoreOriginalNodeNumbering(false) ) {
         // if changes were made, restored numbering and point coordinates are output to file
         ofstream ofs("linear_model-original-node-numbering.txt");
         ofs <<"linear model: restored original node numbers and coordinates:\n";
         for ( size_t i=0U; i<pore_domain.Nodes(); ++i ) {
              Point<dim> xyz(pore_domain.N(i)->Coordinate());
              ofs << pore_domain.N(i)->Idx() <<" "<< xyz[0] <<" "<< xyz[1] <<" "<< xyz[2] <<"\n";
           }
      }



    // check that node coordinates of the linear and quadratic models are the same
    // ---------------------------------------------------------------------------
    cout <<"\n\nmain: checking node numbering of the 2 models for consistency...\n";
    size_t non_matching_nodes(0U);
    auto it2=pore_domain2.NodesBegin();
    for ( auto it=pore_domain.NodesBegin(); it!=pore_domain.NodesEnd(); ++it, ++it2 )
      if ( (*it)->Coordinate() != (*it2)->Coordinate() ) {
           Point<dim> xyz1((*it)->Coordinate()), xyz2((*it2)->Coordinate());
           cerr << (*it)->Idx() <<" "<< xyz1[0] <<" "<< xyz1[1] <<" "<< xyz1[2] <<" ";
           cerr << (*it2)->Idx() <<" "<< xyz2[0] <<" "<< xyz2[1] <<" "<< xyz2[2] <<"\n";
           non_matching_nodes++;
        }
    if ( non_matching_nodes > 0 ) {
         cerr <<"\nmain: "<< non_matching_nodes <<" nodes do not have matching coordinates.\n";
         ofstream ofs("non-matching-nodes.txt");
         ofs <<"node numbers and coordinates from linear model followed by those from quadratic model\n";
         auto it2=pore_domain2.NodesBegin();
         for ( auto it=pore_domain.NodesBegin(); it!=pore_domain.NodesEnd(); ++it, ++it2 )
            if ( (*it)->Coordinate() != (*it2)->Coordinate() ) {
                 Point<dim> xyz1((*it)->Coordinate()), xyz2((*it2)->Coordinate());
                 ofs << (*it)->Idx() <<" "<< xyz1[0] <<" "<< xyz1[1] <<" "<< xyz1[2] <<" ";
                 ofs << (*it2)->Idx() <<" "<< xyz2[0] <<" "<< xyz2[1] <<" "<< xyz2[2] <<"\n";
              }
      }
      
*/


/**
     TODO: allow command-line supply of input file set and runtime options
     TODO: remove functions that are not used
     
     TODO: write pore-radius onto pore walls using quadratic model and obtain statistics on it
     
     TODO: Add support for periodic boundary conditions
     
     TODO: make the model fault tolerant with regard to single grain-boundary points touching the inlet or outlet
     (Q: is it already fault tolerant?)
     
     TESTING RESULTS
     ===============
     Alternative choice of integration rule for quadratic tetrahedron (by default the Gauss-Points are at midside nodes !
       (none of the schemes gave a significant improvement)
    
    FAIL dynamic_cast<IsoparametricQuadraticTetrahedron*>(quadratic_model.FE_Manager().E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON))->QuadratureRules_11Points();

    NO DIFFERENCE TO DEFAULT dynamic_cast<IsoparametricQuadraticTetrahedron*>(quadratic_model.FE_Manager().E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON))->QuadratureRules_24Points();

    NO_DIFFERENCE dynamic_cast<IsoparametricQuadraticTetrahedron*>(quadratic_model.FE_Manager().E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON))->QuadratureRules_5Points();

    NO DIFFERENCE dynamic_cast<IsoparametricQuadraticTetrahedron*>(quadratic_model.FE_Manager().E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON))->QuadratureRules_4Points();

    ZERO IN DIAGONAL dynamic_cast<IsoparametricQuadraticTetrahedron*>(quadratic_model.FE_Manager().E(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON))->QuadratureRules_1Point();

    Quadratic approximation of fluid pressure (works fine with specific settings for SAMG solver
    computePressureLEFT_RIGHT( quadratic_model );
    vtu_quadratic.OutputDataToVTU( string("quadratic_model-2nd-order-fluid-pressure"), string("fluid pressure"), string(pore_space), 0 );

*/


} // end csmp



