//
//  staggered_grid_StokesSolver_main.cpp
//  CSMP_GitHub/applications/pore_scale_modelling/2-stage-Stokes
//
//  Created by Stephan Matthai on 6/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include <iostream>
#include "ANSYS_Model3D.h"
#include "VTU_Interface.h"
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "PDE_Integrator.h"
#include "PDE_IntegratorExperimental.h"
#include "IsoparametricQuadraticTetrahedron.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_DNi_rhsop_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "CSMP_highLevelUtilities.h"
#include "FEM_Data.h"
#include "Box.h"
#include "StatisticalAnalyzer.h"

namespace csmp {

    void displayTitle();

    void assignNoSlipBoundaryConditions( Model<3U>& model,
                                         Region<3U>& stokes_flow_domain,
                                         Region<3U>& grain_edges,
                                         ushort flow_direction );

    void constructVelocityVector( Model<3>& );

    double64 volumeAveragedVelocity( const Model<3U>&, const std::string& flow_domain, double64& velovity_integral );

    double64 equivalentPermeabilityQuadraticFEM( Model<3U>&, const std::string& flow_domain,
                                                 BOX_BOUNDARY inflow_boundary, BOX_BOUNDARY outflow_boundary,
                                                 double64& flux_through_model );
      
    void printAnalysisResults1( const Model<3U>&, double64 porosity, const std::vector<double64>&  k_equiv );

    // TESTING
    /// quadratic fluid pressure computation
    void computePressureLEFT_RIGHT( Model<3>& );
    void nodesToFile( const Model<3U>&, const char* region, SUBDOMAIN_PART );
    void flagReadWriteTest( const Model<3U>&, const char* variable );
    void changeBoundaryFlags( Model<3U>&, const char* boundary, SUBDOMAIN_PART, const char* node_variable, VARIABLE_FLAG );
    void changeBoxBoundaryFlags( Model<3U>&, BOX_BOUNDARY, const char* node_variable, VARIABLE_FLAG );
}


using namespace std;
using namespace csmp;


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
int main()
 {
    displayTitle();
   
    const int dim(3U);
   
    const bool   verbose(true);

// TESTED MODELS
   
   // DONE
//    const string input_model_linear("tubes_all_quadra_linearized");
//    const string input_model_quadratic("tubes_all_quadra");
   
//   const string input_model_linear("spiral_quadra_linearized");
//   const string input_model_quadratic("spiral_quadra");

   //// new bcc
   
//   const string input_model_linear("bcc_0366_5x5x5_raw_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_5x5x5_raw_quadra");
//   const string input_model_linear("bcc_0366_10x10x10_coarse_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_10x10x10_coarse_quadra");
//   const string input_model_linear("bcc_0366_10x10x10_raw_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_10x10x10_raw_quadra");
//   const string input_model_linear("bcc_0366_25x25x25_coarse_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_25x25x25_coarse_quadra");
//   const string input_model_linear("bcc_0366_25x25x25_raw_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_25x25x25_raw_quadra");
//   const string input_model_linear("bcc_0366_50x50x50_coarse_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_50x50x50_coarse_quadra");
//   const string input_model_linear("bcc_0366_50x50x50_raw_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_50x50x50_raw_quadra");
   
//   const string input_model_linear("bcc_0366_quadratic_fine_linarized");
//   const string input_model_quadratic("bcc_0366_quadratic_fine");
//   const string input_model_linear("bcc_0366_basecase_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_basecase_quadra");
//   const string input_model_linear("bcc_0366_coarse_quadra_linearized_new"); == basecase
//   const string input_model_quadratic("bcc_0366_coarse_quadra_new"); == basecase
//    const string input_model_linear("bcc_0366_coarse2_quadra_linearized");
//    const string input_model_quadratic("bcc_0366_coarse2_quadra");
//   const string input_model_linear("bcc_0366_2elements_channel_quadra_linearized");
//   const string input_model_quadratic("bcc_0366_2elements_channel_quadra");
   
  // 5 June 2017, Apoorv's models
//   const string input_model_linear("500Voxel_Berea_Quadratic_linear");
//   const string input_model_quadratic("500Voxel_Berea_Quadratic");
//   const string input_model_linear("MtGambier_500Voxel_Quadratic_Linear");
//   const string input_model_quadratic("MtGambier_500Voxel_Quadratic");
   
   //// rock samples
   // Berea 100 Caroline
//   const string input_model_linear("berea_100x100x100_MP_in_porespace_quadra_linearized");
//   const string input_model_quadratic("berea_100x100x100_MP_in_porespace_quadra");
   
   // Berea subvolume 100 Chloe
   const string input_model_linear("berea_sub100_quadra_linearized");
   const string input_model_quadratic("berea_sub100_quadra");
   // size of seed elements 20
//   const string input_model_linear("berea_sub100_20_quadra_linearized");
//   const string input_model_quadratic("berea_sub100_20_quadra");
   
   // ANLEC 100
//   const string input_model_linear("anlec100_quadra_linearized");
//   const string input_model_quadratic("anlec100_quadra");


	// output variables of initial configuration to VTU
	list<string> outputProps_before;
	outputProps_before.push_back("fluid pressure");
	outputProps_before.push_back("nodal velocity x");
	outputProps_before.push_back("nodal velocity y");
	outputProps_before.push_back("nodal velocity z");
  // testing only
  if ( verbose ) {
       outputProps_before.push_back("node-flag velocity x");
       outputProps_before.push_back("node-flag velocity y");
       outputProps_before.push_back("node-flag velocity z");
    }
	
	// output result variables to VTU
	list<string> outputProps_after;
	outputProps_after.push_back("fluid pressure");
	outputProps_after.push_back("nodal velocity");
	outputProps_after.push_back("nodal velocity x");
	outputProps_after.push_back("nodal velocity y");
	outputProps_after.push_back("nodal velocity z");


    // ===========================================================================================================================
    // A. First linear FEM model (LFEM) to solve the pressure equation
    // ===========================================================================================================================
    const bool irregular_mesh(true),
               binary_file(true),
               use_regions_file(true),
               create_boundaries(true);
   
    // temporary since not needed for quadratic FEM computation
    ANSYS_Model3D* linear_model = new ANSYS_Model3D( input_model_linear.c_str(),
                                                    "staggered_grid_stokes_variables_p.txt",
                                                     irregular_mesh, binary_file, use_regions_file, create_boundaries );

    const size_t LFEM_nodes    = linear_model->Mesh().Nodes();
    const size_t LFEM_elements = linear_model->Mesh().Elements();
 
	  const double64 edge_length = printModelDimensions(*linear_model, true);

    VTU_Interface<dim>  vtu_linear(*linear_model);

    const string pore_space("PORES");
    const string pore_edge("PORE_EDGE");
    const double64 viscosity(1.0e-3);
    const double64 patm(100325.);
//    const double64 hydrostatic_grad(1.e+3 * 9.86); // to get consistent histograms ... 9.80665?
       const double64 hydrostatic_grad(1.e+6 * 9.86); // to get consistent histograms
   
    // SAMG solver and integrator for linear algebraic system
    SAMG_Settings settings;
    SAMG_Solver   solver(&settings);
    // FEM assembly for pressure equation: Laplacian = 0
    PDE_Integrator<dim,Region>       pressure_solver(solver);
    NumIntegral_dNT_dN_dV<dim>       laplacian( linear_model->Database(), "fluid pressure", "fluid pressure" );
    NumIntegral_SetRHS_to_Zero<dim>  zero( linear_model->Database(), "fluid pressure" );
    pressure_solver.Add( &laplacian );
    pressure_solver.Add( &zero );
   
    // flow domain (subregion of the model that must be called PORES)
    Region<dim>&              pore_domain(linear_model->Region(pore_space));
//    FEM_Data<ScalarVariable>  pdata[dim];
    /// pressure map of linear model
    map<Point<3U>,double64>  pressureX_linmodel, pressureY_linmodel, pressureZ_linmodel;
//   Index nodeKey_linmod( linear_model->Database().StorageKey( "fluid pressure") );
   
//linear_model->RestoreOriginalNodeNumbering(verbose);

    // 1. LEFT-RIGHT fluid-pressure computation: pressure gradient in X direction
    // --------------------------------------------------------------------------
    // boundary conditions
    linear_model->Boundary("LEFT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * edge_length) );
    linear_model->Boundary("RIGHT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
   vtu_linear.OutputDataToVTU( string("fluid-pressure-x_pboundaries"), string("fluid pressure"), string("Model"), 0 );
    // computation
    pressure_solver.IntegrateOver( pore_domain );
    pore_domain.ChangePropertyStatus( "fluid pressure", ANY, PERIMETER );
    // output
    printRangeOfVariable( *linear_model, pore_space.c_str(), "fluid pressure" );
    vtu_linear.OutputDataToVTU( string("fluid-pressure-x"), string("fluid pressure"), string("Model"), 0 );
    // output to FEM_Data
//    linear_model->RestoreOriginalNodeNumbering(verbose);
//    linear_model->OutputVariableTo("fluid pressure", pdata[0] );
//    if ( verbose ) {
//         ScalarVariable  tmin, tmax;
//         pdata[0].MinMaxOf( tmin, tmax );
//         cout <<"\nmain: 'fluid pressure' range in FEM_Data container: ";
//         cout << tmin <<" - "<< tmax << endl << endl;
//      }
    ///// Pressure mapping of linear model
    const auto nodesEnd_linmodel(linear_model->Mesh().NodesEnd());
   cout<<"X: pressure under patm"<<endl;
     const csmp::Index p_key (linear_model->Database().StorageKey("fluid pressure"));
    for ( auto nit=linear_model->Mesh().NodesBegin(); nit!=nodesEnd_linmodel; ++nit  ) {
      // TODO use property_name, pressure, for ex. pref_.StorageKey(property_name)
      pressureX_linmodel.insert( make_pair( nit->Coordinate(), nit->Read(p_key)) );
      if (nit->Read(p_key) < patm ) cout<<nit->Coordinate()<<" "<<nit->Read(p_key)<<endl;
    }
    assert( linear_model->Mesh().Nodes() == pressureX_linmodel.size() );
    cerr  <<"\nnumber of nodes in pressure map of linear model, X direction: "<< pressureX_linmodel.size();
 
    // 2. BOTTOM-TOP fluid-pressure computation: pressure gradient in Y direction
    // --------------------------------------------------------------------------
    // boundary conditions
    linear_model->Boundary("BOTTOM").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * edge_length) );
    linear_model->Boundary("TOP").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
   vtu_linear.OutputDataToVTU( string("fluid-pressure-y_pboundaries"), string("fluid pressure"), string("Model"), 0 );
    // computation
    pressure_solver.IntegrateOver( pore_domain );
    pore_domain.ChangePropertyStatus( "fluid pressure", ANY, PERIMETER );
    // output
    printRangeOfVariable( *linear_model, pore_space.c_str(), "fluid pressure" );
    vtu_linear.OutputDataToVTU( string("fluid-pressure-y"), string("fluid pressure"), string("Model"), 0 );
    // output to FEM_Data
//    linear_model->RestoreOriginalNodeNumbering(verbose);
//    linear_model->OutputVariableTo("fluid pressure", pdata[1] );
    ///// Pressure mapping of linear model
   cout<<"Y: pressure under patm"<<endl;
    for ( auto nit=linear_model->Mesh().NodesBegin(); nit!=nodesEnd_linmodel; ++nit  ) {
      // TODO use property_name, pressure, for ex. pref_.StorageKey(property_name)
      pressureY_linmodel.insert( make_pair( nit->Coordinate(), nit->Read(p_key)) );
      if (nit->Read(p_key) < patm ) cout<<nit->Coordinate()<<" "<<nit->Read(p_key)<<endl;
    }
    assert( linear_model->Mesh().Nodes() == pressureY_linmodel.size() );
    cerr  <<"\nnumber of nodes in pressure map of linear model, Y direction: "<< pressureY_linmodel.size();
   
   
    // 3. BACK-FRONT fluid-pressure computation: pressure gradient in Z direction
    // --------------------------------------------------------------------------
    // boundary conditions
    linear_model->Boundary("BACK").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm + hydrostatic_grad * edge_length) );
    linear_model->Boundary("FRONT").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
   vtu_linear.OutputDataToVTU( string("fluid-pressure-z_pboundaries"), string("fluid pressure"), string("Model"), 0 );
    // computation
    pressure_solver.IntegrateOver( pore_domain );
    pore_domain.ChangePropertyStatus( "fluid pressure", ANY, PERIMETER );
    // output
    printRangeOfVariable( *linear_model, pore_space.c_str(), "fluid pressure" );
    vtu_linear.OutputDataToVTU( string("fluid-pressure-z"), string("fluid pressure"), string("Model"), 0 );
    // output to FEM_Data
//    linear_model->RestoreOriginalNodeNumbering(verbose);
//    linear_model->OutputVariableTo("fluid pressure", pdata[2] );
    ///// Pressure mapping of linear model
   cout<<"Z: pressure under patm"<<endl;
    for ( auto nit=linear_model->Mesh().NodesBegin(); nit!=nodesEnd_linmodel; ++nit  ) {
      // TODO use property_name, pressure, for ex. pref_.StorageKey(property_name)
      pressureZ_linmodel.insert( make_pair( nit->Coordinate(), nit->Read(p_key)) );
      if (nit->Read(p_key) < patm ) cout<<nit->Coordinate()<<" "<<nit->Read(p_key)<<endl;
    }
    assert( linear_model->Mesh().Nodes() == pressureZ_linmodel.size() );
   cerr  <<"\nnumber of nodes in pressure map of linear model, Z direction: "<< pressureZ_linmodel.size()<<endl;
   
   cout<<"Linear model, Number of interior nodes= "<<pore_domain.InteriorNodes()<<endl;
    // output of the model to CSMP native format (for future runs)
    string binary_output_file_set_name(linear_model->Name());
    binary_output_file_set_name += "_csmp_LFEM";
    linear_model->OutputToBinaryFile( binary_output_file_set_name.c_str() );
   
    // end of LFEM computations
    delete linear_model;
   
   
   
   
   
   
    // ===========================================================================================================================
    // B. Second model (quadratic FEM) to solve the Stokes-flow equation
    // ===========================================================================================================================
    // output from the Ansys mesh generator via the CSMP interface
    ANSYS_Model3D  quadratic_model( input_model_quadratic.c_str(),
                                   "staggered_grid_stokes_variables_u.txt",
                                    irregular_mesh, binary_file, use_regions_file, create_boundaries );

    recreateBoxBoundaryFlags( quadratic_model );
   
    Region<dim>&  pore_domain2(quadratic_model.Region(pore_space));
    Region<dim>&  grain_edges(quadratic_model.Region(pore_edge));
    assert( quadratic_model.Mesh().Nodes() == LFEM_nodes );
    cout <<"\nmain: LFEM / QFEM element ratio: "<< LFEM_elements / quadratic_model.Mesh().Elements() << endl;


    // 1. QFEM algorithm: mu laplacian u = grad p
    // ------------------------------------------
    PDE_IntegratorExperimental<dim,Region>  velocity_solver(solver);
   
/*
    // settings usually used to compute solution for quadratic FEM
    settings.Set_ncgtyp(5);
    settings.Set_nxtyp(0);
    settings.Set_ndefault(40);
*/

/*  
    // settings for variable-by-variable smoothing
    // DOES NOT CONVERGENCE
    settings.Set_napproach(3); // interpolation seperate for each unknown
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

    velocity_solver.Add( &viscosity_matr_x );
    velocity_solver.Add( &viscosity_matr_y );
    velocity_solver.Add( &viscosity_matr_z );
    velocity_solver.Add( &gradient_x );
    velocity_solver.Add( &gradient_y );
    velocity_solver.Add( &gradient_z );


    // to get velocity histograms for each flow direction
    StatisticalAnalyzer<3U>  modelAnalyzer( quadratic_model );
    HistogramBins            velocityRanges;
    // specify via ascii input file, the x-axis range of the histogram columns, i.e. the bin size
    modelAnalyzer.DefineBins( "pore_velocity_bins.bins", velocityRanges );
    // histogram output
    map<string,pair<HistogramBins,size_t> > velo_histograms;
    // paraview output
    VTU_Interface<dim> vtu_quadratic(quadratic_model);

    // fluid viscosity
    quadratic_model.InputPropertyValue( "viscosity", makeScalar(ANY,viscosity) );
    // model porosity
    const double64  porosity = pore_domain2.Volume() / pow(edge_length,3);
    cout <<"\n\n\nmodel: '"<<  quadratic_model.Name() <<"', porosity(computed): "<< std::fixed << std::setprecision(2) << porosity <<"\n";
    cout << std::scientific << std::setprecision(5);
    // permeability as measured in the 3 coordinate directions
    vector<double64> k_equiv(3,0.);


    // 2. Computation of flow velocities from the XYZ axes-aligned pressure gradients
    // ------------------------------------------------------------------------------
    BOX_BOUNDARY inlet(LEFT), outlet(RIGHT);
const auto nodesEnd_quadrmodel(quadratic_model.Mesh().NodesEnd());
const Index nodeKey_quadrmod( quadratic_model.Database().StorageKey( "fluid pressure") );
   
    for ( ushort flow_direction=X_DIRECTION; flow_direction <= Z_DIRECTION; flow_direction++ )
      {
         // assigning pre-calculated pressure values from linear to quadratic model
//         quadratic_model.RestoreOriginalNodeNumbering(verbose);
//         quadratic_model.InputVariableFrom("fluid pressure", pdata[flow_direction] );
        /// From pressure map to quadratic model
        if (flow_direction==X_DIRECTION) {
          const auto nodesEnd_pressureX_linmodel(pressureX_linmodel.end());
          assert( quadratic_model.Mesh().Nodes() == pressureX_linmodel.size() );
          for ( auto nit=quadratic_model.Mesh().NodesBegin(); nit!=nodesEnd_quadrmodel; ++nit ) {
            auto lin_it( pressureX_linmodel.find( (*nit).Coordinate() ) );
            if ( lin_it != nodesEnd_pressureX_linmodel ) {
                (*nit).Store(nodeKey_quadrmod, makeScalar((*nit).Status(nodeKey_quadrmod) , (*lin_it).second ) );
              // cout<<"(*nit).getData().data["<<nodeKey_quadrmod.dataOffset<<"] = "<<(*nit).getData().data[nodeKey_quadrmod.dataOffset]
              // <<" (*lin_it).second = "<<(*lin_it).second<<endl;
            }
            else {
              (*nit).Coordinate().Out();
              throw csmp::Exception( ERROR, "quadratic_model:",
                                  "node could not be identified; has it been newly created?" );
            }
          }
        } else if (flow_direction==Y_DIRECTION) {
          const auto nodesEnd_pressureY_linmodel(pressureY_linmodel.end());
          assert( quadratic_model.Mesh().Nodes() == pressureY_linmodel.size() );
          for ( auto nit=quadratic_model.Mesh().NodesBegin(); nit!=nodesEnd_quadrmodel; ++nit ) {
            auto lin_it( pressureY_linmodel.find( (*nit).Coordinate() ) );
            if ( lin_it != nodesEnd_pressureY_linmodel ) {
                (*nit).Store( nodeKey_quadrmod, makeScalar((*nit).Status(nodeKey_quadrmod) , (*lin_it).second ) );
            }
            else {
              (*nit).Coordinate().Out();
              throw csmp::Exception( ERROR, "quadratic_model:",
                                    "node could not be identified; has it been newly created?" );
            }
          }
        } else { //(flow_direction==Z_DIRECTION)
          const auto nodesEnd_pressureZ_linmodel(pressureZ_linmodel.end());
          assert( quadratic_model.Mesh().Nodes() == pressureZ_linmodel.size() );
          for ( auto nit=quadratic_model.Mesh().NodesBegin(); nit!=nodesEnd_quadrmodel; ++nit ) {
            auto lin_it( pressureZ_linmodel.find( (*nit).Coordinate() ) );
            if ( lin_it != nodesEnd_pressureZ_linmodel ) {
                (*nit).Store(nodeKey_quadrmod, makeScalar((*nit).Status(nodeKey_quadrmod) , (*lin_it).second ));
            }
            else {
              (*nit).Coordinate().Out();
              throw csmp::Exception( ERROR, "quadratic_model:",
                                    "node could not be identified; has it been newly created?" );
            }
          }
        }

        cout <<"\n\nmain: analysing permeability / flow velocity distribution in "<< parse(static_cast<SPATIAL_DERIVATIVE>(flow_direction));
         cout <<" direction for the given boundary pressure differential...\n";
         printRangeOfVariable( quadratic_model, pore_space.c_str(), "fluid pressure" );
         if ( verbose ) vtu_quadratic.OutputDataToVTU( string("quadratic_model-interpolated-fluid-pressure"),
                                                       string("fluid pressure"),
                                                       string(pore_space), flow_direction+1 );

         // assigning property values and boundary conditions
         assignNoSlipBoundaryConditions( quadratic_model, pore_domain2, grain_edges, flow_direction );

         // testing of boundary conditions
         if ( verbose ) {
              flagToNumber( quadratic_model, "nodal velocity x", "node-flag velocity x" );
              flagToNumber( quadratic_model, "nodal velocity y", "node-flag velocity y" );
              flagToNumber( quadratic_model, "nodal velocity z", "node-flag velocity z" );
              vtu_quadratic.OutputDataToVTU( "input", outputProps_before, pore_space.c_str(), static_cast<long>(flow_direction+1) );
           }
         velocity_solver.IntegrateOver( pore_domain2 );
         constructVelocityVector( quadratic_model );
         printRangeOfVariable( quadratic_model, pore_space.c_str(), "nodal velocity" );
 
         // analysis
         if      ( flow_direction == Y_DIRECTION ) {  inlet = BOTTOM, outlet = TOP; }
         else if ( flow_direction == Z_DIRECTION ) { inlet = BACK, outlet = FRONT; }
         double64 flux_through_model;
         k_equiv[flow_direction] = equivalentPermeabilityQuadraticFEM( quadratic_model, pore_space, inlet, outlet, flux_through_model );
         cout <<"\n\tmain: k_equivalent old ("<< parse(static_cast<SPATIAL_DERIVATIVE>(flow_direction)) <<"): "<< k_equiv[flow_direction] <<" m2.\n";

        //modelAnalyzer.RegionPropertyHistogramsElement( "nodal velocity", velocityRanges, velo_histograms );
        // pore domain for tubes but cell volume for BCC models
        //modelAnalyzer.RegionPropertyHistogramsIntegrationPoint("nodal velocity", velocityRanges, velo_histograms, pore_space, pore_domain2.Volume() );
    // TODO add from the version on Mac        modelAnalyzer.RegionPropertyHistogramsIntegrationPoint("nodal velocity", velocityRanges, velo_histograms, pore_space);
         // writes Maple statlist format with weights; results must be copied and pasted into Maple workbook
         //const bool log10_of_bin_values(true);
        const bool log10_of_bin_values(false);
         modelAnalyzer.OutputRegionPropertyHistogramsMaple( (string("nodal velocity") + to_string(flow_direction+1)).c_str(),
                                                            velocityRanges, velo_histograms, log10_of_bin_values );
   
         vtu_quadratic.OutputDataToVTU( string("velocity"), outputProps_after, string("Model"), static_cast<long>(flow_direction+1) );
        
        // Tubes: extracting a region of interest in VTU format for statistical analysis
/*        if (flow_direction==X_DIRECTION) {
          quadratic_model.FormRegionFrom( "tube_x", "fluid pressure", patm, patm + hydrostatic_grad * edge_length, true );
          vtu_quadratic.OutputDataToVTU( string("velocity"), outputProps_after, string("tube_x"), 0 );
        } else if (flow_direction==Y_DIRECTION) {
          quadratic_model.FormRegionFrom( "tube_y", "fluid pressure", patm, patm + hydrostatic_grad * edge_length, true );
          vtu_quadratic.OutputDataToVTU( string("velocity"), outputProps_after, string("tube_y"), 0 );
        } else {
          quadratic_model.FormRegionFrom( "tube_z", "fluid pressure", patm, patm + hydrostatic_grad * edge_length, true );
          vtu_quadratic.OutputDataToVTU( string("velocity"), outputProps_after, string("tube_z"), 0 );
        }
*/
//        if (flow_direction==Z_DIRECTION) {
//          quadratic_model.PartitionRegionIntoContiguousSubRegions(pore_space.c_str());
//        }
    
      } // end flow_direction loop
   cout<<"Quadratic model, Number of interior nodes= "<<pore_domain2.InteriorNodes()<<endl;
    printAnalysisResults1( quadratic_model, porosity, k_equiv );

    // output of quadratic model to CSMP native format (for future runs)
    string binary_output_file_set_name2(quadratic_model.Name());
    binary_output_file_set_name2 += "_csmp_QFEM";
    quadratic_model.OutputToBinaryFile( binary_output_file_set_name2.c_str() );
  
   
   
    return 0;
   
 } // end main






namespace csmp {


void displayTitle() {
      cout <<"CSMP-Stokes 2-Stage\n";
      cout <<"===================\n";
      cout <<"Stephan Matthai, 15 Dec. 2016.\n";
      cout <<"Sequential pressure (piecewise linear) - velocity (piecewise quadratic) FEM solver for incompressible Stokes lubrication equation.\n";
      cout <<"Note: uses node-matched staggered tetrahedral element meshes as input; each quadratic tetrahedron consists of 10 linear tetrahedra.";
      cout <<"\nUsage Instructions\n\n";
      cout <<"1. input:  mesh-files from ANSYS (*.asc and *.dat), created using CSP-output interface; ";
      cout <<" 2 sets of files containing linear and quadratic meshes, respectively.\n";
      cout <<"2. Input models: must be box-shaped, consisting of contiguous mesh domains. Side boundaries must be labeled LEFT, RIGHT etc., see CSMP User's guide.\n";
      cout <<"\nmodels must contain a region that denotes the pore space called 'PORES' ";
      cout <<" and the intersection curves of the grain boundaries with the box boundaries.";
      cout <<"This line-element region must be called 'PORE_EDGE'.\n\n\n";
      cout.flush();
  }



/**
    composes the 3D vector variable velocity from its components all across the 
    entire domain "Model"
    
    @attention velocity components x, y, z must have been initialised everywhere
    before this function is called.
*/
void constructVelocityVector( Model<3>& mdl )
{
	VectorVariable<3> v;
	const csmp::Index vx_key(mdl.Database().StorageKey("nodal velocity x")),
		                vy_key(mdl.Database().StorageKey("nodal velocity y")),
		                vz_key(mdl.Database().StorageKey("nodal velocity z")),
		                v_key(mdl.Database().StorageKey("nodal velocity"));

	Region<3>& mref = mdl.Region("Model");

	for ( vector<Node<3>* >::const_iterator
        nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++)
    {
       v(0) = (*nit)->Read(vx_key);
       v(1) = (*nit)->Read(vy_key);
       v(2) = (*nit)->Read(vz_key);
       (*nit)->Store(v_key, v);
    }
}



/**
     Assigns 
     
     1. No slip conditions to grain edges
 
     2. Box-boundary parallel slip conditions to model side boundaries and edges
 
     @param flow_direction X=0, Y=1, Z=2
 
     @attention Method does not correct for ill-configured meshes where all the nodes
     of boundary triangles touch grain edges.
     
     @todo SKM distinguish properly between value and flag changes
     
     tested:  22-11-2016 SKM
*/
void assignNoSlipBoundaryConditions( Model<3U>& model,
                                     Region<3U>& stokes_flow_domain,
                                     Region<3U>& grain_edges,
                                     ushort flow_direction )
 {
    // flow direction 0=x
    string inlet("LEFT"), outlet("RIGHT");
    if ( flow_direction == 1 ) { // Y-direction
         inlet="BOTTOM", outlet="TOP";
      }
    else if ( flow_direction == 2 ) { // Z-direction with increasing coordinates
         inlet="BACK", outlet="FRONT";
      }
    assert( flow_direction <=2 );
   
    // --------------------------------------------------------
    // 1. Zapping pre-existing conditions
    // --------------------------------------------------------
    // by default all velocity scalars are flagged ANY
    const double64  slip(0.);
    stokes_flow_domain.InputPropertyValue("nodal velocity x", makeScalar(DIRICH,slip), PERIMETER );
    stokes_flow_domain.InputPropertyValue("nodal velocity y", makeScalar(DIRICH,slip), PERIMETER );
    stokes_flow_domain.InputPropertyValue("nodal velocity z", makeScalar(DIRICH,slip), PERIMETER );

    // opening the inflow and outflow boundaries
    // -----------------------------------------
/*
    model.Boundary(inlet).ChangePropertyStatus( "nodal velocity x", ANY, INTERIOR );
    model.Boundary(inlet).ChangePropertyStatus( "nodal velocity y", ANY, INTERIOR );
    model.Boundary(inlet).ChangePropertyStatus( "nodal velocity z", ANY, INTERIOR );
    model.Boundary(outlet).ChangePropertyStatus( "nodal velocity x", ANY, INTERIOR );
    model.Boundary(outlet).ChangePropertyStatus( "nodal velocity y", ANY, INTERIOR );
    model.Boundary(outlet).ChangePropertyStatus( "nodal velocity z", ANY, INTERIOR );
*/

/* does not work !
    changeBoundaryFlags( model, inlet.c_str(), INTERIOR, "nodal velocity x", ANY );
    changeBoundaryFlags( model, inlet.c_str(), INTERIOR, "nodal velocity y", ANY );
    changeBoundaryFlags( model, inlet.c_str(), INTERIOR, "nodal velocity z", ANY );
    changeBoundaryFlags( model, outlet.c_str(), INTERIOR, "nodal velocity x", ANY );
    changeBoundaryFlags( model, outlet.c_str(), INTERIOR, "nodal velocity y", ANY );
    changeBoundaryFlags( model, outlet.c_str(), INTERIOR, "nodal velocity z", ANY );
*/
    // works
    changeBoxBoundaryFlags( model, parseBoundary(inlet), "nodal velocity x", ANY );
    changeBoxBoundaryFlags( model, parseBoundary(inlet), "nodal velocity y", ANY );
    changeBoxBoundaryFlags( model, parseBoundary(inlet), "nodal velocity z", ANY );
    changeBoxBoundaryFlags( model, parseBoundary(outlet), "nodal velocity x", ANY );
    changeBoxBoundaryFlags( model, parseBoundary(outlet), "nodal velocity y", ANY );
    changeBoxBoundaryFlags( model, parseBoundary(outlet), "nodal velocity z", ANY );
/*  assigns to the entire boundary
    vector<VARIABLE_FLAG> flags(1,ANY);
    model.InputBoundaryFlags( parseBoundary(inlet), "nodal velocity x", flags );
    model.InputBoundaryFlags( parseBoundary(inlet), "nodal velocity y", flags );
    model.InputBoundaryFlags( parseBoundary(inlet), "nodal velocity z", flags );
    model.InputBoundaryFlags( parseBoundary(outlet), "nodal velocity x", flags );
    model.InputBoundaryFlags( parseBoundary(outlet), "nodal velocity y", flags );
    model.InputBoundaryFlags( parseBoundary(outlet), "nodal velocity z", flags );
*/
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
    // setting up the side boundaries to allow boundary-parallel flow
    // these are TOP, BOTTOM, FRONT BACK
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
              if ( edge == EDGE1 or edge == EDGE3 or edge == EDGE9 or edge == EDGE11 ) {
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
              if ( edge == EDGE2 or edge == EDGE4 or edge == EDGE9 or edge == EDGE12 ) {
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
              if ( edge == EDGE5 or edge == EDGE6 or edge == EDGE7 or edge == EDGE8 ) {
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
     Reports phi, k and averages to file
*/
void printAnalysisResults1( const Model<3U>& model, double64 porosity, const vector<double64>&  k_equiv )
 {
    assert( k_equiv.size() == 3U );
    ofstream  ofs( string(model.Name()) + "-simulation-results.txt" );
  
    ofs << "model: '"<<  model.Name() <<"', ";
    ofs << model.Mesh().Nodes() <<" nodes, "<< model.Mesh().Elements() <<" elements\n";
    ofs <<"porosity(computed): "<< porosity <<"\nequivalent permeability (kx,ky,kz):";
    double64 average_k(0.);
    for ( size_t i=0U; i<3U; ++i ) {
         ofs <<" "<< k_equiv[i];
         average_k += k_equiv[i];
      }
    average_k /= static_cast<double64>(3U);
    ofs <<"\naverage k: "<< average_k;
 }



/**
    Computes flux through the boundary of the model using quadratic finite-element integration.
    
    1. we start with velocities on the nodes of the quadratic elements
    
    2. we integrate these using the quadrature rules of the boundary (lower-dimensional) faces
    
    A single unit normal is used for the flux calculation. 
    It is taken from the Box class.
*/
template<size_t dim>
double64 boxBoundaryFluxFEM_Quadratic( const Boundary<dim>& boundary, const csmp::Index& velo_key, double64& boundary_area)
 {
   // todo to find face area and maximum velocity magnitude at boundary
    vector<double64>     nrml;
    Box().UnitNormalTo( parseBoundary(boundary.Name()), dim, nrml );
//cout<<boundary.Name()<<endl;

    VectorVariable<dim>  velo;
   double64             boundary_flux(0.);
   // verification
   double64 maxflux(0.);
  
    for ( typename vector<Face<dim>*>::const_iterator
          it=boundary.ElementsBegin(); it!=boundary.ElementsEnd(); ++it )
      {
         double64 face_flux(0.), face_area(0.);
         for ( size_t i=0U; i<(*it)->IntegrationPoints(); ++i )
           {
              (*it)->PropertyValueAtIntegrationPoint( velo_key, i, velo );
              // project velocity onto boundary normal
              double64 face_normal_flux(0.);
              for ( size_t j=0U; j<dim; ++j )
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

template double64 boxBoundaryFluxFEM_Quadratic( const Boundary<3U>&, const csmp::Index&, double64& );




/**
    Assuming a box-shaped model, function calculates the equivalent permeability parallel to 
    its side boundaries.
*/
double64 equivalentPermeabilityQuadraticFEM( Model<3U>& model, const string& region,
                                             BOX_BOUNDARY inflow_boundary, BOX_BOUNDARY outflow_boundary,
                                             double64& flux_through_model )
 {
    const size_t DIM(3);

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    const Region<DIM>&  flow_domain(model.Region(region));
   
    // 1. calculating the cross-sectional area of the model
    // ---------------------------------------------------
    // (2 options: left-right and top bottom)
    Point<DIM>  xyz_min, xyz_max;
    model.MinMaxCoordinates( xyz_min, xyz_max );
    // assumption is that flow boundaries are opposite to one another
    double64 length = xyz_max[0] - xyz_min[0];
    if      ( inflow_boundary == BOTTOM ) length = xyz_max[1] - xyz_min[1];
    else if ( inflow_boundary == BACK )   length = xyz_max[2] - xyz_min[2];
cout<<"length = "<<length<<endl;
   // = true only for box-shaped model
    double64 xsect_area(length * length);
   
    // 2. calculating the far-field fluid pressure gradient
    // -----------------------------------------------------
    // O.K.
    const csmp::Index pf_key = model.Database().StorageKey("fluid pressure");
    double64 inflow_p(0.), outflow_p(0.);
    size_t   in_nodes(0U), out_nodes(0U);
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
    inflow_p  /= static_cast<double64>(in_nodes);
    outflow_p /= static_cast<double64>(out_nodes);
   
    double64 farfield_pf_gradient = (inflow_p - outflow_p) / length;
   cout<<"inflow p= "<<inflow_p<<" outflow p= "<<outflow_p<<endl;

	  // 3a. fluid throughput through model as determined by velocity projections on model boundaries
    // -------------------------------------------------------------------------------------------
    const csmp::Index vt_key = model.Database().StorageKey("nodal velocity");
   double64 inlet_area(0.), outlet_area(0.);
    const double64  inflow  = fabs(boxBoundaryFluxFEM_Quadratic( model.Boundary(parseBoundary(inflow_boundary)), vt_key, inlet_area ));
    const double64  outflow = fabs(boxBoundaryFluxFEM_Quadratic( model.Boundary(parseBoundary(outflow_boundary)), vt_key, outlet_area ));

    // computing the mean of the in- and outflow
    flux_through_model = (fabs(inflow) + fabs(outflow)) / 2.;
    // making sure that the difference between in- and outflux from the model is less than 1%
//    const double64 flux_discrepancy_tolerance(0.01);
//    if ( fabs(fabs(inflow)-fabs(outflow)) / flux_through_model > flux_discrepancy_tolerance ) {
         cerr <<"\naveraged influx versus averaged outflux: "<< inflow <<" vs. "<< outflow <<" m3/s.";
   cout<<"relative flux error = "<<fabs(fabs(inflow)-fabs(outflow)) / flux_through_model <<endl;
//         csmp_error.notice( ERROR, "equivalentPermeabilityQuadraticFEM:", "flux analysis is inaccurate; check integrity of boundary mesh." );
         // using the maximum flux estimate in this case
//         flux_through_model = std::max( fabs(inflow), fabs(outflow) );
//      }
   
	  cout <<"\nequivalentPermeabilityQuadraticFEM: average flux through inlet and outlet boundaries: "<< flux_through_model;
	  // for bcc
    //cout <<" m3/s, flux per m2: "<< flux_through_model / xsect_area <<" m/s."<< endl;
    //for tubes
    cout <<" m3/s, flux per m2: "<< 2* flux_through_model / (inlet_area+outlet_area) <<" m/s."<< endl;
//    cout <<"\n\testimated error of calculation (%): <= "<< 100. * (fabs(fabs(inflow)-fabs(outflow)) / flux_through_model) << endl << endl;
   
    // recovering fluid viscosity
    double64  mumin, fluid_viscosity;
    flow_domain.MinMaxOf( "viscosity", mumin, fluid_viscosity );
    // asserting that there is only a single viscosity value
    assert( fabs(fluid_viscosity - mumin) <= numeric_limits<double64>::epsilon() );
   
cout<<"tube model: permeability 1 (over inlet/outlet area) = "<<(2*flux_through_model * fluid_viscosity) / ((inlet_area+outlet_area) * farfield_pf_gradient)<<" m2"<<endl;
   
   // 3b. fluid throughput through model, averaged over pore space
   // -------------------------------------------------------------------------------------------
   double64 velo_integral(0.);
   double64 avg_velocity = volumeAveragedVelocity( model, region, velo_integral);
   
   // directly from Darcy's law where v_avg = Darcy's velocity
   cout<<"averaged velocity over pore space, m/s= "<<avg_velocity<<endl;
   cout <<"tube model: permeability 2 (from volume-averaged velocity) = "<< (avg_velocity * fluid_viscosity) / farfield_pf_gradient <<" m2"<<endl;
   cout <<"not for tubes, permeability 3( from volume-averaged velocity/unit cell volume) = "<< (velo_integral * fluid_viscosity) / (farfield_pf_gradient* length* length *length ) <<" m2"<<endl;
   // velocity analysed element by element after interpolation to barycenter
   // TODO: analyse velocity at the element integration points using their weights (rather than at the element barycentre)

   
    // equivalent permeability, old version
	  return (flux_through_model * fluid_viscosity) / (xsect_area * farfield_pf_gradient);

 } // end equivalentPermeabilityQuadraticFEM (using box-model boundaries)






/**
    Computation using velocity magnitude in the averaging.
*/
double64  volumeAveragedVelocity( const Model<3U>& model, const std::string& flow_domain, double64& velo_integral )
 {
   double64 volume_integral(0.);
//double64 tot_velocity_element(0.), tot_volume_element(0.);
     const Region<3U>& model_domain(model.Region(flow_domain));
     const csmp::Index v_key(model.Database().StorageKey("nodal velocity"));
     VectorVariable<3> velocity;

   double64 patm(100325.);
   
   const csmp::Index p_key(model.Database().StorageKey("fluid pressure"));

     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         for ( size_t i=0; i<(*it)->IntegrationPoints(); ++i )
           {
             if ((*it)->PropertyValueAtIntegrationPoint( p_key, i) >= patm)
             {
              (*it)->PropertyValueAtIntegrationPoint( v_key, i, velocity );
              double64 velo_magnitude = velocity.Length();
              double64 det_J((*it)->det_JINV_AtIntegrationPoint(i));
              velo_integral   += velo_magnitude * det_J * (*it)->WeightAtIntegrationPoint(i);
              // test: error
              volume_integral += det_J * (*it)->WeightAtIntegrationPoint(i);
              }
           }
       // diff way to calc
//       tot_volume_element +=(*it)->Volume();
 //      (*it)->Read(v_key,velocity);
 //      tot_velocity_element += velocity.Length() * (*it)->Volume();
       }
   cout<<"volumeAveragedVelocity: velo_integral= "<< velo_integral << " volume = "<<volume_integral<<endl;
// diff way to calc
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
void computePressureLEFT_RIGHT( Model<3>& model )
 {
    Boundary<3>&  inlet(model.Boundary("LEFT"));
    Boundary<3>&  outlet(model.Boundary("RIGHT"));
    const double64 patm(100325.);
    inlet.InputPropertyValue( "fluid pressure", makeScalar(DIRICH,2. * patm) );
    outlet.InputPropertyValue( "fluid pressure", makeScalar(DIRICH, patm) );
   
    // 1. Set-up to compute fluid pressure
    SAMG_Settings settings;
    settings.Set_ncgtyp(5);    // strong connectivity
    settings.Set_nxtyp(0);     // Gauss-Seidl relaxation
    settings.Set_ndefault(40); // standard coarsening strategy
    SAMG_Solver   solver(&settings);
    PDE_IntegratorExperimental<3,Region>  pressure_solver(solver);
    NumIntegral_dNT_dN_dV<3>              laplacian( model.Database(), "fluid pressure", "fluid pressure" );
    NumIntegral_SetRHS_to_Zero<3>         zero( model.Database(), "fluid pressure" );
    pressure_solver.Add( &laplacian );
    pressure_solver.Add( &zero );
   
    Region<3>& pore_domain(model.Region("PORES"));
    pressure_solver.IntegrateOver( pore_domain );
    printRangeOfVariable( model, "PORES", "fluid pressure" );
 }




// TESTING helper for function below
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




/**
    Sets the flags of scalar node variables on specified boundary / subregion thereof
*/
void changeBoundaryFlags( Model<3U>& model, const char* boundary, SUBDOMAIN_PART part, const char* node_variable, VARIABLE_FLAG flag )
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
void changeBoxBoundaryFlags( Model<3U>& model, BOX_BOUNDARY boundary, const char* node_variable, VARIABLE_FLAG flag )
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








/*  CODE FRAGMENTS FOR CHECKING

    // Extract from the first mesh the pressure values at the nodes.
    if ( linear_model.RestoreOriginalNodeNumbering(false) ) {
         // if changes were made, restored numbering and point coordinates are output to file
         ofstream ofs("linear_model-original-node-numbering.txt");
         ofs <<"linear model: restored original node numbers and coordinates:\n";
         for ( size_t i=0U; i<pore_domain.Nodes(); ++i ) {
              Point<dim> xyz(pore_domain.N(i)->Coordinate());
              ofs << pore_domain.N(i)->Idx() <<" "<< xyz[0] <<" "<< xyz[1] <<" "<< xyz[2] <<"\n";
           }
      }



    // checking that the node coordinates of the linear and quadratic models are the same
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

} // end csmp



