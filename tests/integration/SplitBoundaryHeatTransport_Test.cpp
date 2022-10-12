//
//  SplitBoundaryHeatTransport_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 2/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "Model.h"
#include "vsetMakers.h"
#include "SplitBoundaryHeatTransport_Test.h"
#include "VTK_Interface.h"
#include "Index.h"
#include "ANSYS_Model3D.h"
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "LinearSolver.h"
#include "ExplicitTransport.h"

// is being tested
#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp {

SplitBoundaryHeatTransport_Test::SplitBoundaryHeatTransport_Test( const char* test_model,
                                                                  const char* test_variables )
 : test_model_(test_model),
   test_variable_file_(test_variables),
   model_ptr_(nullptr)
 {
 }



SplitBoundaryHeatTransport_Test::~SplitBoundaryHeatTransport_Test()
{
   delete model_ptr_;
}


/*

Model<3U>*  SplitBoundaryHeatTransport_Test::CreateModel( const char* ansys_input_data )
 {
    delete model_ptr_;
    model_ptr_ = new ANSYS_Model3D( ansys_input_data, ansys_input_data, "SplitBoundaryHeatTransport_Test-variables.txt" );
    Point<3U> min_coord, max_coord;
    model_ptr_->MinMaxCoordinates( min_coord, max_coord );
    model_length_ = max_coord[0] - min_coord[0];
    model_height_ = max_coord[1] - min_coord[1];
    model_width_  = max_coord[2] - min_coord[2];
    if ( verbose_ ) printModelDimensions( *model_ptr_ );
    return model_ptr_;
 }


Model<3U>*  SplitBoundaryHeatTransport_Test::CreateTetrahedralModel()
 {
    return CreateModel("BOX40x3x10m");
 }



Model<3U>*  SplitBoundaryHeatTransport_Test::CreateHexahedralModel()
 {
    delete model_ptr_;
    VSet<3U>  vset;
 
    // create test model
    const bool bSkewed(false);
    // is already isoparametric
    test_Create_Hexahedra_VSet( vset, bSkewed ); // only hexahedral elements
    //                                singlePhase_advection-variables.txt
    model_ptr_ = new Model<3U>( vset, test_variable_file_.c_str() );
    Point<3U> min_coord, max_coord;
    model_ptr_->MinMaxCoordinates( min_coord, max_coord );
    model_length_ = max_coord[0] - min_coord[0];
    model_height_ = max_coord[1] - min_coord[1];
    model_width_  = max_coord[2] - min_coord[2];
    if ( verbose_ ) printModelDimensions( *model_ptr_ );
    return model_ptr_;
 }


void SplitBoundaryHeatTransport_Test::run()
 {
    delete model_ptr_;
    model_ptr_ = new ANSYS_Model3D( ansys_input_data, ansys_input_data, "SplitBoundaryHeatTransport_Test-variables.txt" );
    Point<3U> min_coord, max_coord;
    model_ptr_->MinMaxCoordinates( min_coord, max_coord );
    model_length_ = max_coord[0] - min_coord[0];
    model_height_ = max_coord[1] - min_coord[1];
    model_width_  = max_coord[2] - min_coord[2];
    if ( verbose_ ) printModelDimensions( *model_ptr_ );
    return model_ptr_;


    model_ptr_ = CreateTetrahedralModel();
    Region<3U>  model_domain = model_ptr_->Region("Model");
    AssignFlowProperties();
    const bool initialize_flux( true ); // prescibed 'total velocity'
    initializeFiniteVolumeProperties( *model_ptr_, model_domain, initialize_flux );
 
    // test 1: volume, pore volume and prescribed 'total velocity'
    // -----------------------------------------------------------
    Test_initializeFiniteVolumeProperties();

    // test 2: flux balance for computed divergence free 'total velocity' field
    // ------------------------------------------------------------------------
    // computing divergence free 'total velocity' field and 'facet flux'
    DivergenceFreeTotalVelocityField();
    ExplicitTransport<3U>  transport( *model_ptr_, "Model" );

    // 2.1 flux balance in the model interior for computed velocity
    // -------------------------------------------------------------------------
    TestInteriorFluxBalance();

    // 2.2 flux balance at no-flow boundaries
    // -------------------------------------------------------------------------
    TestNoFlowBoundaryFluxBalance();
    _equal( transport.IncomingVolumetricFlow(),
            transport.OutgoingVolumetricFlow(), numeric_limits<double>::epsilon() * transport.IncomingVolumetricFlow() );
 
    // test 3: flow through model with TVD concentration
    // -------------------------------------------------------------------------
    const bool prescribed_velocity(true);
    TestFlowThroughModel( "BOX40x3x10m", prescribed_velocity );
    TestFlowThroughModel( "BOX40x3x10m", false );

    // void test_Create_Prism_VSet(VSet<3U>& vset, bool bSkewed=false );
    

    // void test_Create_Prism_Hexa_VSet(VSet<3U>& vset, bool bSkewed=false );
 }





    /// thickness, permeability, porosity, total velocity
void  SplitBoundaryHeatTransport_Test::AssignFlowProperties()
 {
    assert( model_ptr_ != nullptr );
     // property assigment
    model_ptr_->InputPropertyValue( "thickness",     makeScalar(ANY,1.0) );
    model_ptr_->InputPropertyValue( "permeability",  makeScalar(ANY,1.0e-12) );
    model_ptr_->InputPropertyValue( "conductivity",  makeScalar(ANY,1.0e-9) ); // dyn visc = 1.0e-3
    model_ptr_->InputPropertyValue( "porosity",      makeScalar(ANY,0.25) );
    model_ptr_->InputPropertyValue( "concentration", makeScalar(ANY,0.) );
    model_ptr_->InputPropertyValue( "fluid volume source", makeScalar(ANY,0.) );
    model_ptr_->InputPropertyValue( "nodal fluid volume source", makeScalar(ANY,0.) );
    model_ptr_->InputPropertyValue( "fluid viscosity", makeScalar(ANY,1.0e-3) );
    VTK_Interface<3U> vtk_output;
    vtk_output.OutputDataToVTK( *model_ptr_, "test_output1", "permeability", 0 );
    // call of testee
    // here the facet fluxes as precomputed are used
    VectorVariable<3U> velo(ANY,ANY,ANY, 1., 1., 1. );
    velo.EuclideanNormalize(); // to 1.
    model_ptr_->InputPropertyValue( "total velocity", velo );
    if ( verbose_ ) {
         cout <<"\nSplitBoundaryHeatTransport_Test::AssignFlowProperties: 'total velocity' magnitude: "<< velo.Length() <<"\n";
         printRangeOfVariable( *model_ptr_, "permeability" );
         printRangeOfVariable( *model_ptr_, "porosity" );
      }
//    printRangeOfVariable( *model_ptr_, "total velocity" );
}




void  SplitBoundaryHeatTransport_Test::DivergenceFreeTotalVelocityField( double delta_pf )
 {
    // 1.  Building the steady-state FE Algorithm "fluid_pressure"
    // -----------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings;
    settings.Set_eps(0.);
    SAMG_Solver    solver(&settings);
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
    PDE_Integrator<3U,Element>  fluid_pressure( solver );

    NumIntegral_dNT_op_dN_dV<3U> conductance( model_ptr_->Database(), "conductivity", "fluid pressure",  "fluid pressure" );
    NumIntegral_NT_op_N_dV<3U>   source( model_ptr_->Database(),  "fluid volume source", "fluid pressure" );
  /// @todo influx surface integral: NumIntegral_NT_op_N_dS<2U>   influx( model.Database(), "influx", "fluid pressure" );
    VelocityAndVolumeFlux<3U>    velocity( *model_ptr_,  "conductivity", "porosity", "fluid pressure", false );

    // add PDE_Operators and post-processor to the FE Algorithm
    fluid_pressure.Add( &conductance );
    fluid_pressure.Add( &source );
    //  fluid_pressure.Add( &influx );
    fluid_pressure.AddPostProcess( &velocity );


    // 2. Assign boundary conditions
    // -----------------------------
    // pressure range between 1 bar and (1 bar + delta_pf)
    const double bar(100325.);
    model_ptr_->InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH,bar+delta_pf) );
    model_ptr_->InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,bar) );


    // 3.  Pass the FE algorithm to the Model and reset it after the solution
    // ----------------------------------------------------------------------
    fluid_pressure.IntegrateOver( model_ptr_->Region("Model") );
    model_ptr_->CopyReplace( "velocity", "total velocity" );


    // 4.  Output the initial range of the variables
    // ---------------------------------------------
    printRangeOfVariable( *model_ptr_, "fluid pressure" );
    printRangeOfVariable( *model_ptr_, "velocity" );
    printRangeOfVariable( *model_ptr_, "pore velocity" );
    printRangeOfVariable( *model_ptr_, "volume flux" );


    // 5.  Output the initial conditions to VTK
    // -----------------------------------------
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( *model_ptr_, "fluid-pressure", "fluid pressure", 0, true );
    vtk_output.OutputDataToVTK( *model_ptr_, "velocity",       "velocity",       0, true );
    vtk_output.OutputDataToVTK( *model_ptr_, "volume-flux",    "volume flux",    0, true );
 
 } // end




double SplitBoundaryHeatTransport_Test::FluxMultiplier( const Element<3U>* const eptr, size_t sector, size_t facet ) const
 {
    // TODO: do this only once inside the finite-volume stencil:
    // for each sector, record whether the facet normal is inward or outward pointing
    vector<vector<short> > sign_of_facet( eptr->Facets(), vector<short>(eptr->Nodes(),0) );
    for ( auto i{0}; i<eptr->Nodes(); ++i ) {
         for ( auto j=0U; j<eptr->FV()->FacetsPerSector(i); ++j ) {
               auto s_facet      = eptr->FV()->FacetSurroundingSector( i, j );
               auto inside_node  = eptr->FV()->InsideNode( s_facet );
               sign_of_facet[s_facet][i] = (inside_node==i) ? 1 : -1;
            }
      }
 
    // using the new array to return the facet flux multiplier
    if ( sign_of_facet[facet][sector] == -1 ) return -1.;
    return 1.;
 
 } // end FluxMultiplier

*/

} // end csmp
