//
//  ExplicitTransport_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 2/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "Model.h"
#include "vsetMakers.h"
#include "ExplicitTransport_Test.h"
#include "VTK_Interface.h"
#include "Index.h"
#include "ANSYS_Model3D.h"
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#include "PDE_IntegratorExperimental.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "ExplicitTransport.h"

// is being tested
#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp {

ExplicitTransport_Test::ExplicitTransport_Test( const char* test_model, const char* test_variables )
 : test_model_(test_model),
   test_variable_file_(test_variables)
 {
 }

ExplicitTransport_Test::~ExplicitTransport_Test()
{
   delete model_ptr_;
}


/**
    Any kind, dependent on input data
*/
Model<3U>*  ExplicitTransport_Test::CreateModel( const char* ansys_input_data )
 {
    delete model_ptr_;
    model_ptr_ = new ANSYS_Model3D( ansys_input_data, ansys_input_data, "ExplicitTransport_Test-variables.txt" );
    Point<3U> min_coord, max_coord;
    model_ptr_->MinMaxCoordinates( min_coord, max_coord );
    model_length_ = max_coord[0] - min_coord[0];
    model_height_ = max_coord[1] - min_coord[1];
    model_width_  = max_coord[2] - min_coord[2];
    if ( verbose_ ) printModelDimensions( *model_ptr_ );
    return model_ptr_;
 }


/**
    Using "BOX40x3x10m" as an input dataset.
*/
Model<3U>*  ExplicitTransport_Test::CreateTetrahedralModel()
 {
    return CreateModel("BOX40x3x10m");
 }



/**
    // BUILD VSET MAKER PATCH POLY-ELEMENT MODEL
    // =========================================
*/
Model<3U>*  ExplicitTransport_Test::CreateHexahedralModel()
 {
    delete model_ptr_;
    VSet<3U>  vset;
 
    // create test model
    const bool bSkewed(false);
    const bool isoparametric(true);
    test_Create_Hexahedra_VSet( vset, bSkewed ); // only hexahedral elements
    //                                singlePhase_advection-variables.txt
    model_ptr_ = new Model<3U>( vset, test_variable_file_.c_str(), isoparametric );
    Point<3U> min_coord, max_coord;
    model_ptr_->MinMaxCoordinates( min_coord, max_coord );
    model_length_ = max_coord[0] - min_coord[0];
    model_height_ = max_coord[1] - min_coord[1];
    model_width_  = max_coord[2] - min_coord[2];
    if ( verbose_ ) printModelDimensions( *model_ptr_ );
    return model_ptr_;
 }


/**
    Tests performed:
    
    1. Verify that we have indeed the right porevolumes facet areas and facet normals
    
    use BoxHalves3D - tetra model with 2 regions for tetra testing
    For poly-elemement mesh testing use model b25
 
    2. Flux balance test
 
    3. TVD test
 
    4. Boundary interaction tests
*/
void ExplicitTransport_Test::run()
 {
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
            transport.OutgoingVolumetricFlow(), numeric_limits<double64>::epsilon() * transport.IncomingVolumetricFlow() );
 
    // test 3: flow through model with TVD concentration
    // -------------------------------------------------------------------------
    const bool prescribed_velocity(true);
    TestFlowThroughModel( "BOX40x3x10m", prescribed_velocity );
    TestFlowThroughModel( "BOX40x3x10m", false );

    // void test_Create_Prism_VSet(VSet<3U>& vset, bool bSkewed=false );
    

    // void test_Create_Prism_Hexa_VSet(VSet<3U>& vset, bool bSkewed=false );
 }





    /// thickness, permeability, porosity, total velocity
void  ExplicitTransport_Test::AssignFlowProperties()
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
         cout <<"\nExplicitTransport_Test::AssignFlowProperties: 'total velocity' magnitude: "<< velo.Length() <<"\n";
         printRangeOfVariable( *model_ptr_, "permeability" );
         printRangeOfVariable( *model_ptr_, "porosity" );
      }
//    printRangeOfVariable( *model_ptr_, "total velocity" );
}




/**
    Applies uniform pressures on left and right side of the model,
    and solves for the steady-state fluid pressure distribution in the
    absence of fluid sources and sinks.
*/
void  ExplicitTransport_Test::DivergenceFreeTotalVelocityField( double64 delta_pf )
 {
    // 1.  Building the steady-state FE Algorithm "fluid_pressure"
    // -----------------------------------------------------------
    SAMG_Settings  settings;
    settings.Set_eps(0.);
    SAMG_Solver  samg_solver(&settings);
    PDE_IntegratorExperimental<3U,Region>  fluid_pressure(samg_solver);

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
    const double64 bar(100325.);
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




/**
    Establishes whether a facet normal is outward pointing relative to the specific FV sector
 
    for each facet, there are 2 sectors it contributes fluxes to
    [facets][sectors] = multiplier
 
    trial version for wider use in finite-element finite-volume computations.
 
*/
double64 ExplicitTransport_Test::FluxMultiplier( const Element<3U>* const eptr, size_t sector, size_t facet ) const
 {
    // TODO: do this only once inside the finite-volume stencil:
    // for each sector, record whether the facet normal is inward or outward pointing
    vector<vector<short> > sign_of_facet( eptr->Facets(), vector<short>(eptr->Nodes(),0) );
    for ( size_t i=0U; i<eptr->Nodes(); ++i ) {
         for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(i); ++j ) {
               size_t facet        = eptr->FV()->FacetSurroundingSector( i, j );
               size_t inside_node  = eptr->FV()->InsideNode( facet );
               sign_of_facet[facet][i] = (inside_node==i) ? 1 : -1;
            }
      }
 
// TEST output
/*
cout <<"\nFluxMultiplier: element (facets,sectors) multiplier matrix:\n";
for ( size_t i=0U; i<eptr->Facets(); ++i ) {
     for ( size_t j=0U; j<sign_of_facet[i].size(); ++j ) cout <<"("<< i <<","<< j <<"): "<< sign_of_facet[i][j] <<" ";
     cout << endl;
  }
*/
    // using the new array to return the facet flux multiplier
    if ( sign_of_facet[facet][sector] == -1 ) return -1.;
    return 1.;
 
 } // end FluxMultiplier


// =====================================================================================================================
//
//      TESTING
//
// =====================================================================================================================



/**
    Testing (facet area, facet normal were already tested in FV_Stencil_Test)

    "finite volume", "FV pore volume", "sector volume", "sector pore volume",
    "facet flux", "flux balance"
 
    The flux balance is tested with a prescribed velocity field.
*/
void ExplicitTransport_Test::Test_initializeFiniteVolumeProperties( double64 tolerance_relaxation_factor )
 {
    // output variables
    const csmp::Index fv_key  = model_ptr_->Database().StorageKey("finite volume");
    const csmp::Index pv_key  = model_ptr_->Database().StorageKey("FV pore volume");
    const csmp::Index sv_key  = model_ptr_->Database().StorageKey("sector volume");
    const csmp::Index spv_key = model_ptr_->Database().StorageKey("sector pore volume");
    const csmp::Index phi_key = model_ptr_->Database().StorageKey("porosity");
    const csmp::Index ff_key  = model_ptr_->Database().StorageKey("facet flux");
    const csmp::Index fb_key  = model_ptr_->Database().StorageKey("flux balance");

    // test 1: is the finite volume equal to the element volume
    // --------------------------------------------------------
    const csmp::Index  vol_key(model_ptr_->Database().StorageKey("finite volume"));
    double64           total_volume(0.), total_PV(0.);
    Region<3U>         model_domain(model_ptr_->Region("Model"));
    const double64     model_volume = model_domain.Volume(); // finite element estimate

    for ( vector<Node<3U>*>::iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
         total_volume += (*nit)->Read( vol_key );
         total_PV     += (*nit)->Read( pv_key );
      }
    if ( verbose_ ) cout <<"\nrun: model volume vs. finite volume integrated: "<< model_volume <<" vs "<< total_volume << endl;
    // TODO: test fails for skewed hexahedra, fix:
    _equal( total_volume, model_volume, numeric_limits<double64>::epsilon() * model_volume * tolerance_relaxation_factor );
    // sum of sector volumes
    double64     sector_volume(0.), sector_PV(0.), FE_PV(0.);
    const size_t sector_ip(0U);
    for ( vector<Element<3U>*>::const_iterator it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        FE_PV += (*it)->Volume() * (*it)->Read( phi_key );
        for ( size_t i=0U; i<(*it)->Nodes(); i++ ) {
             sector_volume += (*it)->Read( i, sector_ip, sv_key ); // larger tolerance needed presumable because sectors are hexahedra
             sector_PV     += (*it)->Read( i, sector_ip, spv_key );
          }
      }
    _equal( sector_volume, model_volume, numeric_limits<double64>::epsilon() * model_volume * tolerance_relaxation_factor ); // FV vs. FE
    _equal( sector_PV, FE_PV, numeric_limits<double64>::epsilon() * FE_PV * tolerance_relaxation_factor );       // sector PV vs. FE PV
    _equal( sector_PV, total_PV, numeric_limits<double64>::epsilon() * total_PV * tolerance_relaxation_factor ); // sector PV, vs. FV PV

 
    // test 2: is the flux conserved in the interior of the model (prescribed velocity case)
    //         loop over the element sectors
    // -------------------------------------------------------------------------------------
    const csmp::Index  pf_key(model_ptr_->Database().StorageKey("fluid pressure"));
    const csmp::Index  velo_key(model_ptr_->Database().StorageKey("velocity"));
    const csmp::Index  flux_key(model_ptr_->Database().StorageKey("facet flux"));
    model_domain.UpdateMemberIndexes();
    vector<double64>  flux_balance( model_domain.Nodes(), 0. );
    // NB: establishing the flux balance in an element loop, which must include the perimeter elements,
    //     but avoiding truncated FVs at boundaries
    // double64 sectorFlux( const Element<dim>* const eptr, size_t sector, const csmp::Index& flux_key ); // tested: O.K.
    for ( vector<Element<3U>*>::iterator it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         for ( size_t i=0U; i<(*it)->Nodes(); ++i ) {
              if ( (*it)->N(i)->AtBoundary() == NOT ) {
                   size_t node = (*it)->N(i)->Idx();
                   flux_balance[ node ] += sectorFlux( (*it), i, flux_key );
                }
           }
      }
    auto min_value = *min_element(flux_balance.begin(),flux_balance.end());
    auto max_value = *max_element(flux_balance.begin(),flux_balance.end());
    _equal( fabs(min_value), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );
    _equal( fabs(max_value), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );
    if ( verbose_ )
      cout <<"\nrun: velocity vs. finite volume flux balances (min/max) for total velocity of 1.: "<< min_value <<" to "<< max_value << endl;

    // second version, loop over the element facets
    // --------------------------------------------
    fill( flux_balance.begin(), flux_balance.end(), 0. );
    for ( vector<Element<3U>*>::iterator it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         // Attention: inside/outside is with reference to facet, and NOT to finite-element sector
         //            when outside-inside>1 the opposite sign needs to be applied because facet is between first and last node !
         for ( size_t i=0U; i<(*it)->Facets(); ++i ) {
              const size_t facet_ip(0U);
              // finding the orientation of the facet
              // inside sector
              size_t   node = (*it)->FV()->InsideNode(i);
              double64 sign = FluxMultiplier( (*it), node, i );
              assert( sign != 0. );
              if ( (*it)->N(node)->AtBoundary() == NOT )
                flux_balance[ (*it)->N(node)->Idx() ]  += sign * (*it)->Read( i, facet_ip, flux_key );
              // outside sector
              node = (*it)->FV()->OutsideNode(i);
              sign = FluxMultiplier( (*it), node, i );
              assert( sign != 0. );
              if ( (*it)->N(node)->AtBoundary() == NOT )
                flux_balance[ (*it)->N(node)->Idx() ]  += sign * (*it)->Read( i, facet_ip, flux_key );
            }
      }
    min_value = *min_element(flux_balance.begin(),flux_balance.end());
    max_value = *max_element(flux_balance.begin(),flux_balance.end());
    _equal( fabs(min_value), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );
    _equal( fabs(max_value), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );

    // reporting in case of a failure
    if ( verbose_ && fabs(max_value) >= numeric_limits<double64>::epsilon() ) {
         cerr <<"\nExplicitTransport_Test::run: error in flux balance for elements:\n";
         for ( size_t i=0U; i<flux_balance.size(); ++i ) {
             if ( fabs(flux_balance[i]) >= numeric_limits<double64>::epsilon() ) {
                  cerr <<"\n\tflux balance error: "<< flux_balance[i];
                  model_domain.N(i)->Out();
               }
           }
      }

} // end Test_initializeFiniteVolumeProperties




/**
    Tests the flux balance on intact finite volumes within the domain and its halo.
 
    @attention method assumes that the facet flux values were computed elsewhere.
 
    TODO: use a non-arbitrary definition of the flux balance
*/
void ExplicitTransport_Test::TestInteriorFluxBalance( double64 tolerance_relaxation_factor )
 {
    const Region<3U>   model_domain(model_ptr_->Region("Model"));
    const csmp::Index  pf_key(model_ptr_->Database().StorageKey("fluid pressure"));
    const csmp::Index  flux_key(model_ptr_->Database().StorageKey("facet flux"));
    model_domain.UpdateMemberIndexes();
    vector<double64>  flux_balance( model_domain.Nodes(), 0. );
    // NB: establishing the flux balance in an element loop, including the perimeter elements,
    //     but avoiding truncated FVs at the model boundary
    for ( vector<Element<3U>*>::const_iterator it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         for ( size_t i=0U; i<(*it)->Nodes(); ++i ) {
              if ( (*it)->N(i)->AtBoundary() == NOT ) {
                   size_t node = (*it)->N(i)->Idx();
                   flux_balance[ node ] += sectorFlux( (*it), i, flux_key );
                }
           }
      }
    auto min_value = *min_element(flux_balance.begin(),flux_balance.end());
    auto max_value = *max_element(flux_balance.begin(),flux_balance.end());
    // the relaxation factor addresses deviations from 1 in the final value
    _equal( fabs(min_value), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );
    _equal( fabs(max_value), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );
//    if ( verbose_ )
    cout <<"\nTestInteriorFluxBalance: velocity vs. finite volume flux balances (min/max): "<< min_value <<" to "<< max_value << endl;
 
 } // end InteriorFluxBalanceTest





/**
    // Is the flux conserved along no-flow boundaries
    // ------------------------------------------------------
    // TODO: test whether this applies when lower-dimensional elements are present at the boundary
    // TODO: deal with the case of Neumann boundary conditions (they are however manifest in source terms)
*/
void ExplicitTransport_Test::TestNoFlowBoundaryFluxBalance( double64 tolerance_relaxation_factor )
 {
    const csmp::Index  pf_key(model_ptr_->Database().StorageKey("fluid pressure"));
    const csmp::Index  flux_key(model_ptr_->Database().StorageKey("facet flux"));
    Region<3U>         model_domain(model_ptr_->Region("Model"));
 
    double64  min_val(1.0e30), max_val(-1.0e30);
    for ( vector<Node<3U>*>::const_iterator nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      if ( (*nit)->Status(pf_key) != DIRICH )
        {
           // looping over the parent elements accumulating their flux contributions
           double64 boundary_flux(0.);
           for ( size_t i=0U; i<(*nit)->Parents(); ++i ) {
                const Element<3U>* eptr = (*nit)->Parent(i);
                const size_t       nid  = (*nit)->ParentNodeNumber(i);
                boundary_flux += sectorFlux( eptr, nid, flux_key );
             }
          // recording range
          min_val = min( min_val, boundary_flux );
          max_val = max( max_val, boundary_flux );
          // since the velocity field is divergence free / there are no fluid sources or sinks, there should not be any flow across the boundary
          _equal( fabs(boundary_flux), 0., numeric_limits<double64>::epsilon() * tolerance_relaxation_factor );
        }

//    if ( verbose_ )
    cout <<"\nTestNoFlowBoundaryFluxBalance: velocity vs. finite volume flux balances (min/max): "<< min_val <<" to "<< max_val << endl;

 } // end TestNoFlowBoundaryFluxBalance








/**
    Model 'cube' - 100m long stick shaped tetrahedral element model.
 
    TODO: not clear what to compare results with, compute reference solution using old transport scheme.
*/
void  ExplicitTransport_Test::TestFlowThroughModel( const char* model, bool prescribed_velocity )
 {
    AssignFlowProperties();
 
    // constant velocity field, left-to-right, velocity = 1m/s
    double64  velo_magnitude(1.);
    if ( prescribed_velocity ) {
          VectorVariable<3U>  vc1(ANY,ANY,ANY,1.,0.,0.);
          model_ptr_->InputPropertyValue( "velocity", vc1 ); // NB: transport scheme uses 'velocity'
          // TODO: adjust other parameters so that method performs the same when velocity is calculated internally
      }
    else {
         // left->right pressure gradient and flow (hydrostatic)
         const double64 delta_pf( 9.8 * 1000. * model_length_ );
         DivergenceFreeTotalVelocityField( delta_pf );
      }
    velo_magnitude = printRangeOfVariable( *model_ptr_, "velocity" ); 
    
    // initial and boundary conditions for tracer transport
    // - concentration
    // - constraints at boundary
    // - initial amount of tracer in the system
    model_ptr_->InputPropertyValue( "concentration", makeScalar(ANY,0.) );
    const double64 inlet_concentration(3.);
    model_ptr_->InputBoundaryValue( LEFT, "concentration", makeScalar(DIRICH,inlet_concentration) );
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 0, true );

    // uses the current velocity field, to initialise facet fluxes in construction
    ExplicitTransport<3U>  transport( *model_ptr_, "Model" );
 
     if ( verbose_ ) {
         cout <<"\nExplicitTransport_Test::TestFlowThroughModel: key variable ranges in current model '"<< model <<"'\n";
         printRangeOfVariable( *model_ptr_, "FV pore volume" );
         printRangeOfVariable( *model_ptr_, "concentration" );
      }
 
    // assumes flow and model long axis are aligned with the X-axis
    Point<3U> xyz_min, xyz_max;
    model_ptr_->MinMaxCoordinates( xyz_min, xyz_max ); 
    const double64 model_length(xyz_max[0]-xyz_min[0]), xsect_area((xyz_max[1]-xyz_min[1]) * (xyz_max[2]-xyz_min[2]));
    const double64 time_interval( (model_length/velo_magnitude) / 100. ); // ~10-m travel distance
    double64       duration(0.); // calculated from velocity and model length

    // 0. testing whether inflow and outflow from the model have the expected values
    // -----------------------------------------------------------------------------
    if ( prescribed_velocity ) {
          cout <<"\nrun: prescribed_velocity 'velocity' magnitude: "<< velo_magnitude << endl;
          const double64 expected_volume_flux(xsect_area * velo_magnitude);
          _equal( transport.IncomingVolumetricFlow(), expected_volume_flux, numeric_limits<double64>::epsilon() * expected_volume_flux );
          _equal( transport.OutgoingVolumetricFlow(), expected_volume_flux, numeric_limits<double64>::epsilon() * expected_volume_flux );
      }

    // 1. tracer tranport and conservation tests
    // -----------------------------------------
    // 1.1 getting some tracer into model
    transport.AdvectVariable( time_interval );
    duration += time_interval;
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 1, true );

    // 1.2 switching supply off and transporting more
    model_ptr_->InputBoundaryValue( LEFT, "concentration", makeScalar(DIRICH,0.) );
    transport.AdvectVariable( time_interval );
    duration += time_interval;
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 2, true );

    // 1.3 integrating this initial tracer concentration
    Region<3U>  model_domain = model_ptr_->Region("Model");
    const bool  multiply_with_porosity(true);
    const double64 initial_concentration = model_domain.VolumeIntegral_x_Thickness( "concentration",  multiply_with_porosity );
 
    // 1.4 transporting for trice the time
    transport.AdvectVariable( time_interval * 3. );
    duration += time_interval * 3.;
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 3, true );

    // 1.5 integration final tracer concentration and comparing total amount of tracer
    const double64 final_concentration = model_domain.VolumeIntegral_x_Thickness( "concentration",  multiply_with_porosity );
 
    // testing
    // -------
    // first TVD test
    const bool print_maximum(true);
    const double64 max_concentration = printRangeOfVariable( *model_ptr_, "concentration", print_maximum );
    _test( max_concentration <= inlet_concentration );
    _test( printRangeOfVariable( *model_ptr_, "concentration", !print_maximum ) >= 0. );
    // tracer conservation test
    _equal( initial_concentration, final_concentration, numeric_limits<double64>::epsilon() * initial_concentration );


    // 2. transporting tracer across outflow boundary, verifying that there is no build up
    // -----------------------------------------------------------------------------------
    transport.AdvectVariable( time_interval * 20. );
    duration += time_interval * 20.;
    _test( printRangeOfVariable( *model_ptr_, "concentration", print_maximum ) <= inlet_concentration );
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 4, true );


    // 3. test that the arrival time of tracer is modelled correctly
    // -----------------------------------------------------------------------------------
    // resetting the model
    model_ptr_->InputPropertyValue( "concentration", makeScalar(ANY,0.) );
    model_ptr_->InputBoundaryValue( LEFT, "concentration", makeScalar(DIRICH,inlet_concentration) );
    double64 phi_min, phi_max;
    model_ptr_->MinMaxOf( "porosity", phi_min, phi_max );
    assert( phi_min == phi_max );
    const double64 porosity(phi_max);
    const double64 expected_arrival_time( model_length / (velo_magnitude/porosity) );
    const double64 threshold_value(inlet_concentration * 0.1);

    // tracer should not be there yet
    transport.AdvectVariable( expected_arrival_time * 0.9 );
    _test( !TestForTracerArrival( "RIGHT", threshold_value ) );
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 5, true );

    // tracer should have arrived at expected_arrival_time because the scheme will be diffusive
    transport.AdvectVariable( expected_arrival_time * 0.1 );
    _test( TestForTracerArrival( "RIGHT", threshold_value ) );
    vtk_output.OutputDataToVTK( *model_ptr_, "concentration", "concentration", 6, true );
 
    return;

 } // end TestFlowThroughModel



/**
    Has the front arrived ?
*/
bool  ExplicitTransport_Test::TestForTracerArrival( const char* boundary, double64 threshold_value ) const
 {
     assert( model_ptr_->ContainsBoundary( string(boundary) ) );
     const Boundary<3U>& boundary_domain(model_ptr_->Boundary(string(boundary)));
 
     double64 var_min, var_max;
     boundary_domain.MinMaxOf( "concentration", var_min, var_max );
 
     // if the average of min-max values is above the threshold value, we detect the tracer arrival
     if ( (var_min+var_max)/2. >= threshold_value ) return true;
 
     return false;
 
 } // end TestForTracerArrival



} // end csmp
