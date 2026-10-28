// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "StokesDiscrepancyMeasure_Example.h"

#include "Region.h"
#include "PDE_Integrator.h"
#include "ANSYS_Model3D.h"
#include "PropertyHandle.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "LinearSolver.h"
#include "CSMP_highLevelUtilities.h"
#include "ExtractVectorVariableLength.h"
#include "VTK_Interface.h"
#include "InputDataManager.h"

using namespace std;

namespace csmp {

void StokesDiscrepancyMeasure_Example::Specifications()
{
    SetTitle( "Stokes discrepancy measure" );
    SetDifficulty( 2 );
    SetCategory( "Numerical Methods" );
    AddAuthor( "SKM" );
    AddDescription( "source in: StokesDiscrepancyMeasure_Example.cpp" );
    AddDescription( "Estimates the Taylor-series truncation error of the Darcy "
                    "approximation relative to Stokes flow by computing a synthetic "
                    "second-derivative (Hessian) term from the parabolic scalar field "
                    "and the FEM pressure solution." );
    AddRequirement( "one_sphere_0.45_tetra (CSMP binary files, -configuration.txt), "
                    "example25.txt (variable file)" );
}



/**
    @file StokesDiscrepancyMeasure_Example.cpp

    @brief Estimates the Taylor-series truncation error arising from the
           finite-element approximation of Stokes flow in pore space.

    @section background Background

    In pore-scale flow modelling, the Darcy velocity field can be approximated
    by solving a Poisson equation with a parabolic (Poiseuille-like) scalar
    field psi as the permeability operator.  The resulting velocity field
    satisfies the Darcy equation exactly but deviates from the true Stokes
    solution by a truncation error that depends on the second derivatives of
    the pressure field.

    This example quantifies that discrepancy by computing the dimensionless
    error measure E defined in Lateef's thesis:

    @code
        E = || 2 (grad psi . nabla) grad p  +  psi grad(Laplacian p) ||
            -------------------------------------------------------
                            || grad p ||
    @endcode

    where:
      - psi   is the parabolic function (solution of Laplacian psi = -1
              with psi = 0 on pore walls), whose maximum value equals the
              local pore radius squared / 8.
      - p     is the fluid pressure (solution of div(psi grad p) = 0).
      - The numerator is the synthetic second-derivative (Hessian) term
              that represents the Taylor-series truncation error of the
              Darcy approximation relative to the full Stokes operator.
      - E = 0 everywhere would mean the Darcy approximation is exact.

    @section algorithm Algorithm

    1. Solve Laplacian(psi) = -1 with psi = 0 on pore walls to obtain
       the parabolic function.
    2. Compute grad(psi) and map psi to element centres.
    3. Solve div(psi grad p) = 0 with pressure boundary conditions to
       obtain the fluid pressure.
    4. Compute the Hessian of p as grad(grad p) via two successive
       gradient operations with cell-to-node extrapolation.
    5. Assemble the truncation error vector:
           term1 = 2 * Hessian(p) * grad(psi)
           term2 = psi * grad(Laplacian(p))
           E_vec = term1 + term2
    6. Normalise by || grad p || to obtain the scalar measure E.

    @section requirements Requirements

    - one_sphere_0.45_tetra  (CSMP binary mesh files + -configuration.txt)
    - example25.txt          (variable definition file)

    @author S.K. Matthai
*/
void StokesDiscrepancyMeasure_Example::Run()
{
    constexpr uint32_t dim{ 3U };

    // -----------------------------------------------------------------------
    // Model setup
    // -----------------------------------------------------------------------
    string model_name( "one_sphere_0.45_tetra" );
    cout << "\nPlease enter the name of the input model, "
         << "or press ENTER to use the default 'one_sphere_0.45_tetra': \n";
    cin.ignore();
    string user_input;
    getline( cin, user_input );
    if ( !user_input.empty() ) model_name = user_input;

    const string file_name     = GetExampleFileName( __FILE__ );
    const string variable_file = "example25.txt";
    const string config_file   = model_name;

    CreateWorkingDirectoryAndCopyInputModelFiles( file_name, model_name,
                                                  variable_file, config_file );

    Model<dim>  model( model_name, variable_file );
    printModelDimensions( model, true );

    InputDataManager<dim>  model_configuration;
    model_configuration.ConfigureFromFile( model, model_name.c_str(),
                                           false, true, true, true, false );

    VTK_Interface<dim>  vtk_output;

    const string  channel_region( "Model" );
    Region<dim>&  gref( model.Region( channel_region.c_str() ) );


    // -----------------------------------------------------------------------
    // 1. Parabolic function psi: solve Laplacian(psi) = -1, psi = 0 on walls
    //
    //    The maximum value of psi equals the local pore radius squared / 8,
    //    making it a proxy for the local hydraulic aperture.
    // -----------------------------------------------------------------------
    PropertyHandle<dim>  src( model, "source term",        SCALAR, ELEMENT ); src = 1.;
    PropertyHandle<dim>  pbf( model, "parabolic function", SCALAR, NODE    );

    printRangeOfVariable( model, "source term" );

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_nxtyp(0);
    settings.Set_ncgtyp(5);
    settings.Set_ndefault(40);
    SAMG_Solver                  solver( &settings );
#else
    CSMP_DEFAULT_LINEAR_SOLVER   solver;
#endif
    PDE_Integrator<dim,Element>  parabolic_profile( solver );

    NumIntegral_dNT_dN_dV<dim>      laplacian( model.Database(),
                                               "parabolic function",
                                               "parabolic function" );
    NumIntegral_NT_rhsop_N_dV<dim>  rhs( model.Database(),
                                         "source term",
                                         "parabolic function" );
    parabolic_profile.Add( &laplacian );
    parabolic_profile.Add( &rhs );

    // Dirichlet: psi = 0 on pore walls (all perimeter nodes)
    gref.InputPropertyValue( "parabolic function", makeScalar(DIRICH, 0.), PERIMETER );

    // Release Dirichlet at inlet (FRONT) and outlet (BACK) so that psi
    // is free to take non-zero values there — we cannot use the standard
    // boundary functions here because we must not release the condition
    // at the edges of those faces, only at the interior face nodes.
    for ( auto nit = gref.PerimeterNodesBegin(); nit != gref.NodesEnd(); ++nit )
        if ( (*nit)->AtBoundary() == FRONT || (*nit)->AtBoundary() == BACK )
            (*nit)->Status( pbf.Key(), PLAIN );

    parabolic_profile.IntegrateOver( model, gref );
    printRangeOfVariable( model, channel_region.c_str(), "parabolic function" );
    vtk_output.OutputDataToVTK( model, "parabolic-function", "parabolic function", 1 );


    // -----------------------------------------------------------------------
    // 2. Gradient of the parabolic function
    //
    //    psi is mapped to element barycentres for use as the permeability
    //    operator in the pressure solve (step 3).
    // -----------------------------------------------------------------------
    PropertyHandle<dim>  lapbc( model, "element parabolic function", SCALAR, ELEMENT );

    {
        ScalarVariable sc;
        for ( auto& eit : gref.CellVector() ) {
            eit->PropertyValueAtBaryCenter( pbf.Key(), sc );
            eit->Store( lapbc.Key(), sc );
        }
    }
    printRangeOfVariable( model, channel_region.c_str(), "element parabolic function" );

    PropertyHandle<dim>  grad( model, "grad parabolic function", VECTOR, ELEMENT );
    model.CopyGradientOfProperty_A_To_B( "parabolic function",
                                         "grad parabolic function" );
    printRangeOfVariable( model, channel_region.c_str(), "grad parabolic function" );

    PropertyHandle<dim>  magn( model, "magnitude grad parabolic function",
                               SCALAR, ELEMENT );
    ExtractVectorVariableLength<dim>  grad_magnitude( model.Database(),
                                                      "grad parabolic function",
                                                      "magnitude grad parabolic function" );
    gref.Apply( grad_magnitude );
    printRangeOfVariable( model, channel_region.c_str(),
                          "magnitude grad parabolic function" );


    // -----------------------------------------------------------------------
    // 3. Fluid pressure: solve div(psi grad p) = 0
    //
    //    psi acts as the spatially varying permeability (Poiseuille weight).
    // -----------------------------------------------------------------------
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver                  samg_solver;
#else
    CSMP_DEFAULT_LINEAR_SOLVER   linear_solver;
#endif

#ifdef CSMP_WITH_SAMG_SOLVER
    PDE_Integrator<dim,Element>  steady_state_pressure( samg_solver );
#else
    PDE_Integrator<dim,Element>  steady_state_pressure( linear_solver );
#endif

    NumIntegral_dNT_lhsop_dN_dV<dim>  flow_resistance( model.Database(),
                                                        "element parabolic function",
                                                        "fluid pressure",
                                                        "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<dim>    fsrc( model.Database(),
                                            "fluid volume source",
                                            "fluid pressure" );
    steady_state_pressure.Add( &flow_resistance );
    steady_state_pressure.Add( &fsrc );

    model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN, 0.) );
    steady_state_pressure.IntegrateOver( model, gref );
    printRangeOfVariable( model, channel_region.c_str(), "fluid pressure" );
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );


    // -----------------------------------------------------------------------
    // 4. Hessian of fluid pressure: H(p) = grad(grad p)
    //
    //    Two successive gradient operations with cell-to-node extrapolation
    //    yield the second-derivative (Hessian) tensor at element centres.
    //    This is the key quantity for the truncation error estimate.
    // -----------------------------------------------------------------------
    model.CopyGradientOfProperty_A_To_B( "fluid pressure",
                                         "fluid pressure gradient" );
    model.ExtrapolateCellToNodeProperty( "fluid pressure gradient",
                                         "nodal fluid pressure gradient" );
    model.CopyGradientOfProperty_A_To_B( "nodal fluid pressure gradient",
                                         "fluid pressure gradient2" );
    model.ExtrapolateCellToNodeProperty( "fluid pressure gradient2",
                                         "nodal fluid pressure gradient2" );

    vtk_output.OutputDataToVTK( model, "fluid-pressure-gradient",
                                "fluid pressure gradient", 1 );
    vtk_output.OutputDataToVTK( model, "hessian",
                                "fluid pressure gradient2", 1 );


    // -----------------------------------------------------------------------
    // 5. Taylor-series truncation error measure E
    //
    //    Following Lateef's thesis, the discrepancy between the Darcy
    //    approximation and the full Stokes solution is:
    //
    //        E = || term1 + term2 || / || grad p ||
    //
    //    where:
    //        term1 = 2 * H(p) * grad(psi)          [Hessian contracted with grad psi]
    //        term2 = psi * grad(Laplacian(p))       [psi times gradient of pressure Laplacian]
    //
    //    E = 0 everywhere would indicate that the Darcy approximation is exact.
    //    Large values of E identify regions where higher-order Stokes terms
    //    are significant and the Darcy model introduces truncation error.
    // -----------------------------------------------------------------------

    PropertyHandle<dim>  hessian( model, "fluid pressure gradient2", TENSOR, ELEMENT );
    PropertyHandle<dim>  term1(   model, "term1",                    VECTOR, ELEMENT );

    // --- term1 = 2 * H(p) * grad(psi) ---
    {
        TensorVariable<dim>  ts;
        VectorVariable<dim>  vc;
        for ( auto& eit : gref.CellVector() ) {
            eit->Read( hessian.Key(), ts );
            eit->Read( grad.Key(),    vc );
            eit->Store( term1.Key(), ts * vc * 2.0 );
        }
    }

    // --- Laplacian(p) = trace of H(p) ---
    PropertyHandle<dim>  lap_p( model, "laplacian of pressure", SCALAR, ELEMENT );
    {
        TensorVariable<dim>  ts;
        ScalarVariable       sc;
        for ( auto& eit : gref.CellVector() ) {
            eit->Read( hessian.Key(), ts );
            // Laplacian = sum of diagonal (trace of Hessian)
            double trace = 0.;
            for ( uint32_t i{0U}; i < dim; ++i ) trace += ts(i,i);
            sc = makeScalar( PLAIN, trace );
            eit->Store( lap_p.Key(), sc );
        }
    }
    printRangeOfVariable( model, channel_region.c_str(), "laplacian of pressure" );
    vtk_output.OutputDataToVTK( model, "laplacian-pressure",
                                "laplacian of pressure", 1 );

    // --- term2 = psi * grad(Laplacian(p)) ---
    PropertyHandle<dim>  nlap_pf( model, "nodal laplacian of pressure", SCALAR, NODE );
    model.ExtrapolateCellToNodeProperty( "laplacian of pressure",
                                         "nodal laplacian of pressure" );

    PropertyHandle<dim>  term2( model, "gradient of laplacian of pressure",
                                VECTOR, ELEMENT );
    model.CopyGradientOfProperty_A_To_B( "nodal laplacian of pressure",
                                         "gradient of laplacian of pressure" );
    printRangeOfVariable( model, channel_region.c_str(),
                          "gradient of laplacian of pressure" );
    vtk_output.OutputDataToVTK( model, "grad-laplacian-pressure",
                                "gradient of laplacian of pressure", 1 );

    {
        ScalarVariable       sc;
        VectorVariable<dim>  vc;
        for ( auto& eit : gref.CellVector() ) {
            eit->Read( lapbc.Key(), sc );   // psi at element centre
            eit->Read( term2.Key(), vc );   // grad(Laplacian(p))
            vc *= sc;
            // term2 now holds psi * grad(Laplacian(p)) — overwrites intermediate gradient
            eit->Store( term2.Key(), vc );
        }
    }

    // --- Numerator: term1 + term2 ---
    term1 += term2;

    // --- Normalise by || grad p || to obtain scalar measure E ---
    PropertyHandle<dim>  gradp(  model, "fluid pressure gradient",   VECTOR, ELEMENT );
    PropertyHandle<dim>  E_vec(  model, "Stokes Discrepancy vector", VECTOR, ELEMENT );
    PropertyHandle<dim>  E_mag(  model, "Stokes Discrepancy Measure",SCALAR, ELEMENT );

    {
        ScalarVariable       sc, em;
        VectorVariable<dim>  vc;
        for ( auto& eit : gref.CellVector() ) {
            eit->Read( gradp.Key(), vc );
            const double grad_p_length = vc.Length();

            eit->Read( term1.Key(), vc );
            eit->Store( E_vec.Key(), vc );

            // Guard against division by zero at stagnation points
            if ( grad_p_length > numeric_limits<double>::epsilon() )
                em = makeScalar( PLAIN, vc.Length() / grad_p_length );
            else
                em = makeScalar( PLAIN, 0. );

            eit->Store( E_mag.Key(), em );
        }
    }

    printRangeOfVariable( model, channel_region.c_str(),
                          "Stokes Discrepancy Measure" );
    vtk_output.OutputDataToVTK( model, "Stokes-discrepancy-vector",
                                "Stokes Discrepancy vector", 1 );
    vtk_output.OutputDataToVTK( model, "Stokes-discrepancy-measure",
                                "Stokes Discrepancy Measure", 1 );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nStokesDiscrepancyMeasure_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp

