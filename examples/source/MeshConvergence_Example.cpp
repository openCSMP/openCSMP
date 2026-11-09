// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only


#include "MeshConvergence_Example.h"
#include "exampleHelpers.h"

#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"
#include "meshManagementUtilities.h"

// Mesh construction
#include "VSet.h"
#include "VSetConverter.h"
#include "vsetMakers.h"          // create_2_Triangle_VSet()
#include "ModelTopology.h"

// PDE operators
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "LinearSolver.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#endif

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

using namespace std;

namespace csmp {

// ============================================================================
// Specifications
// ============================================================================

void MeshConvergence_Example::Specifications()
{
    SetTitle( "Mesh convergence study (interactive)" );
    SetDifficulty( 2 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: MeshConvergence_Example.cpp" );
    AddDescription( "Demonstrates FEM mesh convergence using the method of "
                    "manufactured solutions. Solves the Poisson equation on "
                    "a unit square with a known exact solution, measuring the "
                    "L2 error norm at each refinement level for linear and "
                    "quadratic triangles." );
    AddRequirement( "example7-variables.txt (variable definition file)" );
}


// ============================================================================
// Local helpers
// ============================================================================

namespace {

/// Exact solution: p(x,y) = sin(pi*x) * sin(pi*y)
inline double p_exact( double x, double y )
{
    return sin(M_PI*x) * sin(M_PI*y);
}

/// Source term: f = -div(grad p_exact) = 2*pi^2 * sin(pi*x) * sin(pi*y)
inline double source_term( double x, double y )
{
    return 2. * M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y);
}


/**
    Assigns the manufactured source term f(x,y) to all elements
    using their barycentre coordinates.
*/
void assign_source_term( Model<2U>& model )
{
    const csmp::Index f_key =
        model.Database().StorageKey("fluid volume source");

    // Assign source term at element barycentres
    for ( auto& eit : model.Region("Model").CellVector() ) {
        auto bctr = eit->BaryCenter();
        const double f = source_term( bctr[0], bctr[1] );
        eit->Store( f_key, makeScalar(PLAIN, f) );
    }
}


/**
    Applies Dirichlet BCs: p = p_exact on all four boundaries.
    For the manufactured solution with p_exact = sin(pi*x)*sin(pi*y),
    p = 0 on all sides of the unit square.
*/
void apply_boundary_conditions( Model<2U>& model )
{
    // p_exact = sin(pi*x)*sin(pi*y) = 0 on all sides of the unit square.
    // Cannot use CSMP Boundary objects because there are none in the original model.
    // (no lower-dimensional Face elements in the mesh).
    // Instead, flag all perimeter nodes as Dirichlet directly.

    const csmp::Index p_key = model.Database().StorageKey("fluid pressure");
    Region<2U>& gref( model.Region("Model") );

    // All perimeter nodes get p = 0 (Dirichlet)
    // p_exact = 0 on all four sides of the unit square
    for ( auto nit = gref.PerimeterNodesBegin(); nit != gref.NodesEnd(); ++nit )
        (*nit)->Store( p_key, makeScalar(DIRICH, 0.) );
}


/**
    Computes the L2 error norm:
        ||p_h - p_exact||_L2 = sqrt( sum_nodes (p_h - p_exact)^2 / N )

    Also stores the error at each node in 'pressure error' variable
    for VTK visualisation.

    @return L2 error norm
*/
double compute_L2_error( Model<2U>& model )
{
    const csmp::Index p_key   = model.Database().StorageKey("fluid pressure");
    const csmp::Index err_key = model.Database().StorageKey("pressure error");

    double sum_sq{ 0. };
    size_t n_nodes{ 0U };

    for ( const auto& nit : model.Region("Model").NodeVector() ) {
        const double p_h   = nit->Read( p_key );
        const double p_ref = p_exact( nit->x(), nit->y() );
        const double err   = p_h - p_ref;

        // Store pointwise error for visualisation
        nit->Store( err_key, makeScalar(PLAIN, err) );

        sum_sq  += err * err;
        ++n_nodes;
    }

    return sqrt( sum_sq / static_cast<double>(n_nodes) );
}


/**
    Solves the Poisson equation on the given model
    with homogeneous boundary conditions.
    Returns the L2 error norm.
*/
double solve_and_measure( Model<2U>& model )
{
    // Assign coordinates, source term, and BCs
    assign_source_term( model );
    apply_boundary_conditions( model );

    // Build and run the PDE integrator
    EigenSolver                 solver; // needs direct solver with flexible convergence criteria
    PDE_Integrator<2U,Element>  poisson( solver );

    // conductivity = 1 (Laplacian)
    model.InputPropertyValue( "conductivity", makeScalar(PLAIN, 1.) );

    NumIntegral_dNT_dN_dV<2U>      lhs( model.Database(),
                                        "fluid pressure", "fluid pressure" );
    NumIntegral_NT_rhsop_N_dV<2U>  rhs( model.Database(),
                                        "fluid volume source",
                                        "fluid pressure" );
    poisson.Add( &lhs );
    poisson.Add( &rhs );

    poisson.IntegrateOver( model, model.Region("Model") );

    return compute_L2_error( model );
}


/**
    Prints one row of the convergence table.
*/
void print_convergence_row( size_t level,
                            size_t nodes,
                            size_t elements,
                            double h,
                            double L2_error,
                            double rate )
{
    cout << "  " << left
         << setw(8)  << level
         << setw(10) << nodes
         << setw(12) << elements
         << setw(12) << dbl_str(h,5)
         << setw(16) << sci_str(L2_error,4);
    if ( rate > 0. )
        cout << setw(10) << dbl_str(rate,3);
    else
        cout << setw(10) << "---";
    cout << "\n";
}


/// Stores one row of convergence data for deferred printing
struct ConvergenceResult {
    size_t level;
    size_t nodes;
    size_t elements;
    double h;
    double L2;
    double rate;
};

} // anonymous namespace


// ============================================================================
// Run
// ============================================================================

/**
    @file MeshConvergence_Example.cpp

    @brief CSMP for Beginners -- Example 7: Mesh Convergence Study

    @section purpose Purpose

    This interactive example demonstrates FEM mesh convergence using
    the Method of Manufactured Solutions (MMS). We solve the Poisson
    equation on a unit square with a known exact solution, then measure
    how the error decreases as the mesh is refined.

    @section physics Physics

    We solve:
        -div( grad p ) = f    in [0,1] x [0,1]

    with the manufactured exact solution:
        p_exact(x,y) = sin(pi*x) * sin(pi*y)

    The source term is derived by substitution:
        f(x,y) = 2*pi^2 * sin(pi*x) * sin(pi*y)

    Boundary conditions:
        p = 0 on all four sides (Dirichlet)

    @section convergence Expected Convergence Rates

    For the L2 error norm ||p_h - p_exact||:
        Linear triangles  (P1): O(h^2)  -- rate ~2
        Quadratic triangles (P2): O(h^4) -- rate ~4 (or O(h^3) in H1)

    @section method Method of Manufactured Solutions

    The MMS works as follows:
      1. Choose any smooth function p_exact(x,y)
      2. Compute f = -div(grad p_exact) analytically
      3. Solve the PDE with this f as the source term
      4. Compare the numerical solution with p_exact
      5. The error should converge at the expected rate

    This is the standard way to verify FEM implementations.

    @section meshes Mesh Refinement Strategy

    Starting from a 2-triangle unit square (create_2_Triangle_VSet),
    each refinement level applies refineSimplexMesh() which:
      - splits each triangle into 4 sub-triangles
      - adds midside nodes on each edge
      - halves the element size h

    After convergence with linear triangles, the mesh is converted
    to quadratic triangles using the VSetConverter's method
    ConvertLinearToQuadraticTriangles()
    and the convergence study is repeated.

    @section requirements Requirements

    - example7-variables.txt  (variable definition file)
    - No mesh file needed -- mesh built programmatically

    @author S.K. Matthai
*/
void MeshConvergence_Example::Run()
{
    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 7\n"
         << "  Mesh Convergence Study\n\n"
         << "  We solve the Poisson equation with a known exact solution\n"
         << "  and measure how the error decreases as the mesh is refined.\n\n"
         << "  Exact solution: p(x,y) = sin(pi*x) * sin(pi*y)\n"
         << "  Source term:    f(x,y) = 2*pi^2 * sin(pi*x) * sin(pi*y)\n"
         << "  Domain:         unit square [0,1] x [0,1]\n"
         << "  BCs:            p = 0 on all four sides\n\n"
         << "  Steps:\n"
         << "    1.  The method of manufactured solutions\n"
         << "    2.  The unit square mesh and refinement strategy\n"
         << "    3.  Expected convergence rates\n"
         << "    4.  Linear triangle convergence study\n"
         << "    5.  Quadratic triangle convergence study\n"
         << "    6.  Comparison and interpretation\n";
    divider();
    pause( "  [ Press ENTER to begin ]" );

    const string model_name    = "";
    const string variable_file = "CSMP-1phase-variables.txt";
    // get "MeshConvergence_Example"
    const string file_name     = GetExampleFileName( __FILE__ );
    CreateWorkingDirectoryAndCopyInputModelFiles( file_name, model_name, variable_file );


    // ========================================================================
    // Step 1: The method of manufactured solutions
    // ========================================================================
    heading( 1, 6, "The method of manufactured solutions" );

    explain(
"The Method of Manufactured Solutions (MMS) is the standard way\n"
"to verify a FEM implementation.\n"
"\n"
"The idea:\n"
"  1. Choose any smooth function p_exact(x,y)\n"
"  2. Compute f = -div(grad p_exact) analytically\n"
"  3. Solve the PDE with this f as the source term\n"
"  4. Compare the numerical solution with p_exact\n"
"  5. The error should converge at the expected rate\n"
"\n"
"We choose:\n"
"    p_exact(x,y) = sin(pi*x) * sin(pi*y)\n"
"\n"
"Then:\n"
"    -div(grad p_exact) = 2*pi^2 * sin(pi*x) * sin(pi*y) = f\n"
"\n"
"On the unit square, p_exact = 0 on all four sides,\n"
"so the Dirichlet BC is simply p = 0 everywhere on the boundary."
    );
    pause();

    subheading( "Why this exact solution?" );
    explain(
"sin(pi*x)*sin(pi*y) is ideal because:\n"
"  - It is smooth (infinitely differentiable)\n"
"  - It is zero on all sides of the unit square\n"
"  - Its Laplacian is proportional to itself (eigenfunction)\n"
"  - It has a single interior maximum at (0.5, 0.5)\n"
"  - It is the fundamental mode of the Laplacian on [0,1]^2\n"
"\n"
"For smooth solutions, FEM theory predicts:\n"
"  Linear triangles (P1):    L2 error = O(h^2)\n"
"  Quadratic triangles (P2): L2 error = O(h^4)\n"
"\n"
"where h is the characteristic element size."
    );
    pause();

    subheading( "Live verification of the source term" );
    {
        // Print exact solution and source term at a few points
        cout << "  p_exact and source term at selected points:\n\n";
        cout << "  " << left
             << setw(8) << "x"
             << setw(8) << "y"
             << setw(18) << "p_exact(x,y)"
             << setw(18) << "f(x,y)" << "\n";
        divider('-');

        const double pts[][2] = {
            {0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0},
            {0.0, 0.5}, {0.5, 0.5}, {1.0, 0.5},
            {0.0, 1.0}, {0.5, 1.0}, {1.0, 1.0}
        };
        for ( const auto& pt : pts ) {
            const double x = pt[0], y = pt[1];
            cout << "  " << left
                 << setw(8) << dbl_str(x,1)
                 << setw(8) << dbl_str(y,1)
                 << setw(18) << dbl_str(p_exact(x,y),6)
                 << setw(18) << dbl_str(source_term(x,y),6) << "\n";
        }
        divider('-');
        output( "p_exact = 0 on all boundaries (x=0, x=1, y=0, y=1) ✓" );
        output( "Maximum p_exact = " + dbl_str(p_exact(0.5,0.5),6)
                + " at (0.5, 0.5)" );
        output( "f = 2*pi^2 * p_exact  (eigenfunction of Laplacian)" );
    }
    pause();

    keypoint(
"MMS: choose p_exact, compute f = -div(grad p_exact), solve, compare.\n"
"The error must converge at the rate predicted by FEM theory.\n"
"If it does not, there is a bug in the implementation.\n"
"sin(pi*x)*sin(pi*y) is the ideal test function for the unit square."
    );
    pause();


    // ========================================================================
    // Step 2: The unit square mesh and refinement strategy
    // ========================================================================
    heading( 2, 6, "The unit square mesh and refinement strategy" );

    explain(
"We start from the simplest possible mesh: two triangles on the\n"
"unit square (create_2_Triangle_VSet).\n"
"\n"
"    3-----------2\n"
"    |  \\  T1   |\n"
"    |    \\     |\n"
"    |  T0  \\   |\n"
"    0-----------1\n"
"\n"
"Each refinement level applies refineSimplexMesh() which:\n"
"  - Splits each triangle into 4 sub-triangles\n"
"  - Adds a midside node on each edge\n"
"  - Halves the element size h\n"
"\n"
"After N refinements:\n"
"  Elements = 2 * 4^N\n"
"  h ~ 1 / sqrt(2 * 4^N) = 1 / (sqrt(2) * 2^N)"
    );
    pause();

    subheading( "refineSimplexMesh() -- how it works" );
    code(
"// Starting mesh: 2 triangles, 4 nodes\n"
"VSet<2U>  vset;\n"
"create_2_Triangle_VSet( vset );\n"
"\n"
"// Refinement level 1: 8 triangles, 9 nodes\n"
"// parent_cell_indices[new_idx] = old_parent_idx\n"
"vector<size_t> parent_cell_indices = vset.Refine();\n"
"\n"
"// Refinement level 2: 32 triangles, 25 nodes\n"
"parent_cell_indices = vset.Refine();\n"
"\n"
"// Each call halves h and quadruples the element count.\n"
"// Rebuild ModelTopology using parent_cell_indices."
    );
    pause();

    subheading( "Mesh sizes at each refinement level" );
    {
        cout << "  " << left
             << setw(10) << "Level"
             << setw(12) << "Elements"
             << setw(10) << "Nodes"
             << setw(12) << "h (approx)" << "\n";
        divider('-');

        size_t elements = 2;
        size_t nodes    = 4;
        for ( size_t level = 0; level <= 6; ++level ) {
            const double h = 1. / sqrt( static_cast<double>(elements) / 2. );
            cout << "  " << left
                 << setw(10) << level
                 << setw(12) << elements
                 << setw(10) << nodes
                 << setw(12) << dbl_str(h,5) << "\n";
            elements *= 4;
            // nodes: (2^level + 1)^2 for a structured grid
            const size_t n = static_cast<size_t>(pow(2.,level+1)) + 1;
            nodes = n * n;
        }
        divider('-');
    }
    pause();

keypoint(
"vset.Refine() splits each triangle into 4 sub-triangles.\n"
"Each refinement halves h and quadruples the element count.\n"
"The ModelTopology must be rebuilt after each refinement\n"
"using parent_cell_indices returned by vset.Refine():\n"
"  parent_cell_indices[new_idx] = old_parent_idx"
);
    pause();


    // ========================================================================
    // Step 3: Expected convergence rates
    // ========================================================================
    heading( 3, 6, "Expected convergence rates" );

    explain(
"FEM theory predicts the following convergence rates for the\n"
"L2 error norm ||p_h - p_exact||_L2:\n"
"\n"
"  Linear triangles (P1):\n"
"    L2 error = C * h^2\n"
"    Halving h reduces the error by a factor of 4\n"
"    Convergence rate = 2\n"
"\n"
"  Quadratic triangles (P2):\n"
"    L2 error = C * h^4  (superconvergence for smooth solutions)\n"
"    Halving h reduces the error by a factor of 16\n"
"    Convergence rate = 4\n"
"\n"
"The convergence rate is measured as:\n"
"    rate = log(error_i / error_{i+1}) / log(h_i / h_{i+1})\n"
"         = log(error_i / error_{i+1}) / log(2)   (since h halves)\n"
"\n"
"A rate close to the theoretical value confirms the implementation\n"
"is correct and the mesh is in the asymptotic regime."
    );
    pause();

    subheading( "Log-log convergence plot" );
    explain(
"On a log-log plot of error vs h:\n"
"  - A straight line confirms power-law convergence\n"
"  - The slope of the line = convergence rate\n"
"  - P1 slope = 2,  P2 slope = 4\n"
"\n"
"In ParaView, plot the 'pressure error' variable to see where\n"
"the error is largest (near the maximum of sin(pi*x)*sin(pi*y))."
    );
    pause();

    keypoint(
"P1 triangles: rate ~2  (error halves when h halves, then halves again)\n"
"P2 triangles: rate ~4  (error reduces 16x when h halves)\n"
"Rate < expected: implementation error or mesh too coarse\n"
"Rate > expected: superconvergence (possible for special meshes)"
    );
    pause();


        // ====================================================================
        // Step 4: Linear triangle convergence study
        // ====================================================================
        heading( 4, 6, "Linear triangle convergence study" );

        explain(
"We now run the convergence study for linear (P1) triangles.\n"
"Starting from the 2-triangle mesh, we refine several times and\n"
"measure the L2 error at each level.\n"
"\n"
"Expected result: convergence rate ~2 for all levels."
        );
        pause();

        // variables for reuse in Step 5
        size_t n_levels{ 6 };

        subheading( "How many refinement levels?" );
        {
            cout << "\n  Enter number of refinement levels (4-7, recommended 6): ";
            cout.flush();
            cin >> n_levels;
            n_levels = max( size_t(4), min( size_t(7), n_levels ) );
            output( "Running " + to_string(n_levels) + " refinement levels." );
            pause();

            // ------------------------------------------------------------------
            // Build initial 2-triangle mesh
            // ------------------------------------------------------------------
            output( "Building initial 2-triangle mesh..." );

            VSet<2U>       vset_current;
            create_2_Triangle_VSet( vset_current );

            // Build initial topology: one region containing both elements
            ModelTopology  topo_current( true );  // isoparametric
            {
                set<string>    fem_types;
                fem_types.insert( "ISOPARAMETRIC_LINEAR_TRIANGLE" );
                vector<size_t>  elms;
                for ( size_t i{0}; i < vset_current.Elements(); ++i )
                    elms.push_back(i);
                topo_current.AddDomain( "Model", fem_types, elms );
            }

            output( "Initial mesh: " + to_string(vset_current.Elements())
                    + " elements, " + to_string(vset_current.Vertices())
                    + " nodes." );

            // ------------------------------------------------------------------
            // Convergence table header
            // ------------------------------------------------------------------
            blank();
            divider('-');
            cout << "  LINEAR TRIANGLE (P1) CONVERGENCE\n";
            divider('-');
            cout << "  " << left
                 << setw(8)  << "Level"
                 << setw(10) << "Nodes"
                 << setw(12) << "Elements"
                 << setw(12) << "h"
                 << setw(16) << "L2 error"
                 << setw(10) << "Rate" << "\n";
            divider('-');

            vector<double>  h_values;
            vector<double>  L2_errors;
            VTK_Interface<2U>  vtk;

            // ------------------------------------------------------------------
            // Start at level 1: the 2-triangle mesh (level 0) has only 4 corner
            // nodes which all become Dirichlet constraints, leaving zero free DOFs.
            // One refinement gives 8 triangles and 9 nodes (5 interior) which is
            // the coarsest mesh with meaningful degrees of freedom.
            // ------------------------------------------------------------------

            // Perform the first refinement before the loop
            {
                output( "Pre-refining to level 1 (coarsest usable mesh)..." );
                cout.flush();

                const vector<size_t>  parent_cell_indices =
                    vset_current.Refine();

                ModelTopology  topo_refined( true );
                topo_refined.ModelName( topo_current.ModelName().c_str() );

                for ( auto sdom  = topo_current.DomainsBegin();
                           sdom != topo_current.DomainsEnd(); ++sdom )
                {
                    vector<size_t>  new_cell_ids;
                    new_cell_ids.reserve( (*sdom).second.second.size() * 4 );

                    for ( const size_t old_cid : (*sdom).second.second )
                        for ( size_t new_idx{0};
                              new_idx < parent_cell_indices.size(); ++new_idx )
                            if ( parent_cell_indices[new_idx] == old_cid )
                                new_cell_ids.push_back( new_idx );

                    assert( (*sdom).second.first.size() == 1 );
                    topo_refined.AddDomain( (*sdom).first.c_str(),
                                            (*sdom).second.first,
                                            new_cell_ids );
                }

                topo_refined.CheckCellNumbering();
                topo_current = topo_refined;

                output( "Level 1: " + to_string(vset_current.Elements())
                        + " elements, " + to_string(vset_current.Vertices())
                        + " nodes." );
            }

            vector<ConvergenceResult>  p1_results;
            p1_results.reserve( n_levels );

            for ( size_t level = 1; level <= n_levels; ++level )
            {
                // ----------------------------------------------------------
                // Build model
                // ----------------------------------------------------------
                VSet<2U>       vset_model( vset_current );
                ModelTopology  topo_model( topo_current );

                VSetConverter().ConvertElementTypesToOnesUsingLocalCoordinateSystem( vset_model );

                Model<2U>  model( topo_model, vset_model, variable_file.c_str(), true );
                model.CreateProperty( "pressure error", "perr", "Pa" );

                // Solve and measure
                const double L2 = solve_and_measure( model );
                const double h  = 1. / sqrt( static_cast<double>(
                                     model.Mesh().Elements()) / 2. );

                h_values.push_back( h );
                L2_errors.push_back( L2 );

                const double rate = ( level > 1 )
                    ? log( L2_errors[L2_errors.size()-2] / L2_errors.back() )
                    / log( h_values[h_values.size()-2]   / h_values.back()  )
                    : -1.;

                // Store result for deferred printing
                p1_results.push_back({ level,
                                       model.Mesh().Nodes(),
                                       model.Mesh().Elements(),
                                       h, L2, rate });

                // VTK output for first and last level
                if ( level == 1 || level == n_levels ) {
                    vtk.OutputDataToVTK( model, "p1_pressure",
                                         "fluid pressure",
                                         static_cast<long>(level) );
                    vtk.OutputDataToVTK( model, "p1_error",
                                         "pressure error",
                                         static_cast<long>(level) );
                    vtk.OutputDataToVTK( model, "p1_source",
                                         "fluid volume source",
                                         static_cast<long>(level) );
                }

                // Refine for next level
                if ( level < n_levels )
                {
                    const vector<size_t>  parent_cell_indices =
                        vset_current.Refine();

                    ModelTopology  topo_refined( true );
                    topo_refined.ModelName( topo_current.ModelName().c_str() );

                    for ( auto sdom  = topo_current.DomainsBegin();
                               sdom != topo_current.DomainsEnd(); ++sdom )
                    {
                        vector<size_t>  new_cell_ids;
                        new_cell_ids.reserve(
                            (*sdom).second.second.size() * 4 );

                        for ( const size_t old_cid : (*sdom).second.second )
                            for ( size_t new_idx{0};
                                  new_idx < parent_cell_indices.size();
                                  ++new_idx )
                                if ( parent_cell_indices[new_idx] == old_cid )
                                    new_cell_ids.push_back( new_idx );

                        assert( (*sdom).second.first.size() == 1 );
                        topo_refined.AddDomain( (*sdom).first.c_str(),
                                                (*sdom).second.first,
                                                new_cell_ids );
                    }

                    topo_refined.CheckCellNumbering();
                    topo_current = topo_refined;
                }

            } // end for level

            // ------------------------------------------------------------------
            // Print the complete convergence table AFTER all solves are done
            // ------------------------------------------------------------------
            blank();
            divider('-');
            cout << "  LINEAR TRIANGLE (P1) CONVERGENCE\n";
            divider('-');
            cout << "  " << left
                 << setw(8)  << "Level"
                 << setw(10) << "Nodes"
                 << setw(12) << "Elements"
                 << setw(12) << "h"
                 << setw(16) << "L2 error"
                 << setw(10) << "Rate" << "\n";
            divider('-');

            for ( const auto& r : p1_results )
                print_convergence_row( r.level, r.nodes, r.elements,
                                       r.h, r.L2, r.rate );

            divider('-');

            // Check convergence rate
            if ( p1_results.size() >= 2 ) {
                const double last_rate = p1_results.back().rate;
                blank();
                output( "Final P1 convergence rate: "
                        + dbl_str(last_rate,3) );
                if ( last_rate >= 1.8 && last_rate <= 2.2 )
                    output( "PASS: rate ~2 as expected for P1 triangles." );
                else
                    output( "WARNING: rate deviates from expected value of 2." );
            }

            // VTK summary
            blank();
            output( "VTK files written:" );
            output( "  p1_pressure_1.vtk  (coarsest mesh)" );
            output( "  p1_pressure_" + to_string(n_levels)
                    + ".vtk  (finest mesh)" );
            output( "  p1_error_1.vtk     (pointwise error, coarsest)" );
            output( "  p1_error_" + to_string(n_levels)
                    + ".vtk     (pointwise error, finest)" );
            blank();

            // Check convergence rate from last two levels
            if ( L2_errors.size() >= 3 ) {
                const double last_rate =
                    log( L2_errors[L2_errors.size()-2] / L2_errors.back() )
                  / log( h_values[h_values.size()-2]   / h_values.back()  );
                output( "Final P1 convergence rate: "
                        + dbl_str(last_rate,3) );
                if ( last_rate >= 1.8 && last_rate <= 2.2 )
                    output( "PASS: rate ~2 as expected for P1 triangles." );
                else
                    output( "WARNING: rate deviates from expected value of 2." );
            }

        } // end n_levels scope
        pause();

        keypoint(
"P1 triangles converge at rate ~2 for smooth solutions.\n"
"vset.Refine() returns parent_cell_indices[new_idx] = old_parent_idx.\n"
"The topology is rebuilt by finding all new cells per old domain.\n"
"topo.CheckCellNumbering() verifies the rebuilt topology is consistent.\n"
"VTK output shows the error peaks at (0.5, 0.5) -- the solution maximum."
        );
        pause();


// ====================================================================
// Step 5: Quadratic triangle convergence study
// ====================================================================
heading( 5, 6, "Quadratic triangle convergence study" );

        explain(
"We now repeat the convergence study with quadratic (P2) triangles.\n"
"\n"
"Starting from the same refined linear mesh, we convert each\n"
"mesh to quadratic triangles using VSetConverter's:\n"
"    ConvertLinearToQuadraticTriangles( vset )\n"
"\n"
"This adds midside nodes to each triangle edge, converting\n"
"ISOPARAMETRIC_LINEAR_TRIANGLE to ISOPARAMETRIC_QUADRATIC_TRIANGLE.\n"
"\n"
"Expected result: convergence rate ~4 for all levels.\n"
"\n"
"Note: we use fewer refinement levels for P2 because the error\n"
"decreases much faster -- 4 levels is usually sufficient."
        );
        pause();

        subheading( "ConvertLinearToQuadraticTriangles" );
        code(
"// Convert linear triangles to quadratic using VSetConverter's:\n"
"ConvertLinearToQuadraticTriangles( vset );\n"
"\n"
"// This:\n"
"//   - adds a midside node on each edge\n"
"//   - changes element type from LINEAR_TRIANGLE to QUADRATIC_TRIANGLE\n"
"//   - doubles the number of nodes approximately\n"
"//   - does NOT change the number of elements\n"
"\n"
"// The quadratic elements use the same mesh topology\n"
"// but with higher-order interpolation."
        );
        pause();

subheading( "Running P2 convergence study..." );
{
    const size_t n_levels_p2 = min( n_levels, size_t(5) );

    // Rebuild from scratch for P2
    VSet<2U>       vset_p2_current;
    create_2_Triangle_VSet( vset_p2_current );

    ModelTopology  topo_p2_current( true );
    {
        set<string>    fem_types;
        fem_types.insert( "ISOPARAMETRIC_LINEAR_TRIANGLE" );
        vector<size_t>  elms;
        for ( size_t i{0}; i < vset_p2_current.Elements(); ++i )
            elms.push_back(i);
        topo_p2_current.AddDomain( "Model", fem_types, elms );
    }

    vector<double>  h_p2;
    vector<double>  L2_p2;
    VTK_Interface<2U>  vtk;

    vector<ConvergenceResult>  p2_results;
    p2_results.reserve( n_levels_p2 );

    for ( size_t level = 0; level <= n_levels_p2; ++level )
      {
        // Step 1: convert VSet to quadratic
        // This adds midside nodes AND updates element types in the VSet.
        VSet<2U>  vset_q( vset_p2_current );
        VSetConverter().ConvertLinearToQuadraticTriangles( vset_q );

        // Step 2: build quadratic topology from the linear topology.
        // The cell indices are unchanged (same elements, just more nodes).
        // We only need to update the element type strings.
        //
        // Note: ConvertLinearToQuadraticTriangles() already updated
        // the element types in vset_q. Here we update the ModelTopology
        // to match, because Model<2U> is constructed from BOTH.
        ModelTopology  topo_q( true );
        topo_q.ModelName( topo_p2_current.ModelName().c_str() );

        const bool isoparametric = topo_p2_current.IsoparametricFiniteElements();

        for ( auto it  = topo_p2_current.DomainsBegin();
                   it != topo_p2_current.DomainsEnd(); ++it )
        {
            assert( (*it).second.first.size() == 1 );
            const string  linear_etype = *(*it).second.first.begin();

            // Convert linear element type string to quadratic
            const auto q_type = getQuadraticType(
                isoparametric,
                static_cast<int8_t>( parseFiniteElementType(linear_etype) ) );

            const set<string>  q_fem_types{ parseFiniteElementType(q_type) };

            // Cell indices are unchanged -- same elements, just more nodes per element
            topo_q.AddDomain( (*it).first.c_str(),
                              q_fem_types,
                              (*it).second.second );
        }

        // Sanity check
        if ( !topo_q.QuadraticElementMesh() )
            output( "WARNING: quadratic topology not recognised." );

        // Convert to isoparametric local coordinate system
        VSetConverter().ConvertElementTypesToOnesUsingLocalCoordinateSystem( vset_q );

        Model<2U>  model_q( topo_q, vset_q, variable_file.c_str(), true );
        model_q.CreateProperty( "pressure error", "perr", "Pa" );

        const double L2 = solve_and_measure( model_q );
        const double h  = 1. / sqrt( static_cast<double>( model_q.Mesh().Elements()) / 2. );

        h_p2.push_back( h );
        L2_p2.push_back( L2 );

        const double rate = ( level > 0 )
            ? log( L2_p2[L2_p2.size()-2] / L2_p2.back() )
            / log( h_p2[h_p2.size()-2]   / h_p2.back()  )
            : -1.;

        // Store for deferred printing
        p2_results.push_back({ level,
                               model_q.Mesh().Nodes(),
                               model_q.Mesh().Elements(),
                               h, L2, rate });

        if ( level == 0 || level == n_levels_p2 ) {
            vtk.OutputDataToVTK( model_q, "p2_pressure", "fluid pressure", static_cast<long>(level) );
            vtk.OutputDataToVTK( model_q, "p2_error", "pressure error", static_cast<long>(level) );
        }

        // ------------------------------------------------------------------
        // Refine the LINEAR mesh for next level
        // (we always refine the linear mesh, then convert to quadratic)
        // ------------------------------------------------------------------
        if ( level < n_levels_p2 )
        {
            output( "Refining linear mesh to level " + to_string(level+1) + "..." );
            cout.flush();

            const vector<size_t>  parent_cell_indices = vset_p2_current.Refine();

            // Rebuild linear topology using parent cell mapping
            ModelTopology  topo_refined( true );
            topo_refined.ModelName( topo_p2_current.ModelName().c_str() );

            for ( auto sdom  = topo_p2_current.DomainsBegin();
                       sdom != topo_p2_current.DomainsEnd(); ++sdom )
            {
                vector<size_t>  new_cell_ids;
                new_cell_ids.reserve( (*sdom).second.second.size() * 4 );

                for ( const size_t old_cid : (*sdom).second.second )
                    for ( size_t new_idx{0};
                          new_idx < parent_cell_indices.size(); ++new_idx )
                        if ( parent_cell_indices[new_idx] == old_cid )
                            new_cell_ids.push_back( new_idx );

                assert( (*sdom).second.first.size() == 1 );
                topo_refined.AddDomain( (*sdom).first.c_str(),
                                        (*sdom).second.first,
                                        new_cell_ids );
            }

            topo_refined.CheckCellNumbering();
            topo_p2_current = topo_refined;

            output( "Level " + to_string(level+1) + ": "
                    + to_string(vset_p2_current.Elements()) + " elements." );
        }
    }
    divider('-');

    // Print P2 table after all solves
    blank();
    divider('-');
    cout << "  QUADRATIC TRIANGLE (P2) CONVERGENCE\n";
    divider('-');
    cout << "  " << left
         << setw(8)  << "Level"
         << setw(10) << "Nodes"
         << setw(12) << "Elements"
         << setw(12) << "h"
         << setw(16) << "L2 error"
         << setw(10) << "Rate" << "\n";
    divider('-');

    for ( const auto& r : p2_results )
        print_convergence_row( r.level, r.nodes, r.elements,
                               r.h, r.L2, r.rate );

    divider('-');

    if ( p2_results.size() >= 2 ) {
        const double last_rate = p2_results.back().rate;
        blank();
        output( "Final P2 convergence rate: "
                + dbl_str(last_rate,3) );
        if ( last_rate >= 3.5 && last_rate <= 4.5 )
            output( "PASS: rate ~4 as expected for P2 triangles." );
        else
            output( "WARNING: rate deviates from expected value of 4." );
    }
 }
pause();



    // ========================================================================
    // Step 6: Comparison and interpretation
    // ========================================================================
    heading( 6, 6, "Comparison and interpretation" );

    explain(
"The convergence study demonstrates two fundamental FEM concepts:\n"
"\n"
"1. CONVERGENCE: the numerical solution approaches the exact solution\n"
"   as the mesh is refined. This is guaranteed for well-posed problems\n"
"   and consistent, stable FEM formulations.\n"
"\n"
"2. CONVERGENCE RATE: the rate at which the error decreases.\n"
"   This depends on the polynomial order of the elements:\n"
"     P1 (linear):    rate ~2  -- error ~ h^2\n"
"     P2 (quadratic): rate ~4  -- error ~ h^4\n"
"\n"
"Practical implication:\n"
"  To achieve error < 1e-4 starting from error ~ 1e-1:\n"
"    P1: need 3 halvings of h  (8x more elements per dimension)\n"
"    P2: need 1.5 halvings of h (3x more elements per dimension)\n"
"\n"
"P2 elements are more expensive per element but far more efficient\n"
"per unit of accuracy for smooth problems."
    );
    pause();

    subheading( "Efficiency comparison" );
    {
        // Compute approximate cost comparison
        // Assume error_P1 ~ C1 * h^2, error_P2 ~ C2 * h^4
        // For same error epsilon:
        //   P1: h_P1 = sqrt(epsilon/C1)  --> N_P1 ~ 1/h_P1^2
        //   P2: h_P2 = (epsilon/C2)^0.25 --> N_P2 ~ 1/h_P2^2

        cout << "  Cost to achieve L2 error < epsilon:\n\n";
        cout << "  " << left
             << setw(14) << "Target error"
             << setw(20) << "P1 elements needed"
             << setw(20) << "P2 elements needed"
             << setw(12) << "P2/P1 ratio" << "\n";
        divider('-');

        // Use approximate constants from the convergence study
        // C1 ~ 0.1 (typical for P1 on unit square)
        // C2 ~ 0.01 (typical for P2 on unit square)
        constexpr double C1 = 0.1;
        constexpr double C2 = 0.01;

        const double targets[] = { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6 };
        for ( double eps : targets ) {
            const double h_p1 = sqrt( eps / C1 );
            const double h_p2 = pow( eps / C2, 0.25 );
            const double N_p1 = 2. / (h_p1 * h_p1);  // 2 triangles per h^2
            const double N_p2 = 2. / (h_p2 * h_p2);
            const double ratio = N_p2 / N_p1;

            cout << "  " << left
                 << setw(14) << sci_str(eps,0)
                 << setw(20) << to_string(static_cast<size_t>(N_p1))
                 << setw(20) << to_string(static_cast<size_t>(N_p2))
                 << setw(12) << dbl_str(ratio,3) << "\n";
        }
        divider('-');
        output( "P2 requires far fewer elements for the same accuracy." );
        output( "The advantage grows as the target error decreases." );
    }
    pause();

    subheading( "When does high-order NOT help?" );
    explain(
"P2 elements are NOT always better than P1:\n"
"\n"
"  1. NON-SMOOTH SOLUTIONS: if p has discontinuities or sharp gradients\n"
"     (e.g. at material interfaces), the convergence rate drops to ~1\n"
"     for both P1 and P2. Mesh refinement near the singularity\n"
"     (adaptive refinement) is more effective than high order.\n"
"\n"
"  2. HETEROGENEOUS PERMEABILITY: if k(x,y) is discontinuous,\n"
"     the solution gradient is discontinuous at material interfaces.\n"
"     P2 elements do not help here -- use mesh refinement instead.\n"
"\n"
"  3. VERY COARSE MESHES: on the coarsest meshes, P2 may not yet\n"
"     be in the asymptotic regime. The rate may be lower than 4.\n"
"\n"
"For the Jura-slope1 model in Example 6, P1 elements are appropriate\n"
"because the permeability is discontinuous at formation boundaries."
    );
    pause();

    subheading( "Reading the VTK output in ParaView" );
    explain(
"Open the VTK files in ParaView:\n"
"\n"
"  p1_pressure_0.vtk   -- P1 solution on coarsest mesh\n"
"  p1_pressure_N.vtk   -- P1 solution on finest mesh\n"
"  p1_error_0.vtk      -- pointwise error on coarsest mesh\n"
"  p1_error_N.vtk      -- pointwise error on finest mesh\n"
"  p2_pressure_0.vtk   -- P2 solution on coarsest mesh\n"
"  p2_error_0.vtk      -- P2 pointwise error on coarsest mesh\n"
"\n"
"What to look for:\n"
"  1. The pressure field should look like sin(pi*x)*sin(pi*y)\n"
"     -- a smooth dome with maximum at (0.5, 0.5)\n"
"  2. The error field should also peak at (0.5, 0.5)\n"
"     -- where the solution curvature is largest\n"
"  3. The P2 error on the coarsest mesh should be comparable\n"
"     to the P1 error on a much finer mesh\n"
"  4. Use 'Plot Over Line' to extract a 1D profile and compare\n"
"     the numerical and exact solutions directly"
    );
    pause();

    keypoint(
"Mesh convergence confirms the FEM implementation is correct.\n"
"P1 rate ~2, P2 rate ~4 for smooth solutions on the unit square.\n"
"P2 is more efficient per unit of accuracy for smooth problems.\n"
"For non-smooth problems (discontinuous k), use mesh refinement.\n"
"The MMS is the standard verification tool for FEM codes."
    );
    pause();


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 7: Mesh Convergence Study\n\n"
         << "  Steps covered:\n\n"
         << "    1.  The method of manufactured solutions\n"
         << "    2.  The unit square mesh and refinement strategy\n"
         << "    3.  Expected convergence rates\n"
         << "    4.  Linear triangle (P1) convergence study\n"
         << "    5.  Quadratic triangle (P2) convergence study\n"
         << "    6.  Comparison and interpretation\n\n"
         << "  Key concepts:\n\n"
         << "    - MMS: choose p_exact, compute f, solve, compare\n"
         << "    - refineSimplexMesh(): splits each triangle into 4\n"
         << "    - ConvertLinearToQuadraticTriangles(): adds midside nodes\n"
         << "    - P1 convergence rate ~2  (error ~ h^2)\n"
         << "    - P2 convergence rate ~4  (error ~ h^4)\n"
         << "    - L2 error norm: sqrt( sum (p_h - p_exact)^2 / N )\n"
         << "    - Convergence rate = log(e_i/e_{i+1}) / log(h_i/h_{i+1})\n"
         << "    - P2 far more efficient for smooth problems\n"
         << "    - High order does NOT help for non-smooth solutions\n\n"
         << "  Next: Example 8 -- Non-linear iteration.\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nMeshConvergence_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp

