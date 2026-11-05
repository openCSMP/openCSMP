// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "AtmosphericPressureDiffusion_Example.h"
#include "exampleHelpers.h"

#include "Model.h"
#include "Model1D.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyHandle.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"
#include "VTU_Interface.h"
#include "Standard_IO_Handler.h"
#include "RegionMonitor.h"

// PDE operators
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"

// Interrelation for hydraulic conductivity
#include "ConstantFactor.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

using namespace std;

namespace csmp {

// ============================================================================
// Specifications
// ============================================================================

void AtmosphericPressureDiffusion_Example::Specifications()
{
    SetTitle( "Atmospheric pressure diffusion into rock (interactive)" );
    SetDifficulty( 3 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: AtmosphericPressureDiffusion_Example.cpp" );
    AddDescription( "Transient pressure diffusion driven by atmospheric "
                    "pressure fluctuations. Verified against the analytical "
                    "solution in 1D, then applied to a 2D field-scale "
                    "geological cross-section of the Swiss Jura Mountains." );
    AddRequirement( "Jura-slope1 (CSMP binary files, -configuration.txt), "
                    "CSMP_field_scale_mechanics_variables.txt" );
}


// ============================================================================
// Analytical solution helper
// ============================================================================

namespace {

/**
    Analytical solution for pressure diffusion from a step change at z=0.

    p(z,t) = p_low + (p_high - p_low) * erf( z / (2*sqrt(D*t)) )

    @param z      depth below surface [m]
    @param t      time since step change [s]
    @param D      hydraulic diffusivity [m^2/s]
    @param p_high initial pressure (t<0) [Pa]
    @param p_low  surface pressure after step change [Pa]
    @return       pressure at depth z and time t [Pa]
*/
double analytical_pressure( double z, double t,
                             double D,
                             double p_high, double p_low )
{
    if ( t <= 0. ) return p_high;
    const double arg = z / ( 2. * sqrt(D * t) );
    return p_low + ( p_high - p_low ) * erf( arg );
}

} // anonymous namespace


// ============================================================================
// Run
// ============================================================================

/**
    @file AtmosphericPressureDiffusion_Example.cpp

    @brief CSMP for Beginners -- Example 6: Atmospheric Pressure Diffusion

    @section purpose Purpose

    This interactive example demonstrates transient pressure diffusion
    using atmospheric pressure fluctuations as the driving signal.
    It proceeds in two stages:

    Stage 1 -- 1D analytical verification:
      A sudden pressure drop at the surface (from 104000 Pa to 96000 Pa,
      simulating a weather system passing through) diffuses into a
      homogeneous half-space. The CSMP result is compared with the
      analytical solution and both are written to the same VTU file.

    Stage 2 -- 2D field-scale simulation:
      The same pressure drop is applied to the surface of the
      Jura-slope1 model, a 2D geological cross-section of the
      Swiss Jura Mountains containing five contrasting formations.

    @section physics Physics

    The governing equation is the pressure diffusion equation:

        beta_t * dp/dt = div( k/mu * grad p )

    where:
        beta_t = phi*beta_f + Kb  [Pa^-1]  total system compressibility
        k      [m^2]               permeability
        mu     [Pa s]              fluid viscosity

    The hydraulic diffusivity is:
        D = k / (mu * beta_t)     [m^2/s]

    The analytical solution for a step change at the surface is:
        p(z,t) = p_low + (p_high - p_low) * erf( z / (2*sqrt(D*t)) )

    @section model The Jura-slope1 Model

    A 2D cross-section of the Jura Mountains in southern Germany containing:
        REEF_DOLOMITE    -- stiff carbonate, high diffusivity
        IRON_SANDSTONE   -- iron-cemented sandstone
        OPALINUS_CLAY    -- soft clay shale (Mont Terri formation)
        AMALTHEA_CLAY    -- soft clay, very low diffusivity
        RHAET_SANDSTONE  -- Rhaetian sandstone, higher porosity
 
    @section requirements Requirements

    - Jura-slope1  (CSMP binary files + -configuration.txt)
    - CSMP_field_scale_mechanics_variables.txt

    @author S.K. Matthai
*/
void AtmosphericPressureDiffusion_Example::Run()
{
    string example_name  = "AtmosphericPressureDiffusion_Example";
    string model_name    = "Jura-slope1";
    string variable_file = "CSMP_field_scale_mechanics_variables.txt";
    string config_file = model_name;
    //find the name of this current example source file
    string file_name = GetExampleFileName(__FILE__);
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(example_name, model_name, variable_file, config_file);

    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 6\n"
         << "  Atmospheric Pressure Diffusion into Rock\n\n"
         << "  How quickly do weather-driven pressure changes\n"
         << "  penetrate into the subsurface?\n\n"
         << "  Stage 1: 1D verification against analytical solution\n"
         << "  Stage 2: 2D field-scale simulation (Jura Mountains)\n\n"
         << "  Steps:\n"
         << "    1.  The physics: pressure diffusion equation\n"
         << "    2.  The atmospheric pressure signal\n"
         << "    3.  The analytical solution\n"
         << "    4.  Stage 1: 1D CSMP solve and comparison\n"
         << "    5.  Hydraulic diffusivity of the Jura formations\n"
         << "    6.  Stage 2: loading the Jura-slope1 model\n"
         << "    7.  Computing storativity and conductivity\n"
         << "    8.  Setting up the transient PDE_Integrator\n"
         << "    9.  The time stepping loop\n"
         << "   10.  Results and interpretation\n";
    divider();
    pause( "  [ Press ENTER to begin ]" );


    // ========================================================================
    // Step 1: The physics
    // ========================================================================
    heading( 1, 10, "The physics: pressure diffusion equation" );

    explain(
"When atmospheric pressure changes at the surface, the pressure\n"
"signal diffuses into the subsurface through the porous rock.\n"
"\n"
"The governing equation is:\n"
"\n"
"    beta_t * dp/dt = div( k/mu * grad p )\n"
"\n"
"where:\n"
"    beta_t = phi*beta_f + Kb   [Pa^-1]  total system compressibility\n"
"    k                          [m^2]    permeability\n"
"    mu                         [Pa s]   fluid viscosity\n"
"\n"
"The hydraulic diffusivity D = k / (mu * beta_t) [m^2/s]\n"
"controls how fast the signal propagates.\n"
"\n"
"The skin depth (depth at which amplitude decays to 1/e) is:\n"
"\n"
"    delta = sqrt( D * T / pi )\n"
"\n"
"where T is the period of the pressure signal."
    );
    pause();

    subheading( "The CSMP discretisation -- Backward Euler" );
    code(
"// Backward Euler: ([S]/dt + [C]){p}^{t+dt} = [S]/dt * {p}^t + {Q}\n"
"\n"
"// LHS: conductance matrix [C]\n"
"NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(\n"
"    model.Database(), \"conductivity\",\n"
"    \"fluid pressure\", \"fluid pressure\" );\n"
"\n"
"// LHS: storage matrix [S]/dt\n"
"NumIntegral_NT_lhsop_N_dV<dim>  capacitance_lhs(\n"
"    model.Database(), \"storativity\",\n"
"    \"fluid pressure\", \"fluid pressure\" );\n"
"capacitance_lhs.LumpedFormulation(true);\n"
"capacitance_lhs.MultiplyWithTimeIncrement(true);\n"
"\n"
"// RHS: initial conditions [S]/dt * {p}^t\n"
"NumIntegral_NT_rhsop_N_dV<dim>  capacitance_rhs(\n"
"    model.Database(), \"storativity\", \"fluid pressure\" );\n"
"capacitance_rhs.LumpedFormulation(true);\n"
"capacitance_rhs.MultiplyWithTimeIncrement(true);"
    );
    pause();

    keypoint(
"The hydraulic diffusivity related storage term [S]/dt (S ~= beta_t)\n"
" appears on BOTH sides of the time-discretized equation.\n"
"LHS: capacitance_lhs -- the implicit part (unknown p^{t+dt})\n"
"RHS: capacitance_rhs -- the explicit part (known p^t)\n"
"LumpedFormulation(true) uses volume-weighted nodal lumping.\n"
"MultiplyWithTimeIncrement(true) means the operator stores beta_t,\n"
"and the integrator divides by dt automatically."
    );
    pause();


    // ========================================================================
    // Step 2: The atmospheric pressure signal
    // ========================================================================
    heading( 2, 10, "The atmospheric pressure signal" );

    explain(
"Air pressure in Southern Germany varies between approximately\n"
"960 hPa (deep low) and 1040 hPa (strong high) over the year.\n"
"\n"
"We model a sudden pressure drop from a high to a low as transmitted to the pore pressure:\n"
"    p_high = 104000 Pa  (1040 hPa -- strong anticyclone)\n"
"    p_low  =  96000 Pa  ( 960 hPa -- deep cyclone)\n"
"    Delta p = 8000 Pa   (  80 hPa)\n"
"\n"
"This is a step function in time -- the simplest possible signal.\n"
"It allows comparison with the analytical error function solution.\n"
"\n"
"The initial condition is p = p_high everywhere.\n"
"At t=0 the surface pressure drops to p_low (Dirichlet BC).\n"
"We then watch how this pressure drop diffuses into the rock."
    );
    pause();

    subheading( "Physical context" );
    explain(
"This experiment is directly relevant to:\n"
"\n"
"  Borehole monitoring: pressure sensors in boreholes record\n"
"  atmospheric loading. The barometric efficiency B = dp_bh/dp_atm\n"
"  depends on the storativity and Biot coefficient.\n"
"\n"
"  Sealing capabilities of Opalinus clay: atmospheric pressure variations have\n"
"  no influence on the pressure in these tight formations,\n"
"  but they are immediately felt in the fractured dolomite above.\n"
"\n"
    );
    pause();

    keypoint(
"A step pressure change at the surface is the simplest transient BC.\n"
"It has an exact analytical solution (error function).\n"
"This makes it ideal for verifying the CSMP transient solver.\n"
"The Jura slope is a real-world application of the same physics."
    );
    pause();


    // ========================================================================
    // Step 3: The analytical solution
    // ========================================================================
    heading( 3, 10, "The analytical solution" );

    explain(
"For a semi-infinite homogeneous medium with:\n"
"  - initial pressure p_high everywhere\n"
"  - step change to p_low at z=0 at t=0\n"
"\n"
"The exact solution is:\n"
"\n"
"    p(z,t) = p_low + (p_high - p_low) * erf( z / (2*sqrt(D*t)) )\n"
"\n"
"where erf() is the error function.\n"
"\n"
"Key properties:\n"
"  - At z=0: p = p_low for all t > 0  (BC enforced)\n"
"  - At z=inf: p = p_high for all t   (undisturbed)\n"
"  - The 'penetration depth' ~ 2*sqrt(D*t) grows as sqrt(t)\n"
"  - Doubling the time increases penetration by sqrt(2) ~ 41%"
    );
    pause();

    subheading( "Skin depth for different rock types" );
    {
        constexpr double mu      = 1.0e-3;   // Pa s (water at 20 degC)
        constexpr double T_year  = 365.*86400.; // annual period [s]
        constexpr double T_day   = 86400.;      // daily period [s]

        // Rock properties
        struct RockType {
            const char* name;
            double k;        // permeability [m^2]
            double beta_t;   // total compressibility [Pa^-1]
        };

        const RockType rocks[] = {
            { "Amalthea/Opalinus Clay", 1e-20, 5.275e-10 },
            { "Iron Sandstone",         1e-15, 9.5e-11   },
            { "Rhaet Sandstone",        1e-14, 1.9e-10   },
            { "Reef Dolomite",          1e-13, 3.65e-11  }
        };

        cout << "  Hydraulic diffusivity and skin depth:\n\n";
        cout << "  " << left << setw(28) << "Rock"
             << setw(14) << "D (m^2/s)"
             << setw(16) << "delta annual"
             << setw(16) << "delta 1-day" << "\n";
        divider('-');

        for ( const auto& r : rocks ) {
            const double D       = r.k / ( mu * r.beta_t );
            const double delta_y = sqrt( D * T_year / M_PI );
            const double delta_d = sqrt( D * T_day  / M_PI );
            cout << "  " << left << setw(28) << r.name
                 << setw(14) << sci_str(D,2)
                 << setw(16) << (dbl_str(delta_y,2) + " m")
                 << setw(16) << (dbl_str(delta_d,3) + " m") << "\n";
        }
        divider('-');
    }
    pause();

    keypoint(
"The skin depth delta = sqrt(D*T/pi) tells you how deep the signal\n"
"penetrates in one period T. For the annual cycle:\n"
"  Clay: ~cm  -->  atmospheric signal barely enters\n"
"  Sandstone: ~10-40 m  -->  signal penetrates tens of metres\n"
"  Dolomite: ~100s m  -->  signal felt hundreds of metres deep\n"
"This contrast is what makes the Jura cross-section so instructive."
    );
    pause();


        // ====================================================================
        // Step 4: Stage 1 -- 1D CSMP solve and comparison
        // ====================================================================
        heading( 4, 10, "Stage 1: 1D CSMP solve and comparison with analytical solution" );

        explain(
"Before running the complex 2D Jura model, we verify the CSMP\n"
"transient solver against the analytical solution in a simple\n"
"1D column.\n"
"\n"
"We use a single representative hydraulic diffusivity:\n"
"    k       = 1e-14 m^2   (Rhaet Sandstone)\n"
"    beta_t  = 1.9e-10 Pa^-1\n"
"    mu      = 1e-3 Pa s\n"
"    D       = k/(mu*beta_t) = 5.26e-2 m^2/s\n"
"\n"
"The mesh is graded: fine (1 m) near the surface where gradients\n"
"are steep, coarsening to 50 m at depth. This keeps the total\n"
"node count manageable while resolving the sharp front.\n"
"\n"
"Boundary conditions:\n"
"    x=0 (surface, CNR1): p = p_low  = 96000 Pa  (Dirichlet)\n"
"    x=2000 m (CNR2):     no-flux               (natural Neumann)\n"
"\n"
"Initial condition:\n"
"    p = p_high = 104000 Pa everywhere\n"
"\n"
"The analytical solution is:\n"
"    p(x,t) = p_low + (p_high - p_low) * erf( x / (2*sqrt(D*t)) )\n"
"\n"
"Both CSMP and analytical solutions are stored as NODE variables\n"
"and written to the same VTU file for direct comparison in ParaView."
        );
        pause();

        subheading( "Graded mesh with Model1D" );
        code(
"// Model1D with graded mesh:\n"
"// - min_size = 1 m  (near surface, where gradient is steep)\n"
"// - max_size = 50 m (at depth, where solution is smooth)\n"
"// - transition zone = 200 m (gradual coarsening)\n"
"\n"
"const vector<double>  no_splitnodes;\n"
"const Point<1U>  origin(    0. );\n"
"const Point<1U>  destination( 2000. );\n"
"\n"
"Model1D<1U>  model_1d(\n"
"    \"1D_column\",\n"
"    \"CSMP_field_scale_mechanics_variables.txt\",\n"
"    2000.,    // total length [m]\n"
"    1.,       // min element size [m]  -- near surface\n"
"    50.,      // max element size [m]  -- at depth\n"
"    200.,     // width of transition zone [m]\n"
"    nullptr,  // MeshDensity* (nullptr = default linear grading)\n"
"    no_splitnodes,\n"
"    origin,\n"
"    destination );\n"
"\n"
"// Nodes are NOT sorted by x-coordinate in the internal storage.\n"
"// We build a sorted map<depth, node*> for the comparison table."
        );
        pause();

        subheading( "Live 1D model construction and solve" );
        {
            // ------------------------------------------------------------------
            // Physical parameters
            // ------------------------------------------------------------------
            constexpr double k_1d      = 1.0e-14;   // m^2  (Rhaet Sandstone)
            constexpr double beta_t_1d = 1.9e-10;   // Pa^-1
            constexpr double mu_1d     = 1.0e-3;    // Pa s
            constexpr double D_1d      = k_1d / ( mu_1d * beta_t_1d );
            constexpr double p_high    = 104000.;   // Pa
            constexpr double p_low     =  96000.;   // Pa
            constexpr double day       = 86400.;    // s

            output( "D       = " + sci_str(D_1d,3) + " m^2/s" );
            output( "Skin depth (1 day)    = "
                    + dbl_str(sqrt(D_1d*1.*day/M_PI),0)   + " m" );
            output( "Skin depth (10 days)  = "
                    + dbl_str(sqrt(D_1d*10.*day/M_PI),0)  + " m" );
            output( "Skin depth (100 days) = "
                    + dbl_str(sqrt(D_1d*100.*day/M_PI),0) + " m" );

            // ------------------------------------------------------------------
            // Build the graded 1D model
            // ------------------------------------------------------------------
            output( "Building graded 1D model..." );
            cout.flush();

            const vector<double>  no_splitnodes;
            const Point<1U>  origin(      0.    );
            const Point<1U>  destination( 2000. );

            ErfDensity  density;   // default: a=3.0, b=1.0
                                   // density(t) = erf(3*t)
                                   // rapid coarsening away from surface

            Model1D<1U>  model_1d(
                "1D_column",
                "CSMP_field_scale_mechanics_variables.txt",
                2000.,          // total length [m]
                1.,             // dx_min [m]  -- finest element at surface
                50.,            // dx_max [m]  -- coarsest element at depth
                1000.,          // width of transition zone [m]
                &density,       // ErfDensity -- mirrors the erf solution gradient
                no_splitnodes,
                origin,
                destination );

            output( "1D model built." );
            printModelDimensions( model_1d );
            output( "Nodes: " + to_string(model_1d.Mesh().Nodes()) );

            Region<1U>&  gref_1d( model_1d.Region("Model") );

            // ------------------------------------------------------------------
            // Assign material properties
            // ------------------------------------------------------------------
            model_1d.InputPropertyValue( "permeability",
                                         makeScalar(FIELD_DATA, k_1d) );
            model_1d.InputPropertyValue( "storativity",
                                         makeScalar(FIELD_DATA, beta_t_1d) );
            model_1d.InputPropertyValue( "fluid volume source",
                                         makeScalar(PLAIN, 0.) );

            {
                ConstantFactor<1U,divides>  conductivity(
                    model_1d.Database(),
                    "conductivity", "permeability", mu_1d );
                model_1d.Apply( conductivity );
            }
            output( "conductivity = " + sci_str(k_1d/mu_1d,3) + " m^2/(Pa s)" );

            // ------------------------------------------------------------------
            // Create 'analytical pressure' variable
            // ------------------------------------------------------------------
            PropertyHandle<1U>  analytical_p(
                model_1d, "analytical pressure", SCALAR, NODE );
            model_1d.InputPropertyValue( "analytical pressure",
                                         makeScalar(PLAIN, p_high) );

            // ------------------------------------------------------------------
            // Initial and boundary conditions
            // ------------------------------------------------------------------
            model_1d.InputPropertyValue( "fluid pressure",
                                         makeScalar(ANY, p_high) );
            // Surface (CNR1, x=0): Dirichlet p = p_low
            model_1d.InputBoundaryValue( CNR1, "fluid pressure",
                                         makeScalar(DIRICH, p_low) );
            // Base (CNR2, x=2000 m): no BC -- natural Neumann = 0

            output( "Initial p = " + dbl_str(p_high/100.,0) + " hPa everywhere." );
            output( "BC at surface (CNR1): DIRICH p = "
                    + dbl_str(p_low/100.,0) + " hPa." );
            output( "BC at base (CNR2): no-flux (natural Neumann = 0)." );

            // ------------------------------------------------------------------
            // Build the transient PDE_Integrator
            // ------------------------------------------------------------------
            EigenSolver                 solver_1d; // direct solver
            PDE_Integrator<1U,Element>  pressure_1d( solver_1d );

            NumIntegral_dNT_lhsop_dN_dV<1U>  conductance_1d(
                model_1d.Database(), "conductivity",
                "fluid pressure", "fluid pressure" );

            NumIntegral_NT_lhsop_N_dV<1U>  cap_lhs_1d(
                model_1d.Database(), "storativity",
                "fluid pressure", "fluid pressure" );
            cap_lhs_1d.LumpedFormulation(true);
            cap_lhs_1d.MultiplyWithTimeIncrement(true);

            NumIntegral_NT_rhsop_N_dV<1U>  cap_rhs_1d(
                model_1d.Database(), "storativity", "fluid pressure" );
            cap_rhs_1d.LumpedFormulation(true);
            cap_rhs_1d.MultiplyWithTimeIncrement(true);

            NumIntegral_NT_rhsop_N_dV<1U>  source_1d(
                model_1d.Database(), "fluid volume source", "fluid pressure" );
            source_1d.AddAccumulateLater();
            source_1d.LumpedFormulation(true);

            pressure_1d.Add( &conductance_1d );
            pressure_1d.Add( &cap_lhs_1d );
            pressure_1d.Add( &cap_rhs_1d );
            pressure_1d.Add( &source_1d );

            output( "PDE_Integrator built." );

            // ------------------------------------------------------------------
            // Get index keys
            // ------------------------------------------------------------------
            const csmp::Index p_key  =
                model_1d.Database().StorageKey("fluid pressure");
            const csmp::Index ap_key =
                model_1d.Database().StorageKey("analytical pressure");

            // ------------------------------------------------------------------
            // Build a depth-sorted map of node pointers
            // Nodes are NOT stored in x-order internally.
            // map<double,Node<1U>*> sorts by x automatically.
            // ------------------------------------------------------------------
            map<double, Node<1U>*>  nodes_by_depth;
            for ( auto& nit : gref_1d.NodeVector() )
                nodes_by_depth[ nit->x() ] = nit;

            output( "Built depth-sorted node map ("
                    + to_string(nodes_by_depth.size()) + " nodes)." );

            // ------------------------------------------------------------------
            // VTU output setup
            // ------------------------------------------------------------------
            list<string>  output_vars{
                "fluid pressure",
                "analytical pressure"
            };
            VTU_Interface<1U>  vtu_1d( model_1d );

            vtu_1d.OutputDataToVTU( "1D_verification",
                                    output_vars, gref_1d, 0 );
            output( "Written: 1D_verification_0.vtu  (initial state)" );
            output( "\n   Check results and BE PATIENT: we are using very small time-steps and a serial direct solver..." );
            pause();

            // ------------------------------------------------------------------
            // Time loop: solve and compare at t = 1, 10, 100 days
            // ------------------------------------------------------------------
            // Time step: Courant condition for the finest element (1 m)
            // dt_courant = dx_min^2 / (2*D) = 1^2 / (2*5.26e-2) ~ 9.5 s
            // Use dt = 10 s initially, allow geometric growth but cap at
            // dt_max = dx_min^2 / D = 19 s to maintain accuracy near surface
            constexpr double dx_min    = 1.;
            constexpr double dt_start  = 0.5 * dx_min * dx_min / D_1d;
            constexpr double dt_max    = dx_min * dx_min / D_1d;

            output( "dt_start = " + dbl_str(dt_start,2) + " s" );
            output( "dt_max   = " + dbl_str(dt_max,2)   + " s" );
            output( "Steps to 1 day:   ~"
                    + to_string(static_cast<size_t>(day/dt_start)) );

            const double output_times[] = { 1.*day, 10.*day, 100.*day };
            double model_time_1d = 0.;
            double dt_1d         = dt_start;
            size_t step_1d       = 0;
            size_t output_idx    = 0;

            output( "Running 1D transient solve..." );
            blank();

            while ( output_idx < 3 )
            {
                // Solve one step
                pressure_1d.TimeIncrement( 1. / dt_1d );
                model_1d.Apply( pressure_1d );
                model_time_1d += dt_1d;
                step_1d++;

                // Gradually increase dt but cap at dt_max
                dt_1d = min( dt_1d * 1.05, dt_max );

                // Check whether we have reached the next output time
                if ( model_time_1d >= output_times[output_idx] )
                {
                    const double t_days = model_time_1d / day;

                    // ----------------------------------------------------------
                    // Update analytical solution at ALL nodes using their
                    // actual x-coordinate (depth from surface).
                    // Simultaneously build error statistics.
                    // ----------------------------------------------------------
                    double max_error     = 0.;
                    double max_error_pct = 0.;

                    for ( auto& [depth, nptr] : nodes_by_depth ) {
                        const double p_anal = analytical_pressure(
                            depth, model_time_1d, D_1d, p_high, p_low );
                        nptr->Store( ap_key, makeScalar(PLAIN, p_anal) );

                        ScalarVariable sc;
                        nptr->Read( p_key, sc );
                        const double error     = fabs( sc() - p_anal );
                        const double error_pct = 100. * error / (p_high - p_low);
                        max_error     = max( max_error,     error );
                        max_error_pct = max( max_error_pct, error_pct );
                    }

                    // Write VTU with both solutions
                    vtu_1d.OutputDataToVTU( "1D_verification",
                                            output_vars, gref_1d,
                                            static_cast<long>(t_days) );

                    // ----------------------------------------------------------
                    // Print comparison table.
                    // Iterate over depth-sorted map and print every node
                    // up to 50 m (fine mesh), then every 10th node beyond.
                    // This shows the steep gradient near the surface clearly.
                    // ----------------------------------------------------------
                    cout << "\n  t = " << fixed << setprecision(1)
                         << t_days << " days  ("
                         << step_1d << " steps, dt = "
                         << dbl_str(dt_1d,2) << " s)\n\n";

                    cout << "  " << left
                         << setw(12) << "depth (m)"
                         << setw(18) << "CSMP (Pa)"
                         << setw(18) << "Analytical (Pa)"
                         << setw(14) << "error (Pa)"
                         << setw(10) << "error (%)" << "\n";
                    divider('-');

                    size_t print_count = 0;
                    size_t coarse_counter = 0;

                    for ( const auto& [depth, nptr] : nodes_by_depth ) {

                        // Print all nodes up to 50 m (fine mesh region)
                        // then every 10th node beyond 50 m
                        const bool in_fine_region  = ( depth <= 50. );
                        const bool print_this_node = in_fine_region ||
                                                     ( coarse_counter % 10 == 0 );

                        if ( !in_fine_region ) ++coarse_counter;

                        if ( print_this_node ) {
                            ScalarVariable sc;
                            nptr->Read( p_key, sc );
                            const double p_csmp = sc();
                            const double p_anal = analytical_pressure(
                                depth, model_time_1d, D_1d, p_high, p_low );
                            const double error     = fabs( p_csmp - p_anal );
                            const double error_pct = 100.*error/(p_high-p_low);

                            cout << "  " << left
                                 << setw(12) << dbl_str(depth,1)
                                 << setw(18) << dbl_str(p_csmp,1)
                                 << setw(18) << dbl_str(p_anal,1)
                                 << setw(14) << dbl_str(error,2)
                                 << setw(10) << dbl_str(error_pct,3) << "\n";
                            ++print_count;
                        }
                    }
                    divider('-');
                    output( "Rows printed: " + to_string(print_count)
                            + " (all nodes <= 50 m, every 10th beyond)" );
                    output( "Max error (all " + to_string(nodes_by_depth.size())
                            + " nodes): " + dbl_str(max_error,2)
                            + " Pa  = " + dbl_str(max_error_pct,3)
                            + "% of Delta p" );
                    output( "Written: 1D_verification_"
                            + to_string(static_cast<long>(t_days)) + ".vtu" );
                    blank();

                    output_idx++;
                }
            }

            // ------------------------------------------------------------------
            // Summary
            // ------------------------------------------------------------------
            blank();
            divider('-');
            cout << "  1D VERIFICATION SUMMARY\n";
            divider('-');
            cout << "  Rock type:  Rhaet Sandstone\n"
                 << "  k       = " << sci_str(k_1d,2)      << " m^2\n"
                 << "  beta_t  = " << sci_str(beta_t_1d,2) << " Pa^-1\n"
                 << "  D       = " << sci_str(D_1d,3)       << " m^2/s\n"
                 << "  Mesh:     graded 1-50 m, 2000 m total\n"
                 << "  Nodes:    " << nodes_by_depth.size() << "\n\n"
                 << "  Files written:\n"
                 << "    1D_verification_0.vtu    (initial state)\n"
                 << "    1D_verification_1.vtu    (t = 1 day)\n"
                 << "    1D_verification_10.vtu   (t = 10 days)\n"
                 << "    1D_verification_100.vtu  (t = 100 days)\n\n"
                 << "  In ParaView:\n"
                 << "    1. Open all four files as a group\n"
                 << "    2. Apply 'Plot Over Line' filter along x-axis\n"
                 << "    3. Plot 'fluid pressure' and 'analytical pressure'\n"
                 << "       on the same axis\n"
                 << "    4. The two curves should overlap closely\n"
                 << "    5. Largest errors near x=0 (steep gradient)\n\n"
                 << "  Acceptance criterion:\n"
                 << "    Max error < 1% of Delta p = "
                 << dbl_str(0.01*(p_high-p_low),0) << " Pa\n";
            divider('-');

        } // end 1D model scope
        pause();

        keypoint(
"Model1D with graded mesh: fine near surface (1 m), coarse at depth (50 m).\n"
"Nodes are NOT stored in x-order -- build a map<double,Node*> sorted by depth.\n"
"Iterate over the sorted map to update analytical solution and print table.\n"
"Print all nodes in the fine region (<=50 m) to show the steep gradient.\n"
"Max error over ALL nodes is the acceptance criterion, not just printed rows."
        );
        pause();


    // ========================================================================
    // Step 5: Hydraulic diffusivity of the Jura formations
    // ========================================================================
    heading( 5, 10, "Hydraulic diffusivity of Jurassic succession in southern Germany" );

    explain(
"The Jura-slope1 model contains five geological formations\n"
"with very different hydraulic properties:\n"
"\n"
"  REEF_DOLOMITE    -- Reef dolomite (stiff carbonate)\n"
"  IRON_SANDSTONE   -- Iron-cemented sandstone (Dogger beta)\n"
"  OPALINUS_CLAY    -- Opalinus Clay (like in Mont Terri nuclear waste site)\n"
"  AMALTHEA_CLAY    -- Amalthea member, soft clay\n"
"  RHAET_SANDSTONE  -- Rhaetian sandstone (important aquifer)\n"
"\n"
"The contrast in hydraulic diffusivity spans 7 orders of magnitude."
    );
    pause();

    subheading( "Storativity calculation" );
    code(
"// Total system compressibility:\n"
"// beta_t = phi * beta_f + Kb\n"
"//\n"
"// beta_f = 4.5e-10 Pa^-1  (water compressibility)\n"
"// Kb = 1/K_dr             (drained bulk compressibility)\n"
"\n"
"// Amalthea/Opalinus Clay: phi=0.15, K_dr=2.2 GPa\n"
"// beta_t = 0.15 * 4.5e-10 + 4.6e-10 = 5.275e-10 Pa^-1\n"
"model.Region(\"AMALTHEA_CLAY\").InputPropertyValue(\n"
"    \"storativity\", makeScalar(FIELD_DATA, 5.275e-10) );\n"
"model.Region(\"OPALINUS_CLAY\").InputPropertyValue(\n"
"    \"storativity\", makeScalar(FIELD_DATA, 5.275e-10) );\n"
"\n"
"// Iron Sandstone: phi=0.10, K_dr=20 GPa\n"
"// beta_t = 0.10 * 4.5e-10 + 5.0e-11 = 9.5e-11 Pa^-1\n"
"model.Region(\"IRON_SANDSTONE\").InputPropertyValue(\n"
"    \"storativity\", makeScalar(FIELD_DATA, 9.5e-11) );\n"
"\n"
"// Rhaet Sandstone: phi=0.20, K_dr=10 GPa\n"
"// beta_t = 0.20 * 4.5e-10 + 1.0e-10 = 1.9e-10 Pa^-1\n"
"model.Region(\"RHAET_SANDSTONE\").InputPropertyValue(\n"
"    \"storativity\", makeScalar(FIELD_DATA, 1.9e-10) );\n"
"\n"
"// Reef Dolomite: phi=0.05, K_dr=70 GPa\n"
"// beta_t = 0.05 * 4.5e-10 + 1.4e-11 = 3.65e-11 Pa^-1\n"
"model.Region(\"REEF_DOLOMITE\").InputPropertyValue(\n"
"    \"storativity\", makeScalar(FIELD_DATA, 3.65e-11) );"
    );
    pause();

    subheading( "Why FIELD_DATA flag?" );
    explain(
"We use FIELD_DATA instead of PLAIN for the storativity.\n"
"\n"
"FIELD_DATA means: 'this value comes from measured/computed data\n"
"and should not be overwritten by the solver or post-processing.'\n"
"\n"
"This can be used to prevent computations from accidentally overwriting\n"
"storativity values.\n"
"\n"
"It is good practice to flag property values as FIELD_DATA when they represent actual measurements\n"
"to distinguish them from variables computed via constitutive relationships or computations (ANY)."
    );
    pause();

    keypoint(
"Storativity beta_t = phi * beta_f + (1-phi) * Kb combines fluid and rock compressibility.\n"
"Use FIELD_DATA flag for material properties -- prevents accidental overwrite.\n"
"The 7-order-of-magnitude contrast in D across the Jura formations\n"
"means that atmospheric pressure variations penetrates cm in clay but hm in dolomite."
    );
    pause();


    // ========================================================================
    // Step 6: Loading the Jura-slope1 model
    // ========================================================================
    heading( 6, 10, "Stage 2: loading the Jura-slope1 model" );

    explain(
"The Jura-slope1 model is a 2D cross-section of the German Jura\n"
"table-top mountains south of Stuttgart, oriented perpendicular to their slope.\n"
"\n"
"The model represents a slope where the geological formations\n"
"dip to the south away from the valley. The surface boundary (IRREGULAR) follows\n"
"the topography and is where the atmospheric pressure BC is applied.\n"
"\n"
"The model is loaded from CSMP binary files and configured\n"
"using the InputDataManager as in Example 5."
    );
    pause();

    subheading( "Loading and configuring" );
    {
        constexpr uint32_t dim{ 2U };
        cout << "  Loading model 'Jura-slope1'..." << endl;
        Model<dim>  model( model_name, variable_file );

        output( "Model loaded: " + string(model.Name()) );
        printModelDimensions( model );
        model.RegionsOut();
        model.BoundariesOut();
        pause();

        // ====================================================================
        // Step 7: Computing storativity and conductivity
        // ====================================================================
        heading( 7, 10, "Computing storativity and conductivity" );

        explain(
"Two derived properties must be computed before the PDE solve:\n"
"\n"
"  1. Hydraulic conductivity = k / mu\n"
"     Computed from permeability using a ConstantFactor -> Interrelation subclass.\n"
"\n"
"  2. Storativity = phi*beta_f + Kb\n"
"     Assigned per region using the values derived in Step 5."
        );
        pause();

        subheading( "Applying configuration file" );
        {
            const bool region_specifications{false},
                       default_property_values{true},
                       region_property_values{true},
                       box_boundary_conditions{true},
                       region_property_conditions{true},
                       boundary_conditions{true};

            ComputationalSettings settings;
            InputDataManager<dim>  model_configuration;
            model_configuration.ConfigureFromFile(
                model, model.Name(),
                region_specifications,
                default_property_values,
                region_property_values,
                box_boundary_conditions,
                region_property_conditions,
                boundary_conditions,
                settings );

            output( "Configuration applied." );
        }
        pause();

        subheading( "Computing hydraulic conductivity" );
        {
            constexpr double fluid_viscosity = 1.0e-3;  // Pa s (water at 20 degC)
            ConstantFactor<dim,divides>  conductivity(
                model.Database(), "conductivity", "permeability", fluid_viscosity );
            model.Apply( conductivity );
            output( "Hydraulic conductivity = permeability / viscosity:" );
            printRangeOfVariable( model, "conductivity" );
        }
        pause();

        subheading( "Assigning storativity per region" );
        {
            // Total compressibility beta_t = phi*beta_f + Kb [Pa^-1]
            // beta_f = 4.5e-10 Pa^-1 (water at 20 degC)
            model.Region("AMALTHEA_CLAY").InputPropertyValue(
                "storativity", makeScalar(FIELD_DATA, 5.275e-10) );
            model.Region("OPALINUS_CLAY").InputPropertyValue(
                "storativity", makeScalar(FIELD_DATA, 5.275e-10) );
            model.Region("IRON_SANDSTONE").InputPropertyValue(
                "storativity", makeScalar(FIELD_DATA, 9.5e-11) );
            model.Region("RHAET_SANDSTONE").InputPropertyValue(
                "storativity", makeScalar(FIELD_DATA, 1.9e-10) );
            model.Region("REEF_DOLOMITE").InputPropertyValue(
                "storativity", makeScalar(FIELD_DATA, 3.65e-11) );

            output( "Storativity assigned per region:" );
            printRangeOfVariable( model, "storativity" );
        }
        pause();


        // ====================================================================
        // Step 8: Setting up the transient PDE_Integrator
        // ====================================================================
        heading( 8, 10, "Setting up the transient PDE_Integrator" );

        explain(
"The transient pressure solver uses Backward Euler time-stepping.\n"
"This is unconditionally stable -- however, time step size is capped to 14 days.\n"
"This is because accuracy requires dt << L^2/D (diffusion timescale).\n"
"\n"
"We start with a small time step (5 seconds) and gradually\n"
"increase it by a factor of 1.2 each step. This geometric\n"
"progression efficiently covers many orders of magnitude in time."
        );
        pause();

        subheading( "Initial and boundary conditions" );
        code(
"// Initial condition: high pressure everywhere\n"
"constexpr double p_high = 104000.;  // Pa (1040 hPa)\n"
"constexpr double p_low  =  96000.;  // Pa ( 960 hPa)\n"
"\n"
"model.InputPropertyValue( \"fluid pressure\",\n"
"                          makeScalar(ANY, p_high) );\n"
"\n"
"// Dirichlet BC: pressure drop at surface (IRREGULAR boundary)\n"
"model.Boundary(\"IRREGULAR\").InputPropertyValue(\n"
"    \"fluid pressure\", makeScalar(DIRICH, p_low) );\n"
"\n"
"// Note: IRREGULAR is the topographic surface boundary.\n"
"// It is named IRREGULAR because it follows the terrain\n"
"// rather than being a flat box boundary."
        );
        pause();

        subheading( "Building the PDE_Integrator" );
        code(
"PDE_Integrator<dim,Element>  transient_pressure;\n"
"\n"
"// LHS: conductance [C]\n"
"NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(\n"
"    model.Database(), \"conductivity\",\n"
"    \"fluid pressure\", \"fluid pressure\" );\n"
"\n"
"// LHS: storage [S]/dt\n"
"NumIntegral_NT_lhsop_N_dV<dim>  capacitance_lhs(\n"
"    model.Database(), \"storativity\",\n"
"    \"fluid pressure\", \"fluid pressure\" );\n"
"capacitance_lhs.LumpedFormulation(true);\n"
"capacitance_lhs.MultiplyWithTimeIncrement(true);\n"
"\n"
"// RHS: initial conditions [S]/dt * p^t\n"
"NumIntegral_NT_rhsop_N_dV<dim>  capacitance_rhs(\n"
"    model.Database(), \"storativity\", \"fluid pressure\" );\n"
"capacitance_rhs.LumpedFormulation(true);\n"
"capacitance_rhs.MultiplyWithTimeIncrement(true);\n"
"\n"
"// RHS: source term (zero in this example)\n"
"NumIntegral_NT_rhsop_N_dV<dim>  source(\n"
"    model.Database(), \"fluid volume source\", \"fluid pressure\" );\n"
"source.AddAccumulateLater();\n"
"source.LumpedFormulation(true);\n"
"\n"
"transient_pressure.Add( &conductance );\n"
"transient_pressure.Add( &capacitance_lhs );\n"
"transient_pressure.Add( &capacitance_rhs );\n"
"transient_pressure.Add( &source );"
        );
        pause();

        // ====================================================================
        // Declare solver and integrator HERE -- outside any subheading scope
        // so they are visible to Step 9 and Step 10
        // ====================================================================

#ifdef CSMP_WITH_SAMG_SOLVER
        SAMG_Solver  samg_solver;
        PDE_Integrator<dim,Element>  transient_pressure( samg_solver );
#else
        CSMP_DEFAULT_LINEAR_SOLVER   linear_solver;
        PDE_Integrator<dim,Element>  transient_pressure( linear_solver );
#endif

        // PDE operators -- declared here so they outlive the subheading scopes
        NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(
            model.Database(), "conductivity",
            "fluid pressure", "fluid pressure" );

        NumIntegral_NT_lhsop_N_dV<dim>  capacitance_lhs(
            model.Database(), "storativity",
            "fluid pressure", "fluid pressure" );
        capacitance_lhs.LumpedFormulation(true);
        capacitance_lhs.MultiplyWithTimeIncrement(true);

        NumIntegral_NT_rhsop_N_dV<dim>  capacitance_rhs(
            model.Database(), "storativity", "fluid pressure" );
        capacitance_rhs.LumpedFormulation(true);
        capacitance_rhs.MultiplyWithTimeIncrement(true);

        NumIntegral_NT_rhsop_N_dV<dim>  source(
            model.Database(), "fluid volume source", "fluid pressure" );
        source.AddAccumulateLater();
        source.LumpedFormulation(true);

        transient_pressure.Add( &conductance );
        transient_pressure.Add( &capacitance_lhs );
        transient_pressure.Add( &capacitance_rhs );
        transient_pressure.Add( &source );

        output( "PDE_Integrator built with 4 operators." );
        pause();

        keypoint(
"Backward Euler: storage appears on BOTH LHS and RHS.\n"
"LHS capacitance_lhs: implicit term (unknown p^{t+dt})\n"
"RHS capacitance_rhs: explicit term (known p^t)\n"
"MultiplyWithTimeIncrement(true): operator stores total systems compressibility, beta_t.\n"
"                                 integrator divides by dt.\n"
"AddAccumulateLater(): source term is added to RHS after initial conditions."
        );
        pause();


        // ====================================================================
        // Step 9: The time stepping loop
        // ====================================================================
        heading( 9, 10, "The time stepping loop" );

        explain(
"The time stepping loop:\n"
"  1. Sets the time increment on the integrator\n"
"  2. Calls model.Apply() to assemble and solve\n"
"  3. Monitors variable ranges\n"
"  4. Writes output at user-specified frequency\n"
"  5. Increases the time step geometrically\n"
"\n"
"The geometric time step increase (factor 1.2 per step) is a\n"
"common strategy for diffusion problems:\n"
"  - Early times: small steps needed (sharp gradients)\n"
"  - Late times: large steps acceptable (smooth solution)\n"
"  - Total steps: ~log(t_max/t_0) / log(1.2) ~ 200 steps\n"
"    to cover 4 orders of magnitude in time"
        );
        pause();

        subheading( "Time stepping loop structure" );
        code(
"const double day           = 86400.;          // 1 day in seconds\n"
"double       model_time    = 0.;\n"
"double       maxtime       = 100. * 365. * day; // 100 years\n"
"double       time_increment = 5.;              // start: 5 seconds\n"
"size_t       timestep      = 1;\n"
"size_t       save_counter  = 0;\n"
"size_t       save_frequency = 5;              // user-defined\n"
"\n"
"while ( model_time <= maxtime )\n"
"{\n"
"    // 1. Set time increment and solve\n"
"    transient_pressure.TimeIncrement( 1. / time_increment );\n"
"    model.Apply( transient_pressure );\n"
"\n"
"    // 2. Monitor\n"
"    printRangeOfVariable( model, \"fluid pressure\" );\n"
"    monitor.ScalarPropertyRanges( model, model_time );\n"
"\n"
"    // 3. Output at user-specified frequency\n"
"    if ( save_counter == save_frequency ) {\n"
"        vtu_output.OutputDataToVTU(\n"
"            (model_name + \"_output\").c_str(),\n"
"            runtime_variables,\n"
"            model.Region(\"Model\"),\n"
"            static_cast<long>(model_time/day) );\n"
"        save_counter = 0;\n"
"    }\n"
"\n"
"    // 4. Advance time and increase step\n"
"    model_time      += time_increment;\n"
"    time_increment  *= 1.2;   // geometric increase\n"
"    timestep++;\n"
"    save_counter++;\n"
"}"
        );
        pause();

        subheading( "Why model.Apply() not IntegrateOver()?" );
        explain(
"In previous examples we used:\n"
"    integrator.IntegrateOver( model, gref );\n"
"\n"
"Here we use:\n"
"    model.Apply( transient_pressure );\n"
"\n"
"Both work, but model.Apply() is preferred for this simple model:\n"
"  - Only the RHS needs to be rebuilt each step -- the LHS matrix\n"
"    structure is retained (RetainGlobalSolutionMatrix)\n"
"  - This is significantly faster for many time steps\n"
"\n"
"IntegrateOver() is more versatile, see advanced examples."
        );
        pause();

        keypoint(
"Use either model.Apply() or IntegrateOver() for solves. The latter has more features.\n"
"Geometric time step increase: dt *= 1.2 each step.\n"
"This covers 6 orders of magnitude in ~200 steps."
        );
        pause();


        // ====================================================================
        // Step 10: Running the simulation and results
        // ====================================================================
        heading( 10, 10, "Running the simulation and interpreting results" );

        explain(
"We now run the full transient simulation.\n"
"You will be asked how often to save output files.\n"
"\n"
"The simulation runs from t=0 to t=100 years,\n"
"but the interesting physics happens in the first few years:\n"
"  - Days to weeks: signal enters the sandstone layers\n"
"  - Months to years: signal penetrates tens of metres\n"
"  - Decades: signal reaches the dolomite at depth\n"
"  - Centuries: full equilibration in clay layers\n"
"\n"
"Watch the pressure range -- it should start at 104000 Pa\n"
"and gradually approach 96000 Pa as the signal diffuses in."
        );
        pause();

        subheading( "Setting up initial and boundary conditions" );
        {
            constexpr double p_high = 104000.;  // Pa (1040 hPa)
            constexpr double p_low  =  96000.;  // Pa ( 960 hPa)

            model.InputPropertyValue( "fluid pressure", makeScalar(ANY, p_high) );
            model.Boundary("IRREGULAR").InputPropertyValue( "fluid pressure", makeScalar(DIRICH, p_low) );

            output( "Initial pressure = " + dbl_str(p_high/100.,0) + " hPa everywhere." );
            output( "Surface BC: DIRICH p = " + dbl_str(p_low/100.,0) + " hPa on IRREGULAR." );
            printRangeOfVariable( model, "fluid pressure" );
        }
        pause();

        subheading( "How often to save output?" );
        size_t save_frequency{ 5 };
        cout << "\n  Enter save frequency (1=every step, 5=every 5th step,\n"
             << "  recommended: 5): ";
        cout.flush();
        cin >> save_frequency;
        if ( save_frequency < 1 ) save_frequency = 1;
        output( "Saving every " + to_string(save_frequency) + " steps." );
        pause();

        subheading( "Running the transient simulation..." );
        {
            // Output setup
            list<string>  material_props{ "porosity", "permeability", "storativity" };
            list<string>  runtime_variables{ "fluid pressure" };

            VTU_Interface<dim>  vtu_output( model );

            // Write material properties once
            vtu_output.OutputDataToVTU(
                (string(model.Name()) + "_material_properties").c_str(),
                material_props,
                model.Region("Model"), 0 );
            output( "Written: " + model_name + "_material_properties_0.vtu" );

            // Region monitor
            string first_integral_property{ "volume flux" };
            string first_range_property{ "fluid pressure" };
            RegionMonitor<dim>  monitor( model,
                                         first_integral_property,
                                         first_range_property );

            // Time loop parameters
            constexpr double day                = 86400.;            //   1 day in seconds, max-time 100 years
            constexpr double max_time_increment =  14. * day;        //   1 fortnight
            constexpr double maxtime            = 100. * 365. * day; // 100 years
            double       model_time     = 0.;
            double       time_increment = 5.;   // start: 5 seconds
            size_t       timestep       = 1;
            size_t       save_counter   = 0;

            output( "Starting time loop..." );
            blank();

            while ( model_time <= maxtime )
            {
                // 1. Solve - divide the storage (=capacitance) matrices by the time increment as in the Backward-Euler discretisation shown above
                transient_pressure.TimeIncrement( 1. / time_increment );
                model.Apply( transient_pressure );

                // 2. Monitor
                monitor.ScalarPropertyRanges( model, model_time );

                // 4. Screen output at key times
                const double t_days = model_time / day;
                if ( timestep <= 5 ||
                     ( t_days >= 1.   && t_days < 1.1   ) ||
                     ( t_days >= 10.  && t_days < 10.1  ) ||
                     ( t_days >= 100. && t_days < 100.1 ) ||
                     ( t_days >= 365. && t_days < 365.1 ) )
                {
                    cout << "\n  t = " << fixed << setprecision(2)
                         << t_days << " days  (dt = "
                         << dbl_str(time_increment,1) << " s)\n";
                    printRangeOfVariable( model, "fluid pressure" );
                }

                // 5. Save output
                if ( save_counter == save_frequency ) {
                    vtu_output.OutputDataToVTU(
                        (string(model.Name()) + "_output").c_str(),
                        runtime_variables,
                        model.Region("Model"),
                        static_cast<long>(model_time / day) );
                    save_counter = 0;
                    output( "Written output at t = "
                            + dbl_str(model_time/day,1) + " days" );
                }

                // 6. Advance time
                model_time += time_increment;
                if ( time_increment < max_time_increment ) time_increment  *= 1.2;
                timestep++;
                save_counter++;
            }

            // Final output
            vtu_output.OutputDataToVTU(
                (string(model.Name()) + "_output").c_str(),
                runtime_variables,
                model.Region("Model"),
                static_cast<long>(model_time / day) );

            monitor.Out( model.Name() );

            output( "Simulation complete." );
            output( "Final time: " + dbl_str(model_time/(365.*day),1) + " years" );
            output( "Total timesteps: " + to_string(timestep) );
        }
        pause();

        subheading( "Interpreting the results in ParaView" );
        explain(
"Open the output files in ParaView:\n"
"\n"
"  Jura-slope1_material_properties_0.vtu\n"
"    -- shows permeability, porosity, storativity\n"
"    -- note the 7-order-of-magnitude contrast in permeability\n"
"\n"
"  Jura-slope1_output_*.vtu  (time series)\n"
"    -- 'fluid pressure': CSMP transient solution\n"
"\n"
"What to look for:\n"
"  1. The pressure drop enters the dolomite and sandstone layers first\n"
"     (high diffusivity -- signal propagates quickly)\n"
"  2. The clay layers act as barriers\n"
"     (low diffusivity -- signal barely penetrates)\n"
"  3. The dolomite equilibrates rapidly at depth\n"
"\n"
"Use ParaView's 'Plot Over Line' filter to extract a 1D profile\n"
"of 'fluid pressure'."
        );
        pause();

        subheading( "Final variable ranges" );
        printRangeOfVariable( model, "fluid pressure" );
        printRangeOfVariable( model, "storativity" );
        pause();

    } // end model scope


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 6: Atmospheric Pressure Diffusion\n\n"
         << "  Steps covered:\n\n"
         << "    1.  The pressure diffusion equation\n"
         << "    2.  The atmospheric pressure signal (Southern Germany)\n"
         << "    3.  The analytical solution (error function)\n"
         << "    4.  Stage 1: 1D verification table\n"
         << "    5.  Hydraulic diffusivity of the Jura formations\n"
         << "    6.  Stage 2: loading the Jura-slope1 model\n"
         << "    7.  Computing storativity and conductivity\n"
         << "    8.  Setting up the transient PDE_Integrator\n"
         << "    9.  The time stepping loop\n"
         << "   10.  Results and interpretation\n\n"
         << "  Key concepts:\n\n"
         << "    - Backward Euler: storage on BOTH LHS and RHS\n"
         << "    - MultiplyWithTimeIncrement(true): operator stores beta_t\n"
         << "    - model.Apply() reuses matrix structure in time loops\n"
         << "    - Geometric time step: dt *= 1.2 covers 6 orders of magnitude\n"
         << "    - 'analytical pressure' variable enables direct VTU comparison\n"
         << "    - FIELD_DATA flag protects material properties\n"
         << "    - 7-order-of-magnitude contrast in D across Jura formations\n\n"
         << "  Next: Example 7 -- Mesh convergence study\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nAtmosphericPressureDiffusion_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp

