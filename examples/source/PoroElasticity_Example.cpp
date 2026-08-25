/*
 *  PoroElasticity_Example.cpp
 *  opencsmp
 *
 *  Biot poromechanics — coupled displacement and fluid pressure.
 *  Demonstrates block-structured assembly using PDE_Integrator.
 */

#include "PoroElasticity_Example.h"

#include "Region.h"
#include "Model.h"
#include "Boundary.h"
#include "PDE_Integrator.h"
#include "VTU_Interface.h"

// FE operators — mechanical block
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dV.h"

// FE operators — fluid block
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"

// FE operators — Biot coupling (Jacobian cross-terms)
#include "NumJacobianIntegral_BT_m_N_dV.h"
#include "NumJacobianIntegral_N_mT_B_dV.h"

// post-processing
#include "StressesAndStrains.h"

// mesh and model I/O
#include "ANSYS_Model2D.h"

// utilities
#include "CSMP_definitions.h"
#include "CSMP_highLevelUtilities.h"
#include "LinearSolver.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#endif

using namespace std;

namespace csmp {

// ----------------------------------------------------------------------------
//  Specifications
// ----------------------------------------------------------------------------

/**
  @section What this Example Demonstrates

  @subsection 1. Mandel-Cryer Effect

  When a load is applied suddenly to a saturated porous medium, the pore pressure initially rises above the applied stress before dissipating.
  This counterintuitive behaviour — pore pressure exceeding the applied load — is a hallmark of coupled poromechanics and cannot be captured by uncoupled (Terzaghi) consolidation theory.

  To demonstrate this the user would:

  - Apply a vertical compressive stress on the TOP boundary as a Neumann traction
  - Fix the bottom boundary (roller)
  - Set the lateral boundaries to zero flux (no drainage)
  - Observe that interior pore pressure initially exceeds the applied stress

  @subsection 2. Undrained Loading — Hydrostatic to Near-Lithostatic Pore Pressure

  This is the scenario you describe. When loading is fast relative to the drainage timescale (undrained conditions),
  the pore pressure response approaches:
  
  Δp = B · Δσ_mean
  
  where B is the Skempton coefficient:
  
  B = α M / (K_dr + α² M)

  For α = 1 and a stiff skeleton (K_dr >> 1/M), B → 1 and the pore pressure increment equals the mean stress increment
  — approaching lithostatic pressure under full overburden loading.

  Set this up by:
  
  @code
    // overburden stress [Pa] — e.g. 2 km of rock at 2500 kg/m³
    constexpr double sigma_v = 2500. * 9.81 * 2000.;  // ~49 MPa

    // apply as Neumann traction on TOP boundary
    model.InputBoundaryValue( TOP,
        "Neumann stress",
        makeVector( NEUMANN,NEUMANN, 0., -sigma_v ) );

    // undrained: no drainage boundaries — all boundaries zero flux
    // (do NOT apply Dirichlet pressure conditions)

    // very low permeability to simulate undrained response
    model.InputPropertyValue( "Biot modulus", makeScalar( PLAIN, 1.0e-10 ) );
  @endcode
  The expected result is that pore pressure throughout the domain rises to approximately B · σ_v.
  
  @subsection 3. Consolidation — Drainage Over Time
  
  By adding the transient storage term (following the PressureDiffusion_Example pattern) and opening a drainage boundary,
  you can observe:

  - Initial undrained state: p ≈ B · σ_v everywhere
  - Progressive drainage from the open boundary
  - Final drained state: p = 0 (or hydrostatic), full stress carried by skeleton
  
  The timescale of consolidation is controlled by the hydraulic diffusivity:
  
    c_v = k / (μ (1/M + α²/K_dr))

  @subsection 4. Stress Partitioning Between Skeleton and Fluid
  
  By varying the Biot coefficient α between 0 (no coupling) and 1 (full coupling) the user can observe how the load is partitioned
   0	    No coupling — pore pressure unaffected by loading
  0.5	Partial coupling — moderate pore pressure rise
  1.0	Full coupling — maximum pore pressure response.
  
  @subsection 5. Failure Assessment
 
  Combined with the MohrCoulombFailure_Visitor, the user can assess whether the elevated pore pressure triggers shear or tensile failure:
  
  @code
  // after solving the coupled system
  MohrCoulombFailure_Visitor<2U> failure( model, ELEMENT_INTEGRATION_POINT );
  model.Region( "Model" ).Accept( failure );
  @endcode
  
  High pore pressure reduces the effective confining stress, moving the Mohr circle toward the failure envelope. This is directly relevant to:

  - Induced seismicity from fluid injection
  - Slope instability after rainfall
  - Wellbore stability during drilling
  
 */
void PoroElasticity_Example::Specifications()
{
    SetTitle( "Biot poromechanics — coupled displacement and fluid pressure" );
    SetDifficulty( 4 );
    SetCategory( "Simulation of Physical Processes" );
    AddAuthor( "SKM" );
    AddDescription( "2D quasi-static Biot consolidation in plane strain" );
    AddDescription( "Block-structured coupled FE assembly via PDE_Integrator" );
    AddDescription( "Displacement (vector) coupled with fluid pressure (scalar)" );
    AddDescription( "Dirichlet elimination of essential boundary conditions" );
    AddDescription( "Source: PoroElasticity_Example.cpp" );
    AddRequirement( "2D mesh file compatible with ANSYS_Model2D" );
    AddRequirement( "Variable file: PoroElasticity_Example-variables.txt" );
}


// ----------------------------------------------------------------------------
//  Run
// ----------------------------------------------------------------------------

/**
    @brief Assembles and solves the coupled Biot poromechanics system.

    The governing block system is:

        | K    Q  | { u }   { f_u }
        |         |       =
        | Qᵀ   S  | { p }   { f_p }

    Assembly proceeds element by element. Each element contributes to one
    or more blocks depending on which solution variables it carries. The
    PDE_Integrator framework handles the block structure automatically via
    the basic and test operand keys supplied to each operator.

    @par Block assembly via PDE_Integrator

    The key insight is that the PDE_Integrator assembles into the correct
    block of the global matrix based on the (basic, test) operand pair:

    | (basic, test)              | Block assembled |
    |----------------------------|----------------|
    | (displacement, displacement) | K — stiffness |
    | (fluid pressure, fluid pressure) | S — storage |
    | (fluid pressure, displacement) | Q — coupling |
    | (displacement, fluid pressure) | Qᵀ — coupling |

    @par Biot coupling sign convention

    The coupling matrices Q and Qᵀ use `SubtractAccumulate()` because
    the Biot effective stress principle gives:

        σ' = σ - α p I

    so pore pressure reduces the effective compressive stress. The negative
    sign is absorbed into the accumulation mode of the operator.

    @par DOF ordering

    Displacement DOFs are interleaved: for node n, DOFs are [2n, 2n+1]
    for [u_x, u_y]. Pressure DOFs follow after all displacement DOFs in
    the global system. The PDE_Integrator handles this ordering internally
    via the DOF index maps built during EstablishMatrixSetup.
*/
void PoroElasticity_Example::Run()
{
    // =========================================================================
    //  Model parameters
    // =========================================================================

    // elastic properties (plane strain, carbonate-like)
    constexpr double E   = 10000.;   // Young's modulus [Pa]
    constexpr double nu  = 0.25;     // Poisson's ratio [-]

    // Biot poromechanical parameters
    constexpr double alpha    = 1.0;     // Biot coefficient [-]
    constexpr double biot_M   = 1.0e-5; // 1/M — inverse Biot modulus [Pa⁻¹]
                                         // M = Ks / (alpha - phi + phi*Ks/Kf)

    // boundary conditions
    constexpr double u_x_right  =  0.01;  // prescribed horizontal displacement [m]
    constexpr double p_left     =  0.;    // drained left boundary [Pa]
    constexpr double p_right    =  500.;  // applied fluid pressure on right [Pa]

    // =========================================================================
    //  Model construction
    // =========================================================================

    //  Load the mesh and variable definitions.
    //  The variable file must define at minimum:
    //    "displacement"      — vector, node
    //    "fluid pressure"    — scalar, node
    //    "Young's modulus"   — scalar, element
    //    "Poisson's ratio"   — scalar, element
    //    "Biot alpha"        — scalar, element
    //    "Biot modulus"      — scalar, element  (stores 1/M)
    //    "gravity term"      — vector, element  (body force = rho*g)
    //    "fluid volume source" — scalar, element
    //    "stress"            — tensor, ipoint
    //    "strain"            — tensor, ipoint
    //    "mean stress"       — scalar, ipoint

    string model_name;
    cout << "\nEnter mesh name (or ENTER for default 'poroelasticity_tiny'): ";
    cin.ignore();
    getline( cin, model_name );
    if ( model_name.empty() ) model_name = "poroelasticity_tiny";

    const string variable_file = "PoroElasticity_Example-variables.txt";
    const string file_name     = GetExampleFileName( __FILE__ );
    CreateWorkingDirectoryAndCopyInputModelFiles( file_name, model_name,
                                                  variable_file );

    Model<2U> model( model_name, variable_file );
    printModelDimensions( model, true );

    // =========================================================================
    //  Material properties
    // =========================================================================

    // assign uniform material properties to the entire model domain
    model.InputPropertyValue( "Young's modulus",    makeScalar( PLAIN, E       ) );
    model.InputPropertyValue( "Poisson's ratio",    makeScalar( PLAIN, nu      ) );
    model.InputPropertyValue( "Biot alpha",         makeScalar( PLAIN, alpha   ) );
    model.InputPropertyValue( "Biot modulus",       makeScalar( PLAIN, biot_M  ) );

    // zero body force and fluid source — modify here to add gravity or injection
    model.InputPropertyValue( "gravity term",       makeVector( PLAIN,PLAIN, 0.,0. ) );
    model.InputPropertyValue( "fluid volume source",makeScalar( PLAIN, 0.          ) );

    // initialise solution variables to zero
    model.InputPropertyValue( "displacement",   makeVector( PLAIN,PLAIN, 0.,0. ) );
    model.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.          ) );
    model.InputPropertyValue( "mean stress",    makeScalar( PLAIN, 0.          ) );

    // =========================================================================
    //  Boundary conditions
    // =========================================================================

    // --- displacement ---
    // Left boundary: fully fixed (no horizontal or vertical displacement)
    model.InputBoundaryValue( LEFT,
        "displacement", makeVector( DIRICH,DIRICH, 0., 0. ) );

    // Bottom boundary: roller — no vertical displacement, free horizontally
    model.InputBoundaryValue( BOTTOM,
        "displacement", makeVector( PLAIN,DIRICH, 0., 0. ) );

    // Right boundary: prescribed horizontal displacement, free vertically
    model.InputBoundaryValue( RIGHT,
        "displacement", makeVector( DIRICH,PLAIN, u_x_right, 0. ) );

    // --- fluid pressure ---
    // Left boundary: drained (zero pressure)
    model.InputBoundaryValue( LEFT,
        "fluid pressure", makeScalar( DIRICH, p_left ) );

    // Right boundary: applied fluid pressure
    model.InputBoundaryValue( RIGHT,
        "fluid pressure", makeScalar( DIRICH, p_right ) );

    // =========================================================================
    //  Diagnostic output of input state
    // =========================================================================

    printRangeOfVariable( model, "Young's modulus"  );
    printRangeOfVariable( model, "Poisson's ratio"  );
    printRangeOfVariable( model, "Biot alpha"       );
    printRangeOfVariable( model, "Biot modulus"     );
    printRangeOfVariable( model, "displacement"     );
    printRangeOfVariable( model, "fluid pressure"   );

    // =========================================================================
    //  Solver
    // =========================================================================

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_napproach( 2 );
    SAMG_Solver solver( &settings );
    PDE_Integrator<2U,Element> poromechanics( solver );
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    PDE_Integrator<2U,Element> poromechanics( solver );
#endif

    // =========================================================================
    //  FE operators
    // =========================================================================

    // --- Mechanical stiffness K: ∫ Bᵀ D(E,ν) B dV ---
    // Assembles the plane-strain elastic stiffness matrix.
    // (basic=displacement, test=displacement) → K block
    NumIntegral_BT_D_B_dV<2U> lhs_K( model.Database(),
                                       "Young's modulus",
                                       "Poisson's ratio",
                                       "displacement",
                                       "displacement" );
    poromechanics.Add( &lhs_K );

    // --- Storage matrix S: ∫ Nᵀ (1/M) N dV ---
    // Assembles the fluid compressibility / storage matrix.
    // (basic=fluid pressure, test=fluid pressure) → S block
    // Lumped formulation is appropriate when 1/M is element-constant.
    NumIntegral_NT_lhsop_N_dV<2U> lhs_S( model.Database(),
                                           "Biot modulus",
                                           "fluid pressure",
                                           "fluid pressure" );
    lhs_S.LumpedFormulation( true );
    poromechanics.Add( &lhs_S );

    // --- Biot coupling Q: -α ∫ Bᵀ m N dV ---
    // Assembles the upper off-diagonal coupling block.
    // m = [1,1,0]ᵀ extracts volumetric strain from the strain vector.
    // SubtractAccumulate: pore pressure reduces effective compressive stress
    //   σ' = σ - α p I  →  negative sign on coupling term.
    // (basic=fluid pressure, test=displacement) → Q block
    NumJacobianIntegral_BT_m_N_dV<2U> lhs_Q( model.Database(),
                                               "Biot alpha",
                                               "fluid pressure",
                                               "displacement" );
    lhs_Q.SubtractAccumulate();
    poromechanics.Add( &lhs_Q );

    // --- Biot coupling Qᵀ: -α ∫ Nᵀ mᵀ B dV ---
    // Assembles the lower off-diagonal coupling block (transpose of Q).
    // For a symmetric system, Qᵀ = Qᵀ exactly.
    // SubtractAccumulate: same sign convention as Q.
    // (basic=displacement, test=fluid pressure) → Qᵀ block
    NumJacobianIntegral_N_mT_B_dV<2U> lhs_QT( model.Database(),
                                                "Biot alpha",
                                                "displacement",
                                                "fluid pressure" );
    lhs_QT.SubtractAccumulate();
    poromechanics.Add( &lhs_QT );

    // --- Mechanical body force f_u: ∫ Pᵀ f dV ---
    // Integrates the gravitational body force over the element volume.
    // Set "gravity term" = rho*g vector to activate; zero by default.
    // (test=displacement) → f_u vector
    NumIntegral_PT_op_dV<2U> rhs_fu( model.Database(),
                                      "gravity term",
                                      "displacement" );
    poromechanics.Add( &rhs_fu );

    // --- Fluid source f_p: ∫ Nᵀ q N dV ---
    // Integrates the volumetric fluid source over the element volume.
    // Set "fluid volume source" > 0 for injection, < 0 for production.
    // (test=fluid pressure) → f_p vector
    NumIntegral_NT_rhsop_N_dV<2U> rhs_fp( model.Database(),
                                            "fluid volume source",
                                            "fluid pressure" );
    poromechanics.Add( &rhs_fp );

    // =========================================================================
    //  Post-processing: stress and strain invariants
    // =========================================================================

    // StressesAndStrains computes principal stresses (sigma1, sigma2, sigma3),
    // principal strains, mean stress and von Mises stress from the computed
    // displacement field. Results are stored at element integration points.
    constexpr bool output_principal_vectors    { true  };
    constexpr bool extrapolate_results_to_nodes{ false };
    StressesAndStrains<2U> postpro( model,
                                    "Young's modulus",
                                    "Poisson's ratio",
                                    "displacement",
                                    output_principal_vectors,
                                    extrapolate_results_to_nodes );
    poromechanics.AddPostProcess( &postpro );

    // =========================================================================
    //  Solve
    // =========================================================================

    // IntegrateOver assembles the full block system, applies Dirichlet
    // elimination, solves the reduced system, and back-substitutes the
    // solution into the model's nodal variables.
    model.Apply( poromechanics );

    // =========================================================================
    //  Diagnostic output of solution
    // =========================================================================

    printRangeOfVariable( model, "displacement"   );
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "mean stress"    );

    // =========================================================================
    //  VTU output
    // =========================================================================

    VTU_Interface<2U> vtu( model, "PoroElasticity" );
    vtu.OmitZeroInFileName( true );

    const set<string> output_variables {
        "displacement",    // nodal displacement vector [m]
        "fluid pressure",  // nodal fluid pressure [Pa]
        "stress",          // Cauchy stress tensor at integration points [Pa]
        "strain",          // strain tensor at integration points [-]
        "mean stress",     // isotropic stress p = (σ₁+σ₂+σ₃)/3 [Pa]
        "sigma1",          // maximum principal stress [Pa]
        "sigma2",          // minimum principal stress [Pa]
        "strain1",         // maximum principal strain [-]
        "strain2"          // minimum principal strain [-]
    };

    vtu.OutputDataToVTU( "PoroElasticity_result", output_variables, "Model", 0 );
    
 
    // =========================================================================
    // --- interactive exploration ---
    // =========================================================================
 
    cout << "\nEnter overburden stress [Pa] (e.g. 49e6 for 2 km depth): ";
    double sigma_v;
    cin >> sigma_v;

    cout << "\nEnter Biot coefficient alpha [0-1]: ";
    double alpha_user;
    cin >> alpha_user;

    model.InputPropertyValue( "Biot alpha", makeScalar( PLAIN, alpha_user ) );
    model.InputBoundaryValue( TOP,
        "Neumann stress",
        makeVector( NEUMANN,NEUMANN, 0., -sigma_v ) );

    // re-solve
    model.Apply( poromechanics );

    // report Skempton B and expected undrained pressure
    const double K_dr    = E / ( 3. * ( 1. - 2.*nu ) );
    const double biot_M_val = /* read from model */ 1.0e-5;
    const double B_skempton = alpha_user * biot_M_val
                            / ( 1./K_dr + alpha_user*alpha_user * biot_M_val );
    cout << "\nSkempton B = " << B_skempton << "\n";
    cout << "Expected undrained Δp = " << B_skempton * sigma_v / 3. << " Pa\n";
    cout << "Lithostatic pressure  = " << sigma_v << " Pa\n";

    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "mean stress"    );

    cout << "\nPoroElasticity_Example: complete.\n";

    filesystem::current_path( "../../example_inputs/" );

} // end Run

} // csmp
