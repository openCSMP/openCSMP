// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SteadyStateHeatConduction_Example.h"
#include "exampleHelpers.h"

#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"
#include "CSMP_highLevelUtilities.h"
#include "InputDataManager.h"

// PDE integrator
#include "PDE_Integrator.h"

// LHS operator: conductance matrix
#include "NumIntegral_dNT_lhsop_dN_dV.h"

// RHS operators
#include "NumIntegral_NT_rhsop_N_dV.h"   // volumetric source term
#include "NumIntegral_NT_op_N_dS.h"      // boundary flux (Neumann)

// Solver
#include "LinearSolver.h"


using namespace std;

namespace csmp {

// ============================================================================
// Specifications
// ============================================================================

void SteadyStateHeatConduction_Example::Specifications()
{
    SetTitle( "Steady-state heat conduction (interactive)" );
    SetDifficulty( 2 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: SteadyStateHeatConduction_Example.cpp" );
    AddDescription( "Solves div(lambda grad T) = 0 in a layered geological "
                    "model with Dirichlet surface temperature and Neumann "
                    "basal heat flow. Introduces PDE_Integrator, LHS/RHS "
                    "operators, and boundary flux integrals." );
    AddRequirement( "fault_boundary_test (CSMP binary files, "
                    "-configuration.txt), heat_conduction-variables.txt" );
}


// ============================================================================
// Run
// ============================================================================


/**
    @file SteadyStateHeatConduction_Example.cpp

    @brief CSMP for Beginners -- Example 4: Your First PDE Solve

    @section purpose Purpose

    This interactive example introduces the complete workflow for solving
    a partial differential equation with CSMP. The physics is kept simple
    -- steady-state heat conduction -- so that the focus is entirely on
    the CSMP mechanics.

    @section physics Physics

    We solve the steady-state heat equation:

        div( lambda * grad T ) = 0

    where lambda is the thermal conductivity [W/(m K)] and T is
    temperature [K or degC].

    The model is a 3D layered geological cross-section (km scale):

        LAYER_TOP       -- porous shale,       lambda = 1.5 W/(m K)
        LAYER_RESERVOIR -- quartz sandstone,   lambda = 4.0 W/(m K)
        LAYER_BOTTOM    -- tight shale,        lambda = 2.0 W/(m K)
        NORMAL_FAULT    -- fault zone,         lambda = 1.0 W/(m K)

    Boundary conditions:

        TOP    -- Dirichlet: T = 15 degC  (mean surface temperature)
        BOTTOM -- Neumann:   q = 30 mW/m2 (basal heat flow)
        LEFT, RIGHT, FRONT, BACK -- no-flux (natural Neumann = 0)

    A local heat sink (cold fluid injection) is optionally added
    to demonstrate source terms.

    @section operators PDE Operators Used

        NumIntegral_dNT_lhsop_dN_dV -- conductance matrix [L] (LHS)
                                       integral of (grad N)^T * lambda * grad N dV
        NumIntegral_NT_rhsop_N_dV -- source term (RHS)
                                     integral of N^T * q * N dV
        NumIntegral_NT_op_N_dS    -- boundary flux (RHS, on Boundary faces)
                                     integral of N^T * q_n dS

    @section requirements Requirements

    - fault_boundary_test  (CSMP binary files + -configuration.txt)
    - heat_conduction-variables.txt  (variable definition file)

    @author S.K. Matthai
*/
void SteadyStateHeatConduction_Example::Run()
{
    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 4\n"
         << "  Your First PDE Solve: Steady-State Heat Conduction\n\n"
         << "  We will solve:  div( lambda * grad T ) = 0\n\n"
         << "  Model: layered geological cross-section (km scale)\n"
         << "    LAYER_TOP       -- porous shale\n"
         << "    LAYER_RESERVOIR -- quartz sandstone (most conductive)\n"
         << "    LAYER_BOTTOM    -- tight shale\n"
         << "    NORMAL_FAULT    -- fault zone\n\n"
         << "  Boundary conditions:\n"
         << "    TOP    -- Dirichlet: T = 15 degC\n"
         << "    BOTTOM -- Neumann:   q = 30 mW/m2 (basal heat flow)\n"
         << "    Sides  -- no-flux (natural Neumann = 0)\n\n"
         << "  Steps:\n"
         << "    1.  The heat equation and its FEM discretisation\n"
         << "    2.  Loading the model\n"
         << "    3.  Assigning thermal conductivities\n"
         << "    4.  Setting boundary conditions\n"
         << "    5.  Building the PDE_Integrator\n"
         << "    6.  The LHS operator: conductance matrix\n"
         << "    7.  The RHS operator: basal heat flow (Neumann)\n"
         << "    8.  Solving and checking the result\n"
         << "    9.  Adding a heat sink (source term)\n"
         << "   10.  Output and visualisation\n";
    divider();
    pause( "  [ Press ENTER to begin ]" );


    // -----------------------------------------------------------------------
    // Model setup
    // -----------------------------------------------------------------------
    constexpr uint32_t dim{ 3U };

    string model_name( "fault_boundary_test" );
    cout << "\n  Enter model name or press ENTER for default '"
         << model_name << "':\n";
    string user_input;
    getline( cin, user_input );
    if ( !user_input.empty() ) model_name = user_input;

    const string file_name     = GetExampleFileName( __FILE__ );
    const string variable_file = "SteadyStateHeatConduction_Example-variables.txt";

    CreateWorkingDirectoryAndCopyInputModelFiles( file_name, model_name,
                                                  variable_file, model_name );


    // ========================================================================
    // Step 1: The heat equation and its FEM discretisation
    // ========================================================================
    heading( 1, 10, "The heat equation and its FEM discretisation" );

    explain(
"We want to find the temperature field T(x,y,z) that satisfies:\n"
"\n"
"    div( lambda * grad T ) = 0    in the domain\n"
"\n"
"where lambda [W/(m K)] is the thermal conductivity.\n"
"\n"
"This is the steady-state heat equation -- no time derivative,\n"
"no fluid flow, no phase changes. Just conduction.\n"
"\n"
"The FEM weak form multiplies by a test function N and integrates:\n"
"\n"
"    integral( (grad N)^T * lambda * grad T dV ) = integral( N * q_n dS )\n"
"    \\___________________________________/         \\___________________/\n"
"              LHS (conductance matrix)                RHS (boundary flux)"
    );
    pause();

    subheading( "What each term means" );
    explain(
"LHS -- conductance matrix [L]:\n"
"  Assembled by NumIntegral_dNT_lhsop_dN_dV.\n"
"  For each element: L_ij = integral( (grad N_i)^T * lambda * grad N_j dV )\n"
"  This is the FEM equivalent of the Laplacian operator.\n"
"\n"
"RHS -- basal heat flow (Neumann condition):\n"
"  Assembled by NumIntegral_NT_op_dS on the BOTTOM boundary faces.\n"
"  For each face: f_i = integral( N_i * q_n dS )\n"
"  where q_n = 30 mW/m2 is the prescribed heat flux.\n"
"\n"
"Dirichlet condition (TOP boundary):\n"
"  T = 15 degC is enforced by flagging TOP nodes as DIRICH.\n"
"  The assembler eliminates these DOFs from the system."
    );
    pause();

    subheading( "The CSMP operator names" );
    code(
"// LHS: conductance matrix\n"
"NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(\n"
"    model.Database(),\n"
"    \"thermal conductivity\",  // material property (lambda)\n"
"    \"temperature\",           // basic operand (row)\n"
"    \"temperature\" );         // test operand (column)\n"
"\n"
"// RHS: volumetric heat source (optional)\n"
"NumIntegral_NT_rhsop_N_dV<dim>  heat_source(\n"
"    model.Database(),\n"
"    \"heat source\",           // source term variable\n"
"    \"temperature\" );         // test operand\n"
"\n"
"// RHS: boundary heat flux (Neumann)\n"
"NumIntegral_NT_op_dS<dim,Face>  basal_flux(\n"
"    model.Database(),\n"
"    \"basal heat flow\",       // flux variable on boundary faces\n"
"    \"temperature\" );         // test operand"
    );
    pause();

    keypoint(
"The LHS operator integrates over ELEMENT elements.\n"
"The Neumann RHS operator integrates over FACE elements on the boundary.\n"
"The Dirichlet condition is NOT an operator -- it is a flag on the variable.\n"
"This separation of concerns is the key design principle of CSMP."
    );
    pause();


    // ========================================================================
    // Step 2: Loading the model
    // ========================================================================
    heading( 2, 10, "Loading the model" );

    explain(
"The model 'fault_boundary_test' is a 3D layered geological cross-section.\n"
"It contains four volumetric regions and six box boundaries.\n"
"\n"
"We load it from CSMP binary files and a variable definition file\n"
"that specifies which physical variables exist and their valid ranges."
    );
    pause();

    subheading( "Loading now..." );
    cout << "  Loading: " << model_name << " + " << variable_file << "\n";
    cout.flush();

    Model<dim>  model( model_name, variable_file );

    output( "Model loaded successfully." );
    output( "Nodes:    " + to_string(model.Mesh().Nodes())    );
    output( "Elements: " + to_string(model.Mesh().Elements()) );
    output( "Faces:    " + to_string(model.Mesh().Faces())    );
    pause();

    subheading( "Checking regions and boundaries" );
    {
        const string regions[]{ "LAYER_TOP", "LAYER_RESERVOIR",
                                "LAYER_BOTTOM", "NORMAL_FAULT" };
        const string boundaries[]{ "TOP", "BOTTOM",
                                   "LEFT", "RIGHT", "FRONT", "BACK" };

        cout << "  Regions:\n";
        for ( const auto& r : regions ) {
            if ( model.ContainsRegion(r.c_str()) ) {
                Region<dim>& reg( model.Region(r.c_str()) );
                output( "  " + r + ": "
                        + to_string(reg.Nodes()) + " nodes, "
                        + to_string(reg.Cells()) + " elements" );
            } else {
                output( "  WARNING: region '" + r + "' not found!" );
            }
        }

        cout << "\n  Boundaries:\n";
        for ( const auto& b : boundaries ) {
            if ( model.ContainsBoundary(b.c_str()) ) {
                Boundary<dim>& bnd( model.Boundary(b.c_str()) );
                output( "  " + b + ": "
                        + to_string(bnd.Nodes()) + " nodes, "
                        + to_string(bnd.Cells()) + " faces" );
            } else {
                output( "  WARNING: boundary '" + b + "' not found!" );
            }
        }
    }
    pause();

    // Get the main computation domain
    Region<dim>&  gref( model.Region("Model") );

    keypoint(
"Verify with CountainsRegion(const char*) that the expected regions and boundaries exist\n"
"after loading. A missing region means the mesh file does not\n"
"match the variable file or the regions file."
    );
    pause();


    // ========================================================================
    // Step 3: Assigning thermal conductivities
    // ========================================================================
    heading( 3, 10, "Assigning thermal conductivities" );

    explain(
"Thermal conductivity lambda [W/(m K)] is an ELEMENT variable.\n"
"Different geological layers have different conductivities:\n"
"\n"
"  LAYER_TOP       -- porous shale:       1.5 W/(m K)\n"
"  LAYER_RESERVOIR -- quartz sandstone:   4.0 W/(m K)  (most conductive)\n"
"  LAYER_BOTTOM    -- tight shale:        2.0 W/(m K)\n"
"  NORMAL_FAULT    -- fault zone:         1.0 W/(m K)  (least conductive)\n"
"\n"
"We assign these using InputPropertyValue() on each region separately."
    );
    pause();

    subheading( "Assignment pattern" );
    code(
"// Create the variable (if not already in variable file):\n"
"PropertyHandle<dim>  lambda( model, \"thermal conductivity\",\n"
"                             SCALAR, ELEMENT );\n"
"\n"
"// Assign uniform value to each region:\n"
"model.Region(\"LAYER_TOP\").InputPropertyValue(\n"
"    \"thermal conductivity\", makeScalar(PLAIN, 1.5) );\n"
"\n"
"model.Region(\"LAYER_RESERVOIR\").InputPropertyValue(\n"
"    \"thermal conductivity\", makeScalar(PLAIN, 4.0) );\n"
"\n"
"// ... etc for other regions"
    );
    pause();

    subheading( "Live assignment" );
    {
        // NOTE: handle must be static, else variable will be deleted when it goes out of the subheading scope
        static PropertyHandle<dim>  lambda( model, "thermal conductivity", SCALAR, ELEMENT );

        // Physical thermal conductivities [W/(m K)]
        constexpr double lambda_top        = 1.5;  // porous shale
        constexpr double lambda_reservoir  = 4.0;  // quartz sandstone
        constexpr double lambda_bottom     = 2.0;  // tight shale
        constexpr double lambda_fault      = 1.0;  // fault zone

        model.Region("LAYER_TOP").InputPropertyValue(
            "thermal conductivity", makeScalar(PLAIN, lambda_top) );
        model.Region("LAYER_RESERVOIR").InputPropertyValue(
            "thermal conductivity", makeScalar(PLAIN, lambda_reservoir) );
        model.Region("LAYER_BOTTOM").InputPropertyValue(
            "thermal conductivity", makeScalar(PLAIN, lambda_bottom) );
        model.Region("NORMAL_FAULT").InputPropertyValue(
            "thermal conductivity", makeScalar(PLAIN, lambda_fault) );

        output( "Assigned thermal conductivities:" );
        output( "  LAYER_TOP       = " + dbl_str(lambda_top,1)       + " W/(m K)" );
        output( "  LAYER_RESERVOIR = " + dbl_str(lambda_reservoir,1) + " W/(m K)" );
        output( "  LAYER_BOTTOM    = " + dbl_str(lambda_bottom,1)    + " W/(m K)" );
        output( "  NORMAL_FAULT    = " + dbl_str(lambda_fault,1)     + " W/(m K)" );

        printRangeOfVariable( model, "Model", "thermal conductivity" );
    }
    pause();

    // Also assign zero heat source (needed by RHS operator)
    {
        static PropertyHandle<dim>  heat_src( model, "heat source", SCALAR, ELEMENT );
        model.InputPropertyValue( "heat source", makeScalar(PLAIN, 0.) );
        output( "Set 'heat source' = 0 everywhere (no volumetric heat generation)." );
    }
    pause();

    keypoint(
"Assign material properties BEFORE building the PDE_Integrator.\n"
"The operators read the property values during assembly.\n"
"Use InputPropertyValue() on individual regions for heterogeneous models.\n"
"Use model.InputPropertyValue() for uniform values across the whole model."
    );
    pause();


    // ========================================================================
    // Step 4: Setting boundary conditions
    // ========================================================================
    heading( 4, 10, "Setting boundary conditions" );

    explain(
"Two types of boundary conditions are used:\n"
"\n"
"  Dirichlet (TOP):  T = 15 degC -- fixed surface temperature\n"
"    --> flag the TOP boundary nodes as DIRICH\n"
"    --> the assembler eliminates these DOFs from the system\n"
"\n"
"  Neumann (BOTTOM): q_n = 30 mW/m2 -- prescribed basal heat flow\n"
"    --> store the flux value on BOTTOM boundary faces\n"
"    --> flag the faces as NEUMANN\n"
"    --> the boundary integral operator adds this to the RHS\n"
"\n"
"  No-flux (sides):  q_n = 0 -- natural Neumann condition\n"
"    --> this is the DEFAULT -- no action needed\n"
"    --> if no Neumann condition is specified, flux = 0 automatically"
    );
    pause();

    subheading( "Dirichlet condition: fixed surface temperature" );
    code(
"// Set T = 15 degC on all TOP boundary nodes:\n"
"model.Boundary(\"TOP\").InputPropertyValue( TOP, \"temperature\",\n"
"                          makeScalar(DIRICH, 15.) );\n"
"\n"
"// Verify:\n"
"printRangeOfVariable( model, \"TOP\", \"temperature\" );"
    );
    pause();

    subheading( "Neumann condition: basal heat flow" );
    code(
"// The basal heat flow is stored on BOTTOM boundary FACES.\n"
"// The NumIntegral_NT_op_dS operator reads this value and\n"
"// integrates it over each face to get the nodal heat flux.\n"
"\n"
"// 30 mW/m2 = 0.030 W/m2\n"
"constexpr double basal_heat_flow = 0.030;  // W/m2\n"
"\n"
"// Create the variable on FACE placement:\n"
"PropertyHandle<dim>  bhf( model, \"basal heat flow\",\n"
"                          SCALAR, FACE );\n"
"\n"
"// Set value and flag on BOTTOM faces:\n"
"model.Boundary(\"BOTTOM\").InputPropertyValue(\n"
"    \"basal heat flow\",\n"
"    makeScalar(NEUMANN, basal_heat_flow) );"
    );
    pause();

    subheading( "Live assignment of boundary conditions" );
    {
        // Initialise temperature to zero everywhere first
        static PropertyHandle<dim>  temperature( model, "temperature", SCALAR, NODE );
        model.InputPropertyValue( "temperature", makeScalar(PLAIN, 0.) );
        output( "Initialised temperature = 0 everywhere." );

        // Dirichlet: fixed surface temperature (assignment for BOX_BOUNDARY models)
        constexpr double T_surface = 15.0;  // degC
        model.InputBoundaryValue( TOP, "temperature", makeScalar(DIRICH, T_surface) );
        output( "Set DIRICH T = " + dbl_str(T_surface,1) + " degC on TOP boundary." );

        // Count Dirichlet nodes
        const csmp::Index T_key = model.Database().StorageKey("temperature");
        size_t n_dirich{ 0U };
        ScalarVariable sc;
        for ( const auto& nit : gref.NodeVector() ) {
            nit->Read( T_key, sc );
            if ( sc.Flag() == DIRICH ) ++n_dirich;
        }
        output( "Dirichlet nodes (TOP): " + to_string(n_dirich) );

        // Neumann: basal heat flow on BOTTOM faces
        constexpr double basal_heat_flow = 0.030;  // W/m2 = 30 mW/m2
        static PropertyHandle<dim>  bhf( model, "basal heat flow", SCALAR, FACE );
        model.Boundary("BOTTOM").InputPropertyValue( "basal heat flow", makeScalar(NEUMANN, basal_heat_flow) );
        output( "Set NEUMANN q = " + dbl_str(basal_heat_flow*1000.,1) + " mW/m2 on BOTTOM boundary." );

        // Verify
        printRangeOfVariable( model, "TOP",    "temperature" );
        printRangeOfVariable( model, "BOTTOM", "basal heat flow" );
    }
    pause();

    keypoint(
"Dirichlet conditions: flag the NODE variable as DIRICH.\n"
"Neumann conditions:   flag the FACE variable as NEUMANN.\n"
"No-flux (natural Neumann = 0): do nothing -- it is the default.\n"
"The assembler reads these flags automatically during matrix assembly."
    );
    pause();


    // ========================================================================
    // Step 5: Building the PDE_Integrator
    // ========================================================================
    heading( 5, 10, "Building the PDE_Integrator" );

    explain(
"PDE_Integrator is the central object that:\n"
"  1. Assembles the global sparse matrix from element contributions\n"
"  2. Applies boundary conditions\n"
"  3. Calls the solver\n"
"  4. Writes the solution back to the mesh\n"
"\n"
"You build it by:\n"
"  (a) Creating a solver\n"
"  (b) Creating PDE operators (LHS and RHS)\n"
"  (c) Adding operators to the integrator\n"
"  (d) Calling IntegrateOver()"
    );
    pause();

    subheading( "The complete setup pattern" );
    code(
"// (a) Create solver:\n"
"CSMP_DEFAULT_LINEAR_SOLVER   solver;\n"
"PDE_Integrator<dim,Element>  heat_conductor( solver );\n"
"\n"
"// (b) Create LHS operator (conductance matrix):\n"
"NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(\n"
"    model.Database(),\n"
"    \"thermal conductivity\",  // lambda\n"
"    \"temperature\",           // basic operand\n"
"    \"temperature\" );         // test operand\n"
"\n"
"// (c) Create RHS operator (volumetric source):\n"
"NumIntegral_NT_rhsop_N_dV<dim>  heat_source(\n"
"    model.Database(),\n"
"    \"heat source\",\n"
"    \"temperature\" );\n"
"\n"
"// (c) Create boundary flux operator (Neumann):\n"
"NumIntegral_NT_op_N_dS<dim,Face>  basal_flux(\n"
"    model.Database(),\n"
"    \"basal heat flow\",\n"
"    \"temperature\" );\n"
"\n"
"// (d) Add operators:\n"
"heat_conductor.Add( &conductance );\n"
"heat_conductor.Add( &heat_source );\n"
"heat_conductor.AddBoundaryIntegral( &basal_flux );\n"
"\n"
"// (e) Solve:\n"
"heat_conductor.IntegrateOver( model, gref );"
    );
    pause();

    subheading( "Why AddBoundaryIntegral() not Add()?" );
    explain(
"Add() registers operators that integrate over VOLUME elements.\n"
"AddBoundaryIntegral() registers operators that integrate over FACE\n"
"elements on the model boundary.\n"
"\n"
"The integrator automatically finds which boundary faces are shared\n"
"with the computational domain and accumulates their contributions\n"
"into the global RHS vector.\n"
"\n"
"This is why the Neumann condition is so clean in CSMP:\n"
"you just store the flux value on the face and add the operator.\n"
"The integrator handles the rest."
    );
    pause();

    keypoint(
"PDE_Integrator separates the WHAT (operators) from the HOW (assembly).\n"
"You describe the physics via operators.\n"
"The integrator handles matrix assembly, BC application, and solving.\n"
"Add() = volume integral.  AddBoundaryIntegral() = surface integral."
    );
    pause();


    // ========================================================================
    // Step 6: The LHS operator -- conductance matrix
    // ========================================================================
    heading( 6, 10, "The LHS operator: conductance matrix" );

    explain(
"NumIntegral_dNT_lhsop_dN_dV assembles the conductance matrix:\n"
"\n"
"    L_ij = integral( (grad N_i)^T * lambda * grad N_j dV )\n"
"\n"
"This is the FEM discretisation of the Laplacian operator.\n"
"For heat conduction, lambda is the thermal conductivity.\n"
"For pressure diffusion, lambda would be the permeability / viscosity.\n"
"For elasticity, lambda would be the stiffness tensor.\n"
"\n"
"The same operator works for any scalar diffusion problem."
    );
    pause();

    subheading( "Operator arguments" );
    code(
"NumIntegral_dNT_lhs_dN_dV<dim>  conductance(\n"
"    model.Database(),\n"
"    \"thermal conductivity\",  // material property (lambda)\n"
"    \"temperature\",           // basic operand: the unknown\n"
"    \"temperature\" );         // test operand: the test function\n"
"\n"
"// For a coupled system (e.g. thermo-mechanical):\n"
"// basic operand = row variable (what we solve for)\n"
"// test operand  = column variable (what we test against)\n"
"// For a scalar problem, both are the same variable."
    );
    pause();

    subheading( "What happens inside the operator" );
    explain(
"For each element, the operator:\n"
"  1. Reads lambda from the element (ELEMENT variable)\n"
"  2. Evaluates grad N at each Gauss integration point\n"
"  3. Computes L_ij = sum_gp w_gp * (grad N_i)^T * lambda * grad N_j * detJ\n"
"  4. Stores the result in the local DenseMatrix LHS\n"
"\n"
"AssignToGlobal() then scatters LHS into the global CompressedRowMatrix\n"
"using the DOF_indexes_ mapping (which accounts for Dirichlet elimination)."
    );
    pause();

    keypoint(
"NumIntegral_dNT_lhsop_dN_dV is the most commonly used LHS operator in CSMP.\n"
"It works for any scalar diffusion problem by changing the material property.\n"
"The basic and test operands are the same for scalar problems.\n"
"For coupled systems, they can be different variables."
    );
    pause();


    // ========================================================================
    // Step 7: The RHS operator -- basal heat flow (Neumann)
    // ========================================================================
    heading( 7, 10, "The RHS operator: basal heat flow (Neumann)" );

    explain(
"NumIntegral_NT_op_N_dS assembles the boundary flux vector:\n"
"\n"
"    f_i = integral( N_i * q_n dS )\n"
"\n"
"where q_n is the prescribed normal heat flux [W/m2].\n"
"\n"
"This operator integrates over FACE elements, not volume elements.\n"
"It is added with AddBoundaryIntegral(), not Add().\n"
"\n"
"The operator only accumulates contributions from faces where the\n"
"material operand ('basal heat flow') is flagged NEUMANN.\n"
"Faces with other flags (PLAIN, ANY) are skipped."
    );
    pause();

    subheading( "Operator arguments" );
    code(
"NumIntegral_NT_op_N_dS<dim,Face>  basal_flux(\n"
"    model.Database(),\n"
"    \"basal heat flow\",  // material operand: the flux value\n"
"    \"temperature\" );   // test operand: the test function\n"
"\n"
"// The operator checks: if face.Status(\"basal heat flow\") == NEUMANN\n"
"//   --> accumulate N_i * q_n * dS into the RHS vector\n"
"// else\n"
"//   --> skip this face (no contribution)\n"
"\n"
"// This is why we flagged the BOTTOM faces as NEUMANN in step 4."
    );
    pause();

    subheading( "Physical meaning of 30 mW/m2" );
    explain(
"The basal heat flow of 30 mW/m2 is a typical value for stable\n"
"continental crust. It represents the heat conducted upward from\n"
"the mantle through the lithosphere.\n"
"\n"
"In a 1D model with uniform conductivity lambda = 2 W/(m K) and\n"
"basal heat flow q = 30 mW/m2, the temperature gradient is:\n"
"\n"
"    dT/dz = q / lambda = 0.030 / 2.0 = 0.015 K/m = 15 K/km\n"
"\n"
"So a 3 km deep model would have a temperature difference of ~45 K\n"
"between top and bottom (plus the 15 degC surface temperature)."
    );
    pause();

    keypoint(
"NumIntegral_NT_op_N_dS integrates over FACE elements.\n"
"It only accumulates where the material operand flag == NEUMANN.\n"
"The no-flux condition on the sides is automatic -- no operator needed.\n"
"Physical insight: 30 mW/m2 gives ~15 K/km geothermal gradient."
    );
    pause();


    // ========================================================================
    // Step 8: Solving and checking the result
    // ========================================================================
    heading( 8, 10, "Solving and checking the result" );

    explain(
"Now we assemble and solve the system. The call to IntegrateOver() does:\n"
"\n"
"  1. EstablishMatrixSetup() -- size the global matrix\n"
"  2. Accumulate()           -- assemble LHS from volume elements\n"
"  3. AccumulateBoundaryIntegrals() -- assemble RHS from boundary faces\n"
"  4. AssignEssentialConditions()   -- apply Dirichlet conditions\n"
"  5. Solve()                -- call the linear solver\n"
"  6. OutputResults()        -- write solution back to nodes\n"
"  7. PostProcess()          -- optional post-processing\n"
"\n"
"After the call, 'temperature' at every node contains the solution."
    );
    pause();

    subheading( "Solving now..." );
    {
#ifdef CSMP_WITH_SAMG_SOLVER
        SAMG_Settings  settings;
        settings.Set_nxtyp(0);
        settings.Set_ncgtyp(5);
        settings.Set_ndefault(40);
        settings.Set_eps(0.);  // solve to round-off
        SAMG_Solver                  solver( &settings );
#else
        CSMP_DEFAULT_LINEAR_SOLVER   solver;
#endif
        PDE_Integrator<dim,Element>  heat_conductor( solver );

        // LHS: conductance matrix
        NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(
            model.Database(),
            "thermal conductivity",
            "temperature",
            "temperature" );

        // RHS: volumetric heat source (zero in this case)
        NumIntegral_NT_rhsop_N_dV<dim>  heat_source(
            model.Database(),
            "heat source",
            "temperature" );

        // RHS: basal heat flow (Neumann on BOTTOM)
        NumIntegral_NT_op_N_dS<dim>  basal_flux(
            model.Database(),
            "basal heat flow",
            "temperature" );

        heat_conductor.Add( &conductance );
        heat_conductor.Add( &heat_source );
        heat_conductor.AddBoundaryIntegral( &basal_flux );

        output( "Assembling and solving..." );
        cout.flush();

        heat_conductor.IntegrateOver( model, gref );

        output( "Solve complete." );
    }
    pause();

    subheading( "Checking the result" );
    {
        output( "Temperature range in the model:" );
        printRangeOfVariable( model, "Model", "temperature" );

        output( "Temperature at TOP boundary (should be ~15 degC):" );
        printRangeOfVariable( model, "TOP", "temperature" );

        output( "Temperature at BOTTOM boundary:" );
        printRangeOfVariable( model, "BOTTOM", "temperature" );

        // Compute expected bottom temperature for a simple 1D estimate
        // Using average conductivity and model depth
        Point<dim> xyz_min, xyz_max;
        model.MinMaxCoordinates( xyz_min, xyz_max );
        const double depth_m = xyz_max[1] - xyz_min[1];  // y = vertical
        const double lambda_avg = 2.0;  // approximate average [W/(m K)]
        const double q_basal    = 0.030; // [W/m2]
        const double dT_expected = q_basal * depth_m / lambda_avg;
        const double T_bottom_expected = 15.0 + dT_expected;

        output( "Model depth: " + dbl_str(depth_m/1000.,2) + " km" );
        output( "Expected bottom temperature (1D estimate): "
                + dbl_str(T_bottom_expected,1) + " degC" );
        output( "  (actual will differ due to layering and fault)" );
    }
    pause();

    keypoint(
"IntegrateOver() does everything: assemble, apply BCs, solve, output.\n"
"Always check the result with printRangeOfVariable().\n"
"The TOP boundary should show T = 15 degC (Dirichlet enforced).\n"
"The BOTTOM temperature depends on the conductivity structure."
    );
    pause();


    // ========================================================================
    // Step 9: Adding a heat sink (source term)
    // ========================================================================
    heading( 9, 10, "Adding a heat sink (source term)" );

    explain(
"A volumetric heat source or sink is added via NumIntegral_NT_rhsop_N_dV.\n"
"This operator assembles:\n"
"\n"
"    f_i = integral( N_i * q_vol * N_j dV )\n"
"\n"
"where q_vol [W/m3] is the volumetric heat generation rate.\n"
"\n"
"A negative value represents a heat sink (cooling).\n"
"This could model cold fluid injection into the reservoir,\n"
"or radioactive heat generation in the basement.\n"
"\n"
"We add a heat sink in LAYER_RESERVOIR to simulate cold fluid injection."
    );
    pause();

    subheading( "Adding a heat sink in the reservoir" );
    code(
"// Set a heat sink in the reservoir layer:\n"
"// -1e-4 W/m3 = -0.1 mW/m3 (modest cooling)\n"
"model.Region(\"LAYER_RESERVOIR\").InputPropertyValue(\n"
"    \"heat source\", makeScalar(PLAIN, -1e-4) );\n"
"\n"
"// Re-solve with the source term active:\n"
"heat_conductor.Reset( false );  // false = keep operators\n"
"heat_conductor.IntegrateOver( model, gref );"
    );
    pause();

    subheading( "Live demonstration -- heat sink in reservoir" );
    {
        // Add heat sink in reservoir
        constexpr double heat_sink = -1e-4;  // W/m3 (cooling)
        model.Region("LAYER_RESERVOIR").InputPropertyValue(
            "heat source", makeScalar(PLAIN, heat_sink) );
        output( "Set heat sink = " + sci_str(heat_sink,1)
                + " W/m3 in LAYER_RESERVOIR." );

        // Re-apply boundary conditions (Reset clears them)
        constexpr double T_surface      = 15.0;
        constexpr double basal_heat_flow = 0.030;
        model.InputBoundaryValue( TOP, "temperature",
                                  makeScalar(DIRICH, T_surface) );
        model.Boundary("BOTTOM").InputPropertyValue(
            "basal heat flow", makeScalar(NEUMANN, basal_heat_flow) );

#ifdef CSMP_WITH_SAMG_SOLVER
        SAMG_Settings  settings2;
        settings2.Set_eps(0.);
        SAMG_Solver                  solver2( &settings2 );
#else
        CSMP_DEFAULT_LINEAR_SOLVER   solver2;
#endif
        PDE_Integrator<dim,Element>  heat_conductor2( solver2 );

        NumIntegral_dNT_lhsop_dN_dV<dim>  conductance2(
            model.Database(),
            "thermal conductivity", "temperature", "temperature" );
        NumIntegral_NT_rhsop_N_dV<dim>  heat_source2(
            model.Database(),
            "heat source", "temperature" );
        NumIntegral_NT_op_N_dS<dim>  basal_flux2(
            model.Database(),
            "basal heat flow", "temperature" );

        heat_conductor2.Add( &conductance2 );
        heat_conductor2.Add( &heat_source2 );
        heat_conductor2.AddBoundaryIntegral( &basal_flux2 );

        output( "Re-solving with heat sink..." );
        cout.flush();
        heat_conductor2.IntegrateOver( model, gref );
        output( "Solve complete." );

        output( "Temperature range with heat sink:" );
        printRangeOfVariable( model, "Model", "temperature" );
        output( "Temperature in LAYER_RESERVOIR:" );
        printRangeOfVariable( model, "LAYER_RESERVOIR", "temperature" );
        output( "  (should be cooler than without the sink)" );

        // Reset heat source to zero for clean state
        model.InputPropertyValue( "heat source", makeScalar(PLAIN, 0.) );
    }
    pause();

    keypoint(
"NumIntegral_NT_rhsop_N_dV adds a volumetric source/sink term.\n"
"Negative values = heat sink (cooling), positive = heat source (warming).\n"
"The operator is already registered -- just change the variable value\n"
"and re-solve. No need to rebuild the integrator."
    );
    pause();


    // ========================================================================
    // Step 10: Output and visualisation
    // ========================================================================
    heading( 10, 10, "Output and visualisation" );

    explain(
"The solution is now stored at every node as 'temperature'.\n"
"We output it to VTK files for visualisation in ParaView.\n"
"\n"
"We also compute the temperature gradient (heat flux direction)\n"
"and extrapolate it to nodes for smooth visualisation.\n"
"\n"
"Finally, we print a summary of the simulation."
    );
    pause();

    subheading( "Computing derived quantities" );
    code(
"// Temperature gradient at element barycentres:\n"
"model.CopyGradientOfProperty_A_To_B(\n"
"    \"temperature\",           // source: NODE\n"
"    \"temperature gradient\" ); // target: ELEMENT\n"
"\n"
"// Extrapolate gradient to nodes for smooth VTK output:\n"
"model.ExtrapolateCellToNodeProperty(\n"
"    \"temperature gradient\",\n"
"    \"nodal temperature gradient\" );\n"
"\n"
"// Heat flux = -lambda * grad T:\n"
"// (computed in a post-processing loop)"
    );
    pause();

    subheading( "VTK output" );
    code(
"VTK_Interface<dim>  vtk;\n"
"\n"
"// Write temperature field:\n"
"vtk.OutputDataToVTK( model, \"temperature\",\n"
"                     \"temperature\", 1 );\n"
"\n"
"// Write thermal conductivity:\n"
"vtk.OutputDataToVTK( model, \"thermal-conductivity\",\n"
"                     \"thermal conductivity\", 1 );\n"
"\n"
"// Write temperature gradient:\n"
"vtk.OutputDataToVTK( model, \"temperature-gradient\",\n"
"                     \"temperature gradient\", 1 );"
    );
    pause();

    subheading( "Live output" );
    {
        // Re-solve without heat sink for clean final result
        constexpr double T_surface       = 15.0;
        constexpr double basal_heat_flow = 0.030;

        model.InputPropertyValue( "temperature",  makeScalar(PLAIN,  0.) );
        model.InputPropertyValue( "heat source",  makeScalar(PLAIN,  0.) );
        model.InputBoundaryValue( TOP, "temperature",
                                  makeScalar(DIRICH, T_surface) );
        model.Boundary("BOTTOM").InputPropertyValue(
            "basal heat flow", makeScalar(NEUMANN, basal_heat_flow) );

#ifdef CSMP_WITH_SAMG_SOLVER
        SAMG_Settings  settings3;
        settings3.Set_eps(0.);
        SAMG_Solver                  solver3( &settings3 );
#else
        CSMP_DEFAULT_LINEAR_SOLVER   solver3;
#endif
        PDE_Integrator<dim,Element>  heat_conductor3( solver3 );

        NumIntegral_dNT_lhsop_dN_dV<dim>  conductance3(
            model.Database(),
            "thermal conductivity", "temperature", "temperature" );
        NumIntegral_NT_rhsop_N_dV<dim>  heat_source3(
            model.Database(),
            "heat source", "temperature" );
        NumIntegral_NT_op_N_dS<dim>  basal_flux3(
            model.Database(),
            "basal heat flow", "temperature" );

        heat_conductor3.Add( &conductance3 );
        heat_conductor3.Add( &heat_source3 );
        heat_conductor3.AddBoundaryIntegral( &basal_flux3 );

        output( "Final solve (no heat sink)..." );
        cout.flush();
        heat_conductor3.IntegrateOver( model, gref );
        output( "Solve complete." );

        // Compute temperature gradient
        PropertyHandle<dim>  T_grad( model, "temperature gradient",
                                     VECTOR, ELEMENT );
        PropertyHandle<dim>  T_grad_n( model, "nodal temperature gradient",
                                       VECTOR, NODE );
        model.CopyGradientOfProperty_A_To_B(
            "temperature", "temperature gradient" );
        model.ExtrapolateCellToNodeProperty(
            "temperature gradient", "nodal temperature gradient" );
        output( "Computed temperature gradient." );

        // Compute heat flux magnitude: q = lambda * |grad T|
        // and store as element scalar for visualisation
        PropertyHandle<dim>  heat_flux_mag( model, "heat flux magnitude",
                                            SCALAR, ELEMENT );
        {
            const csmp::Index lambda_key =
                model.Database().StorageKey("thermal conductivity");
            const csmp::Index grad_key   =
                model.Database().StorageKey("temperature gradient");
            ScalarVariable       lambda_sc;
            VectorVariable<dim>  grad_T;

            for ( auto& eit : gref.CellVector() ) {
                eit->Read( lambda_key, lambda_sc );
                eit->Read( grad_key,   grad_T    );
                const double q_mag = lambda_sc() * grad_T.Length();
                eit->Store( heat_flux_mag.Key(),
                            makeScalar(PLAIN, q_mag) );
            }
        }
        output( "Computed heat flux magnitude." );
        printRangeOfVariable( model, "Model", "heat flux magnitude" );

        // VTK output
        VTK_Interface<dim>  vtk;
        vtk.OutputDataToVTK( model, "temperature",
                             "temperature", 1 );
        vtk.OutputDataToVTK( model, "thermal-conductivity",
                             "thermal conductivity", 1 );
        vtk.OutputDataToVTK( model, "temperature-gradient",
                             "temperature gradient", 1 );
        vtk.OutputDataToVTK( model, "heat-flux-magnitude",
                             "heat flux magnitude", 1 );

        output( "Written VTK files:" );
        output( "  temperature_1.vtu" );
        output( "  thermal-conductivity_1.vtu" );
        output( "  temperature-gradient_1.vtu" );
        output( "  heat-flux-magnitude_1.vtu" );
        output( "Open these in ParaView to visualise the results." );
    }
    pause();

    subheading( "Simulation summary" );
    {
        blank();
        divider('-');
        cout << "  SIMULATION SUMMARY\n";
        divider('-');
        cout << "  Problem:    Steady-state heat conduction\n"
             << "  Equation:   div( lambda * grad T ) = 0\n"
             << "  Model:      " << model_name << "\n"
             << "  Nodes:      " << model.Mesh().Nodes()    << "\n"
             << "  Elements:   " << model.Mesh().Elements() << "\n\n"
             << "  Boundary conditions:\n"
             << "    TOP    -- Dirichlet: T = 15 degC\n"
             << "    BOTTOM -- Neumann:   q = 30 mW/m2\n"
             << "    Sides  -- No-flux (natural Neumann = 0)\n\n"
             << "  Thermal conductivities [W/(m K)]:\n"
             << "    LAYER_TOP       = 1.5\n"
             << "    LAYER_RESERVOIR = 4.0\n"
             << "    LAYER_BOTTOM    = 2.0\n"
             << "    NORMAL_FAULT    = 1.0\n\n"
             << "  Temperature range:\n";
        printRangeOfVariable( model, "Model", "temperature" );
        divider('-');
    }
    pause();

    keypoint(
"The complete PDE solve workflow in CSMP:\n"
"  1. Load model\n"
"  2. Assign material properties (per region)\n"
"  3. Set boundary conditions (Dirichlet = DIRICH flag, Neumann = NEUMANN flag)\n"
"  4. Create solver and PDE_Integrator\n"
"  5. Create and add LHS operators (Add)\n"
"  6. Create and add RHS operators (Add or AddBoundaryIntegral)\n"
"  7. Call IntegrateOver()\n"
"  8. Check result with printRangeOfVariable()\n"
"  9. Compute derived quantities (gradient, flux)\n"
" 10. Output to VTK"
    );
    pause();


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 4: Your First PDE Solve\n\n"
         << "  Steps covered:\n\n"
         << "    1.  The heat equation and its FEM discretisation\n"
         << "    2.  Loading the model\n"
         << "    3.  Assigning thermal conductivities\n"
         << "    4.  Setting boundary conditions\n"
         << "    5.  Building the PDE_Integrator\n"
         << "    6.  The LHS operator: conductance matrix\n"
         << "    7.  The RHS operator: basal heat flow (Neumann)\n"
         << "    8.  Solving and checking the result\n"
         << "    9.  Adding a heat sink (source term)\n"
         << "   10.  Output and visualisation\n\n"
         << "  Next: Example 5 -- InputDataManager and the PropertyDatabase\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nSteadyStateHeatConduction_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp
