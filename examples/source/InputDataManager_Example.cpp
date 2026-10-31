// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "InputDataManager_Example.h"
#include "exampleHelpers.h"

#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyHandle.h"
#include "PropertyDatabase.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"
#include "VTK_Interface.h"
#include "CSMP_highLevelUtilities.h"

#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dS.h"
#include "LinearSolver.h"

using namespace std;

namespace csmp {

// ============================================================================
// Specifications
// ============================================================================

void InputDataManager_Example::Specifications()
{
    SetTitle( "Model configuration with InputDataManager (interactive)" );
    SetDifficulty( 2 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: InputDataManager_Example.cpp" );
    AddDescription( "Interactive introduction to block-structured "
                    "configuration files, the PropertyDatabase, "
                    "variable ranges, whitespace rules, and "
                    "ComputationalSettings." );
    AddRequirement( "fault_boundary_test (CSMP binary files, "
                    "-configuration.txt), heat_conduction-variables.txt" );
}


// ============================================================================
// Local helper: write a configuration file
// ============================================================================

namespace {

/// Writes the example configuration file to disk.
/// Returns the filename written.
string WriteConfigurationFile( const string& model_name )
{
    const string fname = model_name + "-configuration.txt";
    ofstream ofs( fname );

    // Title line (not parsed -- just documentation)
    ofs << "fault_boundary_test configuration -- heat conduction example created by InputDataManager_Example.cpp\n";

    // -----------------------------------------------------------------------
    // Block 2: default property values for the whole model
    // (Block 1 is omitted -- we do not create regions from property ranges)
    // -----------------------------------------------------------------------
    ofs << "\n";  // blank line ends header line Block (Block 1 is empty here)

    ofs << "# Block 2: default property values for the entire model\n";
    ofs << "# Format: variable_name<TAB>value\n";
    ofs << "# NOTE: tokens MUST be separated by TAB characters.\n";
    ofs << "#       Spaces within variable names are allowed.\n";
    ofs << "#       Spaces between tokens will cause a read error.\n";
    ofs << "thermal conductivity\t2.0\n";
    ofs << "heat source\t0.0\n";
    ofs << "temperature\t0.0\n";
    ofs << "basal heat flow\t0.0\n";

    ofs << "\n";  // blank line ends Block 2

    // -----------------------------------------------------------------------
    // Block 3: regional property values
    // Format: region_name<TAB>part<TAB>variable_name<TAB>value
    // part = complete | interior | perimeter
    // -----------------------------------------------------------------------
    ofs << "# Block 3: regional property values\n";
    ofs << "# Format: region_name<TAB>part<TAB>variable_name<TAB>value\n";
    ofs << "# 'complete' applies to all nodes/elements in the region.\n";
    ofs << "# 'interior' applies only to interior nodes/elements.\n";
    ofs << "# 'perimeter' applies only to perimeter nodes/elements.\n";
    ofs << "LAYER_TOP\tcomplete\tthermal conductivity\t1.5\n";
    ofs << "LAYER_RESERVOIR\tcomplete\tthermal conductivity\t4.0\n";
    ofs << "LAYER_BOTTOM\tcomplete\tthermal conductivity\t2.0\n";
    ofs << "NORMAL_FAULT\tcomplete\tthermal conductivity\t1.0\n";

    ofs << "\n";  // blank line ends Block 3

    // -----------------------------------------------------------------------
    // Block 4: boundary conditions for box-shaped models
    // Format: BOX_BOUNDARY<TAB>flag<TAB>variable_name<TAB>value(s)
    // BOX_BOUNDARY = TOP | BOTTOM | LEFT | RIGHT | FRONT | BACK
    // flag = Dirichlet | Neumann | Plain
    // -----------------------------------------------------------------------
    ofs << "# Block 4: boundary conditions for box-shaped model\n";
    ofs << "# Format: BOX_BOUNDARY<TAB>flag<TAB>variable_name<TAB>value(s)\n";
    ofs << "# Dirichlet: fixed value (DOF eliminated from system)\n";
    ofs << "# Neumann:   flux prescribed (added to RHS)\n";
    ofs << "# Plain:     free (no special treatment)\n";
    ofs << "TOP\tDirichlet\ttemperature\t15.0 15.0 15.0 15.0\n";

    ofs << "\n";  // blank line ends Block 4

    // -----------------------------------------------------------------------
    // Block 5: essential conditions for regions
    // Format: region_name<TAB>part<TAB>variable_name<TAB>flag
    // -----------------------------------------------------------------------
    ofs << "# Block 5: essential conditions for regions\n";
    ofs << "# Format: region_name<TAB>part<TAB>variable_name<TAB>flag\n";
    ofs << "# (empty in this example -- no regional Dirichlet conditions)\n";

    ofs << "\n";  // blank line ends Block 5

    // -----------------------------------------------------------------------
    // Block 6: boundary conditions for arbitrary-shaped boundaries
    // Format: BOUNDARY_NAME<TAB>part<TAB>flag<TAB>variable_name<TAB>value
    // -----------------------------------------------------------------------
    ofs << "# Block 6: boundary conditions for arbitrary-shaped boundaries\n";
    ofs << "# Format: BOUNDARY_NAME<TAB>part<TAB>flag<TAB>variable_name<TAB>value\n";
    ofs << "# Used for boundaries that are NOT box boundaries.\n";
    ofs << "# The boundary name must contain the string 'BOUNDARY'.\n";
    ofs << "BOTTOM\tcomplete\tNeumann\tbasal heat flow\t0.030\n";

    ofs << "\n";  // blank line ends Block 6

    // -----------------------------------------------------------------------
    // ComputationalSettings (optional -- appended after Block 6)
    // -----------------------------------------------------------------------
    ofs << "# ComputationalSettings\n";
    ofs << "# Format: keyword<TAB>value\n";
    ofs << "# This block is read only if a ComputationalSettings object\n";
    ofs << "# is passed to ConfigureFromFile().\n";
    ofs << "time stepping\tconservative\n";
    ofs << "duration\t1.0e6\n";
    ofs << "output time\t1.0e5\n";
    ofs << "output time\t5.0e5\n";
    ofs << "output time\t1.0e6\n";

    // -----------------------------------------------------------------------
    // End-of-file convention:
    // Anything after 5 consecutive blank lines is treated as comments
    // and is completely ignored by the parser.
    // -----------------------------------------------------------------------
    ofs << "\n\n\n\n\n";  // 5 blank lines = end of parsed content

    ofs << "=============================================================\n";
    ofs << "COMMENTS IGNORED BY THE SIMULATOR\n";
    ofs << "=============================================================\n";
    ofs << "\n";
    ofs << "This section can contain any text:\n";
    ofs << "  - calibration notes\n";
    ofs << "  - references\n";
    ofs << "  - alternative parameter values (commented out)\n";
    ofs << "  - observations made during the simulation\n";
    ofs << "\n";
    ofs << "# Alternative thermal conductivities (not used):\n";
    ofs << "# LAYER_TOP\tcomplete\tthermal conductivity\t1.2\n";
    ofs << "# LAYER_RESERVOIR\tcomplete\tthermal conductivity\t3.5\n";

    ofs.close();
    return fname;
}

} // anonymous namespace


// ============================================================================
// Run
// ============================================================================

/**
    @file InputDataManager_Example.cpp

    @brief CSMP for Beginners -- Example 5: Model Configuration

    @section purpose Purpose

    This interactive example introduces the InputDataManager -- the
    standard way to configure CSMP models from plain text files without
    recompiling the simulation code.

    @section concepts Concepts Covered

      1.  Why configuration files exist
      2.  The variable definition file (*-variables.txt)
      3.  The PropertyDatabase and physically meaningful ranges
      4.  The configuration file structure -- six blocks
      5.  Block 1: creating regions from property ranges
      6.  Block 2: default property values for the whole model
      7.  Block 3: regional property values
      8.  Block 4: boundary conditions for box-shaped models
      9.  Block 5: essential conditions for regions
     10.  Block 6: boundary conditions for arbitrary-shaped models
     11.  Whitespace rules -- tabs only between tokens
     12.  Comments and the end-of-file convention
     13.  ComputationalSettings
     14.  Calling ConfigureFromFile()
     15.  Verifying the configuration

    @section model The Example Model

    We use 'fault_boundary_test' -- the same layered geological model
    as in Example 4. We write a configuration file for it step by step,
    then load and verify it.

    @section requirements Requirements

    - fault_boundary_test  (CSMP binary files + -configuration.txt)
    - heat_conduction-variables.txt  (variable definition file)

    @author S.K. Matthai
*/

void InputDataManager_Example::Run()
{
    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 5\n"
         << "  Model Configuration with InputDataManager\n\n"
         << "  You will learn how to configure a CSMP model from a\n"
         << "  plain text file -- without recompiling the simulation.\n\n"
         << "  Steps:\n"
         << "    1.  Why configuration files exist\n"
         << "    2.  The variable definition file\n"
         << "    3.  The PropertyDatabase and variable ranges\n"
         << "    4.  The configuration file structure\n"
         << "    5.  Block 2: default property values\n"
         << "    6.  Block 3: regional property values\n"
         << "    7.  Block 4: box boundary conditions\n"
         << "    8.  Block 5: essential conditions for regions\n"
         << "    9.  Block 6: arbitrary boundary conditions\n"
         << "   10.  Whitespace rules -- tabs only\n"
         << "   11.  Comments and end-of-file convention\n"
         << "   12.  ComputationalSettings\n"
         << "   13.  Calling ConfigureFromFile()\n"
         << "   14.  Verifying the configuration\n"
         << "   15.  Running the heat conduction solve\n";
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
    // Step 1: Why configuration files exist
    // ========================================================================
    heading( 1, 15, "Why configuration files exist" );

    explain(
"In Example 4, we assigned material properties and boundary conditions\n"
"directly in the C++ source code:\n"
"\n"
"    model.Region(\"LAYER_TOP\").InputPropertyValue(\n"
"        \"thermal conductivity\", makeScalar(PLAIN, 1.5) );\n"
"\n"
"This works, but it has a serious drawback:\n"
"  - Every time you change a parameter, you must recompile.\n"
"  - Recompilation takes minutes for large projects.\n"
"  - You cannot run parameter studies without code changes.\n"
"  - Collaborators cannot change parameters without C++ knowledge.\n"
"\n"
"The InputDataManager makes this more versatile by reading parameter values from\n"
"a plain text file at runtime. You change the file, re-run the\n"
"simulation -- no recompilation needed."
    );
    pause();

    subheading( "The CSMP input file set" );
    explain(
"A typical CSMP simulation reads from several text files:\n"
"\n"
"  model-variables.txt       -- which physical variables exist\n"
"                               and their valid ranges\n"
"  model-configuration.txt   -- material properties, initial\n"
"                               conditions, boundary conditions, output time steps\n"
"  model-regions.txt         -- (optional) which regions to keep\n"
"\n"
"The variable file is read when the Model is constructed.\n"
"The configuration file is read by InputDataManager::ConfigureFromFile().\n"
"\n"
"This separation means:\n"
"  - The variable file defines WHAT variables exist.\n"
"  - The configuration file defines WHAT VALUES they have."
    );
    pause();

    keypoint(
"Configuration files let you change simulation parameters\n"
"without recompiling. This is essential for parameter studies,\n"
"calibration, and collaboration with non-programmers.\n"
"The variable file and configuration file work together."
    );
    pause();


    // ========================================================================
    // Step 2: The variable definition file
    // ========================================================================
    heading( 2, 15, "The variable definition file" );

    explain(
"The variable definition file (*-variables.txt) tells CSMP:\n"
"  - the name of each physical variable\n"
"  - its mathematical notation\n"
"  - its physical unit\n"
"  - its type (SCALAR, VECTOR, TENSOR, ARRAY, FLAGGEDARRAY)\n"
"  - its placement (NODE, ELEMENT, FACE, ...)\n"
"  - its physically valid range [vmin, vmax]\n"
"\n"
"The Model constructor reads this file and initialises the\n"
"PropertyDatabase with all these specifications."
    );
    pause();

    subheading( "Variable file format" );
    code(
"# Each line defines one variable.\n"
"# Columns separated by TABS:\n"
"# name<TAB>notation<TAB>unit<TAB>type<TAB>vmin<TAB>vmax<TAB>placement\n"
"\n"
"# Example lines from heat_conduction-variables.txt:\n"
"temperature<TAB>T<TAB>degC<TAB>SCALAR<TAB>-273<TAB>1000<TAB>NODE\n"
"thermal conductivity<TAB>lambda<TAB>W/(m K)<TAB>SCALAR<TAB>0<TAB>100<TAB>ELEMENT\n"
"heat source<TAB>q<TAB>W/m3<TAB>SCALAR<TAB>-1e6<TAB>1e6<TAB>ELEMENT\n"
"basal heat flow<TAB>q_b<TAB>W/m2<TAB>SCALAR<TAB>0<TAB>1<TAB>FACE\n"
"\n"
"# Note: variable names CAN contain spaces.\n"
"# This is why TAB is the only valid token separator."
    );
    pause();

    subheading( "Live demonstration -- inspect the variable file" );
    {
        cout << "  Loading model to inspect its PropertyDatabase...\n";
        cout.flush();
        Model<dim>  model_check( model_name, variable_file );
        output( "Variables defined in '" + variable_file + "':" );
        model_check.Database().FlushToScreen();
    }
    pause();

    keypoint(
"The variable file defines the PARAMETERS of the simulation.\n"
"Every variable used in the configuration file MUST be defined here.\n"
"The valid range [vmin, vmax] is checked after configuration --\n"
"values outside the range trigger a warning."
    );
    pause();


    // ========================================================================
    // Step 3: The PropertyDatabase and variable ranges
    // ========================================================================
    heading( 3, 15, "The PropertyDatabase and variable ranges" );

    explain(
"The PropertyDatabase stores the specification of every variable.\n"
"While ConfigureFromFile() runs, the database is queried, checking\n"
"whether all variable values are within their physically valid ranges.\n"
"\n"
"This is a sometimes painful but essential sanity check:\n"
"  - A negative thermal conductivity would be caught immediately.\n"
"  - A temperature below absolute zero would be flagged.\n"
"  - A permeability above the physical maximum would be reported.\n"
"\n"
"It is your task to set tight ranges for your problem. The check is done automatically with IsWithinRange() calls."
    );
    pause();

    subheading( "Range checking" );
    code(
"// Check all variables in the whole model:\n"
"model.IsWithinRange( \"thermal conductivity\" );\n"
"\n"
"// Check in a specific region:\n"
"model.IsWithinRange( \"thermal conductivity\", \"LAYER_TOP\" );\n"
"\n"
"// The range is defined in the variable file:\n"
"// thermal conductivity: vmin=1, vmax=10 W/(m K)\n"
"// If any element has lambda < 1 or lambda > 10, IsWithinRange returns false."
    );
    pause();

    subheading( "Why ranges matter" );
    explain(
"Physical variable ranges serve three purposes:\n"
"\n"
"1. DOCUMENTATION: they record the physically meaningful bounds\n"
"   for each variable, making the simulation self-documenting.\n"
"\n"
"2. VERIFICATION: after reading the configuration file, CSMP checks\n"
"   that all assigned values are within the declared ranges.\n"
"   This catches typos like '1.5e-3' instead of '1.5' for\n"
"   thermal conductivity.\n"
"\n"
"3. SIMULATION VERIFICATION: of Interrelations, Visitors, PDE solves etc.\n"
"   The correct functioning of a complex simulation is almost impossible to verify, yet\n"
"   making sure that each model component receives and computes parameter values\n"
"   in ranges for which its correct operation has been verified,\n"
"   is an important steps toward this goal.\n"
"\n"
"The range check is NOT a hard error -- it is a WARNING.\n"
"The simulation continues, but the user is alerted."
"Warnings are time-stamped and logged into a file by ErrorHandler."
    );
    pause();

    keypoint(
"Always define physically meaningful ranges in the variable file.\n"
"Call var::IsWithinRange() to verify parameter input and output values in own code.\n"
"A range violation is a warning, not a crash -- but it usually\n"
"indicates a typo or unit error raising the possibility for follow-up errors."
    );
    pause();


    // ========================================================================
    // Step 4: The configuration file structure
    // ========================================================================
    heading( 4, 15, "The configuration file structure" );

    explain(
"The configuration file is divided into up to six blocks,\n"
"each separated by a SINGLE blank line.\n"
"\n"
"  Line 1:   Title (not parsed -- documentation only)\n"
"  Block 1:  Create regions from property ranges (optional)\n"
"  Block 2:  Default property values for the whole model\n"
"  Block 3:  Regional property values\n"
"  Block 4:  Boundary conditions for box-shaped models\n"
"  Block 5:  Essential conditions for regions\n"
"  Block 6:  Boundary conditions for arbitrary-shaped boundaries\n"
"  Optional: ComputationalSettings\n"
"\n"
"After 5 consecutive blank lines, everything is treated as\n"
"free-form comments and is completely ignored by the parser."
    );
    pause();

    subheading( "Block structure overview" );
    code(
"Title line -- not parsed\n"
"\n"
"# Block 1: region creation from property ranges (optional)\n"
"region_name<TAB>property_name<TAB>min_value max_value\n"
"\n"
"# Block 2: default values for whole model\n"
"variable_name<TAB>value\n"
"\n"
"# Block 3: regional values\n"
"region_name<TAB>part<TAB>variable_name<TAB>value\n"
"\n"
"# Block 4: box boundary conditions\n"
"BOX_BOUNDARY<TAB>flag<TAB>variable_name<TAB>value(s)\n"
"\n"
"# Block 5: essential conditions for regions\n"
"region_name<TAB>part<TAB>variable_name<TAB>flag\n"
"\n"
"# Block 6: arbitrary boundary conditions\n"
"BOUNDARY_NAME<TAB>part<TAB>flag<TAB>variable_name<TAB>value\n"
"\n"
"# ComputationalSettings\n"
"keyword<TAB>value\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"Everything after 5 blank lines is ignored."
    );
    pause();

    keypoint(
"Blocks are separated by EXACTLY ONE blank line.\n"
"The parser counts blank lines to know which block it is in.\n"
"An extra blank line inside a block will confuse the parser.\n"
"After 5 blank lines: free-form comments, completely ignored."
    );
    pause();


    // ========================================================================
    // Step 5: Block 2 -- default property values
    // ========================================================================
    heading( 5, 15, "Block 2: default property values" );

    explain(
"Block 2 assigns uniform values to the ENTIRE model.\n"
"These are the model 'default' values that apply everywhere\n"
"unless overridden by Block 3 (regional values).\n"
"\n"
"Format:\n"
"    variable_name<TAB>value\n"
"\n"
"The variable name can contain spaces -- this is why TAB\n"
"is the only valid separator between tokens."
    );
    pause();

    subheading( "Example Block 2" );
    code(
"# Block 2: default property values for the entire model\n"
"thermal conductivity<TAB>2.0\n"
"heat source<TAB>0.0\n"
"temperature<TAB>0.0\n"
"basal heat flow<TAB>0.0\n"
"\n"
"# This sets lambda = 2.0 W/(m K) everywhere.\n"
"# Block 3 will then override this for specific regions."
    );
    pause();

    subheading( "What happens if Block 2 is skipped?" );
    explain(
"If you skip Block 2 (pass block2=false to ConfigureFromFile),\n"
"variables retain whatever value they had from the Model constructor.\n"
"For CSMP variables this is NaN (not a number).\n"
"\n"
"NaN values get propagated through calculations and will cause the solver to produce NaN results.\n"
"Always provide meaningful values in Block 2 for input variables\n"
"to PDE operators and calculations."
    );
    pause();

    keypoint(
"Block 2 sets 'background' values for the whole model.\n"
"Block 3 then overrides specific regions.\n"
"Always provide defaults for ALL variables used in computations.\n"
"Unset variables default to NaN -- which will corrupt the solution."
    );
    pause();


    // ========================================================================
    // Step 6: Block 3 -- regional property values
    // ========================================================================
    heading( 6, 15, "Block 3: regional property values" );

    explain(
"Block 3 assigns values to specific model subdomains.\n"
"It overrides the defaults set in Block 2.\n"
"\n"
"Format:\n"
"    region_name<TAB>part<TAB>variable_name<TAB>value\n"
"\n"
"'part' controls which nodes/elements are affected:\n"
"  complete  -- all nodes/elements in the region\n"
"  interior  -- only interior nodes/elements\n"
"  perimeter -- only perimeter nodes/elements\n"
"\n"
"IMPORTANT: perimeter nodes are SHARED between adjacent regions.\n"
"The LAST assignment wins for shared perimeter nodes.\n"
"Order your Block 3 entries carefully."
    );
    pause();

    subheading( "Example Block 3" );
    code(
"# Block 3: regional thermal conductivities\n"
"LAYER_TOP<TAB>complete<TAB>thermal conductivity<TAB>1.5\n"
"LAYER_RESERVOIR<TAB>complete<TAB>thermal conductivity<TAB>4.0\n"
"LAYER_BOTTOM<TAB>complete<TAB>thermal conductivity<TAB>2.0\n"
"NORMAL_FAULT<TAB>complete<TAB>thermal conductivity<TAB>1.0\n"
"\n"
"# Multiple properties for the same region:\n"
"LAYER_RESERVOIR<TAB>complete<TAB>heat source<TAB>-1e-4\n"
"\n"
"# Using 'interior' to avoid overwriting shared boundary nodes:\n"
"LAYER_TOP<TAB>interior<TAB>temperature<TAB>50.0\n"
"# (perimeter nodes of LAYER_TOP are shared with LAYER_RESERVOIR)"
    );
    pause();

    subheading( "The perimeter node sharing issue" );
    explain(
"Consider two adjacent regions A and B sharing a boundary:\n"
"\n"
"    A | B\n"
"    -----\n"
"    shared perimeter nodes\n"
"\n"
"If Block 3 contains:\n"
"    A  complete  temperature  100\n"
"    B  complete  temperature  200\n"
"\n"
"The shared nodes will have temperature = 200 (B wins, assigned last).\n"
"\n"
"To avoid this, use 'interior' for one of the regions:\n"
"    A  interior  temperature  100\n"
"    B  complete  temperature  200\n"
"\n"
"Now A's interior nodes get 100, B's nodes (including shared) get 200."
    );
    pause();

    keypoint(
"Block 3 overrides Block 2 for specific regions.\n"
"'complete' affects all nodes including shared perimeter nodes.\n"
"'interior' affects only non-shared nodes.\n"
"The LAST assignment wins for shared perimeter nodes.\n"
"Use 'interior' to prevent unintended overwriting."
    );
    pause();


    // ========================================================================
    // Step 7: Block 4 -- box boundary conditions
    // ========================================================================
    heading( 7, 15, "Block 4: box boundary conditions" );

    explain(
"Block 4 assigns boundary conditions to the six faces of a\n"
"box-shaped model using the BOX_BOUNDARY identifiers:\n"
"  TOP, BOTTOM, LEFT, RIGHT, FRONT, BACK\n"
"\n"
"Format:\n"
"    BOX_BOUNDARY<TAB>flag<TAB>variable_name<TAB>value(s)\n"
"\n"
"flag can be:\n"
"  Dirichlet (or DIRICH)  -- fixed value\n"
"  Neumann                -- flux prescribed\n"
"  Plain                  -- free (no special treatment)\n"
"\n"
"For scalar variables: one value.\n"
"For vector variables: dim values (one per component).\n"
"For 2D models: 2 values per boundary side (corners).\n"
"For 3D models: 4 values per boundary face (corners).\n"
"  Linear interpolation is applied between the corner values,\n"
"  but you can also compute entire fields on boundaries."
    );
    pause();

    subheading( "Example Block 4" );
    code(
"# Block 4: box boundary conditions\n"
"TOP<TAB>Dirichlet<TAB>temperature<TAB>15.0\n"
"\n"
"# Linear variation across the bottom (4 corner values in 3D):\n"
"BOTTOM<TAB>Dirichlet<TAB>temperature<TAB>60.0 60.0 65.0 65.0\n"
"\n"
"# No-flux on sides (natural Neumann = 0, no entry needed):\n"
"# LEFT, RIGHT, FRONT, BACK are automatically no-flux\n"
"\n"
"# Displacement boundary condition (vector, 3 components):\n"
"BOTTOM<TAB>Dirichlet<TAB>displacement<TAB>0.0 0.0 0.0"
    );
    pause();

    subheading( "Block 4 vs Block 6" );
    explain(
"Block 4 uses BOX_BOUNDARY identifiers (TOP, BOTTOM, etc.) to assign boundary conditions.\n"
"It only works for box-shaped models.\n"
"\n"
"Block 6 uses the actual CSMP Boundary object names containing the string \"boundary\".\n"
"It works for any model shape, including irregular boundaries.\n"
"\n"
"For the fault_boundary_test model:\n"
"  - TOP is a BOX_BOUNDARY --> use Block 4 for Dirichlet T\n"
"  - BOTTOM is also a BOX_BOUNDARY, but we use Block 6 for\n"
"    the Neumann basal heat flow because it is a CSMP Boundary object"
    );
    pause();

    keypoint(
"Block 4 is for box-shaped models only (TOP, BOTTOM, LEFT, etc.).\n"
"Block 6 is for arbitrary-shaped CSMP Boundary objects.\n"
"For 3D box models, 4 corner values enable linear interpolation.\n"
"No-flux (Neumann = 0) is the default in FEM computations -- no entry needed."
    );
    pause();


    // ========================================================================
    // Step 8: Block 5 -- essential conditions for regions
    // ========================================================================
    heading( 8, 15, "Block 5: essential conditions for regions" );

    explain(
"Block 5 assigns VARIABLE FLAGS to specific regions.\n"
"It does NOT assign values -- only flags.\n"
"\n"
"Format:\n"
"    region_name<TAB>part<TAB>variable_name<TAB>flag\n"
"\n"
"This is used when:\n"
"  - You want to fix a variable in a region (DIRICH)\n"
"    without using a box boundary\n"
"  - You want to mark a region as a Neumann source\n"
"  - You want to release a condition set in Block 4\n"
"\n"
"Typical use: fix pressure in a well region,\n"
"or fix temperature in a heat source region."
    );
    pause();

    subheading( "Example Block 5" );
    code(
"# Block 5: essential conditions for regions\n"
"\n"
"# Fix temperature in the fault zone (e.g. geothermal spring):\n"
"NORMAL_FAULT<TAB>complete<TAB>temperature<TAB>Dirichlet\n"
"\n"
"# Fix pressure in a well region:\n"
"WELL<TAB>complete<TAB>fluid pressure<TAB>Dirichlet\n"
"\n"
"# Release Dirichlet condition at inlet/outlet nodes:\n"
"# (set in Block 4 but released here for specific regions)\n"
"INLET<TAB>interior<TAB>fluid pressure<TAB>Plain"
    );
    pause();

    subheading( "Block 5 vs Block 4" );
    explain(
"Block 4 assigns BOTH a flag AND a value to boundary nodes.\n"
"Block 5 assigns ONLY a flag to region nodes.\n"
"\n"
"Typical workflow:\n"
"  1. Block 3: assign the desired value to the region\n"
"  2. Block 5: flag those nodes as DIRICH\n"
"\n"
"This two-step approach gives full control:\n"
"  - The value comes from Block 3 (or Block 2)\n"
"  - The flag comes from Block 5\n"
"  - The assembler sees DIRICH and eliminates those DOFs"
    );
    pause();

    keypoint(
"Block 5 assigns FLAGS only -- not values.\n"
"Use Block 3 to set the value, Block 5 to set the flag.\n"
"This two-step approach gives full control over BCs in regions.\n"
"Typical use: fix pressure in wells, temperature in heat sources, prevent override."
    );
    pause();


    // ========================================================================
    // Step 9: Block 6 -- arbitrary boundary conditions
    // ========================================================================
    heading( 9, 15, "Block 6: arbitrary boundary conditions" );

    explain(
"Block 6 assigns boundary conditions to CSMP Boundary objects.\n"
"These are named boundaries that can have any shape.\n"
"\n"
"Format:\n"
"    BOUNDARY_NAME<TAB>part<TAB>flag<TAB>variable_name<TAB>value\n"
"\n"
"The boundary name must be a valid CSMP Boundary object name.\n"
"It does NOT need to contain the string 'BOUNDARY' --\n"
"any named boundary works.\n"
"\n"
"This is the most flexible block:\n"
"  - Works for box-shaped AND irregular models\n"
"  - Can assign to complete, interior, or perimeter of boundary\n"
"  - Supports Dirichlet, Neumann, and Plain flags"
    );
    pause();

    subheading( "Example Block 6" );
    code(
"# Block 6: arbitrary boundary conditions\n"
"# Format: BOUNDARY_NAME<TAB>part<TAB>flag<TAB>variable_name<TAB>value\n"
"\n"
"# Basal heat flow on BOTTOM boundary (Neumann):\n"
"BOTTOM<TAB>complete<TAB>Neumann<TAB>basal heat flow<TAB>0.030\n"
"\n"
"# Fixed pressure on an irregular inlet boundary:\n"
"BOUNDARY_INLET<TAB>complete<TAB>Dirichlet<TAB>fluid pressure<TAB>1.0e7\n"
"\n"
"# Fixed temperature on a fault boundary:\n"
"BOUNDARY_FAULT_TOP<TAB>complete<TAB>Dirichlet<TAB>temperature<TAB>80.0"
    );
    pause();

    keypoint(
"Block 6 is the most flexible boundary condition block.\n"
"It works for any named CSMP Boundary object and boundary condition.\n"
"Use it for Neumann conditions (heat flux, fluid flux)\n"
"and for Dirichlet conditions on irregular boundaries.\n"
"The boundary name must match exactly (case-sensitive)."
    );
    pause();


    // ========================================================================
    // Step 10: Whitespace rules -- tabs only
    // ========================================================================
    heading( 10, 15, "Whitespace rules -- tabs only between tokens" );

    explain(
"This is the most common source of errors in configuration files.\n"
"\n"
"CSMP variable names can contain spaces:\n"
"    'thermal conductivity'  -- two words, one variable name\n"
"    'fluid pressure'        -- two words, one variable name\n"
"    'brooks corey parameter' -- three words, one variable name\n"
"\n"
"Because of this, the parser CANNOT use spaces as token separators.\n"
"Instead, it uses TAB characters (\\t) as the ONLY valid separator.\n"
"\n"
"RULE: Every token in a configuration line must be separated\n"
"      by exactly one TAB character. No spaces between tokens."
    );
    pause();

    subheading( "Correct vs incorrect formatting" );
    code(
"# CORRECT: tokens separated by single or multiple TAB(s) (shown as --> here)\n"
"thermal conductivity-->2.0\n"
"LAYER_TOP-->complete-->thermal conductivity-->1.5\n"
"TOP-->Dirichlet-->temperature-->15.0\n"
"\n"
"# WRONG: spaces between tokens (will cause read error)\n"
"thermal conductivity  2.0          # WRONG: spaces, not tab\n"
"LAYER_TOP complete thermal conductivity 1.5  # WRONG\n"
"\n"
"# WRONG: extra spaces at end of line\n"
"thermal conductivity-->2.0   # WRONG: trailing spaces\n"
"\n"
"# CORRECT: spaces WITHIN a variable name are fine\n"
"brooks corey parameter-->2.0  # 'brooks corey parameter' is one token"
    );
    pause();

    subheading( "How to check your file" );
    explain(
"In most text editors, you can show invisible characters\n"
"(tabs appear as arrows, spaces as dots).\n"
"\n"
"In vim:    :set list\n"
"In Emacs:  M-x whitespace-mode\n"
"In VS Code: View > Render Whitespace\n"
"\n"
"A common mistake is to use spaces instead of tabs when\n"
"copying lines from a document or email (or automatic TAB to space conversion).\n"
"Always verify that separators are TAB characters,\n"
"not sequences of spaces that look like tabs."
    );
    pause();

    keypoint(
"TAB is the ONLY valid token separator in configuration files.\n"
"Spaces are allowed WITHIN variable names but NOT between tokens.\n"
"Trailing spaces at the end of a line will cause read errors.\n"
"Use your editor's 'show whitespace' and 'zap gremlins' features to verify.\n"
"Watch out when using 'rich text' editors and transferring files between platforms."
    );
    pause();


    // ========================================================================
    // Step 11: Comments and end-of-file convention
    // ========================================================================
    heading( 11, 15, "Comments and end-of-file convention" );

    explain(
"Comments can be inserted anywhere in the configuration file.\n"
"A comment line starts with a # character.\n"
"\n"
"    # This is a comment -- the parser ignores this line\n"
"\n"
"Comments are useful for:\n"
"  - documenting why a particular value was chosen\n"
"  - temporarily disabling a line without deleting it\n"
"  - recording calibration history\n"
"\n"
"End-of-file convention:\n"
"  After 5 consecutive blank lines, the parser stops reading.\n"
"  Everything after the 5th blank line is completely ignored.\n"
"  This allows you to append free-form notes, observations,\n"
"  and alternative parameter sets to the file."
    );
    pause();

    subheading( "The end-of-file section" );
    code(
"# ... last block of configuration ...\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"=============================================================\n"
"COMMENTS IGNORED BY THE SIMULATOR\n"
"=============================================================\n"
"\n"
"This section can contain anything:\n"
"  - calibration notes\n"
"  - references to papers\n"
"  - alternative parameter values (commented out)\n"
"  - observations made during the simulation\n"
"\n"
"# Alternative values not currently used:\n"
"# LAYER_TOP<TAB>complete<TAB>thermal conductivity<TAB>1.2\n"
"\n"
"WellManager::MaximumHistoricProductionRates:\n"
"  DM2X: 8.29e-03\n"
"  R420: 5.06e-03"
    );
    pause();

    keypoint(
"Start comment lines with #.\n"
"After 5 consecutive blank lines: everything is ignored.\n"
"Use the end-of-file section for notes, history, and alternatives.\n"
"This makes the configuration file self-documenting."
    );
    pause();


    // ========================================================================
    // Step 12: ComputationalSettings
    // ========================================================================
    heading( 12, 15, "ComputationalSettings" );

    explain(
"ComputationalSettings is an optional block appended after Block 6.\n"
"It stores runtime parameters for time-dependent simulations:\n"
"  - time stepping strategy\n"
"  - simulation duration\n"
"  - output times\n"
"\n"
"It is only read if a ComputationalSettings object is passed\n"
"to the ConfigureFromFile() overload that accepts one.\n"
"\n"
"For steady-state simulations (like our heat conduction example),\n"
"ComputationalSettings is not needed."
    );
    pause();

    subheading( "ComputationalSettings format" );
    code(
"# ComputationalSettings block (after Block 6)\n"
"# Format: keyword<TAB>value\n"
"\n"
"time stepping<TAB>conservative    # or: aggressive\n"
"duration<TAB>1.0e6               # total simulation time [s]\n"
"output time<TAB>1.0e5            # write output at this time [s]\n"
"output time<TAB>5.0e5            # multiple output times allowed\n"
"output time<TAB>1.0e6"
    );
    pause();

    subheading( "Time stepping strategies" );
    explain(
"'conservative' -- small global time steps, more stable, slower\n"
"'aggressive'   -- large global time steps, faster, may diverge\n"
"'asynchronous' -- adaptive locally optimaltime steps, fast in heterogeneous models\n"
"\n"
"The actual time step size is computed by the transport scheme used\n"
"(e.g. Courant number criterion or other). The strategy keyword modifies\n"
"a safety factor applied to this computed step size."
    );
    pause();

    subheading( "ConfigureFromFile() with ComputationalSettings" );
    code(
"ComputationalSettings  run_settings;\n"
"\n"
"model_configuration.ConfigureFromFile(\n"
"    model,\n"
"    model_name.c_str(),\n"
"    false,   // block1: region creation from property ranges\n"
"    true,    // block2: default property values\n"
"    true,    // block3: regional property values\n"
"    true,    // block4: box boundary conditions\n"
"    false,   // block5: essential conditions for regions\n"
"    true,    // block6: arbitrary boundary conditions\n"
"    run_settings );  // reads ComputationalSettings block\n"
"\n"
"// Access the settings:\n"
"double duration = run_settings.Duration();\n"
"double dt_first = run_settings.FirstOutputTime();"
    );
    pause();

    keypoint(
"ComputationalSettings stores time stepping parameters.\n"
"It is only needed for transient simulations.\n"
"Pass a ComputationalSettings object to ConfigureFromFile()\n"
"to read the time stepping block from the configuration file."
    );
    pause();


    // ========================================================================
    // Step 13: Calling ConfigureFromFile()
    // ========================================================================
    heading( 13, 15, "Calling ConfigureFromFile()" );

    explain(
"ConfigureFromFile() reads the configuration file and applies\n"
"all the specified values and flags to the model.\n"
"\n"
"The boolean arguments control which blocks are read:\n"
"  block1 -- region creation from property ranges\n"
"  block2 -- default property values\n"
"  block3 -- regional property values\n"
"  block4 -- box boundary conditions\n"
"  block5 -- essential conditions for regions\n"
"  block6 -- arbitrary boundary conditions\n"
"\n"
"You can skip blocks by passing false.\n"
"Skipped blocks are treated as empty -- the parser advances\n"
"past the blank line that would end that block."
    );
    pause();

    subheading( "The two main overloads" );
    code(
"InputDataManager<dim>  config;\n"
"\n"
"// Simple overload (reads all blocks 1-5, box-shaped model):\n"
"config.ConfigureFromFile( model, model_name.c_str() );\n"
"\n"
"// Full overload (choose which blocks to read):\n"
"config.ConfigureFromFile(\n"
"    model,\n"
"    model_name.c_str(),\n"
"    false,   // block1\n"
"    true,    // block2\n"
"    true,    // block3\n"
"    true,    // block4\n"
"    false,   // block5\n"
"    true );  // block6\n"
"\n"
"// With ComputationalSettings:\n"
"ComputationalSettings  settings;\n"
"config.ConfigureFromFile(\n"
"    model, model_name.c_str(),\n"
"    false, true, true, true, false, true,\n"
"    settings );"
    );
    pause();

    subheading( "Writing and loading the configuration file" );
    {
        output( "Writing configuration file..." );
        const string config_fname = WriteConfigurationFile( model_name );
        output( "Written: " + config_fname );

        // Show the file contents
        blank();
        divider('.');
        cout << "  Contents of " << config_fname << ":\n";
        divider('.');
        ifstream ifs( config_fname );
        string line;
        int line_num{ 0 };
        while ( getline( ifs, line ) && line_num < 40 ) {
            // Replace tabs with visible marker for display
            string display_line = line;
            size_t pos{ 0U };
            while ( (pos = display_line.find('\t', pos)) != string::npos ) {
                display_line.replace( pos, 1, "<TAB>" );
                pos += 5;
            }
            cout << "  " << setw(3) << ++line_num << "  " << display_line << "\n";
        }
        if ( line_num == 40 )
            cout << "  ... (truncated for display)\n";
        divider('.');
    }
    pause();

    keypoint(
"ConfigureFromFile() reads the file and applies all values.\n"
"Use the boolean arguments to skip blocks you do not need.\n"
"The file name is constructed as: model_name + '-configuration.txt'\n"
"unless you use ConfigureFromSpecificFile() with a custom name."
    );
    pause();


    // ========================================================================
    // Step 14: Verifying the configuration
    // ========================================================================
    heading( 14, 15, "Verifying the configuration" );

    explain(
"After ConfigureFromFile() runs, always verify the result:\n"
"\n"
"  1. Check variable ranges with model.IsWithinRange()\n"
"  2. Print variable ranges with printRangeOfVariable()\n"
"  3. Check boundary conditions using ::IsWithinRange() or other functionality \n"
"\n"
"These checks catch the most common configuration errors:\n"
"  - Wrong units (e.g. mW/m2 instead of W/m2)\n"
"  - Wrong sign (e.g. positive heat sink)\n"
"  - Missing values (NaN from uninitialised variables)\n"
"  - Wrong region name (typo -- region not found, value not set)"
    );
    pause();

    subheading( "Live verification" );
    {
        cout << "  Loading model and applying configuration...\n";
        cout.flush();

        Model<dim>  model( model_name, variable_file );
        Region<dim>& gref( model.Region("Model") );
        WriteConfigurationFile(model_name);

        // making sure the necessary parameters are there
        // NOTE: when defining parameters, name and notation must be unique
        if ( !model.Database().IsDefined("thermal conductivity") )
          model.CreateProperty("thermal conductivity", "K_t", "W m-2", SCALAR, ELEMENT );
        if ( !model.Database().IsDefined("heat source") )
          model.CreateProperty("heat source", "Q", "W m-3", SCALAR, ELEMENT );
        if ( !model.Database().IsDefined("basal heat flow") )
          model.CreateProperty("basal heat flow", "HFU", "W m-2", SCALAR, FACE );
        if ( !model.Database().IsDefined("temperature") )
          model.CreateProperty("temperature", "T", "oC", SCALAR, NODE );

        InputDataManager<dim>  config;
        config.ConfigureFromFile(
            model,
            model_name.c_str(),
            false,   // block1: no region creation
            true,    // block2: default values
            true,    // block3: regional values
            true,    // block4: box BCs
            false,   // block5: no regional flags
            true );  // block6: arbitrary BCs

        output( "Configuration applied." );
        blank();

        // Verify thermal conductivities
        output( "Thermal conductivity ranges by region:" );
        printRangeOfVariable( model, "LAYER_TOP",       "thermal conductivity" );
        printRangeOfVariable( model, "LAYER_RESERVOIR", "thermal conductivity" );
        printRangeOfVariable( model, "LAYER_BOTTOM",    "thermal conductivity" );
        printRangeOfVariable( model, "NORMAL_FAULT",    "thermal conductivity" );

        // Verify boundary conditions
        output( "Temperature at TOP boundary (should be 15 degC):" );
        printRangeOfVariable( model, "TOP", "temperature" );

        output( "Basal heat flow at BOTTOM boundary (should be 0.030 W/m2):" );
        printRangeOfVariable( model, "BOTTOM", "basal heat flow" );

        // Count Dirichlet nodes
        const csmp::Index T_key = model.Database().StorageKey("temperature");
        size_t n_dirich{ 0U };
        ScalarVariable sc;
        for ( const auto& nit : gref.NodeVector() ) {
            nit->Read( T_key, sc );
            if ( sc.Flag() == DIRICH ) ++n_dirich;
        }
        output( "Dirichlet temperature nodes: " + to_string(n_dirich) );

        // Range check
        const bool lambda_ok = model.IsWithinRange( "thermal conductivity" );
        output( "thermal conductivity within range: "
                + string(lambda_ok ? "YES" : "NO -- check configuration!") );

        const bool T_ok = model.IsWithinRange( "temperature" );
        output( "temperature within range: "
                + string(T_ok ? "YES" : "NO -- check configuration!") );
    }
    pause();

    keypoint(
"Always verify after ConfigureFromFile():\n"
"  printRangeOfVariable() -- check values are correct\n"
"  model.IsWithinRange()  -- check values are physical\n"
"  count flagged nodes    -- check BCs were applied\n"
"A missing region name in Block 3 is silent -- the value is not set."
    );
    pause();


    // ========================================================================
    // Step 15: Running the heat conduction solve
    // ========================================================================
    heading( 15, 15, "Running the heat conduction solve" );

    explain(
"Now we put it all together:\n"
"  1. Load the model\n"
"  2. Apply the configuration file\n"
"  3. Verify the configuration\n"
"  4. Solve the heat equation\n"
"  5. Output results\n"
"\n"
"This is the complete workflow for a CSMP simulation\n"
"driven by a configuration file."
    );
    pause();

    subheading( "Complete workflow" );
    code(
"// 1. Load model\n"
"Model<dim>  model( model_name, variable_file );\n"
"Region<dim>& gref( model.Region(\"Model\") );\n"
"\n"
"// 2. Apply configuration\n"
"InputDataManager<dim>  config;\n"
"config.ConfigureFromFile(\n"
"    model, model_name.c_str(),\n"
"    false, true, true, true, false, true );\n"
"\n"
"// 3. Verify\n"
"printRangeOfVariable( model, \"thermal conductivity\" );\n"
"\n"
"// 4. Solve\n"
"CSMP_DEFAULT_LINEAR_SOLVER   solver;\n"
"PDE_Integrator<dim,Element>  heat_conductor( solver );\n"
"NumIntegral_dNT_dN_dV<dim>   conductance( model.Database(),\n"
"    \"thermal conductivity\", \"temperature\", \"temperature\" );\n"
"NumIntegral_NT_rhsop_N_dV<dim> heat_source( model.Database(),\n"
"    \"heat source\", \"temperature\" );\n"
"NumIntegral_NT_op_dS<dim,Face> basal_flux( model.Database(),\n"
"    \"basal heat flow\", \"temperature\" );\n"
"heat_conductor.Add( &conductance );\n"
"heat_conductor.Add( &heat_source );\n"
"heat_conductor.AddBoundaryIntegral( &basal_flux );\n"
"heat_conductor.IntegrateOver( model, gref );\n"
"\n"
"// 5. Output\n"
"VTK_Interface<dim>  vtk;\n"
"vtk.OutputDataToVTK( model, \"temperature\", \"temperature\", 1 );"
    );
    pause();

    subheading( "Live solve" );
    {
        cout << "  Loading model and applying configuration...\n";
        cout.flush();

        Model<dim>  model( model_name, variable_file );
        Region<dim>& gref( model.Region("Model") );
        WriteConfigurationFile(model_name);
        
        // making sure the necessary parameters are there
        // NOTE: when defining parameters, name and notation must be unique
        if ( !model.Database().IsDefined("thermal conductivity") )
          model.CreateProperty("thermal conductivity", "K_t", "W m-2", SCALAR, ELEMENT );
        if ( !model.Database().IsDefined("heat source") )
          model.CreateProperty("heat source", "Q", "W m-3", SCALAR, ELEMENT );
        if ( !model.Database().IsDefined("basal heat flow") )
          model.CreateProperty("basal heat flow", "HFU", "W m-2", SCALAR, FACE );
        if ( !model.Database().IsDefined("temperature") )
          model.CreateProperty("temperature", "T", "oC", SCALAR, NODE );

        InputDataManager<dim>  config;
        config.ConfigureFromFile(
            model, model_name.c_str(),
            false, true, true, true, false, true );

        output( "Configuration applied. Solving..." );
        cout.flush();

#ifdef CSMP_WITH_SAMG_SOLVER
        SAMG_Settings  settings;
        settings.Set_eps(0.);
        SAMG_Solver                  solver( &settings );
#else
        CSMP_DEFAULT_LINEAR_SOLVER   solver;
#endif
        PDE_Integrator<dim,Element>  heat_conductor( solver );

        NumIntegral_dNT_lhsop_dN_dV<dim>  conductance(
            model.Database(),
            "thermal conductivity", "temperature", "temperature" );
        NumIntegral_NT_rhsop_N_dV<dim>  heat_source(
            model.Database(),
            "heat source", "temperature" );
        NumIntegral_NT_op_N_dS<dim>  basal_flux(
            model.Database(),
            "basal heat flow", "temperature" );

        heat_conductor.Add( &conductance );
        heat_conductor.Add( &heat_source );
        heat_conductor.AddBoundaryIntegral( &basal_flux );
        heat_conductor.IntegrateOver( model, gref );

        output( "Solve complete." );
        output( "Temperature range:" );
        printRangeOfVariable( model, "Model", "temperature" );
        output( "TOP boundary (should be ~15 degC):" );
        printRangeOfVariable( model, "TOP", "temperature" );
        output( "BOTTOM boundary:" );
        printRangeOfVariable( model, "BOTTOM", "temperature" );

        VTK_Interface<dim>  vtk;
        vtk.OutputDataToVTK( model, "temperature",          "temperature",          1 );
        vtk.OutputDataToVTK( model, "thermal-conductivity", "thermal conductivity", 1 );
        output( "Written: temperature_1.vtu" );
        output( "Written: thermal-conductivity_1.vtu" );
        output( "Open in ParaView to visualise the layered temperature field." );
    }
    pause();

    keypoint(
"The complete workflow is: load --> configure --> verify --> solve --> output.\n"
"ConfigureFromFile() replaces all the InputPropertyValue() calls\n"
"that were in Example 4's source code.\n"
"Changing parameters now requires only editing the text file."
    );
    pause();


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 5: Model Configuration\n\n"
         << "  Steps covered:\n\n"
         << "    1.  Why configuration files exist\n"
         << "    2.  The variable definition file\n"
         << "    3.  The PropertyDatabase and variable ranges\n"
         << "    4.  The configuration file structure\n"
         << "    5.  Block 2: default property values\n"
         << "    6.  Block 3: regional property values\n"
         << "    7.  Block 4: box boundary conditions\n"
         << "    8.  Block 5: essential conditions for regions\n"
         << "    9.  Block 6: arbitrary boundary conditions\n"
         << "   10.  Whitespace rules -- tabs only\n"
         << "   11.  Comments and end-of-file convention\n"
         << "   12.  ComputationalSettings\n"
         << "   13.  Calling ConfigureFromFile()\n"
         << "   14.  Verifying the configuration\n"
         << "   15.  Running the heat conduction solve\n\n"
         << "  Configuration file rules summary:\n\n"
         << "    - Tokens separated by TAB only (never spaces)\n"
         << "    - Spaces allowed WITHIN variable names\n"
         << "    - No trailing spaces at end of lines\n"
         << "    - Blocks separated by exactly ONE blank line\n"
         << "    - Comment lines start with #\n"
         << "    - After 5 blank lines: everything ignored\n\n"
         << "  Next: Example 6 -- Transient Simulation\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nInputDataManager_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp

