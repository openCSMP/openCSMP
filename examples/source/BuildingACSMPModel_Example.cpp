// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/**
    @file BuildingACSMPModel_Example.cpp

    @brief CSMP for Beginners -- Example 2: Building and Querying a CSMP Model

    @section purpose Purpose

    This interactive example introduces the Model<dim> class -- the central
    object in every CSMP simulation. You will learn how to:

      1.  Load a model from CSMP binary files
      2.  Query mesh dimensions and topology
      3.  Access the PropertyDatabase
      4.  Work with Region, Boundary, and SplitBoundary objects
      5.  Read and write scalar, vector, and tensor variables
      6.  Use PropertyHandle to create and manage variables at runtime
      7.  Set boundary conditions
      8.  Iterate over nodes and elements
      9.  Output results to VTK

    @section model The Example Model

    The example uses 'one_sphere_0.45_tetra' -- a 3D pore-space mesh
    consisting of tetrahedral elements surrounding a sphere.
    It has:
      - a volumetric region  "Model"
      - box boundaries       LEFT, RIGHT, BOTTOM, TOP, FRONT, BACK
      - no split boundaries

    @section requirements Requirements

    - one_sphere_0.45_tetra  (CSMP binary files + -configuration.txt)
    - example25.txt          (variable definition file)

    @author S.K. Matthai
*/

#include "BuildingACSMPModel_Example.h"
#include "exampleHelpers.h"

#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"
#include "CSMP_highLevelUtilities.h"
#include "InputDataManager.h"

using namespace std;

namespace csmp {

// ============================================================================
// Specifications
// ============================================================================

void BuildingACSMPModel_Example::Specifications()
{
    SetTitle( "Building a CSMP Model (interactive)" );
    SetDifficulty( 1 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: BuildingACSMPModel_Example.cpp" );
    AddDescription( "Interactive introduction to Model<dim>: loading, "
                    "querying, reading and writing variables, "
                    "boundary conditions, and VTK output." );
    AddRequirement( "one_sphere_0.45_tetra (CSMP binary files, "
                    "-configuration.txt), example25.txt" );
}


// ============================================================================
// Run
// ============================================================================

void BuildingACSMPModel_Example::Run()
{
    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 2\n"
         << "  Building and Querying a CSMP Model\n\n"
         << "  You will load a real 3D pore-space mesh and learn\n"
         << "  how to work with it step by step.\n\n"
         << "  Steps:\n"
         << "    1.  Loading a model\n"
         << "    2.  Querying mesh dimensions\n"
         << "    3.  The PropertyDatabase\n"
         << "    4.  Regions and boundaries\n"
         << "    5.  Reading and writing variables\n"
         << "    6.  PropertyHandle\n"
         << "    7.  Setting boundary conditions\n"
         << "    8.  Iterating over nodes and elements\n"
         << "    9.  Output to VTK\n";
    divider();
    pause( "  [ Press ENTER to begin ]" );


    // -----------------------------------------------------------------------
    // Model setup (not a numbered step -- happens before the tutorial)
    // -----------------------------------------------------------------------
    constexpr uint32_t dim{ 3U };

    string model_name( "one_sphere_0.45_tetra" );
    cout << "\n  Enter model name or press ENTER for default '"
         << model_name << "':\n";
    string user_input;
    getline( cin, user_input );
    if ( !user_input.empty() ) model_name = user_input;

    const string file_name     = GetExampleFileName( __FILE__ );
    const string variable_file = "example25.txt";

    CreateWorkingDirectoryAndCopyInputModelFiles( file_name, model_name,
                                                  variable_file, model_name );


    // ========================================================================
    // Step 1: Loading a model
    // ========================================================================
    heading( 1, 9, "Loading a model" );

    explain(
"A CSMP Model is loaded from two sources:\n"
"  (a) Binary mesh files  -- store the mesh geometry and connectivity\n"
"  (b) A variable text file -- defines which physical variables exist\n"
"\n"
"The binary files are written by CSMP's mesh import tools.\n"
"The variable file is a plain text file you write yourself."
    );
    pause();

    subheading( "The two constructors you will use most" );
    code(
"// (a) Load mesh + variables from binary files:\n"
"Model<dim>  model( model_name, variable_file );\n"
"//                ^             ^\n"
"//                mesh files    variable definitions\n"
"\n"
"// (b) Load mesh only (no extra variables):\n"
"Model<dim>  model( model_name );"
    );
    pause();

    subheading( "Loading the model now..." );
    cout << "  Loading: " << model_name << " + " << variable_file << "\n";
    cout.flush();

    Model<dim>  model( model_name, variable_file );

    output( "Model loaded successfully." );
    output( "Name: " + string(model.Name()) );
    pause();

    keypoint(
"The Model constructor reads the mesh, builds the finite-element\n"
"connectivity, and initialises storage for all variables defined\n"
"in the variable file. This is the first line of every simulation."
    );
    pause();


    // ========================================================================
    // Step 2: Querying mesh dimensions
    // ========================================================================
    heading( 2, 9, "Querying mesh dimensions" );

    explain(
"Once the model is loaded, you can query its size and geometry.\n"
"These queries are used to check that the mesh is what you expect\n"
"before starting a computation."
    );
    pause();

    subheading( "Key query functions" );
    code(
"model.Mesh().Nodes();     // total number of nodes\n"
"model.Mesh().Elements();  // total number of elements\n"
"model.Mesh().Faces();     // total number of faces\n"
"\n"
"// Bounding box:\n"
"Point<dim> xyz_min, xyz_max;\n"
"model.MinMaxCoordinates( xyz_min, xyz_max );\n"
"\n"
"// Convenience function (also prints to screen):\n"
"printModelDimensions( model, true );\n\n"
"Model subdomains come as Regions, Boundaries and SplitBoundaries:\n"
"model.RegionsOut();       // unique (space exclusive-not overlapping) vs non-unique regions (e.g., 'Model')\n"
"model.BoundariesOut();    // lower dimensional ModelSubDomain objects consisting of Face objects\n"
"model.SplitBoundariesOut(); // ModelSubDomain objects connecting node-matched mesh patches via InterFaces\n"
    );
    pause();

    subheading( "Live query" );
    output( "Nodes:    " + to_string( model.Mesh().Nodes()    ) );
    output( "Elements: " + to_string( model.Mesh().Elements() ) );
    output( "Faces:    " + to_string( model.Mesh().Faces()    ) );

    {
        Point<dim> xyz_min, xyz_max;
        model.MinMaxCoordinates( xyz_min, xyz_max );
        output( "Bounding box min: ("
                + dbl_str(xyz_min[0]) + ", "
                + dbl_str(xyz_min[1]) + ", "
                + dbl_str(xyz_min[2]) + ")" );
        output( "Bounding box max: ("
                + dbl_str(xyz_max[0]) + ", "
                + dbl_str(xyz_max[1]) + ", "
                + dbl_str(xyz_max[2]) + ")" );
    }
    pause();

    subheading( "printModelDimensions" );
    printModelDimensions( model, true );
    pause();

    subheading( "regions, boundaries, their elements and volume/areas" );
    model.RegionsOut();
    model.BoundariesOut();

    keypoint(
"Always query the mesh dimensions after loading to verify that\n"
"the correct mesh was loaded and has the expected size.\n"
"A mesh with 0 nodes or elements indicates a loading error."
    );
    pause();


    // ========================================================================
    // Step 3: The PropertyDatabase
    // ========================================================================
    heading( 3, 9, "The PropertyDatabase" );

    explain(
"The PropertyDatabase stores the specification of every physical\n"
"variable in the model: its name, type, placement, and valid range.\n"
"\n"
"You access it via model.Database().\n"
"It is used to obtain the csmp::Index needed to read and write variables."
    );
    pause();

    subheading( "Key database operations" );
    code(
"// Check whether a variable exists:\n"
"model.Database().IsDefined( \"fluid pressure\" );\n"
"\n"
"// Get the Index for a variable (needed for Read/Store):\n"
"const csmp::Index key = model.Database().StorageKey( \"fluid pressure\" );\n"
"\n"
"// Print all variables to screen:\n"
"model.Database().FlushToScreen();\n"
"\n"
"// Create a new variable at runtime:\n"
"model.CreateProperty( \"my variable\", \"mv\", \"Pa\",\n"
"                      SCALAR, NODE, 1, -1e30, 1e30 );"
    );
    pause();

    subheading( "Live query -- variables in this model" );
    output( "IsDefined(\"fluid pressure\")    = "
            + string( model.Database().IsDefined("fluid pressure")
                      ? "true" : "false" ) );
    output( "IsDefined(\"nonexistent var\")   = "
            + string( model.Database().IsDefined("nonexistent var")
                      ? "true" : "false" ) );
    pause();

    subheading( "All variables defined in this model" );
    model.Database().FlushToScreen();
    pause();

    keypoint(
"The PropertyDatabase is the single source of truth for variable\n"
"specifications. Always use model.Database().IsDefined() to check\n"
"whether a variable exists before trying to access it."
    );
    pause();


    // ========================================================================
    // Step 4: Regions and boundaries
    // ========================================================================
    heading( 4, 9, "Regions and boundaries" );

    explain(
"A CSMP model is divided into named subdomains:\n"
"\n"
"  Region        -- volumetric subdomain (elements)\n"
"  Boundary      -- surface subdomain (faces)\n"
"  SplitBoundary -- internal interface (InterFaces)\n"
"\n"
"Every model has at least one Region called 'Model' which contains\n"
"all elements. Box-shaped models also have boundaries named\n"
"LEFT, RIGHT, BOTTOM, TOP, FRONT, BACK."
    );
    pause();

    subheading( "Accessing regions and boundaries" );
    code(
"// Get a reference to the entire model domain:\n"
"Region<dim>&  gref( model.Region( \"Model\" ) );\n"
"\n"
"// Get a reference to a boundary:\n"
"Boundary<dim>&  left( model.Boundary( \"LEFT\" ) );\n"
"\n"
"// Query sizes:\n"
"gref.Nodes();    // nodes in this region\n"
"gref.Cells();    // elements in this region\n"
"\n"
"// List all regions:\n"
"model.Regions();     // number of regions\n"
"model.Boundaries();  // number of boundaries"
    );
    pause();

    subheading( "Live query -- subdomains in this model" );
    {
        Region<dim>&  gref( model.Region( "Model" ) );
        output( "Region 'Model': nodes="
                + to_string(gref.Nodes())
                + "  elements=" + to_string(gref.Cells()) );
        output( "Total regions:    " + to_string(model.Regions())    );
        output( "Total boundaries: " + to_string(model.Boundaries()) );
    }
    pause();

    subheading( "Listing all boundaries" );
    code(
"for ( auto bit = model.BoundariesBegin();\n"
"      bit != model.BoundariesEnd(); ++bit )\n"
"    cout << bit->first << \"\\n\";"
    );
    for ( auto bit = model.BoundariesBegin();
          bit != model.BoundariesEnd(); ++bit )
        output( bit->first );
    pause();

    keypoint(
"Always store a reference to the region you are working in:\n"
"    Region<dim>& gref( model.Region( \"Model\" ) );\n"
"Then pass 'gref' to PDE_Integrator and post-processing loops.\n"
"This avoids repeated string lookups inside tight loops."
    );
    pause();


    // ========================================================================
    // Step 5: Reading and writing variables
    // ========================================================================
    heading( 5, 9, "Reading and writing variables" );

    explain(
"Variables are read and written via the csmp::Index obtained from\n"
"the PropertyDatabase. The Index encodes the variable's type,\n"
"placement, and storage location.\n"
"\n"
"Read() and Store() are called on individual Node or Element objects."
    );
    pause();

    subheading( "The read/write pattern" );
    code(
"// 1. Get the Index once (outside the loop):\n"
"const csmp::Index key = model.Database().StorageKey( \"fluid pressure\" );\n"
"\n"
"// 2. Read from a node or element:\n"
"ScalarVariable sc;\n"
"node->Read( key, sc );     // sc now holds the value and flag\n"
"double value = sc();       // extract the double value\n"
"VARIABLE_FLAG f = sc.Flag(); // extract the flag\n"
"\n"
"// 3. Write to a node or element:\n"
"sc = makeScalar( PLAIN, 42.0 );\n"
"node->Store( key, sc );"
    );
    pause();

    subheading( "Setting a uniform value across the whole model" );
    code(
"// Fastest way -- sets every node/element at once:\n"
"model.InputPropertyValue( \"fluid pressure\", makeScalar(PLAIN, 0.) );\n"
"\n"
"// Check the range:\n"
"printRangeOfVariable( model, \"Model\", \"fluid pressure\" );"
    );
    pause();

    subheading( "Live demonstration -- set and read fluid pressure" );
    {
        model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN, 0.) );
        output( "Set 'fluid pressure' = 0 everywhere." );
        printRangeOfVariable( model, "Model", "fluid pressure" );
    }
    pause();

    subheading( "Reading from individual nodes" );
    code(
"const csmp::Index key = model.Database().StorageKey( \"fluid pressure\" );\n"
"ScalarVariable sc;\n"
"Region<dim>& gref( model.Region( \"Model\" ) );\n"
"\n"
"// Read from the first node:\n"
"auto& nit = gref.NodeVector();\n"
"nit[0]->Read( key, sc );\n"
"cout << \"Node 0: value=\" << sc() << \"  flag=\" << sc.Flag() << \"\\n\";"
    );
    {
        const csmp::Index key = model.Database().StorageKey( "fluid pressure" );
        ScalarVariable sc;
        Region<dim>& gref( model.Region( "Model" ) );
        gref.NodeVector()[0]->Read( key, sc );
        output( "Node 0: value=" + dbl_str(sc())
                + "  flag=" + string(parseStatus(sc.Flag())) );
    }
    pause();

    keypoint(
"Get the csmp::Index ONCE outside the loop -- not inside.\n"
"StorageKey() involves a map lookup which is slow if called\n"
"millions of times inside an element loop."
    );
    pause();


    // ========================================================================
    // Step 6: PropertyHandle
    // ========================================================================
    heading( 6, 9, "PropertyHandle" );

    explain(
"PropertyHandle<dim> is a Maple/Matlab way to create and manage\n"
"distributed variables at runtime in CSMP simulation code.\n"
"\n"
"It:\n"
"  - creates the variable in the database if it does not already exist\n"
"  - provides a Key() method to get the csmp::Index\n"
"  - automatically deletes the variable when it goes out of scope\n"
"    (if it was created by the handle)\n"
"  - can be assigned a uniform value with operator="
    );
    pause();

    subheading( "Creating a PropertyHandle" );
    code(
"// Creates 'source term' as a SCALAR ELEMENT variable:\n"
"PropertyHandle<dim>  src( model, \"source term\", SCALAR, ELEMENT );\n"
"\n"
"// Assign a uniform value:\n"
"src = 1.;\n"
"\n"
"// Get the Index for use in Read/Store:\n"
"const csmp::Index& key = src.Key();\n"
"\n"
"// Creates a NODE variable:\n"
"PropertyHandle<dim>  pressure( model, \"fluid pressure\", SCALAR, NODE );"
    );
    pause();

    subheading( "Live demonstration" );
    {
        PropertyHandle<dim>  src( model, "source term", SCALAR, ELEMENT );
        src = 1.;
        output( "Created 'source term' as SCALAR ELEMENT." );
        output( "Set uniform value = 1." );
        printRangeOfVariable( model, "Model", "source term" );

        PropertyHandle<dim>  coord_x( model, "x coordinate", SCALAR, NODE );
        model.AssignNodeCoordinatesTo( "x coordinate", 'x' );
        output( "Created 'x coordinate' and filled with node x-coordinates." );
        printRangeOfVariable( model, "Model", "x coordinate" );
    }
    // coord_x and src are destroyed here -- variables removed from model
    output( "PropertyHandles went out of scope -- variables removed." );
    pause();

    subheading( "PropertyHandle vs CreateProperty" );
    code(
"// PropertyHandle (recommended):\n"
"PropertyHandle<dim>  src( model, \"source term\", SCALAR, ELEMENT );\n"
"// -- creates if new, deletes on destruction, provides Key()\n"
"\n"
"// CreateProperty (lower level):\n"
"csmp::Index key = model.CreateProperty( \"source term\", \"q\", \"1/s\",\n"
"                                        SCALAR, ELEMENT );\n"
"// -- you must call model.DeleteProperty() manually"
    );
    pause();

    keypoint(
"Always use PropertyHandle for variables that are local to a\n"
"simulation step. It prevents memory leaks and makes the code\n"
"self-documenting: the variable exists exactly as long as the handle."
    );
    pause();


    // ========================================================================
    // Step 7: Setting boundary conditions
    // ========================================================================
    heading( 7, 9, "Setting boundary conditions" );

    explain(
"Boundary conditions in CSMP are set by changing the VARIABLE_FLAG\n"
"of a variable at specific nodes or elements.\n"
"\n"
"  DIRICH  -- Dirichlet: the value is fixed\n"
"  NEUMANN -- Neumann: the flux is specified\n"
"  PLAIN   -- free: no special treatment\n"
"\n"
"The PDE_Integrator reads these flags and applies the conditions\n"
"automatically during assembly."
    );
    pause();

    subheading( "Setting conditions on a whole boundary" );
    code(
"// Set p = 1 (Dirichlet) on the FRONT boundary:\n"
"model.InputBoundaryValue( FRONT, \"fluid pressure\",\n"
"                          makeScalar(DIRICH, 1.) );\n"
"\n"
"// Set p = 0 (Dirichlet) on the BACK boundary:\n"
"model.InputBoundaryValue( BACK, \"fluid pressure\",\n"
"                          makeScalar(DIRICH, 0.) );"
    );
    pause();

    subheading( "Setting conditions on a named region perimeter" );
    code(
"Region<dim>& gref( model.Region( \"Model\" ) );\n"
"\n"
"// Set p = 0 on all perimeter nodes:\n"
"gref.InputPropertyValue( \"fluid pressure\",\n"
"                         makeScalar(DIRICH, 0.), PERIMETER );\n"
"\n"
"// Then release the condition at FRONT and BACK:\n"
"for ( auto nit = gref.PerimeterNodesBegin();\n"
"      nit != gref.NodesEnd(); ++nit )\n"
"    if ( (*nit)->AtBoundary() == FRONT ||\n"
"         (*nit)->AtBoundary() == BACK  )\n"
"        (*nit)->Status( pressure_key, PLAIN );"
    );
    pause();

    subheading( "Live demonstration" );
    {
        PropertyHandle<dim>  pressure( model, "fluid pressure", SCALAR, NODE );

        // Set p = 0 everywhere first
        model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN, 0.) );

        // Apply Dirichlet conditions
        model.InputBoundaryValue( FRONT, "fluid pressure",
                                  makeScalar(DIRICH, 1.) );
        model.InputBoundaryValue( BACK,  "fluid pressure",
                                  makeScalar(DIRICH, 0.) );

        output( "Set DIRICH p=1 on FRONT, DIRICH p=0 on BACK." );

        // Count Dirichlet nodes
        const csmp::Index key = model.Database().StorageKey( "fluid pressure" );
        size_t n_dirich{ 0U };
        ScalarVariable sc;
        Region<dim>& gref( model.Region( "Model" ) );
        for ( const auto& nit : gref.NodeVector() ) {
            nit->Read( key, sc );
            if ( sc.Flag() == DIRICH ) ++n_dirich;
        }
        output( "Dirichlet nodes: " + to_string(n_dirich) );
    }
    pause();

    keypoint(
"Boundary conditions are just variable flags.\n"
"DIRICH tells the assembler: 'do not solve for this DOF -- use\n"
"the stored value instead.' PLAIN means 'solve for this DOF.'\n"
"Setting the wrong flag is the most common simulation error."
    );
    pause();


    // ========================================================================
    // Step 8: Iterating over nodes and elements
    // ========================================================================
    heading( 8, 9, "Iterating over nodes and elements" );

    explain(
"Post-processing in CSMP is done by looping over the NodeVector()\n"
"or CellVector() of a Region and calling Read() or Store() on each\n"
"node or element.\n"
"\n"
"Always use range-based for with 'auto&' to avoid copying."
    );
    pause();

    subheading( "Element loop pattern" );
    code(
"{\n"
"    ScalarVariable sc;\n"
"    const csmp::Index key = model.Database().StorageKey( \"source term\" );\n"
"    Region<dim>& gref( model.Region( \"Model\" ) );\n"
"\n"
"    for ( auto& eit : gref.CellVector() ) {\n"
"        eit->Read( key, sc );          // read current value\n"
"        sc = makeScalar( PLAIN, 2. ); // modify\n"
"        eit->Store( key, sc );         // write back\n"
"    }\n"
"}"
    );
    pause();

    subheading( "Node loop pattern" );
    code(
"{\n"
"    ScalarVariable sc;\n"
"    const csmp::Index key = model.Database().StorageKey( \"fluid pressure\" );\n"
"    Region<dim>& gref( model.Region( \"Model\" ) );\n"
"\n"
"    for ( const auto& nit : gref.NodeVector() ) {\n"
"        nit->Read( key, sc );          // read only\n"
"        // ... process sc() ...\n"
"    }\n"
"}"
    );
    pause();

    subheading( "Live demonstration -- compute mean pressure" );
    {
        // Set a non-uniform pressure for demonstration
        model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN, 0.) );
        model.InputBoundaryValue( FRONT, "fluid pressure",
                                  makeScalar(PLAIN, 1.) );
        model.InputBoundaryValue( BACK,  "fluid pressure",
                                  makeScalar(PLAIN, 0.) );

        const csmp::Index key = model.Database().StorageKey( "fluid pressure" );
        Region<dim>& gref( model.Region( "Model" ) );

        double sum{ 0. };
        size_t count{ 0U };
        ScalarVariable sc;

        for ( const auto& nit : gref.NodeVector() ) {
            nit->Read( key, sc );
            sum += sc();
            ++count;
        }

        const double mean = ( count > 0 ) ? sum / static_cast<double>(count)
                                          : 0.;
        output( "Nodes iterated: " + to_string(count) );
        output( "Mean pressure:  " + dbl_str(mean) );
    }
    pause();

    subheading( "Perimeter node loop -- boundary condition release" );
    code(
"// ModelSubDomain objects have INTERIOR and PERIMETER cells and nodes; latter have faces on their boundary.\n"
"// PerimeterNodesBegin() starts at the first perimeter node.\n"
"// NodesEnd() is the end of ALL nodes (interior + perimeter).\n"
"// So this loop covers only perimeter nodes:\n"
"\n"
"for ( auto nit = gref.PerimeterNodesBegin();\n"
"      nit != gref.NodesEnd(); ++nit )\n"
"    if ( (*nit)->AtBoundary() == FRONT )\n"
"        (*nit)->Status( key, ANY );"
    );
    pause();

    keypoint(
"Get the csmp::Index ONCE before the loop.\n"
"Use 'const auto&' in the loop when you only read.\n"
"Use 'auto&' when you also write.\n"
"Use scoped blocks {} to limit the lifetime of variables involved."
    );
    pause();


    // ========================================================================
    // Step 9: Output to VTK
    // ========================================================================
    heading( 9, 9, "Output to VTK" );

    explain(
"CSMP can write results to VTK/VTU(XML) files for visualisation in ParaView\n"
"or VisIt. The VTK_Interface and VTU_Interface classes handle this.\n"
"\n"
"Each call writes one variable to one .vtu file.\n"
"The file name is constructed from the label you supply."
    );
    pause();

    subheading( "VTK output pattern" );
    code(
"VTK_Interface<dim>  vtk;\n"
"\n"
"// Write a scalar variable:\n"
"vtk.OutputDataToVTK( model, \"pressure\", \"fluid pressure\", 1 );\n"
"//                          ^label        ^variable name     ^step\n"
"\n"
"// Write a vector variable:\n"
"vtk.OutputDataToVTK( model, \"velocity\", \"fluid velocity\", 1 );\n"
"\n"
"// The output file will be named: pressure_1.vtk"
    );
    pause();

    subheading( "Live demonstration" );
    {
        // Create a simple scalar field to visualise
        PropertyHandle<dim>  coord_x( model, "x coordinate", SCALAR, NODE );
        model.AssignNodeCoordinatesTo( "x coordinate", 'x' );

        VTK_Interface<dim>  vtk;
        vtk.OutputDataToVTK( model, "x-coordinate", "x coordinate", 1 );
        output( "Written: x-coordinate_1.vtu" );
        output( "Open this file in ParaView to see the x-coordinate field." );
    }
    pause();

    keypoint(
"VTK/VTU output is the primary way to visualise CSMP results.\n"
"The label becomes the file name prefix.\n"
"The step number (last argument) is appended to the file name\n"
"and is used for time series output."
    );
    pause();


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 2: Building and Querying a CSMP Model\n\n"
         << "  Steps covered:\n\n"
         << "    1.  Loading a model\n"
         << "    2.  Querying mesh dimensions\n"
         << "    3.  The PropertyDatabase\n"
         << "    4.  Regions and boundaries\n"
         << "    5.  Reading and writing variables\n"
         << "    6.  PropertyHandle\n"
         << "    7.  Setting boundary conditions\n"
         << "    8.  Iterating over nodes and elements\n"
         << "    9.  Output to VTK\n\n"
         << "  Next: Example 3 -- Your First PDE Solve\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nBuildingACSMPModel_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp

