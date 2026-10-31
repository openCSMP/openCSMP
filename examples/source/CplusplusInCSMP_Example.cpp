
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CplusplusInCSMP_Example.h"
#include "exampleHelpers.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "CSMP_definitions.h"

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <functional>
#include <limits>

using namespace std;

namespace csmp {

// ============================================================================
void CplusplusInCSMP_Example::Specifications()
{
    SetTitle( "C++ in CSMP (interactive)" );
    SetDifficulty( 1 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: CplusplusInCSMP_Example.cpp" );
    AddDescription( "If you already know some C++, then move on to the C++ examples" );
    AddDescription( "Interactive introduction to the C++ language features "
                    "used throughout CSMP simulation code. "
                    "No mesh or file input required." );
    AddRequirement( "none" );
}



// ============================================================================
// Run
// ============================================================================

void CplusplusInCSMP_Example::Run()
{
    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 1\n"
         << "  C++ Concepts Used in CSMP Simulation Code\n\n"
         << "  This interactive tutorial introduces the C++ language\n"
         << "  features you will encounter when writing CSMP simulations.\n\n"
         << "  There are 15 short concepts to work through.\n"
         << "  Press ENTER after each step to move forward.\n";
    divider();
    pause( "  [ Press ENTER to begin ]" );


    // ========================================================================
    // Concept 1: Namespaces
    // ========================================================================
    {
        heading( 1, 15, "Namespaces" );

        explain(
"All CSMP classes and functions live in the 'csmp' namespace.\n"
"The C++ standard library lives in 'std'.\n"
"Namespaces prevent name collisions between libraries."
        );
        pause();

        subheading( "Two ways to use CSMP names" );
        code(
"// (a) Fully qualified -- always clear where the name comes from:\n"
"csmp::ScalarVariable  sc;\n"
"\n"
"// (b) Bring everything in -- less typing:\n"
"using namespace csmp;\n"
"ScalarVariable  sc;\n"
"\n"
"// In simulation files you will typically see both:\n"
"using namespace std;\n"
"using namespace csmp;"
        );
        pause();

        subheading( "Anonymous namespace -- file-local helpers" );
        explain(
"A namespace with no name makes its contents invisible outside\n"
"the current file. CSMP uses this for internal helper functions."
        );
        code(
"namespace {\n"
"    void helper() { ... }   // only visible in this .cpp file\n"
"} // anonymous namespace"
        );
        pause();

        keypoint(
"Use 'csmp::' when you want to be explicit about where a name\n"
"comes from. Use 'using namespace csmp;' in .cpp files to reduce\n"
"verbosity, but never in header files."
        );
        pause();
    }


    // ========================================================================
    // Concept 2: Template parameters and spatial dimension
    // ========================================================================
    {
        heading( 2, 15, "Template parameters -- spatial dimension" );

        explain(
"CSMP classes are templated on 'dim' (1, 2, or 3) and often on\n"
"'CELLTYPE' (Element, Face, InterFace, ...).\n"
"The compiler generates separate optimised code for each dimension."
        );
        pause();

        subheading( "Template class declaration" );
        code(
"template<uint32_t dim, template<uint32_t> class CELLTYPE = Element>\n"
"class PDE_Integrator { ... };\n"
"\n"
"// Instantiate for 2D and 3D:\n"
"PDE_Integrator<2U, Element>  integrator_2d;\n"
"PDE_Integrator<3U, Element>  integrator_3d;"
        );
        pause();

        subheading( "Always declare dim as constexpr" );
        code(
"constexpr uint32_t dim{ 3U };   // compile-time constant\n"
"                                 // NOT: const size_t dim(3U);"
        );

        constexpr uint32_t dim{ 3U };
        output( "dim = " + to_string(dim) );
        pause();

        subheading( "if constexpr -- compile-time branching" );
        code(
"if constexpr ( dim == 3U )\n"
"    // 3D branch -- compiled\n"
"else if constexpr ( dim == 2U )\n"
"    // 2D branch -- discarded at compile time\n"
"else\n"
"    // 1D branch -- discarded at compile time"
        );

        if constexpr ( dim == 3U )
            output( "dim == 3U: 3D branch compiled, others discarded" );

        pause();

        keypoint(
"'constexpr uint32_t dim{ 3U };' is the correct way to declare\n"
"the spatial dimension. It is required for template arguments.\n"
"'const size_t dim(3U)' looks similar but is subtly different."
        );
        pause();
    }


    // ========================================================================
    // Concept 3: constexpr and const
    // ========================================================================
    {
        heading( 3, 15, "constexpr and const" );

        explain(
"'constexpr' means the value is known at compile time.\n"
"'const'     means the value cannot change at runtime.\n"
"Both prevent accidental modification, but constexpr is stronger."
        );
        pause();

        subheading( "When to use each" );
        code(
"constexpr uint32_t dim{ 3U };          // compile-time: use for templates\n"
"const     string   model_name(\"mesh\"); // runtime: fixed after init\n"
"const     string   channel_region( \"Model\" );\n"
"\n"
"// constexpr is required for template arguments:\n"
"Model<dim>  model( model_name, \"vars.txt\" );  // OK\n"
"\n"
"// const prevents accidental modification:\n"
"channel_region = \"other\";  // COMPILE ERROR"
        );
        pause();

        subheading( "Live demonstration" );
        constexpr double pi{ 3.14159265358979 };
        const     string greeting{ "Hello from CSMP" };
        output( "constexpr double pi      = " + to_string(pi) );
        output( "const string     greeting = " + greeting );
        pause();

        keypoint(
"Prefer constexpr for template arguments and array sizes.\n"
"Prefer const for values that are fixed after initialisation.\n"
"Never use #define for numeric constants in C++ -- use constexpr."
        );
        pause();
    }


    // ========================================================================
    // Concept 4: auto
    // ========================================================================
    {
        heading( 4, 15, "auto -- type deduction" );

        explain(
"'auto' asks the compiler to deduce the type from the initialiser.\n"
"It reduces verbosity without sacrificing type safety.\n"
"The compiler still knows the exact type -- you just do not have to write it."
        );
        pause();

        subheading( "Basic usage" );
        code(
"auto x     = 3.14;                  // double\n"
"auto name  = string(\"pressure\");    // std::string\n"
"auto count = 42U;                   // unsigned int"
        );

        auto x     = 3.14;
        auto name  = string( "fluid pressure" );
        auto count = 42U;
        output( "auto x     = " + to_string(x)    + "  (double)"       );
        output( "auto name  = " + name             + "  (string)"       );
        output( "auto count = " + to_string(count) + "  (unsigned int)" );
        pause();

        subheading( "auto in CSMP loops" );
        code(
"// 'auto&' = reference -- no copy made (important for performance)\n"
"for ( auto& eit : gref.CellVector() )\n"
"    eit->Store( key, value );\n"
"\n"
"// 'const auto&' = read-only reference\n"
"for ( const auto& name : variable_names )\n"
"    cout << name << \"\\n\";"
        );
        pause();

        keypoint(
"Use 'auto&' in range-based for loops to avoid copying large objects.\n"
"Use 'const auto&' when you only need to read, not modify."
        );
        pause();
    }


    // ========================================================================
    // Concept 5: std::string
    // ========================================================================
    {
        heading( 5, 15, "std::string" );

        explain(
"std::string is used throughout CSMP for:\n"
"  - model names:    string model_name(\"one_sphere_0.45_tetra\");\n"
"  - variable names: string var(\"fluid pressure\");\n"
"  - file paths:     string file = model_name + \"-variables.txt\";\n"
"  - region names:   const string channel_region(\"Model\");"
        );
        pause();

        subheading( "Key operations" );
        code(
"string model_name( \"one_sphere_0.45_tetra\" );\n"
"string var_file   = model_name + \"-variables.txt\";  // concatenation\n"
"\n"
"model_name.length();   // number of characters\n"
"model_name.empty();    // true if length == 0\n"
"\n"
"// Checking for substrings -- used in CSMP region name detection:\n"
"if ( region_name.find(\"BOUNDARY\") != string::npos )\n"
"    cout << \"This is a boundary region\\n\";\n"
"\n"
"// Checking for empty user input:\n"
"if ( !user_input.empty() ) model_name = user_input;"
        );
        pause();

        subheading( "Live demonstration" );
        string model_name( "one_sphere_0.45_tetra" );
        string var_file   = model_name + "-variables.txt";
        output( "model_name          = " + model_name );
        output( "var_file            = " + var_file   );
        output( "model_name.length() = " + to_string(model_name.length()) );

        const string region_name( "BOUNDARY_LEFT" );
        if ( region_name.find("BOUNDARY") != string::npos )
            output( "'" + region_name + "' contains 'BOUNDARY'" );
        pause();

        keypoint(
"Variable names in CSMP are plain strings.\n"
"They must match the variable file exactly -- including spaces\n"
"and capitalisation. A typo will cause a runtime error."
        );
        pause();
    }


    // ========================================================================
    // Concept 6: std::vector
    // ========================================================================
    {
        heading( 6, 15, "std::vector" );

        explain(
"std::vector<T> is a dynamically sized array.\n"
"CSMP uses it for node lists, element lists, and solution vectors."
        );
        pause();

        subheading( "CSMP uses of vector" );
        code(
"vector<Node<dim>*>    node_list;    // node pointers\n"
"vector<Element<dim>*> cell_list;    // element pointers\n"
"vector<double>        solution;     // rh_, x_ in PDE_Integrator\n"
"vector<size_t>        dof_indexes;  // DOF_indexes_"
        );
        pause();

        subheading( "Key operations" );
        code(
"pressures.push_back( 5.1 );   // append one element\n"
"pressures.size();              // number of elements\n"
"pressures[2];                  // access by index (no bounds check)\n"
"pressures.at(2);               // access by index (with bounds check)\n"
"pressures.front();             // first element\n"
"pressures.back();              // last element\n"
"\n"
"// Pre-allocate to avoid repeated reallocation:\n"
"solution.reserve( 1000000 );"
        );
        pause();

        subheading( "Live demonstration" );
        vector<double> pressures{ 1.0, 2.5, 3.7, 0.8, 4.2 };
        pressures.push_back( 5.1 );
        output( "pressures.size()  = " + to_string(pressures.size())  );
        output( "pressures.front() = " + to_string(pressures.front()) );
        output( "pressures.back()  = " + to_string(pressures.back())  );
        output( "pressures[2]      = " + to_string(pressures[2])      );

        vector<double> solution;
        solution.reserve( 1000000 );
        output( "solution.capacity() after reserve(1e6) = "
                + to_string(solution.capacity()) );
        pause();

        keypoint(
"Always call reserve() before filling large vectors.\n"
"Without it, the vector may reallocate and copy its data\n"
"many times as it grows -- very slow for million-node meshes."
        );
        pause();
    }


    // ========================================================================
    // Concept 7: Range-based for loops
    // ========================================================================
    {
        heading( 7, 15, "Range-based for loops" );

        explain(
"The range-based for loop is the idiomatic way to iterate in CSMP.\n"
"It replaces the old iterator pattern and is cleaner and safer."
        );
        pause();

        subheading( "Old iterator style vs range-based for" );
        code(
"// OLD -- verbose and error-prone:\n"
"for ( auto eit = gref.CellsBegin(); eit != gref.CellsEnd(); ++eit )\n"
"    (*eit)->Store( key, value );\n"
"\n"
"// NEW -- clean and idiomatic:\n"
"for ( auto& eit : gref.CellVector() )\n"
"    eit->Store( key, value );"
        );
        pause();

        subheading( "Over nodes" );
        code(
"for ( auto& nit : gref.NodeVector() )\n"
"    nit->Store( key, value );"
        );
        pause();

        subheading( "Over a map" );
        code(
"// e.g. lhs_operators_ in PDE_Integrator:\n"
"for ( const auto& [name, op] : lhs_operators_ )\n"
"    cout << name << \"\\n\";"
        );
        pause();

        subheading( "Live demonstration" );
        const vector<string> variables{
            "fluid pressure",
            "parabolic function",
            "fluid pressure gradient",
            "Stokes Discrepancy Measure"
        };
        cout << "  Variables in this simulation:\n";
        for ( const auto& name : variables )
            output( name );
        pause();

        keypoint(
"'auto&' takes a reference -- no copy is made.\n"
"'const auto&' is read-only -- use it when you do not modify the element.\n"
"Always prefer range-based for over manual iterator loops in new code."
        );
        pause();
    }


    // ========================================================================
    // Concept 8: Structured bindings
    // ========================================================================
    {
        heading( 8, 15, "Structured bindings" );

        explain(
"Structured bindings (C++17) unpack pairs and tuples cleanly.\n"
"CSMP uses std::pair extensively in LocalVariableStorage and\n"
"localVariableDispatch."
        );
        pause();

        subheading( "Without structured bindings -- verbose" );
        code(
"const auto offsets    = localVariableDispatch::containerOffset(...);\n"
"const auto dataOffset = offsets.first;\n"
"const auto flagOffset = offsets.second;"
        );
        pause();

        subheading( "With structured bindings -- clean" );
        code(
"const auto [dataOffset, flagOffset] =\n"
"    localVariableDispatch::containerOffset( storeePtr, idx );\n"
"\n"
"const auto [nipCycles1, nipCycles2] =\n"
"    localVariableDispatch::containerIPCycles( storeePtr, idx );"
        );
        pause();

        subheading( "Live demonstration" );
        auto make_range = []( const string& var ) -> pair<double,double> {
            if ( var == "fluid pressure"    ) return { 0.0, 1.0  };
            if ( var == "parabolic function") return { 0.0, 0.25 };
            return { 0.0, 0.0 };
        };

        for ( const string& name : { "fluid pressure", "parabolic function" } ) {
            const auto [vmin, vmax] = make_range( name );
            output( name + "  range: ["
                    + to_string(vmin) + ", " + to_string(vmax) + "]" );
        }
        pause();

        keypoint(
"Structured bindings make pair-returning CSMP functions much easier\n"
"to read. Use them whenever a function returns a std::pair."
        );
        pause();
    }


    // ========================================================================
    // Concept 9: CSMP variable types and enum class
    // ========================================================================
    {
        heading( 9, 15, "CSMP variable types and enum class" );

        explain(
"CSMP uses strongly typed enumerations (enum class) to prevent\n"
"accidental mixing of unrelated integer values."
        );
        pause();

        subheading( "The three key enumerations" );
        code(
"VARIABLE_TYPE:  SCALAR, VECTOR, TENSOR, ARRAY, FLAGGEDARRAY\n"
"\n"
"PLACEMENT:      NODE, ELEMENT, FACE, INTER_FACE,\n"
"                ELEMENT_INTEGRATION_POINT, ...\n"
"\n"
"VARIABLE_FLAG:  PLAIN, DIRICH, NEUMANN, ROBIN, ANY, ..."
        );
        pause();

        subheading( "Creating variables" );
        code(
"ScalarVariable  sc( PLAIN,   1.0 );   // value + flag\n"
"ScalarVariable  bc( DIRICH,  0.0 );   // Dirichlet boundary condition\n"
"ScalarVariable  nc( NEUMANN, 0.0 );   // Neumann boundary condition\n"
"\n"
"// Factory functions (preferred):\n"
"makeScalar( DIRICH, 0. )              // returns ScalarVariable\n"
"makeVector( ANY, PLAIN, ROBIN, 1.,2.,3. )  // returns VectorVariable<3>"
        );
        pause();

        subheading( "Live demonstration -- ScalarVariable" );
        ScalarVariable  sc_plain(  PLAIN,   1.0 );
        ScalarVariable  sc_dirich( DIRICH,  0.0 );
        ScalarVariable  sc_neumann(NEUMANN, 0.0 );
        output( string("PLAIN   flag=") + parseStatus(sc_plain.Flag())
                + "  value=" + to_string(sc_plain()) );
        output( string("DIRICH  flag=") + parseStatus(sc_dirich.Flag())
                + "  value=" + to_string(sc_dirich()) );
        output( string("NEUMANN flag=") + parseStatus(sc_neumann.Flag())
                + "  value=" + to_string(sc_neumann()) );
        pause();

        subheading( "Live demonstration -- VectorVariable<3>" );
        constexpr uint32_t dim{ 3U };
        VectorVariable<dim>  vc( ANY, ANY, ANY, 1.0, 2.0, 3.0 );
        for ( uint32_t i{0U}; i < dim; ++i )
            output( "  component[" + to_string(i) + "] = "
                    + to_string(vc(i))
                    + "  flag = " + string(parseStatus(vc.Flag(i))) );
        pause();

        keypoint(
"The flag on a variable tells the solver how to treat it:\n"
"  DIRICH  = fixed value (Dirichlet boundary condition)\n"
"  NEUMANN = flux condition\n"
"  PLAIN / ANY = free (no special treatment)\n"
"Setting the wrong flag is a common source of simulation errors."
        );
        pause();
    }


    // ========================================================================
    // Concept 10: Scoped blocks
    // ========================================================================
    {
        heading( 10, 15, "Scoped blocks {}" );

        explain(
"In long simulation functions, scoped blocks {} limit the lifetime\n"
"of temporary variables to the block where they are needed.\n"
"This keeps the code readable and prevents accidental reuse."
        );
        pause();

        subheading( "Without scoped blocks -- risky" );
        code(
"TensorVariable<dim> ts;   // lives for the entire function\n"
"VectorVariable<dim> vc;   // easy to accidentally reuse below\n"
"\n"
"for ( auto& eit : gref.CellVector() ) {\n"
"    eit->Read( hessian.Key(), ts );\n"
"    eit->Read( grad.Key(),    vc );\n"
"    eit->Store( term1.Key(), ts * vc * 2.0 );\n"
"}\n"
"\n"
"// ts and vc still exist here -- could be confused with new ts, vc below"
        );
        pause();

        subheading( "With scoped blocks -- safe" );
        code(
"{\n"
"    TensorVariable<dim> ts;   // only exists inside this block\n"
"    VectorVariable<dim> vc;\n"
"    for ( auto& eit : gref.CellVector() ) {\n"
"        eit->Read( hessian.Key(), ts );\n"
"        eit->Read( grad.Key(),    vc );\n"
"        eit->Store( term1.Key(), ts * vc * 2.0 );\n"
"    }\n"
"}  // ts and vc destroyed here\n"
"\n"
"// Cannot accidentally use ts or vc below this point"
        );
        pause();

        subheading( "Live demonstration" );
        {
            constexpr uint32_t dim{ 3U };
            TensorVariable<dim>  ts_local;
            VectorVariable<dim>  vc_local;
            output( "Inside block: ts_local and vc_local exist" );
            (void)ts_local; (void)vc_local;
        }
        output( "Outside block: ts_local and vc_local are destroyed" );
        pause();

        keypoint(
"In CSMP examples, every step that processes element data uses\n"
"a scoped block. This is the standard pattern -- follow it."
        );
        pause();
    }


    // ========================================================================
    // Concept 11: Lambda functions
    // ========================================================================
    {
        heading( 11, 15, "Lambda functions" );

        explain(
"Lambdas are anonymous functions defined inline.\n"
"In CSMP they are used for concise per-element operations\n"
"and local helpers that do not need a name."
        );
        pause();

        subheading( "Syntax" );
        code(
"auto f = [capture]( parameters ) -> return_type { body };\n"
"\n"
"// Capture by value (copy of local variable):\n"
"auto scale = [factor]( double x ) { return x * factor; };\n"
"\n"
"// Capture by reference (no copy):\n"
"auto store = [&result]( double v ) { result.push_back(v); };\n"
"\n"
"// Capture everything by reference:\n"
"auto process = [&]( double v ) { total += v; count++; };"
        );
        pause();

        subheading( "In CSMP post-processing" );
        code(
"auto normalise = [p_min, p_max]( double val ) -> double {\n"
"    const double range = p_max - p_min;\n"
"    return ( range > numeric_limits<double>::epsilon() )\n"
"           ? ( val - p_min ) / range\n"
"           : 0.;\n"
"};\n"
"\n"
"for ( auto& eit : gref.CellVector() ) {\n"
"    eit->Read( pressure_key, sc );\n"
"    sc = makeScalar( PLAIN, normalise(sc()) );\n"
"    eit->Store( normalised_key, sc );\n"
"}"
        );
        pause();

        subheading( "Live demonstration" );
        vector<double> p{ 0.2, 0.8, 0.5, 1.0, 0.0, 0.6 };
        const double p_max = *max_element( p.begin(), p.end() );
        const double p_min = *min_element( p.begin(), p.end() );

        auto normalise = [p_min, p_max]( double val ) -> double {
            const double range = p_max - p_min;
            return ( range > numeric_limits<double>::epsilon() )
                   ? ( val - p_min ) / range : 0.;
        };

        cout << "  Original:   ";
        for ( double v : p )
            cout << fixed << setprecision(1) << setw(5) << v;
        cout << "\n  Normalised: ";
        for ( double v : p )
            cout << fixed << setprecision(2) << setw(5) << normalise(v);
        cout << "\n";
        pause();

        keypoint(
"Lambdas are ideal for short per-element operations that do not\n"
"need a named function. Keep them short -- if the body is more\n"
"than 5 lines, consider a named function instead."
        );
        pause();
    }


    // ========================================================================
    // Concept 12: Smart pointers
    // ========================================================================
    {
        heading( 12, 15, "Smart pointers" );

        explain(
"Smart pointers manage heap-allocated objects automatically.\n"
"They delete the object when it goes out of scope -- no manual delete."
        );
        pause();

        subheading( "Raw pointer -- manual management" );
        code(
"// Raw pointer -- YOU must call delete:\n"
"Solver* solver = new SAMG_Solver();  // allocate\n"
"// ... use solver ...\n"
"delete solver;                        // must not forget this!\n"
"\n"
"// If an exception is thrown before delete, memory leaks."
        );
        pause();

        subheading( "unique_ptr -- automatic management" );
        code(
"// unique_ptr -- deleted automatically:\n"
"auto solver = make_unique<SAMG_Solver>( &settings );\n"
"PDE_Integrator<dim,Element>  integrator( *solver );\n"
"// solver deleted here automatically -- even if an exception is thrown"
        );
        pause();

        subheading( "In CSMP examples -- prefer stack allocation" );
        code(
"// Stack allocation -- simplest, no delete needed:\n"
"SAMG_Solver                  solver( &settings );\n"
"PDE_Integrator<dim,Element>  integrator( solver );\n"
"\n"
"// Use unique_ptr only when the solver must outlive its creation scope."
        );
        pause();

        subheading( "Live demonstration" );
        auto value = make_unique<double>( 42.0 );
        output( "make_unique<double>(42.0)  *value = " + to_string(*value) );
        output( "Automatically deleted when it goes out of scope." );
        pause();

        keypoint(
"Prefer stack allocation for solvers in CSMP examples.\n"
"Use unique_ptr when ownership must be transferred or the object\n"
"must outlive its creation scope.\n"
"Never use raw new/delete in new CSMP code."
        );
        pause();
    }


    // ========================================================================
    // Concept 13: [[nodiscard]]
    // ========================================================================
    {
        heading( 13, 15, "[[nodiscard]]" );

        explain(
"[[nodiscard]] is a C++17 attribute that causes a compiler WARNING\n"
"if the return value of a function is silently ignored.\n"
"CSMP uses it on functions where ignoring the result is a bug."
        );
        pause();

        subheading( "CSMP functions marked [[nodiscard]]" );
        code(
"[[nodiscard]] bool   IsDefined( const char* ) const;\n"
"[[nodiscard]] size_t Nodes()    const noexcept;\n"
"[[nodiscard]] bool   EstablishMatrixSetup( ... );\n"
"[[nodiscard]] bool   CheckCellNumbering() const;"
        );
        pause();

        subheading( "Correct usage" );
        code(
"// WARNING -- return value ignored:\n"
"integrator.EstablishMatrixSetup( domain );\n"
"\n"
"// CORRECT -- check the return value:\n"
"if ( integrator.EstablishMatrixSetup( domain ) )\n"
"    cout << \"Matrix setup established for the first time\\n\";\n"
"else\n"
"    cout << \"Reusing existing matrix setup\\n\";"
        );
        pause();

        keypoint(
"If the compiler warns about a discarded return value, always\n"
"check it -- it usually indicates a logic error in your simulation."
        );
        pause();
    }


    // ========================================================================
    // Concept 14: static_assert
    // ========================================================================
    {
        heading( 14, 15, "static_assert" );

        explain(
"static_assert checks a condition at compile time.\n"
"If the condition is false, the compiler reports an error with\n"
"your custom message -- before the program even runs."
        );
        pause();

        subheading( "Usage in CSMP" );
        code(
"// Check that dim is valid:\n"
"static_assert( dim >= 1U && dim <= 3U,\n"
"               \"dim must be 1, 2, or 3\" );\n"
"\n"
"// Check that CELLTYPE is not a Region:\n"
"static_assert( !is_same<CELLTYPE<dim>, Region<dim>>::value,\n"
"               \"PDE_Integrator: CELLTYPE must be Element or Face\" );\n"
"\n"
"// Check variable type compatibility:\n"
"static_assert( TypeMatchesVariablePlacement<STOREE,place>::value,\n"
"               \"STOREE type does not match PLACEMENT\" );"
        );
        pause();

        subheading( "Live demonstration" );
        constexpr uint32_t dim{ 3U };
        static_assert( dim >= 1U && dim <= 3U, "dim must be 1, 2, or 3" );
        output( "static_assert( dim >= 1U && dim <= 3U ) -- passed at compile time" );
        output( "If dim were 0 or 4, this would be a COMPILE ERROR, not a runtime crash." );
        pause();

        keypoint(
"Use static_assert to document assumptions about template parameters.\n"
"The error appears at compile time with a clear message -- much better\n"
"than a mysterious crash at runtime."
        );
        pause();
    }


    // ========================================================================
    // Concept 15: A complete simulation skeleton
    // ========================================================================
    {
        heading( 15, 15, "A complete simulation skeleton" );

        explain(
"Every CSMP simulation follows this structure.\n"
"The letters (a)-(j) correspond to the concepts you have just learned."
        );
        pause();

        subheading( "(a) Compile-time constant and (b) model setup" );
        code(
"constexpr uint32_t dim{ 3U };                    // Concept 2, 3\n"
"const string model_name( \"my_mesh\" );            // Concept 5\n"
"Model<dim>   model( model_name, \"vars.txt\" );"
        );
        pause();

        subheading( "(c) Region and (d) property handles" );
        code(
"const string  domain_name( \"Model\" );            // Concept 5\n"
"Region<dim>&  gref( model.Region( domain_name.c_str() ) );\n"
"\n"
"PropertyHandle<dim> pressure( model, \"fluid pressure\",\n"
"                              SCALAR, NODE );     // Concept 9\n"
"PropertyHandle<dim> source(   model, \"source term\",\n"
"                              SCALAR, ELEMENT );\n"
"source = 1.;"
        );
        pause();

        subheading( "(e) Solver and PDE integrator" );
        code(
"CSMP_DEFAULT_LINEAR_SOLVER   solver;              // Concept 12\n"
"PDE_Integrator<dim,Element>  integrator( solver );\n"
"\n"
"NumIntegral_dNT_dN_dV<dim>     lhs( model.Database(),\n"
"                                    \"fluid pressure\",\n"
"                                    \"fluid pressure\" );\n"
"NumIntegral_NT_rhsop_N_dV<dim> rhs( model.Database(),\n"
"                                    \"source term\",\n"
"                                    \"fluid pressure\" );\n"
"integrator.Add( &lhs );\n"
"integrator.Add( &rhs );"
        );
        pause();

        subheading( "(f) Boundary conditions and (g) solve" );
        code(
"gref.InputPropertyValue( \"fluid pressure\",\n"
"                         makeScalar(DIRICH, 0.), PERIMETER ); // Concept 9\n"
"\n"
"integrator.IntegrateOver( model, gref );          // Concept 13"
        );
        pause();

        subheading( "(h) Post-process and (i) output" );
        code(
"{                                                 // Concept 10\n"
"    ScalarVariable sc;\n"
"    auto to_percent = []( double v ) { return v * 100.; }; // Concept 11\n"
"    for ( auto& eit : gref.CellVector() ) {       // Concept 7\n"
"        eit->Read( pressure.Key(), sc );\n"
"        // ... process sc ...\n"
"    }\n"
"}\n"
"\n"
"VTK_Interface<dim> vtk;\n"
"vtk.OutputDataToVTK( model, \"pressure\", \"fluid pressure\", 1 );"
        );
        pause();
        keypoint(
"This skeleton is the starting point for every CSMP simulation.\n"
"Example 2 will build it step by step with a real mesh."
        );
        pause();
    }


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 1: C++ in CSMP\n\n"
         << "  Concepts covered:\n\n"
         << "     1.  Namespaces\n"
         << "     2.  Template parameters and spatial dimension\n"
         << "     3.  constexpr and const\n"
         << "     4.  auto -- type deduction\n"
         << "     5.  std::string\n"
         << "     6.  std::vector\n"
         << "     7.  Range-based for loops\n"
         << "     8.  Structured bindings\n"
         << "     9.  CSMP variable types and enum class\n"
         << "    10.  Scoped blocks\n"
         << "    11.  Lambda functions\n"
         << "    12.  Smart pointers\n"
         << "    13.  [[nodiscard]]\n"
         << "    14.  static_assert\n"
         << "    15.  A complete simulation skeleton\n\n"
         << "  Next: Example 2 -- Building and Querying a CSMP Model\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

} // Run()

} // csmp
