// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/**
    @file CSMP_PhysicalVariables_Example.cpp

    @brief CSMP for Beginners -- Example 3: Physical Variables

    @section purpose Purpose

    This interactive example introduces the five CSMP variable types
    and explains how they are used in simulation code.

    @section concepts Concepts Covered

      1.  Why variables have flags
      2.  ScalarVariable -- value + flag
      3.  VectorVariable -- dim components, each with a flag
      4.  TensorVariable -- dim x dim matrix, diagonal flags
      5.  ArrayVariable  -- runtime-sized array, single flag
      6.  FlaggedArrayVariable -- runtime-sized array, per-element flags
      7.  Variable placements -- NODE, ELEMENT, ELEMENT_INTEGRATION_POINT
      8.  Interpolation -- node to element, node to integration point
      9.  Extrapolation -- element to node, integration point to node
     10.  Assignment to DenseMatrix

    @section requirements Requirements

    - one_sphere_0.45_tetra  (CSMP binary files + -configuration.txt)
    - example25.txt          (variable definition file)

    @author S.K. Matthai
*/

#include "CSMP_PhysicalVariables_Example.h"
#include "exampleHelpers.h"

#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "DenseMatrix.h"
#include "VTK_Interface.h"
#include "CSMP_highLevelUtilities.h"
#include "InputDataManager.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

using namespace std;

namespace csmp {

// ============================================================================
// Specifications
// ============================================================================

void CSMP_PhysicalVariables_Example::Specifications()
{
    SetTitle( "CSMP Physical Variables (interactive)" );
    SetDifficulty( 1 );
    SetCategory( "CSMP for Beginners" );
    AddAuthor( "SKM" );
    AddDescription( "source in: CSMP_PhysicalVariables_Example.cpp" );
    AddDescription( "Interactive introduction to CSMP variable types, "
                    "flags, placements, interpolation, extrapolation, "
                    "and DenseMatrix assignment." );
    AddRequirement( "one_sphere_0.45_tetra (CSMP binary files, "
                    "-configuration.txt), example25.txt" );
}


// ============================================================================
// Run
// ============================================================================

void CSMP_PhysicalVariables_Example::Run()
{
    // -----------------------------------------------------------------------
    // Welcome screen
    // -----------------------------------------------------------------------
    blank();
    divider();
    cout << "\n"
         << "  CSMP for Beginners -- Example 3\n"
         << "  Physical Variables in CSMP\n\n"
         << "  You will learn about the five CSMP variable types,\n"
         << "  why they carry flags, how they are placed on the mesh,\n"
         << "  how to interpolate and extrapolate between placements,\n"
         << "  and how to assemble them into a DenseMatrix.\n\n"
         << "  Steps:\n"
         << "    1.  Why variables have flags\n"
         << "    2.  ScalarVariable\n"
         << "    3.  VectorVariable\n"
         << "    4.  TensorVariable\n"
         << "    5.  ArrayVariable\n"
         << "    6.  FlaggedArrayVariable\n"
         << "    7.  Variable placements\n"
         << "    8.  Interpolation\n"
         << "    9.  Extrapolation\n"
         << "   10.  Assignment to DenseMatrix\n";
    divider();
    pause( "  [ Press ENTER to begin ]" );


    // -----------------------------------------------------------------------
    // Load model (needed for steps 7-10)
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

    cout << "\n  Loading model..." << endl;
    Model<dim>  model( model_name, variable_file );
    Region<dim>& gref( model.Region( "Model" ) );
    output( "Model loaded: " + to_string(model.Mesh().Nodes())
            + " nodes, " + to_string(model.Mesh().Elements()) + " elements." );
    pause();


    // ========================================================================
    // Step 1: Why variables have flags
    // ========================================================================
    heading( 1, 10, "Why variables have flags" );

    explain(
"In a finite-element simulation, the same variable (e.g. 'fluid pressure')\n"
"plays different roles at different locations:\n"
"\n"
"  - At most nodes:    it is an UNKNOWN to be solved for  (flag = PLAIN)\n"
"  - At inlet nodes:   it is a FIXED VALUE                (flag = DIRICH)\n"
"  - At outlet nodes:  it is a FIXED VALUE                (flag = DIRICH)\n"
"  - At wall nodes:    it may be a FLUX condition         (flag = NEUMANN)\n"
"\n"
"The flag is stored WITH the variable value -- not with the node.\n"
"This is essential because a node may have multiple variables\n"
"(pressure, temperature, concentration) each with different flags."
    );
    pause();

    subheading( "The VARIABLE_FLAG enumeration" );
    code(
"enum VARIABLE_FLAG {\n"
"    ANY,          // no special treatment (default after construction)\n"
"    PLAIN,        // free variable -- solver will compute its value\n"
"    DIRICH,       // Dirichlet: value is fixed, DOF is eliminated\n"
"    NEUMANN,      // Neumann: flux is specified on boundary\n"
"    ROBIN,        // Robin: mixed condition\n"
"    INIT_COND,    // initial condition (transient problems)\n"
"    INIT_GUESS,   // initial guess for iterative solver\n"
"    PERIODIC,     // periodic boundary condition\n"
"    FIELD_DATA,   // measured/imported data -- do not overwrite\n"
"    CONSTANT_FLUX // constant flux condition\n"
"};"
    );
    pause();

    subheading( "How the PDE assembler uses flags" );
    explain(
"When PDE_Integrator assembles the global matrix, it checks the flag\n"
"of each variable at each node:\n"
"\n"
"  DIRICH  --> eliminate this DOF from the system (row/column removed)\n"
"  PLAIN   --> include this DOF in the system (row/column kept)\n"
"  NEUMANN --> add flux contribution to the right-hand side vector\n"
"\n"
"This is why setting the wrong flag or conflicting combination of flags is a common simulation error:\n"
"without DIRICH constraints on the computation domain the solver will not find a unique solution,\n"
"or the system may be singular."
    );
    pause();

    keypoint(
"Every CSMP variable carries a flag that tells the assembler\n"
"how to treat that degree of freedom of the solution variable.\n"
"DIRICH = fixed value (eliminated from system).\n"
"PLAIN  = free variable (included in system).\n"
"Setting the wrong flag = wrong physics."
    );
    pause();


    // ========================================================================
    // Step 2: ScalarVariable
    // ========================================================================
    heading( 2, 10, "ScalarVariable" );

    explain(
"ScalarVariable stores a single double value and a VARIABLE_FLAG.\n"
"It behaves like a double for arithmetic but carries the flag\n"
"for finite-element treatment.\n"
"\n"
"Storage: 8 bytes (double) + 1 byte (flag) = 9 bytes per variable."
    );
    pause();

    subheading( "Construction and access" );
    code(
"ScalarVariable  sc;                    // default: NaN, flag=ANY\n"
"ScalarVariable  sc( PLAIN,   1.0 );    // value=1.0, flag=PLAIN\n"
"ScalarVariable  sc( DIRICH,  0.0 );    // value=0.0, flag=DIRICH\n"
"\n"
"double value = sc();                   // extract value\n"
"VARIABLE_FLAG f = sc.Flag();           // extract flag\n"
"sc.Flag() = NEUMANN;                   // change flag\n"
"sc = 42.0;                             // change value, flag unchanged\n"
"\n"
"// Factory function (preferred in function arguments):\n"
"makeScalar( DIRICH, 0. )               // returns ScalarVariable"
    );
    pause();

    subheading( "Arithmetic -- behaves like double" );
    code(
"ScalarVariable  a( PLAIN, 3.0 );\n"
"ScalarVariable  b( PLAIN, 2.0 );\n"
"\n"
"auto c = a + b;    // c = 5.0, flag from a\n"
"auto d = a * 2.0;  // d = 6.0\n"
"a += 1.0;          // a = 4.0\n"
"a /= b;            // a = 2.0\n"
"\n"
"// Range check:\n"
"a.IsWithinRange( 0., 10. );  // true\n"
"a.Has_NaN_Values();          // false"
    );
    pause();

    subheading( "Live demonstration" );
    {
        ScalarVariable  a( PLAIN,  3.0 );
        ScalarVariable  b( PLAIN,  2.0 );
        ScalarVariable  bc( DIRICH, 0.0 );

        output( "a = " + dbl_str(a()) + "  flag=" + string(parseStatus(a.Flag())) );
        output( "b = " + dbl_str(b()) + "  flag=" + string(parseStatus(b.Flag())) );
        output( "bc = " + dbl_str(bc()) + "  flag=" + string(parseStatus(bc.Flag())) );

        auto c = a + b;
        output( "a + b = " + dbl_str(c()) + "  flag=" + string(parseStatus(c.Flag())) );

        a *= 2.0;
        output( "a *= 2.0  -->  a = " + dbl_str(a()) );

        output( "a.IsWithinRange(0,10) = "
                + string(a.IsWithinRange(0.,10.) ? "true" : "false") );
        output( "bc.IsWithinRange(0,10) = "
                + string(bc.IsWithinRange(0.,10.) ? "true" : "false") );
    }
    pause();

    keypoint(
"ScalarVariable is the most common CSMP variable type.\n"
"Use makeScalar(flag, value) when passing to functions.\n"
"The flag travels with the value -- never separate them."
    );
    pause();


    // ========================================================================
    // Step 3: VectorVariable
    // ========================================================================
    heading( 3, 10, "VectorVariable<dim>" );

    explain(
"VectorVariable<dim> stores dim double values, each with its own flag.\n"
"This models physical vectors like velocity, displacement, or gradient.\n"
"\n"
"Each component can have a different flag -- for example, the x-component\n"
"of displacement may be fixed (DIRICH) while y and z are free (PLAIN)."
    );
    pause();

    subheading( "Construction and access" );
    code(
"// All components same flag and value:\n"
"VectorVariable<3U>  vc( PLAIN, 0.0 );\n"
"\n"
"// Independent flags and values:\n"
"VectorVariable<3U>  vc( DIRICH, PLAIN, PLAIN,\n"
"                        0.0,    1.0,   2.0 );\n"
"\n"
"double vx = vc(0);          // read component 0\n"
"vc(1) = 5.0;                // write component 1\n"
"VARIABLE_FLAG f = vc.Flag(0); // flag of component 0\n"
"vc.Flag(1) = DIRICH;          // change flag of component 1\n"
"\n"
"double len = vc.Length();   // Euclidean length\n"
"vc.EuclideanNormalize();    // normalise to unit length\n"
"\n"
"// Factory function:\n"
"makeVector( PLAIN, PLAIN, PLAIN, 1., 2., 3. )"
    );
    pause();

    subheading( "Vector operations" );
    code(
"VectorVariable<3U>  a( PLAIN, PLAIN, PLAIN, 1., 0., 0. );\n"
"VectorVariable<3U>  b( PLAIN, PLAIN, PLAIN, 0., 1., 0. );\n"
"\n"
"double dot   = a & b;       // dot product (operator&)\n"
"auto   cross = a % b;       // cross product (operator%)\n"
"auto   proj  = a.ProjectOnto(b); // projection of a onto b\n"
"\n"
"// Arithmetic:\n"
"auto c = a + b;             // element-wise addition\n"
"auto d = a * 2.0;           // scale all components\n"
"a += b;                     // in-place addition"
    );
    pause();

    subheading( "Live demonstration" );
    {
        VectorVariable<dim>  a( PLAIN, PLAIN, PLAIN, 1., 0., 0. );
        VectorVariable<dim>  b( PLAIN, PLAIN, PLAIN, 0., 1., 0. );
        VectorVariable<dim>  bc( DIRICH, PLAIN, PLAIN, 0., 1., 2. );

        output( "a = (" + dbl_str(a(0)) + ", " + dbl_str(a(1)) + ", " + dbl_str(a(2)) + ")" );
        output( "b = (" + dbl_str(b(0)) + ", " + dbl_str(b(1)) + ", " + dbl_str(b(2)) + ")" );
        output( "bc flags: (" + string(parseStatus(bc.Flag(0))) + ", "
                              + string(parseStatus(bc.Flag(1))) + ", "
                              + string(parseStatus(bc.Flag(2))) + ")" );

        const double dot = a & b;
        output( "a . b (dot product) = " + dbl_str(dot) );

        auto cross = a % b;
        output( "a x b (cross product) = ("
                + dbl_str(cross(0)) + ", "
                + dbl_str(cross(1)) + ", "
                + dbl_str(cross(2)) + ")" );

        output( "a.Length() = " + dbl_str(a.Length()) );

        auto c = a + b;
        output( "a + b = ("
                + dbl_str(c(0)) + ", "
                + dbl_str(c(1)) + ", "
                + dbl_str(c(2)) + ")" );
    }
    pause();

    keypoint(
"Each component of VectorVariable has its own flag.\n"
"This is essential for vector Dirichlet conditions:\n"
"e.g. fix x-displacement (DIRICH) but free y and z (PLAIN).\n"
"Use operator& for dot product, operator% for cross product."
    );
    pause();


    // ========================================================================
    // Step 4: TensorVariable
    // ========================================================================
    heading( 4, 10, "TensorVariable<dim>" );

    explain(
"TensorVariable<dim> stores a dim x dim matrix.\n"
"Only the diagonal elements carry flags -- the off-diagonal elements\n"
"are always treated as free.\n"
"\n"
"Used for: permeability, stress, strain, diffusivity tensors."
    );
    pause();

    subheading( "Construction and access" );
    code(
"// Isotropic diagonal tensor (all diagonal = val):\n"
"TensorVariable<3U>  ts( PLAIN, 1.0 );\n"
"\n"
"// Full tensor with independent diagonal flags:\n"
"TensorVariable<3U>  ts( PLAIN, PLAIN, PLAIN,\n"
"                        1., 0., 0.,\n"
"                        0., 1., 0.,\n"
"                        0., 0., 1. );  // identity\n"
"\n"
"double v = ts(i,j);         // read element (i,j)\n"
"ts(1,2) = 0.5;              // write element (1,2)\n"
"VARIABLE_FLAG f = ts.Flag(0); // flag of diagonal element 0\n"
"\n"
"double tr  = ts.Trace();    // sum of diagonal\n"
"double det = ts.Determinant();\n"
"auto   inv = ts.Inverse();\n"
"auto   tr2 = ts.Transposed();"
    );
    pause();

    subheading( "Three distinct squaring operations" );
    code(
"// 1. Hadamard (element-wise): T_ij = T_ij^2\n"
"ts.HadamardSquared();\n"
"\n"
"// 2. Matrix product: result_ij = sum_k T_ik * T_kj\n"
"ts.MatrixSquared();\n"
"\n"
"// 3. Double contraction (Frobenius norm squared): T:T = sum_ij T_ij^2\n"
"double frob = ts.DoubleContraction();\n"
"\n"
"// Matrix-vector multiplication:\n"
"VectorVariable<3U>  result = ts * vc;"
    );
    pause();

    subheading( "Live demonstration" );
    {
        // Anisotropic permeability tensor
        TensorVariable<dim>  K( PLAIN, PLAIN, PLAIN,
                                1e-12, 0.,    0.,
                                0.,    1e-13, 0.,
                                0.,    0.,    1e-14 );

        output( "Anisotropic permeability tensor K:" );
        for ( uint32_t i{0U}; i < dim; ++i ) {
            ostringstream row;
            row << "  row " << i << ": ";
            for ( uint32_t j{0U}; j < dim; ++j )
                row << scientific << setprecision(2) << setw(12) << K(i,j) << " ";
            output( row.str() );
        }

        // Trace: sum of diagonal = 1e-12 + 1e-13 + 1e-14 = 1.11e-12
        output( "Trace(K)       = " + sci_str(K.Trace(), 4) );

        // Determinant: product of diagonal = 1e-12 * 1e-13 * 1e-14 = 1e-39
        output( "Determinant(K) = " + sci_str(K.Determinant(), 4) );

        // Matrix-vector product: K * grad_p
        // grad_p = (1,0,0) so result = (K_00*1, K_11*0, K_22*0) = (1e-12, 0, 0)
        VectorVariable<dim>  grad_p( PLAIN, PLAIN, PLAIN, 1., 0., 0. );
        auto velocity = K * grad_p;
        output( "K * (1,0,0) = ("
                + sci_str(velocity(0),3) + ", "
                + sci_str(velocity(1),3) + ", "
                + sci_str(velocity(2),3) + ")" );

        // Use a more interesting gradient to show all components
        VectorVariable<dim>  grad_p2( PLAIN, PLAIN, PLAIN, 1., 1., 1. );
        auto velocity2 = K * grad_p2;
        output( "K * (1,1,1) = ("
                + sci_str(velocity2(0),3) + ", "
                + sci_str(velocity2(1),3) + ", "
                + sci_str(velocity2(2),3) + ")" );
        output( "  (shows anisotropy: x-flux >> y-flux >> z-flux)" );

        // Eigenvalues of a diagonal tensor = its diagonal entries
        VectorVariable<dim>  eigenvals;
        if ( K.EigenValues( eigenvals ) )
            output( "Eigenvalues: ("
                    + sci_str(eigenvals(0),3) + ", "
                    + sci_str(eigenvals(1),3) + ", "
                    + sci_str(eigenvals(2),3) + ")" );
        output( "  (eigenvalues of diagonal tensor = diagonal entries)" );

        // Demonstrate Trace() on identity tensor
        TensorVariable<dim>  I( PLAIN, 1.0 );  // identity
        output( "Trace(Identity) = " + dbl_str(I.Trace(), 1)
                + "  (= dim = " + to_string(dim) + ")" );

        // Demonstrate DoubleContraction on identity
        output( "Identity:DoubleContraction = " + dbl_str(I.DoubleContraction(), 1)
                + "  (= number of non-zero entries = " + to_string(dim) + ")" );

        // Demonstrate Inverse
        auto K_inv = K.Inverse();
        output( "K_inv diagonal: ("
                + sci_str(K_inv(0,0),3) + ", "
                + sci_str(K_inv(1,1),3) + ", "
                + sci_str(K_inv(2,2),3) + ")" );
        output( "  (inverse of diagonal tensor = 1/diagonal)" );

        // Verify K * K_inv = Identity
        auto KKinv = K * K_inv;
        output( "K * K_inv diagonal (should be 1,1,1): ("
                + dbl_str(KKinv(0,0),6) + ", "
                + dbl_str(KKinv(1,1),6) + ", "
                + dbl_str(KKinv(2,2),6) + ")" );
    }
    pause();

    keypoint(
"TensorVariable flags are on the DIAGONAL only.\n"
"Use Trace() for the Laplacian (sum of diagonal = div(grad)).\n"
"Use operator* for matrix-vector multiplication.\n"
"Three squaring operations exist -- choose the right one:\n"
"  HadamardSquared = element-wise,  MatrixSquared = matrix product,\n"
"  DoubleContraction = scalar T:T."
    );
    pause();


    // ========================================================================
    // Step 5: ArrayVariable
    // ========================================================================
    heading( 5, 10, "ArrayVariable" );

    explain(
"ArrayVariable stores an arbitrary number of doubles with a SINGLE flag\n"
"shared by all elements. The size is determined at runtime.\n"
"\n"
"Used for: multi-component concentrations, phase saturations,\n"
"          spectral data, or any variable-length physical quantity\n"
"          where all components have the same boundary treatment."
    );
    pause();

    subheading( "Construction and access" );
    code(
"ArrayVariable  av;                        // empty\n"
"ArrayVariable  av( 4, 0., PLAIN );        // 4 elements = 0, flag=PLAIN\n"
"ArrayVariable  av( {1., 2., 3., 4.}, ANY ); // initialiser list\n"
"\n"
"// Size from PropertyDatabase:\n"
"ArrayVariable  av( \"element array\", model.Database() );\n"
"\n"
"double v = av[i];    // read element i (no bounds check)\n"
"av(i) = 5.0;         // write element i (with assert in debug)\n"
"av.Size();           // number of elements\n"
"av.Flag();           // the single flag\n"
"av.Flag() = DIRICH;  // change the flag\n"
"\n"
"// Range-based for (C++11):\n"
"for ( double v : av ) cout << v << \" \";"
    );
    pause();

    subheading( "Arithmetic" );
    code(
"ArrayVariable  a( 4, 1., PLAIN );\n"
"ArrayVariable  b( 4, 2., PLAIN );\n"
"\n"
"auto c = a + b;      // element-wise: {3,3,3,3}\n"
"auto d = a * 2.0;    // scale: {2,2,2,2}\n"
"a += b;              // in-place: a = {3,3,3,3}\n"
"auto e = a.Pow(2.0); // element-wise power: {9,9,9,9}\n"
"\n"
"double mn, mx;\n"
"a.MinMax( mn, mx );  // min and max values"
    );
    pause();

    subheading( "Live demonstration" );
    {
        ArrayVariable  av( 5, 0., PLAIN );
        for ( uint32_t i{0U}; i < av.Size(); ++i )
            av(i) = static_cast<double>(i+1) * 1.5;

        cout << "  av = { ";
        for ( uint32_t i{0U}; i < av.Size(); ++i )
            cout << dbl_str(av[i],1) << " ";
        cout << "}  flag=" << parseStatus(av.Flag()) << "\n";

        auto av2 = av * 2.0;
        cout << "  av*2 = { ";
        for ( uint32_t i{0U}; i < av2.Size(); ++i )
            cout << dbl_str(av2[i],1) << " ";
        cout << "}\n";

        double mn, mx;
        av.MinMax( mn, mx );
        output( "MinMax: [" + dbl_str(mn,1) + ", " + dbl_str(mx,1) + "]" );

        auto av3 = av.Pow( 2.0 );
        cout << "  av^2 = { ";
        for ( uint32_t i{0U}; i < av3.Size(); ++i )
            cout << dbl_str(av3[i],2) << " ";
        cout << "}\n";
    }
    pause();

    keypoint(
"ArrayVariable has ONE flag for ALL elements.\n"
"Use it when all components have the same boundary treatment.\n"
"If components need independent flags, like in chemical reaction calculations with GEMS (PSI), use FlaggedArrayVariable.\n"
"Size is set at construction and can be changed with Resize()."
    );
    pause();


    // ========================================================================
    // Step 6: FlaggedArrayVariable
    // ========================================================================
    heading( 6, 10, "FlaggedArrayVariable" );

    explain(
"FlaggedArrayVariable stores an arbitrary number of doubles,\n"
"each with its OWN independent VARIABLE_FLAG.\n"
"\n"
"Used for: multi-species reactive transport where some species\n"
"          concentrations are fixed (DIRICH) and others are free (PLAIN),\n"
"          or any variable-length quantity with mixed boundary conditions."
    );
    pause();

    subheading( "Construction and access" );
    code(
"FlaggedArrayVariable  fav( 4, 0., PLAIN );  // 4 elements, all PLAIN\n"
"\n"
"double v = fav[i];           // read element i\n"
"fav(i) = 5.0;                // write element i\n"
"VARIABLE_FLAG f = fav.Flag(i); // flag of element i\n"
"fav.Flag(i) = DIRICH;          // change flag of element i\n"
"fav.Flag(i, NEUMANN);          // alternative setter\n"
"\n"
"// Key difference from ArrayVariable:\n"
"// Each element has its OWN flag:\n"
"fav.Flag(0) = DIRICH;   // element 0 is fixed\n"
"fav.Flag(1) = PLAIN;    // element 1 is free\n"
"fav.Flag(2) = PLAIN;    // element 2 is free\n"
"fav.Flag(3) = NEUMANN;  // element 3 has flux condition"
    );
    pause();

    subheading( "Comparison with ArrayVariable" );
    code(
"// ArrayVariable: ONE flag for all\n"
"ArrayVariable  av( 4, 1., DIRICH );  // all DIRICH\n"
"av.Flag();                            // returns DIRICH for all\n"
"\n"
"// FlaggedArrayVariable: INDEPENDENT flags\n"
"FlaggedArrayVariable  fav( 4, 1., PLAIN );\n"
"fav.Flag(0) = DIRICH;   // only element 0 is fixed\n"
"fav.Flag(1) = PLAIN;    // element 1 is free\n"
"fav.Flag(0);            // returns DIRICH\n"
"fav.Flag(1);            // returns PLAIN"
    );
    pause();

    subheading( "Live demonstration" );
    {
        FlaggedArrayVariable  fav( 4, 0., PLAIN );
        for ( uint32_t i{0U}; i < fav.Size(); ++i )
            fav(i) = static_cast<double>(i) * 0.25;

        // Fix first and last elements as boundary conditions
        fav.Flag(0) = DIRICH;
        fav(0) = 1.0;
        fav.Flag(3) = DIRICH;
        fav(3) = 0.0;

        cout << "  FlaggedArrayVariable (4 species concentrations):\n";
        for ( uint32_t i{0U}; i < fav.Size(); ++i )
            output( "  species[" + to_string(i) + "] = "
                    + dbl_str(fav[i],2)
                    + "  flag=" + string(parseStatus(fav.Flag(i))) );

        // Arithmetic preserves flags
        auto fav2 = fav * 2.0;
        cout << "  After *2.0:\n";
        for ( uint32_t i{0U}; i < fav2.Size(); ++i )
            output( "  species[" + to_string(i) + "] = "
                    + dbl_str(fav2[i],2)
                    + "  flag=" + string(parseStatus(fav2.Flag(i))) );
    }
    pause();

    keypoint(
"FlaggedArrayVariable gives each element its own flag.\n"
"Use it for multi-species systems where some species are fixed\n"
"(e.g. pH-buffered) and others are free to react.\n"
"Arithmetic operations preserve the flags of the left operand.\n"
"ARRAY type have runtime sizes, and vector norms can be computed using CSMP_mathUtilities"
    );
    pause();


    // ========================================================================
    // Step 7: Variable placements
    // ========================================================================
    heading( 7, 10, "Variable placements" );

    explain(
"Every CSMP variable is associated with a PLACEMENT that determines\n"
"where in the mesh it is stored:\n"
"\n"
"  NODE                      -- at mesh nodes (vertices)\n"
"  ELEMENT                   -- at element barycentres\n"
"  ELEMENT_INTEGRATION_POINT -- at Gauss quadrature points\n"
"  FACE                      -- at face barycentres\n"
"  INTER_FACE                -- at internal interface barycentres\n"
"  ...\n"
"\n"
"The placement determines:\n"
"  - how many storage slots are allocated\n"
"  - which interpolation/extrapolation operations are valid\n"
"  - how the PDE assembler accesses the variable"
    );
    pause();

    subheading( "Placement in PropertyHandle" );
    code(
"// NODE variable -- one value per mesh node:\n"
"PropertyHandle<dim>  pressure( model, \"fluid pressure\",\n"
"                               SCALAR, NODE );\n"
"\n"
"// ELEMENT variable -- one value per element:\n"
"PropertyHandle<dim>  permeability( model, \"permeability\",\n"
"                                   SCALAR, ELEMENT );\n"
"\n"
"// ELEMENT_INTEGRATION_POINT -- one value per Gauss point:\n"
"PropertyHandle<dim>  ip_perm( model, \"ip permeability\",\n"
"                              SCALAR, ELEMENT_INTEGRATION_POINT );"
    );
    pause();

    subheading( "Why placement matters" );
    explain(
"NODE variables are used for:\n"
"  - primary unknowns (pressure, temperature, displacement)\n"
"  - boundary conditions\n"
"  - to obtain continuous fields during post-processing for visualisation\n"
"\n"
"ELEMENT variables are used for:\n"
"  - material properties (permeability, Young's modulus) on a mesh conforming with material domains\n"
"  - element-averaged results (stress, strain)\n"
"\n"
"ELEMENT_INTEGRATION_POINT variables are used for:\n"
"  - exact finite-element integration\n"
"  - fine-grained evaluation of failure\n"
"  - spatially varying material properties in QFEM\n"
"  - history variables in nonlinear problems\n"
"  - high-accuracy flux computations"
    );
    pause();

    subheading( "Live demonstration -- create variables at different placements" );
    {
        PropertyHandle<dim>  node_var(    model, "demo node scalar",
                                          SCALAR, NODE );
        PropertyHandle<dim>  elem_var(    model, "demo element scalar",
                                          SCALAR, ELEMENT );
        PropertyHandle<dim>  node_vec(    model, "demo node vector",
                                          VECTOR, NODE );
        PropertyHandle<dim>  elem_tensor( model, "demo element tensor",
                                          TENSOR, ELEMENT );

        node_var    = 1.0;
        elem_var    = 2.0;
        node_vec    = makeVector( PLAIN, PLAIN, PLAIN, 1., 2., 3. );
        elem_tensor = makeTensor( PLAIN, PLAIN, PLAIN,
                                  1., 0., 0.,
                                  0., 1., 0.,
                                  0., 0., 1. );

        output( "NODE scalar:    " + to_string(model.Mesh().Nodes())
                + " values stored" );
        output( "ELEMENT scalar: " + to_string(model.Mesh().Elements())
                + " values stored" );
        output( "NODE vector:    " + to_string(model.Mesh().Nodes())
                + " x " + to_string(dim) + " = "
                + to_string(model.Mesh().Nodes()*dim) + " values stored" );
        output( "ELEMENT tensor: " + to_string(model.Mesh().Elements())
                + " x " + to_string(dim*dim) + " = "
                + to_string(model.Mesh().Elements()*dim*dim) + " values stored" );
    }
    pause();

    keypoint(
"Choose the placement that matches the physics:\n"
"  NODE for FEM/FVM-primary variables / unknowns and boundary conditions.\n"
"  ELEMENT for material properties and averaged results.\n"
"  ELEMENT_INTEGRATION_POINT for spatially varying QFEM properties.\n"
"Wrong placement = wrong number of storage slots = crash or wrong results."
    );
    pause();


    // ========================================================================
    // Step 8: Interpolation
    // ========================================================================
    heading( 8, 10, "Interpolation" );

    explain(
"Interpolation moves a variable from a finer placement to a coarser one\n"
"using the finite-element shape functions.\n"
"\n"
"  Node --> Element:           weighted average of node values\n"
"  Node --> Integration point: FEM shape function evaluation at Gauss point\n"
"  Integration point --> Element: weighted average of IP values\n"
"\n"
"These operations are available on both Model and ModelSubDomain (Region)."
    );
    pause();

    subheading( "Interpolation methods" );
    code(
"// Node --> Element (arithmetic average):\n"
"model.InterpolateNodeToCellProperty(\n"
"    \"fluid pressure\",          // source: NODE\n"
"    \"element pressure\" );      // target: ELEMENT\n"
"\n"
"// Node --> Integration point (FEM shape functions):\n"
"model.InterpolateNodeToIntegrationPointProperty(\n"
"    \"fluid pressure\",          // source: NODE\n"
"    \"ip pressure\" );           // target: ELEMENT_INTEGRATION_POINT\n"
"\n"
"// Integration point --> Element (weighted average):\n"
"model.InterpolateIntegrationPointToCellProperty(\n"
"    \"ip pressure\",             // source: ELEMENT_INTEGRATION_POINT\n"
"    \"element pressure\" );      // target: ELEMENT\n"
"\n"
"// Same operations on a Region:\n"
"gref.InterpolateNodeToCellProperty( \"fluid pressure\",\n"
"                                    \"element pressure\" );"
    );
    pause();

    subheading( "Live demonstration -- node to element interpolation" );
    {
        // Set a linear pressure field: p = x-coordinate
        PropertyHandle<dim>  pressure( model, "fluid pressure", SCALAR, NODE );
        model.AssignNodeCoordinatesTo( "fluid pressure", 'x' );
        output( "Set 'fluid pressure' = x-coordinate at nodes." );
        printRangeOfVariable( model, "Model", "fluid pressure" );

        // Interpolate to elements
        PropertyHandle<dim>  elem_p( model, "element pressure", SCALAR, ELEMENT );
        gref.InterpolateNodeToCellProperty( "fluid pressure", "element pressure" );
        output( "Interpolated to elements:" );
        printRangeOfVariable( model, "Model", "element pressure" );

        // The ranges should be similar (linear field interpolates exactly)
        output( "For a linear field, node and element ranges should match." );
    }
    pause();

    keypoint(
"Interpolation uses FEM shape functions -- it is exact for fields\n"
"that lie in the FEM approximation space (e.g. linear fields on\n"
"linear elements). For nonlinear fields it introduces approximation error.\n"
"Always interpolate from the finer to the coarser placement."
    );
    pause();


    // ========================================================================
    // Step 9: Extrapolation
    // ========================================================================
    heading( 9, 10, "Extrapolation" );

    explain(
"Extrapolation moves a variable from a coarser placement to a finer one.\n"
"This is the reverse of interpolation and involves averaging.\n"
"\n"
"  Element --> Node:           distance- or volume-weighted average\n"
"  Integration point --> Node: extrapolation + averaging\n"
"  Element --> Integration point: copy element value to all its IPs\n"
"\n"
"Extrapolation is used in post-processing to visualise element results\n"
"(e.g. stress, strain) at nodes for smooth VTK output."
    );
    pause();

    subheading( "Extrapolation methods" );
    code(
"// Element --> Node (distance-weighted average, default):\n"
"model.ExtrapolateCellToNodeProperty(\n"
"    \"element pressure\",        // source: ELEMENT\n"
"    \"nodal pressure\" );        // target: NODE\n"
"\n"
"// Element --> Node (volume-weighted average):\n"
"model.ExtrapolateCellToNodeProperty(\n"
"    \"element pressure\",\n"
"    \"nodal pressure\",\n"
"    false );                    // false = volume-weighted\n"
"\n"
"// Integration point --> Node:\n"
"model.ExtrapolateIntegrationPointToNodeProperty(\n"
"    \"ip pressure\",             // source: ELEMENT_INTEGRATION_POINT\n"
"    \"nodal ip pressure\" );     // target: NODE\n"
"\n"
"// Element --> Integration point (copy):\n"
"gref.ExtrapolateCellToIntegrationPointProperty(\n"
"    \"element pressure\",        // source: ELEMENT\n"
"    \"ip pressure\" );           // target: ELEMENT_INTEGRATION_POINT\n"
"\n"
"// Same operations on a Region:\n"
"gref.ExtrapolateCellToNodeProperty( \"element pressure\",\n"
"                                    \"nodal pressure\" );"
    );
    pause();

    subheading( "Live demonstration -- element to node extrapolation" );
    {
        // Use the element pressure from step 8
        PropertyHandle<dim>  elem_p( model, "element pressure", SCALAR, ELEMENT );
        gref.InterpolateNodeToCellProperty( "fluid pressure", "element pressure" );
        PropertyHandle<dim>  nodal_p( model, "nodal pressure", SCALAR, NODE );

        // Distance-weighted extrapolation
        gref.ExtrapolateCellToNodeProperty( "element pressure",
                                            "nodal pressure",
                                            true );  // true = distance-weighted
        output( "Extrapolated 'element pressure' to nodes (distance-weighted):" );
        printRangeOfVariable( model, "Model", "nodal pressure" );

        // Compare with original node values
        output( "Original 'fluid pressure' at nodes:" );
        printRangeOfVariable( model, "Model", "fluid pressure" );

        output( "The ranges should be similar -- small differences arise\n"
                "  >>  from the averaging in the extrapolation step." );
    }
    pause();

    subheading( "Interpolation vs Extrapolation -- summary" );
    code(
"// INTERPOLATION (fine --> coarse, uses FEM shape functions):\n"
"Node  --> Element:            InterpolateNodeToCellProperty\n"
"Node  --> IntegrationPoint:   InterpolateNodeToIntegrationPointProperty\n"
"IP    --> Element:            InterpolateIntegrationPointToCellProperty\n"
"\n"
"// EXTRAPOLATION (coarse --> fine, uses averaging):\n"
"Element --> Node:             ExtrapolateCellToNodeProperty\n"
"Element --> IntegrationPoint: ExtrapolateCellToIntegrationPointProperty\n"
"IP      --> Node:             ExtrapolateIntegrationPointToNodeProperty"
    );
    pause();

    keypoint(
"Extrapolation introduces smoothing -- element discontinuities\n"
"are averaged at shared nodes. This can be intentional for visualisation\n"
"but may hide physical discontinuities at material interfaces.\n"
"For accurate flux computation, keep variables at their natural placement."
"Note: in CSMP, C++ class member function names start with capitals, and "
"global functions with small letters."
    );
    pause();


    // ========================================================================
    // Step 10: Assignment to DenseMatrix
    // ========================================================================
    heading( 10, 10, "Assignment to DenseMatrix" );

    explain(
"DenseMatrix is CSMP's local compile-time-sized matrix type used inside finite-element\n"
"integration operators (MathOperatorLHS, MathOperatorRHS).\n"
"\n"
"When you write a custom PDE operator, you compute element integrals\n"
"and store them in a DenseMatrix before they are assembled into the\n"
"global sparse solution matrix by AssignToGlobal().\n"
"\n"
"CSMP variables are assigned to DenseMatrix entries using:\n"
"  - direct element access DM(i,j)\n"
"  - ScalarVariable values via sc()\n"
"  - VectorVariable components via vc(i)\n"
"  - TensorVariable entries via ts(i,j)"
    );
    pause();

    subheading( "DenseMatrix interface" );
    code(
"// Construction:\n"
"DenseMatrix  LHS( n_rows, n_cols );   // n x m matrix, initialised to 0\n"
"DenseMatrix  RHS( n_rows, 1 );        // column vector\n"
"\n"
"// Element access:\n"
"LHS(i,j) = value;                     // write\n"
"double v = LHS(i,j);                  // read\n"
"\n"
"// Size queries:\n"
"LHS.Rows();   // number of rows\n"
"LHS.Cols();   // number of columns\n"
"\n"
"// Zeroing:\n"
"LHS.Zero();                           // set all entries to 0\n"
"\n"
"// Matrix-vector product:\n"
"DenseMatrix  result = LHS * RHS;"
    );
    pause();

    subheading( "Assigning CSMP variables to DenseMatrix" );
    code(
"// Scalar --> single matrix entry:\n"
"ScalarVariable  sc( PLAIN, 3.14 );\n"
"LHS(0,0) = sc();                      // extract double with ()\n"
"\n"
"// Vector --> column of matrix:\n"
"VectorVariable<3U>  vc( PLAIN, 1., 2., 3. );\n"
"for ( uint32_t i{0U}; i < dim; ++i )\n"
"    RHS(i,0) = vc(i);                 // fill column vector\n"
"\n"
"// Tensor --> submatrix block:\n"
"TensorVariable<3U>  ts( PLAIN, 1.0 ); // identity\n"
"for ( uint32_t i{0U}; i < dim; ++i )\n"
"    for ( uint32_t j{0U}; j < dim; ++j )\n"
"        LHS(i,j) = ts(i,j);           // fill 3x3 block\n"
"\n"
"// Array --> column vector:\n"
"ArrayVariable  av( dim, 0., PLAIN );\n"
"for ( uint32_t i{0U}; i < av.Size(); ++i )\n"
"    RHS(i,0) = av[i];"
    );
    pause();

    subheading( "Live demonstration -- build a local stiffness matrix" );
    {
        // Build a simple 3x3 local stiffness matrix from a tensor
        // (this is what NumIntegral_dNT_dN_dV does internally)
        DenseMatrix<dim>  K_local( dim, dim );
        K_local.Zero();

        // Anisotropic permeability tensor
        TensorVariable<dim>  K_perm( PLAIN, PLAIN, PLAIN,
                                     1e-12, 0.,    0.,
                                     0.,    5e-13, 0.,
                                     0.,    0.,    2e-13 );

        // Fill local matrix from tensor
        for ( uint32_t i{0U}; i < dim; ++i )
            for ( uint32_t j{0U}; j < dim; ++j )
                K_local(i,j) = K_perm(i,j);

        output( "Local stiffness matrix K_local (from permeability tensor):" );
        for ( uint32_t i{0U}; i < dim; ++i ) {
            ostringstream row;
            row << "  row " << i << ": ";
            for ( uint32_t j{0U}; j < dim; ++j )
                row << scientific << setprecision(2) << setw(12) << K_local(i,j) << " ";
            output( row.str() );
        }

        // Build a local RHS vector from a gradient vector
        DenseMatrix<dim>  f_local( dim, 1 );
        VectorVariable<dim>  grad_p( PLAIN, PLAIN, PLAIN, 1e4, 0., 0. );
        for ( uint32_t i{0U}; i < dim; ++i )
            f_local(i,0) = grad_p(i);

        output( "Local RHS vector f_local (from pressure gradient):" );
        for ( uint32_t i{0U}; i < dim; ++i )
            output( "  f[" + to_string(i) + "] = "
                    + dbl_str(f_local(i,0), 2) );

        // Matrix-vector product: velocity = K * grad_p
        DenseMatrix  velocity_dm = K_local * f_local;
        output( "K_local * f_local (Darcy velocity):" );
        for ( uint32_t i{0U}; i < dim; ++i )
            output( "  v[" + to_string(i) + "] = "
                    + dbl_str(velocity_dm(i,0), 6) );

        // Verify against direct tensor-vector product
        auto velocity_vc = K_perm * grad_p;
        output( "Direct K_perm * grad_p (should match):" );
        for ( uint32_t i{0U}; i < dim; ++i )
            output( "  v[" + to_string(i) + "] = "
                    + dbl_str(velocity_vc(i), 6) );
    }
    pause();

    subheading( "Where DenseMatrix appears in CSMP operator code" );
    code(
"// Inside a custom MathOperatorLHS:\n"
"void MyOperator::ComputeContribution( const Element<dim>& e )\n"
"{\n"
"    // LHS is the local element matrix (member of MathOperatorLHS)\n"
"    LHS.Zero();\n"
"\n"
"    // Read material property:\n"
"    TensorVariable<dim>  K;\n"
"    e.Read( permeability_key_, K );\n"
"\n"
"    // Compute element contribution (simplified):\n"
"    for ( uint32_t i{0U}; i < e.Nodes(); ++i )\n"
"        for ( uint32_t j{0U}; j < e.Nodes(); ++j )\n"
"            for ( uint32_t k{0U}; k < dim; ++k )\n"
"                LHS(i,j) += dN(i,k) * K(k,k) * dN(j,k) * detJ;\n"
"    // LHS is then assembled into the global matrix by AssignToGlobal()\n"
"}"
    );
    pause();

    keypoint(
"DenseMatrix is the bridge between CSMP variables and the global\n"
"sparse matrix. In a custom operator:\n"
"  1. Read material properties into ScalarVariable / TensorVariable.\n"
"  2. Compute the element contribution into DenseMatrix LHS.\n"
"  3. AssignToGlobal() scatters LHS into the global CompressedRowMatrix.\n"
"You never access the global matrix directly -- always via DenseMatrix."
    );
    pause();


    // ========================================================================
    // Closing screen
    // ========================================================================
    blank();
    divider();
    cout << "\n"
         << "  You have completed Example 3: Physical Variables in CSMP\n\n"
         << "  Steps covered:\n\n"
         << "    1.  Why variables have flags\n"
         << "    2.  ScalarVariable -- value + flag\n"
         << "    3.  VectorVariable -- dim components, each with a flag\n"
         << "    4.  TensorVariable -- dim x dim matrix, diagonal flags\n"
         << "    5.  ArrayVariable  -- runtime-sized, single flag\n"
         << "    6.  FlaggedArrayVariable -- runtime-sized, per-element flags\n"
         << "    7.  Variable placements -- NODE, ELEMENT, IP\n"
         << "    8.  Interpolation -- fine to coarse using FEM shape functions\n"
         << "    9.  Extrapolation -- coarse to fine using averaging\n"
         << "   10.  Assignment to DenseMatrix\n\n"
         << "  Next: Example 4 -- Your First PDE Solve\n\n";
    divider();
    pause( "  [ Press ENTER to exit ]" );

    // Restore working directory
    error_code ec;
    filesystem::current_path( "../../example_inputs/", ec );
    if ( ec )
        cerr << "\nCSMP_PhysicalVariables_Example: could not restore "
             << "working directory: " << ec.message() << "\n";

} // Run()

} // csmp

