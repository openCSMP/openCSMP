// NumIntegral_dNT_op_dN_NT_v_dN_dV_Test.cpp
#include "NumIntegral_dNT_op_dN_NT_v_dN_dV_Test.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::NumIntegral_dNT_op_dN_NT_v_dN_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nNumIntegral_dNT_op_dN_NT_v_dN_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "velocity",       makeVector( PLAIN, PLAIN, 0.0, 0.0 ) );
 }


NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::~NumIntegral_dNT_op_dN_NT_v_dN_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::run()
 {
    dispersionOnlyTest();
    advectionOnlyTest();
    combinedTest();
    nonSymmetryTest();
    zeroVelocityTest();
 }


// =============================================================================
// dispersionOnlyTest — v=0, result must equal pure stiffness matrix
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::dispersionOnlyTest()
 {
    NumIntegral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
                                                     "diffusivity",
                                                     "velocity",
                                                     "fluid pressure",
                                                     "fluid pressure" );

    setElementScalar( { 1.0, 1.0, 1.0, 1.0 }, "diffusivity" );
    setElementVector( { 0.0, 0.0 },            "velocity"    );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, adv_disp );

    if ( verbose_ )
      {
        cout << "\ndispersionOnlyTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // With v=0 and uniform D=1 the result equals the pure stiffness matrix
    // K_jk = (g_j . g_k) / (4A), s = 125000/250000 = 0.5
    //
    // Diagonal:
    //   A_00 = A_22 = 2*s = 1.0  (2 elements, each contributing s)
    //   A_11 = A_33 = 2*s = 1.0
    //   A_44 = 4*s2 = 4.0        (4 elements, each contributing s2=1.0)
    //
    // Wait — with uniform D=1 across all 4 elements:
    //   K_00 = s+s = 2s = 1.0  (e0 and e3 each contribute s=0.5)
    //   K_11 = s+s = 2s = 1.0  (e0 and e1)
    //   K_22 = s+s = 2s = 1.0  (e1 and e2)
    //   K_33 = s+s = 2s = 1.0  (e2 and e3)
    //   K_44 = s2+s2+s2+s2 = 4.0 (all 4 elements, s2=1.0)
    //   K_04 = K_40 = -(s+s) = -1.0
    //   K_14 = K_41 = -(s+s) = -1.0
    //   K_24 = K_42 = -(s+s) = -1.0
    //   K_34 = K_43 = -(s+s) = -1.0
    //   All boundary-boundary pairs = 0
    // -----------------------------------------------------------------------

    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    // diagonal
    _equal( sm.At(0,0),  2.0*s, t );  // = 1.0
    _equal( sm.At(1,1),  2.0*s, t );  // = 1.0
    _equal( sm.At(2,2),  2.0*s, t );  // = 1.0
    _equal( sm.At(3,3),  2.0*s, t );  // = 1.0
    _equal( sm.At(4,4),  4.0,   t );  // = 4.0

    // off-diagonal
    _equal( sm.At(0,4), -2.0*s, t );  // = -1.0
    _equal( sm.At(1,4), -2.0*s, t );
    _equal( sm.At(2,4), -2.0*s, t );
    _equal( sm.At(3,4), -2.0*s, t );

    _equal( sm.At(4,0), -2.0*s, t );  // symmetric
    _equal( sm.At(4,1), -2.0*s, t );
    _equal( sm.At(4,2), -2.0*s, t );
    _equal( sm.At(4,3), -2.0*s, t );
 }




// =============================================================================
// advectionOnlyTest — D=0, result must equal pure advection matrix
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::advectionOnlyTest()
 {
    NumIntegral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
                                                     "diffusivity",
                                                     "velocity",
                                                     "fluid pressure",
                                                     "fluid pressure" );

    setElementScalar( { 0.0, 0.0, 0.0, 0.0 }, "diffusivity" );
    setElementVector( { 1.0, 0.0 },            "velocity"    );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, adv_disp );

    if ( verbose_ )
      {
        cout << "\nadvectionOnlyTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // With D=0 and v=(1,0) uniform the result equals the pure advection matrix
    // V_jk = (v . g_k) / 6 = g_k^x / 6
    // q = 250/6
    //
    // Global assembly (from derivation above):
    //   V_00 = -2q,  V_01 = -q,   V_03 = q,   V_04 = 2q
    //   V_10 = -q,   V_11 = -2q,  V_12 = q,   V_14 = 2q
    //   V_21 = -q,   V_22 = 2q,   V_23 = q,   V_24 = -2q
    //   V_30 = -q,   V_32 = q,    V_33 = 2q,  V_34 = -2q
    //   V_40 = -2q,  V_41 = -2q,  V_42 = 2q,  V_43 = 2q,  V_44 = 0
    // -----------------------------------------------------------------------

    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    // row 0
    _equal( sm.At(0,0), -2.0*q, t );
    _equal( sm.At(0,1), -1.0*q, t );
    _equal( sm.At(0,3),  1.0*q, t );
    _equal( sm.At(0,4),  2.0*q, t );

    // row 1
    _equal( sm.At(1,0), -1.0*q, t );
    _equal( sm.At(1,1), -2.0*q, t );
    _equal( sm.At(1,2),  1.0*q, t );
    _equal( sm.At(1,4),  2.0*q, t );

    // row 2
    _equal( sm.At(2,1), -1.0*q, t );
    _equal( sm.At(2,2),  2.0*q, t );
    _equal( sm.At(2,3),  1.0*q, t );
    _equal( sm.At(2,4), -2.0*q, t );

    // row 3
    _equal( sm.At(3,0), -1.0*q, t );
    _equal( sm.At(3,2),  1.0*q, t );
    _equal( sm.At(3,3),  2.0*q, t );
    _equal( sm.At(3,4), -2.0*q, t );

    // row 4
    _equal( sm.At(4,0), -2.0*q, t );
    _equal( sm.At(4,1), -2.0*q, t );
    _equal( sm.At(4,2),  2.0*q, t );
    _equal( sm.At(4,3),  2.0*q, t );

    // V_44 = 0 — entry may not exist in sparse structure
    if ( sm.HasEntry(4,4) )
      _equal( sm.At(4,4), 0.0, t );
 }
 
 
 
 

// =============================================================================
// combinedTest — D=1, v=(1,0), full advection-dispersion matrix
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::combinedTest()
 {
    NumIntegral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
                                                     "diffusivity",
                                                     "velocity",
                                                     "fluid pressure",
                                                     "fluid pressure" );

    setElementScalar( { 1.0, 1.0, 1.0, 1.0 }, "diffusivity" );
    setElementVector( { 1.0, 0.0 },            "velocity"    );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, adv_disp );

    if ( verbose_ )
      {
        cout << "\ncombinedTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // A_jk = K_jk + V_jk
    // With uniform D=1: K values use s=0.5, each boundary node in 2 elements
    //   K_00 = K_11 = K_22 = K_33 = 2s = 1.0
    //   K_44 = 4.0
    //   K_j4 = K_4j = -2s = -1.0
    //
    // Advection with v=(1,0): q = 250/6
    //   V values as derived above
    //
    // Combined:
    //   A_00 = 2s - 2q = 1.0 - 500/6
    //   A_11 = 2s - 2q = 1.0 - 500/6
    //   A_22 = 2s + 2q = 1.0 + 500/6
    //   A_33 = 2s + 2q = 1.0 + 500/6
    //   A_44 = 4.0 + 0 = 4.0
    //   A_04 = -2s + 2q = -1.0 + 500/6
    //   A_14 = -2s + 2q = -1.0 + 500/6
    //   A_24 = -2s - 2q = -1.0 - 500/6
    //   A_34 = -2s - 2q = -1.0 - 500/6
    //   A_40 = -2s - 2q = -1.0 - 500/6  (note: A_04 != A_40)
    //   A_41 = -2s - 2q = -1.0 - 500/6
    //   A_42 = -2s + 2q = -1.0 + 500/6
    //   A_43 = -2s + 2q = -1.0 + 500/6
    //   A_01 = 0 - q = -q
    //   A_03 = 0 + q =  q
    //   A_10 = 0 - q = -q
    //   A_12 = 0 + q =  q
    //   A_21 = 0 - q = -q
    //   A_23 = 0 + q =  q
    //   A_30 = 0 - q = -q
    //   A_32 = 0 + q =  q
    // -----------------------------------------------------------------------

    constexpr double s{ 0.5 };
    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    // diagonal
    _equal( sm.At(0,0), 2.0*s - 2.0*q, t );
    _equal( sm.At(1,1), 2.0*s - 2.0*q, t );
    _equal( sm.At(2,2), 2.0*s + 2.0*q, t );
    _equal( sm.At(3,3), 2.0*s + 2.0*q, t );
    _equal( sm.At(4,4), 4.0,            t );

    // row 0
    _equal( sm.At(0,1), -q,              t );
    _equal( sm.At(0,3),  q,              t );
    _equal( sm.At(0,4), -2.0*s + 2.0*q, t );

    // row 1
    _equal( sm.At(1,0), -q,              t );
    _equal( sm.At(1,2),  q,              t );
    _equal( sm.At(1,4), -2.0*s + 2.0*q, t );

    // row 2
    _equal( sm.At(2,1), -q,              t );
    _equal( sm.At(2,3),  q,              t );
    _equal( sm.At(2,4), -2.0*s - 2.0*q, t );

    // row 3
    _equal( sm.At(3,0), -q,              t );
    _equal( sm.At(3,2),  q,              t );
    _equal( sm.At(3,4), -2.0*s - 2.0*q, t );

    // row 4 — note asymmetry: A_40 != A_04
    _equal( sm.At(4,0), -2.0*s - 2.0*q, t );
    _equal( sm.At(4,1), -2.0*s - 2.0*q, t );
    _equal( sm.At(4,2), -2.0*s + 2.0*q, t );
    _equal( sm.At(4,3), -2.0*s + 2.0*q, t );
 }





// =============================================================================
// nonSymmetryTest — A_jk != A_kj for non-zero v
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::nonSymmetryTest()
 {
    NumIntegral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
                                                     "diffusivity",
                                                     "velocity",
                                                     "fluid pressure",
                                                     "fluid pressure" );

    setElementScalar( { 1.0, 1.0, 1.0, 1.0 }, "diffusivity" );
    setElementVector( { 1.0, 0.0 },            "velocity"    );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, adv_disp );

    if ( verbose_ )
      cout << "\nnonSymmetryTest: checking A_04 != A_40\n";

    // A_04 = -2.0 + 500/6 = -2.0 + 83.33 = 81.33... / ... wait
    // let me recompute in actual values:
    // A_04 = K_04 + V_04 = -1.0 + 500/6 = -1.0 + 83.333 = 82.333... no
    //
    // With uniform D=1 across all 4 elements:
    //   K_04 = -(s+s) = -1.0  (e0 and e3 each contribute -s=-0.5)
    //   V_04 = 2q = 500/6
    //   A_04 = -1.0 + 500/6
    //
    //   K_40 = K_04 = -1.0  (stiffness is symmetric)
    //   V_40 = -2q = -500/6  (advection is NOT symmetric)
    //   A_40 = -1.0 - 500/6
    //
    // So A_04 != A_40 whenever q != 0

    constexpr double s{ 0.5 };
    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    // A_04 = K_04 + V_04 = -2s + 2q = -1.0 + 500/6
    // A_40 = K_40 + V_40 = -2s - 2q = -1.0 - 500/6
    _equal( sm.At(0,4), -2.0*s + 2.0*q, t );
    _equal( sm.At(4,0), -2.0*s - 2.0*q, t );

    // difference = 4q
    const double diff = sm.At(0,4) - sm.At(4,0);
    _equal( diff, 4.0*q, t );
 }
 
 

// =============================================================================
// zeroVelocityTest — v=0 -> result must be symmetric
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::zeroVelocityTest()
 {
    NumIntegral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
                                                     "diffusivity",
                                                     "velocity",
                                                     "fluid pressure",
                                                     "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );
    setElementVector( { 0.0, 0.0 },            "velocity"    );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, adv_disp );

    if ( verbose_ )
      {
        cout << "\nzeroVelocityTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // with v=0 the result equals the pure stiffness matrix which is symmetric
    // check all known connected pairs
    const std::vector<std::pair<size_t,size_t>> pairs{
        {0U,4U}, {1U,4U}, {2U,4U}, {3U,4U}
    };

    for ( const auto& [i,j] : pairs )
      _equal( sm.At(i,j), sm.At(j,i), tol_ );

    // also verify against known stiffness values from NumIntegral_dNT_lhsop_dN_dV_Test
    // with non-uniform D={1,2,3,4}
    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    _equal( sm.At(0,0),  5.0*s, t );
    _equal( sm.At(1,1),  3.0*s, t );
    _equal( sm.At(2,2),  5.0*s, t );
    _equal( sm.At(3,3),  7.0*s, t );
    _equal( sm.At(4,4), 20.0*s, t );

    _equal( sm.At(0,4), -5.0*s, t );
    _equal( sm.At(1,4), -3.0*s, t );
    _equal( sm.At(2,4), -5.0*s, t );
    _equal( sm.At(3,4), -7.0*s, t );
 }


// =============================================================================
// Helper methods
// =============================================================================

void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::setElementScalar(
    const vector<double>& val, const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR,
          "NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::setElementScalar",
          var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }


void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::setElementVector(
    const array<double,2>& v, const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        VectorVariable<2U> vc;
        vc(0) = v[0];
        vc(1) = v[1];
        (*it)->Store( key, vc );
      }
 }


void NumIntegral_dNT_op_dN_NT_v_dN_dV_Test::calculateGlobalMatrix(
    SparseMatrix& sm, MathOperatorLHS<2U>& oper )
 {
    Region<2>& domain( sg_->Region("Model") );
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        oper.GetOperands(        *(*it) );
        oper.ComputeContribution(*(*it) );
        oper.AssignToGlobal(     *(*it), sm );
      }
 }

} // namespace csmp
