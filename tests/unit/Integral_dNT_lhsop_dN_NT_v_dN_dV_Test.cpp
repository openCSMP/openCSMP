// Integral_dNT_lhsop_dN_NT_v_dN_dV_Test.cpp
#include "Integral_dNT_lhsop_dN_NT_v_dN_dV_Test.h"

#include "Integral_dNT_op_dN_NT_v_dN_dV.h"
#include "NumIntegral_DNT_op_DN_NT_v_DN_dV.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::Integral_dNT_lhsop_dN_NT_v_dN_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nIntegral_dNT_lhsop_dN_NT_v_dN_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "velocity",       makeVector( PLAIN, PLAIN, 0.0, 0.0 ) );
 }


Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::~Integral_dNT_lhsop_dN_NT_v_dN_dV_Test()
 { delete sg_; }


// =============================================================================
// Test runner
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::run()
 {
    dispersionOnlyTest();
    advectionOnlyTest();
    combinedTest();
    nonSymmetryTest();
    zeroVelocityTest();
    matchesNumericalTest();
 }


// =============================================================================
// dispersionOnlyTest — v=0, result must equal pure stiffness matrix
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::dispersionOnlyTest()
 {
    Integral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
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
        cout << "\nIntegral_DNT_op_DN_NT_v_DN_dV dispersionOnlyTest:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // With v=0 and uniform D=1: pure stiffness matrix
    // s=0.5, each boundary node in 2 elements
    //   K_jj = 2s = 1.0  for boundary nodes
    //   K_44 = 4.0
    //   K_j4 = K_4j = -2s = -1.0
    // -----------------------------------------------------------------------
    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    _equal( sm.At(0,0),  2.0*s, t );
    _equal( sm.At(1,1),  2.0*s, t );
    _equal( sm.At(2,2),  2.0*s, t );
    _equal( sm.At(3,3),  2.0*s, t );
    _equal( sm.At(4,4),  4.0,   t );

    _equal( sm.At(0,4), -2.0*s, t );
    _equal( sm.At(1,4), -2.0*s, t );
    _equal( sm.At(2,4), -2.0*s, t );
    _equal( sm.At(3,4), -2.0*s, t );

    _equal( sm.At(4,0), -2.0*s, t );
    _equal( sm.At(4,1), -2.0*s, t );
    _equal( sm.At(4,2), -2.0*s, t );
    _equal( sm.At(4,3), -2.0*s, t );
 }


// =============================================================================
// advectionOnlyTest — D=0, result must equal pure advection matrix
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::advectionOnlyTest()
 {
    Integral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
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
        cout << "\nIntegral_DNT_op_DN_NT_v_DN_dV advectionOnlyTest:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // With D=0 and v=(1,0): pure advection matrix
    // q = 250/6
    // -----------------------------------------------------------------------
    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    _equal( sm.At(0,0), -2.0*q, t );
    _equal( sm.At(0,1), -1.0*q, t );
    _equal( sm.At(0,3),  1.0*q, t );
    _equal( sm.At(0,4),  2.0*q, t );

    _equal( sm.At(1,0), -1.0*q, t );
    _equal( sm.At(1,1), -2.0*q, t );
    _equal( sm.At(1,2),  1.0*q, t );
    _equal( sm.At(1,4),  2.0*q, t );

    _equal( sm.At(2,1), -1.0*q, t );
    _equal( sm.At(2,2),  2.0*q, t );
    _equal( sm.At(2,3),  1.0*q, t );
    _equal( sm.At(2,4), -2.0*q, t );

    _equal( sm.At(3,0), -1.0*q, t );
    _equal( sm.At(3,2),  1.0*q, t );
    _equal( sm.At(3,3),  2.0*q, t );
    _equal( sm.At(3,4), -2.0*q, t );

    _equal( sm.At(4,0), -2.0*q, t );
    _equal( sm.At(4,1), -2.0*q, t );
    _equal( sm.At(4,2),  2.0*q, t );
    _equal( sm.At(4,3),  2.0*q, t );

    if ( sm.HasEntry(4,4) )
      _equal( sm.At(4,4), 0.0, t );
 }


// =============================================================================
// combinedTest — D=1, v=(1,0), full advection-dispersion matrix
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::combinedTest()
 {
    Integral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
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
        cout << "\nIntegral_DNT_op_DN_NT_v_DN_dV combinedTest:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // A_jk = K_jk + V_jk, uniform D=1, v=(1,0)
    // s=0.5, q=250/6
    // -----------------------------------------------------------------------
    constexpr double s{ 0.5 };
    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    _equal( sm.At(0,0), 2.0*s - 2.0*q, t );
    _equal( sm.At(1,1), 2.0*s - 2.0*q, t );
    _equal( sm.At(2,2), 2.0*s + 2.0*q, t );
    _equal( sm.At(3,3), 2.0*s + 2.0*q, t );
    _equal( sm.At(4,4), 4.0,            t );

    _equal( sm.At(0,1), -q,              t );
    _equal( sm.At(0,3),  q,              t );
    _equal( sm.At(0,4), -2.0*s + 2.0*q, t );

    _equal( sm.At(1,0), -q,              t );
    _equal( sm.At(1,2),  q,              t );
    _equal( sm.At(1,4), -2.0*s + 2.0*q, t );

    _equal( sm.At(2,1), -q,              t );
    _equal( sm.At(2,3),  q,              t );
    _equal( sm.At(2,4), -2.0*s - 2.0*q, t );

    _equal( sm.At(3,0), -q,              t );
    _equal( sm.At(3,2),  q,              t );
    _equal( sm.At(3,4), -2.0*s - 2.0*q, t );

    _equal( sm.At(4,0), -2.0*s - 2.0*q, t );
    _equal( sm.At(4,1), -2.0*s - 2.0*q, t );
    _equal( sm.At(4,2), -2.0*s + 2.0*q, t );
    _equal( sm.At(4,3), -2.0*s + 2.0*q, t );
 }


// =============================================================================
// nonSymmetryTest — A_jk != A_kj for non-zero v
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::nonSymmetryTest()
 {
    Integral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
                                                  "diffusivity",
                                                  "velocity",
                                                  "fluid pressure",
                                                  "fluid pressure" );

    setElementScalar( { 1.0, 1.0, 1.0, 1.0 }, "diffusivity" );
    setElementVector( { 1.0, 0.0 },            "velocity"    );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, adv_disp );

    constexpr double s{ 0.5 };
    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    _equal( sm.At(0,4), -2.0*s + 2.0*q, t );  // A_04 = -1.0 + 500/6
    _equal( sm.At(4,0), -2.0*s - 2.0*q, t );  // A_40 = -1.0 - 500/6

    const double diff = sm.At(0,4) - sm.At(4,0);
    _equal( diff, 4.0*q, t );
 }


// =============================================================================
// zeroVelocityTest — v=0 -> symmetric result matching stiffness matrix
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::zeroVelocityTest()
 {
    Integral_dNT_op_dN_NT_v_dN_dV<2U> adv_disp( sg_->Database(),
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
        cout << "\nIntegral_DNT_op_DN_NT_v_DN_dV zeroVelocityTest:\n";
        sm.Out();
      }

    // symmetry check
    const std::vector<std::pair<size_t,size_t>> pairs{
        {0U,4U},{1U,4U},{2U,4U},{3U,4U}
    };
    for ( const auto& [i,j] : pairs )
      _equal( sm.At(i,j), sm.At(j,i), tol_ );

    // values must match pure stiffness matrix with E={1,2,3,4}
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
// matchesNumericalTest — analytical result must match numerical counterpart
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::matchesNumericalTest()
 {
    // build a second model with ISOPARAMETRIC_LINEAR_TRIANGLE for the
    // numerical operator
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container_num;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container_num );
    mesh_container_num.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container_num.EstablishElementConnectivity2D();

    Model<2U> sg_num( mesh_container_num, "CSMP-1phase-variables.txt" );
    sg_num.Region("Model").RenumberNodes();
    sg_num.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_num.InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
    sg_num.InputPropertyValue( "velocity",       makeVector( PLAIN, PLAIN, 0.0, 0.0 ) );

    // set identical operands on both models
    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );
    setElementVector( { 1.0, 0.5 },            "velocity"    );

    // set same values on numerical model
    {
      const csmp::Index key_d( sg_num.Database().StorageKey( "diffusivity" ) );
      const csmp::Index key_v( sg_num.Database().StorageKey( "velocity"    ) );
      Region<2>& domain( sg_num.Region("Model") );
      const std::vector<double> dvals{ 1.0, 2.0, 3.0, 4.0 };
      uint32_t i{0U};
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          (*it)->Store( key_d, ScalarVariable( PLAIN, dvals[i++] ) );
          VectorVariable<2U> vc;  vc(0) = 1.0;  vc(1) = 0.5;
          (*it)->Store( key_v, vc );
        }
    }

    // analytical operator
    Integral_dNT_op_dN_NT_v_dN_dV<2U> analytical( sg_->Database(),
                                                    "diffusivity",
                                                    "velocity",
                                                    "fluid pressure",
                                                    "fluid pressure" );

    // numerical operator
    NumIntegral_dNT_op_dN_NT_v_dN_dV<2U> numerical( sg_num.Database(),
                                                      "diffusivity",
                                                      "velocity",
                                                      "fluid pressure",
                                                      "fluid pressure" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm_anal( n ), sm_num( n );

    // assemble analytical
    {
      Region<2>& domain( sg_->Region("Model") );
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          analytical.GetOperands(        *(*it) );
          analytical.ComputeContribution(*(*it) );
          analytical.AssignToGlobal(     *(*it), sm_anal );
        }
    }

    // assemble numerical
    {
      Region<2>& domain( sg_num.Region("Model") );
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          numerical.GetOperands(        *(*it) );
          numerical.ComputeContribution(*(*it) );
          numerical.AssignToGlobal(     *(*it), sm_num );
        }
    }

    if ( verbose_ )
      {
        cout << "\nmatchesNumericalTest: Analytical:\n";
        sm_anal.Out();
        cout << "\nmatchesNumericalTest: Numerical:\n";
        sm_num.Out();
      }

    // all entries must match to within tolerance
    const std::vector<std::vector<size_t>> connected{
        { 0U, 1U, 3U, 4U },
        { 0U, 1U, 2U, 4U },
        { 1U, 2U, 3U, 4U },
        { 0U, 2U, 3U, 4U },
        { 0U, 1U, 2U, 3U, 4U }
    };

    constexpr double q{ 250.0 / 6.0 };
    const double t{ tol_ * q };

    for ( size_t i{0U}; i < n; ++i )
      for ( const size_t j : connected[i] )
        if ( sm_anal.HasEntry(i,j) && sm_num.HasEntry(i,j) )
          _equal( sm_anal.At(i,j), sm_num.At(i,j), t );
 }


// =============================================================================
// Helper methods
// =============================================================================

void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::setElementScalar(
    const vector<double>& val, const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR,
          "Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::setElementScalar",
          var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }


void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::setElementVector(
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


void Integral_dNT_lhsop_dN_NT_v_dN_dV_Test::calculateGlobalMatrix(
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

