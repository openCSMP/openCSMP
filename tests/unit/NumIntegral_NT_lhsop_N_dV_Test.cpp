// NumIntegral_NT_lhsop_N_dV_Test.cpp
#include "NumIntegral_NT_lhsop_N_dV_Test.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

NumIntegral_NT_lhsop_N_dV_Test::NumIntegral_NT_lhsop_N_dV_Test( bool verbose )
  : verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nNumIntegral_NT_lhsop_N_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure",            makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "diffusivity",               makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "permeability",              makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "nodal fluid volume source", makeScalar( PLAIN, 1.0 ) );
 }


NumIntegral_NT_lhsop_N_dV_Test::~NumIntegral_NT_lhsop_N_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void NumIntegral_NT_lhsop_N_dV_Test::run()
 {
    consistentTest();
    lumpedTest();
    rowSumTest();
    symmetryTest();
 }


// =============================================================================
// consistentTest — verifies consistent mass matrix against analytical solution
// =============================================================================

void NumIntegral_NT_lhsop_N_dV_Test::consistentTest()
 {
    NumIntegral_NT_lhsop_N_dV<2U> integral( sg_->Database(),
                                             "diffusivity",
                                             "fluid pressure",
                                             "fluid pressure" );

    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, integral );

    if ( verbose_ )
      {
        cout << "\nconsistentTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // Mesh: 4 right triangles, A = 62500, A/12 = 5208.333
    // Connectivity (0-indexed): e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}
    // Diffusivity:               E0=1,       E1=2,       E2=3,       E3=4
    //
    // For a linear triangle with constant element operand E:
    //   M_jj^(e) = E * (A/6)
    //   M_jk^(e) = E * (A/12)   j != k
    //
    // Node 0: in e0(E=1) and e3(E=4)
    //   M_00 = (A/6)*(1+4)          = 5*A/6
    // Node 1: in e0(E=1) and e1(E=2)
    //   M_11 = (A/6)*(1+2)          = 3*A/6
    // Node 2: in e1(E=2) and e2(E=3)
    //   M_22 = (A/6)*(2+3)          = 5*A/6
    // Node 3: in e2(E=3) and e3(E=4)
    //   M_33 = (A/6)*(3+4)          = 7*A/6
    // Node 4: in all 4 elements
    //   M_44 = (A/6)*(1+2+3+4)      = 10*A/6
    //
    // Off-diagonal (shared edge):
    //   M_01 = (A/12)*E0*1          = 1*A/12   (e0)
    //   M_03 = (A/12)*E3*1          = 4*A/12   (e3)
    //   M_04 = (A/12)*(E0+E3)       = 5*A/12   (e0+e3)
    //   M_12 = (A/12)*E1*1          = 2*A/12   (e1)
    //   M_14 = (A/12)*(E0+E1)       = 3*A/12   (e0+e1)
    //   M_23 = (A/12)*E2*1          = 3*A/12   (e2)
    //   M_24 = (A/12)*(E1+E2)       = 5*A/12   (e1+e2)
    //   M_34 = (A/12)*(E2+E3)       = 7*A/12   (e2+e3)
    // -----------------------------------------------------------------------

    constexpr double A{  62500.0 };
    constexpr double A6{  A / 6.0  };
    constexpr double A12{ A / 12.0 };
    const double t{ tol_ * A12 };

    // diagonal
    _equal( sm.At(0,0),  5.0*A6,  t );
    _equal( sm.At(1,1),  3.0*A6,  t );
    _equal( sm.At(2,2),  5.0*A6,  t );
    _equal( sm.At(3,3),  7.0*A6,  t );
    _equal( sm.At(4,4), 10.0*A6,  t );

    // off-diagonal — connected pairs only
    _equal( sm.At(0,1),  1.0*A12, t );
    _equal( sm.At(0,3),  4.0*A12, t );
    _equal( sm.At(0,4),  5.0*A12, t );

    _equal( sm.At(1,0),  1.0*A12, t );  // symmetric
    _equal( sm.At(1,2),  2.0*A12, t );
    _equal( sm.At(1,4),  3.0*A12, t );

    _equal( sm.At(2,1),  2.0*A12, t );  // symmetric
    _equal( sm.At(2,3),  3.0*A12, t );
    _equal( sm.At(2,4),  5.0*A12, t );

    _equal( sm.At(3,0),  4.0*A12, t );  // symmetric
    _equal( sm.At(3,2),  3.0*A12, t );  // symmetric
    _equal( sm.At(3,4),  7.0*A12, t );

    _equal( sm.At(4,0),  5.0*A12, t );  // symmetric
    _equal( sm.At(4,1),  3.0*A12, t );  // symmetric
    _equal( sm.At(4,2),  5.0*A12, t );  // symmetric
    _equal( sm.At(4,3),  7.0*A12, t );  // symmetric
 }


// =============================================================================
// lumpedTest — lumped diagonal must equal row sums of consistent matrix
// =============================================================================

void NumIntegral_NT_lhsop_N_dV_Test::lumpedTest()
 {
    NumIntegral_NT_lhsop_N_dV<2U> lumped( sg_->Database(),
                                           "diffusivity",
                                           "fluid pressure",
                                           "fluid pressure" );
    lumped.LumpedFormulation( true );

    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, lumped );

    if ( verbose_ )
      {
        cout << "\nlumpedTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Analytical reference: lumped diagonal = row sum of consistent matrix
    // Node 0: M_00 + M_01 + M_03 + M_04 = 10*A/12 + 1*A/12 + 4*A/12 + 5*A/12 = 20*A/12
    // Node 1: M_11 + M_10 + M_12 + M_14 =  6*A/12 + 1*A/12 + 2*A/12 + 3*A/12 = 12*A/12
    // Node 2: M_22 + M_21 + M_23 + M_24 = 10*A/12 + 2*A/12 + 3*A/12 + 5*A/12 = 20*A/12
    // Node 3: M_33 + M_30 + M_32 + M_34 = 14*A/12 + 4*A/12 + 3*A/12 + 7*A/12 = 28*A/12
    // Node 4: M_44 + M_40 + M_41 + M_42 + M_43
    //       = 20*A/12 + 5*A/12 + 3*A/12 + 5*A/12 + 7*A/12 = 40*A/12
    // -----------------------------------------------------------------------

    constexpr double A{   62500.0 };
    constexpr double A12{ A / 12.0 };
    const double t{ tol_ * A12 };

    _equal( sm.At(0,0), 20.0*A12, t );
    _equal( sm.At(1,1), 12.0*A12, t );
    _equal( sm.At(2,2), 20.0*A12, t );
    _equal( sm.At(3,3), 28.0*A12, t );
    _equal( sm.At(4,4), 40.0*A12, t );

    // the lumped matrix must have exactly 5 non-zero entries (diagonal only)
    _equal( static_cast<double>( sm.Entries() ), 5.0, 0.0 );
 }




// =============================================================================
// rowSumTest — lumped diagonal must equal row sums of consistent matrix
// =============================================================================

void NumIntegral_NT_lhsop_N_dV_Test::rowSumTest()
 {
    NumIntegral_NT_lhsop_N_dV<2U> lumped( sg_->Database(),
                                           "diffusivity",
                                           "fluid pressure",
                                           "fluid pressure" );
    lumped.LumpedFormulation( true );

    NumIntegral_NT_lhsop_N_dV<2U> consistent( sg_->Database(),
                                               "diffusivity",
                                               "fluid pressure",
                                               "fluid pressure" );
    consistent.LumpedFormulation( false );

    // set operands explicitly — rowSumTest must be self-contained
    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm_lumped( n ), sm_consistent( n );
    calculateGlobalMatrix( sm_lumped,     lumped     );
    calculateGlobalMatrix( sm_consistent, consistent );

    if ( verbose_ )
      {
        cout << "\nrowSumTest: Lumped:\n";
        sm_lumped.Out();
        cout << "\nrowSumTest: Consistent:\n";
        sm_consistent.Out();
      }

    // known connected pairs from mesh topology — avoids EntryExists entirely
    // Connectivity: e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}
    // Node 0 connects to: 0,1,3,4
    // Node 1 connects to: 0,1,2,4
    // Node 2 connects to: 1,2,3,4
    // Node 3 connects to: 0,2,3,4
    // Node 4 connects to: 0,1,2,3,4
    const std::vector<std::vector<size_t>> connected{
        { 0U, 1U, 3U, 4U },        // node 0
        { 0U, 1U, 2U, 4U },        // node 1
        { 1U, 2U, 3U, 4U },        // node 2
        { 0U, 2U, 3U, 4U },        // node 3
        { 0U, 1U, 2U, 3U, 4U }     // node 4
    };

    for ( size_t i{0U}; i < n; ++i )
      {
        // row sum using only known connected entries
        double row_sum{0.0};
        for ( const size_t j : connected[i] )
          row_sum += sm_consistent.At(i,j);

        _equal( sm_lumped.At(i,i), row_sum, tol_ * row_sum );
      }
 }



// =============================================================================
// symmetryTest — consistent mass matrix must be symmetric
// =============================================================================

void NumIntegral_NT_lhsop_N_dV_Test::symmetryTest()
 {
    NumIntegral_NT_lhsop_N_dV<2U> integral( sg_->Database(),
                                             "diffusivity",
                                             "fluid pressure",
                                             "fluid pressure" );

    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, integral );

    // check symmetry only for known connected off-diagonal pairs
    const std::vector<std::pair<size_t,size_t>> pairs{
        {0U,1U}, {0U,3U}, {0U,4U},
        {1U,2U}, {1U,4U},
        {2U,3U}, {2U,4U},
        {3U,4U}
    };

    for ( const auto& [i,j] : pairs )
      _equal( sm.At(i,j), sm.At(j,i), tol_ );
 }




// =============================================================================
// Helper methods
// =============================================================================

void NumIntegral_NT_lhsop_N_dV_Test::setNodeVariable( const vector<double>& var,
                                                       const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( var.size() != domain.Nodes() )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV_Test::setNodeVariable",
                             var_name, "Variable size does not match number of nodes." );

    uint32_t i{0U};
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it )
      (*it)->Store( key, makeScalar( PLAIN, var[i++] ) );
 }



void NumIntegral_NT_lhsop_N_dV_Test::setElementVariable( const vector<double>& var,
                                                          const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( var.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV_Test::setElementVariable",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, var[i++] ) );
 }



void NumIntegral_NT_lhsop_N_dV_Test::calculateGlobalMatrix( SparseMatrix& sm,
                                                             MathOperatorLHS<2U>& oper )
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

