// Integral_var_NT_lhsop_N_dV_Test.cpp
#include "Integral_var_NT_lhsop_N_dV_Test.h"
#include "Integral_var_NT_lhsop_N_dV.h"
#include "TRIANGLE_Interface.h"
#include "SparseMatrix.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

Integral_var_NT_lhsop_N_dV_Test::Integral_var_NT_lhsop_N_dV_Test( bool verbose )
  : verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nIntegral_var_NT_lhsop_N_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-2phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "permeability",              makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "diffusivity",               makeScalar( PLAIN, 1.0 ) ); // formerly conductivity (element)
    sg_->InputPropertyValue( "fluid pressure",            makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "nodal fluid volume source", makeScalar( PLAIN, 1.0 ) ); // formerly mobility (node)
 }


Integral_var_NT_lhsop_N_dV_Test::~Integral_var_NT_lhsop_N_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void Integral_var_NT_lhsop_N_dV_Test::run()
 {
    valueTest();
    compareConsistentTest();
    compareLumpedTest();
    lumpedTest();
    rowSumTest();
 }


// =============================================================================
// valueTest — verifies consistent matrix with nodal prefactor against
//             analytical solution
// =============================================================================

void Integral_var_NT_lhsop_N_dV_Test::valueTest()
 {
    Integral_var_NT_lhsop_N_dV<2U> integral( sg_->Database(),
                                              "diffusivity",
                                              "fluid pressure",
                                              "fluid pressure",
                                              "nodal fluid volume source" );

    setElementVariable( { 1.0, 2.0, 3.0, 4.0 },      "diffusivity"              );
    setNodeVariable(    { 1.0, 2.0, 3.0, 4.0, 5.0 }, "nodal fluid volume source" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, integral );

    if ( verbose_ )
      {
        cout << "\nvalueTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // Mesh: 4 right triangles, A = 62500, c = A/60 = 1041.667
    // Connectivity (0-indexed): e0={0,1,4}, e1={1,2,4}, e2={2,3,4}, e3={3,0,4}
    // Diffusivity:               E0=1,       E1=2,       E2=3,       E3=4
    // Nodal prefactor op:        {1,2,3,4,5} for nodes {0,1,2,3,4}
    //
    // Formula for linear triangle with nodal prefactor interpolation:
    //   M_jj^(e) = E * c * (6*op_j + 2*op_k + 2*op_l)
    //   M_jk^(e) = E * c * (2*op_j + 2*op_k + op_l)   j!=k, l=remaining node
    // -----------------------------------------------------------------------

    constexpr double c{ 62500.0 / 60.0 };  // = 1041.667
    const double t{ tol_ * c };

    // node 0: e0={0,1,4} E=1, e3={3,0,4} E=4
    const double M00 = c*( 1.0*(6*1+2*2+2*5) + 4.0*(6*1+2*4+2*5) );
    const double M01 = c*( 1.0*(2*1+2*2+1*5) );                      // e0, l=4
    const double M03 = c*( 4.0*(2*1+2*4+1*5) );                      // e3, l=4 — note j=3,k=0
    const double M04 = c*( 1.0*(2*1+2*5+1*2) + 4.0*(2*1+2*5+1*4) );  // e0 l=1, e3 l=3

    // node 1: e0={0,1,4} E=1, e1={1,2,4} E=2
    const double M11 = c*( 1.0*(6*2+2*1+2*5) + 2.0*(6*2+2*3+2*5) );
    const double M12 = c*( 2.0*(2*2+2*3+1*5) );                      // e1, l=4
    const double M14 = c*( 1.0*(2*2+2*5+1*1) + 2.0*(2*2+2*5+1*3) ); // e0 l=0, e1 l=2

    // node 2: e1={1,2,4} E=2, e2={2,3,4} E=3
    const double M22 = c*( 2.0*(6*3+2*2+2*5) + 3.0*(6*3+2*4+2*5) );
    const double M23 = c*( 3.0*(2*3+2*4+1*5) );                      // e2, l=4
    const double M24 = c*( 2.0*(2*3+2*5+1*2) + 3.0*(2*3+2*5+1*4) ); // e1 l=1, e2 l=3

    // node 3: e2={2,3,4} E=3, e3={3,0,4} E=4
    const double M33 = c*( 3.0*(6*4+2*3+2*5) + 4.0*(6*4+2*1+2*5) );
    const double M34 = c*( 3.0*(2*4+2*5+1*3) + 4.0*(2*4+2*5+1*1) ); // e2 l=2, e3 l=0

    // node 4: all 4 elements
    const double M44 = c*( 1.0*(6*5+2*1+2*2) + 2.0*(6*5+2*2+2*3)
                          + 3.0*(6*5+2*3+2*4) + 4.0*(6*5+2*4+2*1) );

    // diagonal
    _equal( sm.At(0,0), M00, t );
    _equal( sm.At(1,1), M11, t );
    _equal( sm.At(2,2), M22, t );
    _equal( sm.At(3,3), M33, t );
    _equal( sm.At(4,4), M44, t );

    // off-diagonal — connected pairs only
    _equal( sm.At(0,1), M01, t );
    _equal( sm.At(0,3), M03, t );
    _equal( sm.At(0,4), M04, t );

    _equal( sm.At(1,0), M01, t );  // symmetric
    _equal( sm.At(1,2), M12, t );
    _equal( sm.At(1,4), M14, t );

    _equal( sm.At(2,1), M12, t );  // symmetric
    _equal( sm.At(2,3), M23, t );
    _equal( sm.At(2,4), M24, t );

    _equal( sm.At(3,0), M03, t );  // symmetric
    _equal( sm.At(3,2), M23, t );  // symmetric
    _equal( sm.At(3,4), M34, t );

    _equal( sm.At(4,0), M04, t );  // symmetric
    _equal( sm.At(4,1), M14, t );  // symmetric
    _equal( sm.At(4,2), M24, t );  // symmetric
    _equal( sm.At(4,3), M34, t );  // symmetric
 }


// =============================================================================
// compareConsistentTest / compareLumpedTest
// =============================================================================

void Integral_var_NT_lhsop_N_dV_Test::compareConsistentTest()
 { compareTest( false ); }

void Integral_var_NT_lhsop_N_dV_Test::compareLumpedTest()
 { compareTest( true ); }

/**
 * When the nodal prefactor is uniform (all 1.0), the variable-coefficient
 * operator must produce the same result as the simple (constant) operator.
 */
void Integral_var_NT_lhsop_N_dV_Test::compareTest( bool lumped )
 {
    Integral_var_NT_lhsop_N_dV<2U> integral( sg_->Database(),
                                              "diffusivity",
                                              "fluid pressure",
                                              "fluid pressure",
                                              "nodal fluid volume source" );
    integral.LumpedFormulation( lumped );

    Integral_NT_lhsop_N_dV<2U> simple( sg_->Database(),
                                           "diffusivity",
                                           "fluid pressure",
                                           "fluid pressure" );
    simple.LumpedFormulation( lumped );

    // uniform nodal prefactor — both operators must agree
    setNodeVariable(    { 1.0, 1.0, 1.0, 1.0, 1.0 }, "nodal fluid volume source" );
    setElementVariable( { 1.0, 2.0, 3.0, 4.0 },       "diffusivity"              );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm_int( n ), sm_simp( n );
    calculateGlobalMatrix( sm_int,  integral );
    calculateGlobalMatrix( sm_simp, simple   );

    if ( verbose_ )
      {
        cout << "\ncompareTest (" << (lumped?"lumped":"consistent")
             << "): sm_int:\n";
        sm_int.Out();
        cout << "\ncompareTest (" << (lumped?"lumped":"consistent")
             << "): sm_simp:\n";
        sm_simp.Out();
      }

    // compareTest — replace full loop
    const std::vector<std::vector<size_t>> connected{
        { 0U, 1U, 3U, 4U },
        { 0U, 1U, 2U, 4U },
        { 1U, 2U, 3U, 4U },
        { 0U, 2U, 3U, 4U },
        { 0U, 1U, 2U, 3U, 4U }
    };

    for ( size_t i{0U}; i < n; ++i )
      for ( const size_t j : connected[i] )
        if ( !lumped || i == j )
          _equal( sm_int.At(i,j), sm_simp.At(i,j), tol_ );

 }


// =============================================================================
// lumpedTest — lumped diagonal must equal row sums of consistent matrix
// =============================================================================

void Integral_var_NT_lhsop_N_dV_Test::lumpedTest()
 {
    Integral_var_NT_lhsop_N_dV<2U> lumped( sg_->Database(),
                                            "diffusivity",
                                            "fluid pressure",
                                            "fluid pressure",
                                            "nodal fluid volume source" );
    lumped.LumpedFormulation( true );

    Integral_var_NT_lhsop_N_dV<2U> consistent( sg_->Database(),
                                                "diffusivity",
                                                "fluid pressure",
                                                "fluid pressure",
                                                "nodal fluid volume source" );
    consistent.LumpedFormulation( false );

    setNodeVariable(    { 1.0, 1.0, 2.0, 2.0, 5.0 }, "nodal fluid volume source" );
    setElementVariable( { 1.0, 2.0, 3.0, 4.0 },       "diffusivity"              );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm_lumped( n ), sm_consistent( n );
    calculateGlobalMatrix( sm_lumped,     lumped     );
    calculateGlobalMatrix( sm_consistent, consistent );

    if ( verbose_ )
      {
        cout << "\nlumpedTest: Lumped:\n";
        sm_lumped.Out();
        cout << "\nlumpedTest: Consistent:\n";
        sm_consistent.Out();
      }

    // lumped diagonal must equal row sums of consistent matrix
    const std::vector<std::vector<size_t>> connected{
        { 0U, 1U, 3U, 4U },
        { 0U, 1U, 2U, 4U },
        { 1U, 2U, 3U, 4U },
        { 0U, 2U, 3U, 4U },
        { 0U, 1U, 2U, 3U, 4U }
    };

    for ( size_t i{0U}; i < n; ++i )
      {
        double row_sum{0.0};
        for ( const size_t j : connected[i] )
          row_sum += sm_consistent.At(i,j);
        _equal( sm_lumped.At(i,i), row_sum, tol_ * row_sum );
        // there are no off-diagonal terms
      }
}



// =============================================================================
// rowSumTest — lumped diagonal must equal row sums of consistent matrix
//              for varying nodal prefactor
// =============================================================================

void Integral_var_NT_lhsop_N_dV_Test::rowSumTest()
 {
    Integral_var_NT_lhsop_N_dV<2U> lumped( sg_->Database(),
                                            "diffusivity",
                                            "fluid pressure",
                                            "fluid pressure",
                                            "nodal fluid volume source" );
    lumped.LumpedFormulation( true );

    Integral_var_NT_lhsop_N_dV<2U> consistent( sg_->Database(),
                                                "diffusivity",
                                                "fluid pressure",
                                                "fluid pressure",
                                                "nodal fluid volume source" );
    consistent.LumpedFormulation( false );

    // non-uniform nodal prefactor
    setNodeVariable(    { 1.0, 2.0, 3.0, 4.0, 5.0 }, "nodal fluid volume source" );
    setElementVariable( { 1.0, 2.0, 3.0, 4.0 },       "diffusivity"              );

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
        cout << "\nrowSumTest: Row sums:\n";
      }

    // lumpedTest and rowSumTest — use known connected pairs
    const std::vector<std::vector<size_t>> connected{
        { 0U, 1U, 3U, 4U },
        { 0U, 1U, 2U, 4U },
        { 1U, 2U, 3U, 4U },
        { 0U, 2U, 3U, 4U },
        { 0U, 1U, 2U, 3U, 4U }
    };

    for ( size_t i{0U}; i < n; ++i )
      {
        double row_sum{0.0};
        for ( const size_t j : connected[i] )
          row_sum += sm_consistent.At(i,j);
        _equal( sm_lumped.At(i,i), row_sum, tol_ * row_sum );
        // no off-diagonal terms
        if ( verbose_ ) cout << "  row " << i << ": " << row_sum << "\n";
      }
 }


// =============================================================================
// Helper methods
// =============================================================================

void Integral_var_NT_lhsop_N_dV_Test::setNodeVariable( const vector<double>& var,
                                                       const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( var.size() != domain.Nodes() )
      throw csmp::Exception( ERROR, "Integral_var_NT_lhsop_N_dV_Test::setNodeVariable",
                             var_name, "Variable size does not match number of nodes." );

    uint32_t i{0U};
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it )
      (*it)->Store( key, makeScalar( PLAIN, var[i++] ) );
 }


void Integral_var_NT_lhsop_N_dV_Test::showNodeVariable( const char* var_name ) const
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );
    uint32_t i{0U};
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it )
      cout << "Node " << i++ << ": " << (*it)->Read(key) << "\n";
 }


void Integral_var_NT_lhsop_N_dV_Test::setElementVariable( const vector<double>& var,
                                                           const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( var.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "Integral_var_NT_lhsop_N_dV_Test::setElementVariable",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, var[i++] ) );
 }


void Integral_var_NT_lhsop_N_dV_Test::calculateGlobalMatrix( SparseMatrix& sm,
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

