// NumIntegral_dNT_rhsop_dN_dV_Test.cpp
#include "NumIntegral_dNT_rhsop_dN_dV_Test.h"

#include "NumIntegral_dNT_rhsop_dN_dV.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

NumIntegral_dNT_rhsop_dN_dV_Test::NumIntegral_dNT_rhsop_dN_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nNumIntegral_dNT_rhsop_dN_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "porosity",       makeScalar( PLAIN, 1.0 ) );
 }


NumIntegral_dNT_rhsop_dN_dV_Test::~NumIntegral_dNT_rhsop_dN_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void NumIntegral_dNT_rhsop_dN_dV_Test::run()
 {
    scalarTest();
    multiplierTest();
    uniformTest();
 }


// =============================================================================
// scalarTest — RHS = K * u for non-uniform nodal values
// =============================================================================

void NumIntegral_dNT_rhsop_dN_dV_Test::scalarTest()
 {
    NumIntegral_dNT_rhsop_dN_dV<2U> rhs_op( sg_->Database(),
                                              "diffusivity",
                                              "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity"    );
    setNodeVariable(  { 1.0, 2.0, 3.0, 4.0, 5.0 }, "fluid pressure" );

    std::vector<double> rhs;
    calculateRHS( rhs, rhs_op );

    if ( verbose_ )
      {
        cout << "\nscalarTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // RHS[j] = sum_k K_jk * u_k
    //
    // Stiffness matrix (s=0.5, E={1,2,3,4}):
    //   K_00=5s,  K_04=-5s  (only non-zero entries for node 0)
    //   K_11=3s,  K_14=-3s
    //   K_22=5s,  K_24=-5s
    //   K_33=7s,  K_34=-7s
    //   K_44=20s, K_40=-5s, K_41=-3s, K_42=-5s, K_43=-7s
    //
    // Nodal values u = {1,2,3,4,5}:
    //   RHS[0] = K_00*u0 + K_04*u4 = 5s*1  + (-5s)*5  = 5s  - 25s = -20s = -10
    //   RHS[1] = K_11*u1 + K_14*u4 = 3s*2  + (-3s)*5  = 6s  - 15s =  -9s = -4.5
    //   RHS[2] = K_22*u2 + K_24*u4 = 5s*3  + (-5s)*5  = 15s - 25s = -10s = -5
    //   RHS[3] = K_33*u3 + K_34*u4 = 7s*4  + (-7s)*5  = 28s - 35s =  -7s = -3.5
    //   RHS[4] = K_44*u4 + K_40*u0 + K_41*u1 + K_42*u2 + K_43*u3
    //          = 20s*5 + (-5s)*1 + (-3s)*2 + (-5s)*3 + (-7s)*4
    //          = 100s - 5s - 6s - 15s - 28s = 46s = 23
    // -----------------------------------------------------------------------

    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    _equal( rhs[0], -20.0*s, t );
    _equal( rhs[1],  -9.0*s, t );
    _equal( rhs[2], -10.0*s, t );
    _equal( rhs[3],  -7.0*s, t );
    _equal( rhs[4],  46.0*s, t );
 }


// =============================================================================
// multiplierTest — RHS = m * K * u where m is a scalar element multiplier
// =============================================================================

void NumIntegral_dNT_rhsop_dN_dV_Test::multiplierTest()
 {
    NumIntegral_dNT_rhsop_dN_dV<2U> rhs_op( sg_->Database(),
                                              "porosity",
                                              "diffusivity",
                                              "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );
    setElementScalar( { 2.0, 2.0, 2.0, 2.0 }, "porosity"    );
    setNodeVariable(  { 1.0, 2.0, 3.0, 4.0, 5.0 }, "fluid pressure" );

    std::vector<double> rhs;
    calculateRHS( rhs, rhs_op );

    if ( verbose_ )
      {
        cout << "\nmultiplierTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // With uniform multiplier m=2 the result is simply 2 * scalarTest result:
    //   RHS[j] = m * sum_k K_jk * u_k
    //
    // However since m is per-element and K_jk has contributions from two
    // elements per boundary node, we must apply m per element:
    //
    // For node 0 (in e0 with m=2 and e3 with m=2):
    //   RHS[0] = m*(K_00^e0 + K_00^e3)*u0 + m*(K_04^e0 + K_04^e3)*u4
    //          = 2*(-20s) = -20  (uniform m=2 simply doubles)
    //
    // Since m=2 uniformly across all elements:
    //   RHS[j] = 2 * scalarTest_RHS[j]
    // -----------------------------------------------------------------------

    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    _equal( rhs[0], 2.0*(-20.0*s), t );
    _equal( rhs[1], 2.0*( -9.0*s), t );
    _equal( rhs[2], 2.0*(-10.0*s), t );
    _equal( rhs[3], 2.0*( -7.0*s), t );
    _equal( rhs[4], 2.0*( 46.0*s), t );
 }



// =============================================================================
// uniformTest — uniform nodal values -> RHS must be zero (row sum property)
// =============================================================================

void NumIntegral_dNT_rhsop_dN_dV_Test::uniformTest()
 {
    NumIntegral_dNT_rhsop_dN_dV<2U> rhs_op( sg_->Database(),
                                              "diffusivity",
                                              "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity"     );
    setNodeVariable(  { 1.0, 1.0, 1.0, 1.0, 1.0 }, "fluid pressure" );

    std::vector<double> rhs;
    calculateRHS( rhs, rhs_op );

    if ( verbose_ )
      {
        cout << "\nuniformTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // With uniform u=1, RHS[j] = sum_k K_jk * 1 = row sum of K = 0
    // This is the discrete conservation property of the stiffness matrix
    // -----------------------------------------------------------------------

    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    for ( size_t i{0U}; i < rhs.size(); ++i )
      _equal( rhs[i], 0.0, t );
 }




// =============================================================================
// Helper methods
// =============================================================================

void NumIntegral_dNT_rhsop_dN_dV_Test::setNodeVariable( const vector<double>& val,
                                                         const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Nodes() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV_Test::setNodeVariable",
                             var_name, "Variable size does not match number of nodes." );

    uint32_t i{0U};
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it )
      (*it)->Store( key, makeScalar( PLAIN, val[i++] ) );
 }


void NumIntegral_dNT_rhsop_dN_dV_Test::setElementScalar( const vector<double>& val,
                                                          const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV_Test::setElementScalar",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }


void NumIntegral_dNT_rhsop_dN_dV_Test::calculateRHS( std::vector<double>& rhs,
                                                      MathOperatorRHS<2U>& oper )
 {
    Region<2>& domain( sg_->Region("Model") );
    const size_t n{ domain.Nodes() };

    rhs.assign( n, 0.0 );

    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        oper.GetOperands(        *(*it) );
        oper.ComputeContribution(*(*it) );
        oper.AssignToGlobal( *(*it), rhs );
      }
 }

} // namespace csmp


