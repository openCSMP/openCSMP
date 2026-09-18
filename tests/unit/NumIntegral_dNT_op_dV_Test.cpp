// NumIntegral_dNT_op_dV_Test.cpp
#include "NumIntegral_dNT_op_dV_Test.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"
#include "VectorVariable.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

NumIntegral_dNT_op_dV_Test::NumIntegral_dNT_op_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nNumIntegral_dNT_op_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );

    // create a vector property for the gravity term
    sg_->CreateProperty( "gravity term", "GT", "Pa/m", VECTOR, ELEMENT );
    sg_->InputPropertyValue( "gravity term", makeVector( PLAIN, PLAIN, 0.0, 0.0 ) );
 }


NumIntegral_dNT_op_dV_Test::~NumIntegral_dNT_op_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void NumIntegral_dNT_op_dV_Test::run()
 {
    xDirectionTest();
    yDirectionTest();
    diagonalTest();
    conservationTest();
 }


// =============================================================================
// xDirectionTest — v=(1,0), isolates x-gradient components
// =============================================================================

void NumIntegral_dNT_op_dV_Test::xDirectionTest()
 {
    NumIntegral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                        "gravity term",
                                        "fluid pressure" );

    // uniform v=(1,0) across all elements
    setElementVector( {1.0, 0.0}, "gravity term" );

    std::vector<double> rhs;
    calculateRHS( rhs, grad_op );

    if ( verbose_ )
      {
        cout << "\nxDirectionTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values for v=(1,0)
    // -----------------------------------------------------------------------
    // RHS[j] = sum_{e containing j} (g_j^x) / 2
    //
    // g_j^x values:
    //   e0: g0^x=-250, g1^x=-250, g4^x=500
    //   e1: g1^x=-250, g2^x=250,  g4^x=0
    //   e2: g2^x=250,  g3^x=250,  g4^x=-500
    //   e3: g3^x=250,  g0^x=-250, g4^x=0
    //
    // Node 0: e0(-250) + e3(-250) = -500  -> /2 = -250
    // Node 1: e0(-250) + e1(-250) = -500  -> /2 = -250
    // Node 2: e1(250)  + e2(250)  =  500  -> /2 =  250
    // Node 3: e2(250)  + e3(250)  =  500  -> /2 =  250
    // Node 4: e0(500)  + e1(0) + e2(-500) + e3(0) = 0 -> /2 = 0
    // -----------------------------------------------------------------------

    const double t{ tol_ * 250.0 };

    _equal( rhs[0], -250.0, t );
    _equal( rhs[1], -250.0, t );
    _equal( rhs[2],  250.0, t );
    _equal( rhs[3],  250.0, t );
    _equal( rhs[4],    0.0, t );
 }


// =============================================================================
// yDirectionTest — v=(0,1), isolates y-gradient components
// =============================================================================

void NumIntegral_dNT_op_dV_Test::yDirectionTest()
 {
    NumIntegral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                        "gravity term",
                                        "fluid pressure" );

    // uniform v=(0,1) across all elements
    setElementVector( {0.0, 1.0}, "gravity term" );

    std::vector<double> rhs;
    calculateRHS( rhs, grad_op );

    if ( verbose_ )
      {
        cout << "\nyDirectionTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values for v=(0,1)
    // -----------------------------------------------------------------------
    // RHS[j] = sum_{e containing j} (g_j^y) / 2
    //
    // g_j^y values:
    //   e0: g0^y=250,  g1^y=-250, g4^y=0
    //   e1: g1^y=-250, g2^y=-250, g4^y=500
    //   e2: g2^y=-250, g3^y=250,  g4^y=0
    //   e3: g3^y=250,  g0^y=250,  g4^y=-500
    //
    // Node 0: e0(250)  + e3(250)  =  500  -> /2 =  250
    // Node 1: e0(-250) + e1(-250) = -500  -> /2 = -250
    // Node 2: e1(-250) + e2(-250) = -500  -> /2 = -250
    // Node 3: e2(250)  + e3(250)  =  500  -> /2 =  250
    // Node 4: e0(0) + e1(500) + e2(0) + e3(-500) = 0 -> /2 = 0
    // -----------------------------------------------------------------------

    const double t{ tol_ * 250.0 };

    _equal( rhs[0],  250.0, t );
    _equal( rhs[1], -250.0, t );
    _equal( rhs[2], -250.0, t );
    _equal( rhs[3],  250.0, t );
    _equal( rhs[4],    0.0, t );
 }


// =============================================================================
// diagonalTest — v=(1,1), combined x+y
// =============================================================================

void NumIntegral_dNT_op_dV_Test::diagonalTest()
 {
    NumIntegral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                        "gravity term",
                                        "fluid pressure" );

    // uniform v=(1,1) across all elements
    setElementVector( {1.0, 1.0}, "gravity term" );

    std::vector<double> rhs;
    calculateRHS( rhs, grad_op );

    if ( verbose_ )
      {
        cout << "\ndiagonalTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values for v=(1,1)
    // -----------------------------------------------------------------------
    // RHS[j] = xDirectionTest[j] + yDirectionTest[j]
    //
    // Node 0: -250 + 250  =    0
    // Node 1: -250 + -250 = -500
    // Node 2:  250 + -250 =    0
    // Node 3:  250 + 250  =  500
    // Node 4:    0 + 0    =    0
    // -----------------------------------------------------------------------

    const double t{ tol_ * 250.0 };

    _equal( rhs[0],    0.0, t );
    _equal( rhs[1], -500.0, t );
    _equal( rhs[2],    0.0, t );
    _equal( rhs[3],  500.0, t );
    _equal( rhs[4],    0.0, t );
 }


// =============================================================================
// conservationTest — sum of RHS entries must be zero for any uniform v
// =============================================================================

void NumIntegral_dNT_op_dV_Test::conservationTest()
 {
    NumIntegral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                        "gravity term",
                                        "fluid pressure" );

    // test with several different uniform vector fields
    const std::vector<std::array<double,2>> test_vectors{
        { 1.0,  0.0 },
        { 0.0,  1.0 },
        { 1.0,  1.0 },
        { 3.0, -2.0 },
        {-1.0,  4.0 }
    };

    for ( const auto& v : test_vectors )
      {
        setElementVector( v, "gravity term" );

        std::vector<double> rhs;
        calculateRHS( rhs, grad_op );

        // sum of all RHS entries must be zero — discrete divergence theorem
        double total{0.0};
        for ( const double r : rhs ) total += r;

        if ( verbose_ )
          cout << "\nconservationTest: v=(" << v[0] << "," << v[1]
               << ") -> sum = " << total << "\n";

        _equal( total, 0.0, tol_ * 250.0 );
      }
 }


// =============================================================================
// Helper methods
// =============================================================================

void NumIntegral_dNT_op_dV_Test::setElementVector( const std::array<double,2>& v,
                                                    const char* var_name )
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


void NumIntegral_dNT_op_dV_Test::calculateRHS( std::vector<double>& rhs,
                                                MathOperatorRHS<2U>& oper )
 {
    Region<2>& domain( sg_->Region("Model") );
    const size_t n{ domain.Nodes() };

    rhs.assign( n, 0.0 );

    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        oper.GetOperands(        *(*it) );
        oper.ComputeContribution(*(*it) );
        oper.AssignToGlobal(     *(*it), rhs );
      }
 }

} // namespace csmp

