// Integral_dNT_rhsop_dV_Test.cpp
#include "Integral_dNT_rhsop_dV_Test.h"

#include "Integral_dNT_op_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "TRIANGLE_Interface.h"
#include "FiniteElement.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

Integral_dNT_rhsop_dV_Test::Integral_dNT_rhsop_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nIntegral_dNT_rhsop_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_->CreateProperty( "gravity term", "GT", "Pa/m", VECTOR, ELEMENT );
    sg_->InputPropertyValue( "gravity term", makeVector( PLAIN, PLAIN, 0.0, 0.0 ) );
 }


Integral_dNT_rhsop_dV_Test::~Integral_dNT_rhsop_dV_Test()
 { delete sg_; }


// =============================================================================
// Test runner
// =============================================================================

void Integral_dNT_rhsop_dV_Test::run()
 {
    xDirectionTest();
    yDirectionTest();
    diagonalTest();
    conservationTest();
    matchesNumericalTest();
 }


// =============================================================================
// xDirectionTest — v=(1,0)
// =============================================================================

void Integral_dNT_rhsop_dV_Test::xDirectionTest()
 {
    Integral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                     "gravity term",
                                     "fluid pressure" );

    setElementVector( {1.0, 0.0}, "gravity term" );

    std::vector<double> rhs;
    calculateRHS( rhs, grad_op );

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_op_dV xDirectionTest: RHS vector:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Identical reference values to NumIntegral_dNT_rhsop_dV_Test::xDirectionTest
    // RHS[j] = sum_{e containing j} g_j^x / 2
    //   Node 0: e0(-250) + e3(-250) = -500 -> /2 = -250
    //   Node 1: e0(-250) + e1(-250) = -500 -> /2 = -250
    //   Node 2: e1(250)  + e2(250)  =  500 -> /2 =  250
    //   Node 3: e2(250)  + e3(250)  =  500 -> /2 =  250
    //   Node 4: e0(500)+e1(0)+e2(-500)+e3(0) = 0 -> /2 = 0
    // -----------------------------------------------------------------------

    const double t{ tol_ * 250.0 };

    _equal( rhs[0], -250.0, t );
    _equal( rhs[1], -250.0, t );
    _equal( rhs[2],  250.0, t );
    _equal( rhs[3],  250.0, t );
    _equal( rhs[4],    0.0, t );
 }


// =============================================================================
// yDirectionTest — v=(0,1)
// =============================================================================

void Integral_dNT_rhsop_dV_Test::yDirectionTest()
 {
    Integral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                     "gravity term",
                                     "fluid pressure" );

    setElementVector( {0.0, 1.0}, "gravity term" );

    std::vector<double> rhs;
    calculateRHS( rhs, grad_op );

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_op_dV yDirectionTest: RHS vector:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    const double t{ tol_ * 250.0 };

    _equal( rhs[0],  250.0, t );
    _equal( rhs[1], -250.0, t );
    _equal( rhs[2], -250.0, t );
    _equal( rhs[3],  250.0, t );
    _equal( rhs[4],    0.0, t );
 }


// =============================================================================
// diagonalTest — v=(1,1)
// =============================================================================

void Integral_dNT_rhsop_dV_Test::diagonalTest()
 {
    Integral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                     "gravity term",
                                     "fluid pressure" );

    setElementVector( {1.0, 1.0}, "gravity term" );

    std::vector<double> rhs;
    calculateRHS( rhs, grad_op );

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_op_dV diagonalTest: RHS vector:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    const double t{ tol_ * 250.0 };

    _equal( rhs[0],    0.0, t );
    _equal( rhs[1], -500.0, t );
    _equal( rhs[2],    0.0, t );
    _equal( rhs[3],  500.0, t );
    _equal( rhs[4],    0.0, t );
 }


// =============================================================================
// conservationTest — sum of RHS must be zero for any uniform v
// =============================================================================

void Integral_dNT_rhsop_dV_Test::conservationTest()
 {
    Integral_dNT_op_dV<2U> grad_op( sg_->Database(),
                                     "gravity term",
                                     "fluid pressure" );

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

        double total{0.0};
        for ( const double r : rhs ) total += r;

        if ( verbose_ )
          cout << "\nIntegral_dNT_op_dV conservationTest: v=("
               << v[0] << "," << v[1] << ") -> sum = " << total << "\n";

        _equal( total, 0.0, tol_ * 250.0 );
      }
 }


// =============================================================================
// matchesNumericalTest — analytical result must match numerical counterpart
// =============================================================================

void Integral_dNT_rhsop_dV_Test::matchesNumericalTest()
 {
    // build a second model with ISOPARAMETRIC_LINEAR_TRIANGLE
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container_num;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container_num );
    mesh_container_num.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container_num.EstablishElementConnectivity2D();

    Model<2U> sg_num( mesh_container_num, "CSMP-1phase-variables.txt" );
    sg_num.Region("Model").RenumberNodes();
    sg_num.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_num.CreateProperty( "gravity term", "GT", "Pa/m", VECTOR, ELEMENT );
    sg_num.InputPropertyValue( "gravity term", makeVector( PLAIN, PLAIN, 0.0, 0.0 ) );

    // set identical vector operand on both models
    const std::array<double,2> v{ 1.0, 0.5 };
    setElementVector( v, "gravity term" );

    {
      const csmp::Index key( sg_num.Database().StorageKey( "gravity term" ) );
      Region<2>& domain( sg_num.Region("Model") );
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          VectorVariable<2U> vc;
          vc(0) = v[0];  vc(1) = v[1];
          (*it)->Store( key, vc );
        }
    }

    // analytical operator
    Integral_dNT_op_dV<2U> analytical( sg_->Database(),
                                        "gravity term",
                                        "fluid pressure" );

    // numerical operator
    NumIntegral_dNT_op_dV<2U> numerical( sg_num.Database(),
                                          "gravity term",
                                          "fluid pressure" );

    std::vector<double> rhs_anal( sg_->Region("Model").Nodes(),   0.0 );
    std::vector<double> rhs_num(  sg_num.Region("Model").Nodes(), 0.0 );

    {
      Region<2>& domain( sg_->Region("Model") );
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          analytical.GetOperands(        *(*it) );
          analytical.ComputeContribution(*(*it) );
          analytical.AssignToGlobal(     *(*it), rhs_anal );
        }
    }

    {
      Region<2>& domain( sg_num.Region("Model") );
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          numerical.GetOperands(        *(*it) );
          numerical.ComputeContribution(*(*it) );
          numerical.AssignToGlobal(     *(*it), rhs_num );
        }
    }

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_op_dV matchesNumericalTest:\n";
        for ( size_t i{0U}; i < rhs_anal.size(); ++i )
          cout << "  analytical[" << i << "] = " << rhs_anal[i]
               << "  numerical[" << i << "] = " << rhs_num[i] << "\n";
      }

    const double t{ tol_ * 250.0 };

    assert( rhs_anal.size() == rhs_num.size() );
    for ( size_t i{0U}; i < rhs_anal.size(); ++i )
      _equal( rhs_anal[i], rhs_num[i], t );
 }


// =============================================================================
// Helper methods
// =============================================================================

void Integral_dNT_rhsop_dV_Test::setElementVector( const array<double,2>& v,
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


void Integral_dNT_rhsop_dV_Test::calculateRHS( vector<double>& rhs,
                                               MathOperatorRHS<2U>& oper )
 {
    Region<2>& domain( sg_->Region("Model") );
    rhs.assign( domain.Nodes(), 0.0 );

    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        oper.GetOperands(        *(*it) );
        oper.ComputeContribution(*(*it) );
        oper.AssignToGlobal(     *(*it), rhs );
      }
 }

} // namespace csmp

