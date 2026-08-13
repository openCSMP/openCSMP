// Integral_dNT_rhsop_dN_dV_Test.cpp
#include "Integral_dNT_rhsop_dN_dV_Test.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

Integral_dNT_rhsop_dN_dV_Test::Integral_dNT_rhsop_dN_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nIntegral_dNT_rhsop_dN_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "porosity",       makeScalar( PLAIN, 1.0 ) );
 }


Integral_dNT_rhsop_dN_dV_Test::~Integral_dNT_rhsop_dN_dV_Test()
 { delete sg_; }


// =============================================================================
// Test runner
// =============================================================================

void Integral_dNT_rhsop_dN_dV_Test::run()
 {
    scalarTest();
    uniformTest();
    matchesNumericalTest(); 
 }


// =============================================================================
// scalarTest — RHS = K * u for non-uniform nodal values
// =============================================================================

void Integral_dNT_rhsop_dN_dV_Test::scalarTest()
 {
    Integral_dNT_rhsop_dN_dV<2U> rhs_op( sg_->Database(),
                                          "diffusivity",
                                          "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity"     );
    setNodeVariable(  { 1.0, 2.0, 3.0, 4.0, 5.0 }, "fluid pressure" );

    std::vector<double> rhs;
    calculateRHS( rhs, rhs_op );

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_rhsop_dN_dV scalarTest: RHS vector:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Identical reference values to NumIntegral_dNT_rhsop_dN_dV_Test::scalarTest
    // RHS[j] = sum_k K_jk * u_k, s=0.5, E={1,2,3,4}, u={1,2,3,4,5}
    //   RHS[0] = K_00*1 + K_04*5 = 5s*1  + (-5s)*5  = -20s = -10
    //   RHS[1] = K_11*2 + K_14*5 = 3s*2  + (-3s)*5  =  -9s = -4.5
    //   RHS[2] = K_22*3 + K_24*5 = 5s*3  + (-5s)*5  = -10s = -5
    //   RHS[3] = K_33*4 + K_34*5 = 7s*4  + (-7s)*5  =  -7s = -3.5
    //   RHS[4] = K_44*5 + K_40*1 + K_41*2 + K_42*3 + K_43*4
    //          = 20s*5 + (-5s)*1 + (-3s)*2 + (-5s)*3 + (-7s)*4 = 46s = 23
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
// uniformTest — uniform nodal values -> RHS must be zero (row sum property)
// =============================================================================

void Integral_dNT_rhsop_dN_dV_Test::uniformTest()
 {
    Integral_dNT_rhsop_dN_dV<2U> rhs_op( sg_->Database(),
                                          "diffusivity",
                                          "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity"     );
    setNodeVariable(  { 1.0, 1.0, 1.0, 1.0, 1.0 }, "fluid pressure" );

    std::vector<double> rhs;
    calculateRHS( rhs, rhs_op );

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_rhsop_dN_dV uniformTest: RHS vector:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    for ( size_t i{0U}; i < rhs.size(); ++i )
      _equal( rhs[i], 0.0, t );
 }



// =============================================================================
// matchesNumericalTest — analytical and numerical operators must agree
//                        exactly on a linear simplex mesh
// =============================================================================

void Integral_dNT_rhsop_dN_dV_Test::matchesNumericalTest()
 {
    // build a single model with ISOPARAMETRIC_LINEAR_TRIANGLE
    // — for linear simplex elements analytical and numerical integration
    // are mathematically identical, so both operators must agree exactly
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    Model<2U> sg_iso( mesh_container, "CSMP-1phase-variables.txt" );
    sg_iso.Region("Model").RenumberNodes();
    sg_iso.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_iso.InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );

    // set identical operands on the isoparametric model
    const std::vector<double> dvals{ 1.0, 2.0, 3.0, 4.0 };
    const std::vector<double> uvals{ 1.0, 2.0, 3.0, 4.0, 5.0 };

    {
      const csmp::Index key_d( sg_iso.Database().StorageKey( "diffusivity"    ) );
      const csmp::Index key_u( sg_iso.Database().StorageKey( "fluid pressure" ) );
      Region<2>& domain( sg_iso.Region("Model") );
      uint32_t i{0U};
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        (*it)->Store( key_d, ScalarVariable( PLAIN, dvals[i++] ) );
      i = 0U;
      for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it )
        (*it)->Store( key_u, makeScalar( PLAIN, uvals[i++] ) );
    }

    // also set on the analytical model (sg_) for consistency
    setElementScalar( dvals, "diffusivity"    );
    setNodeVariable(  uvals, "fluid pressure" );

    // analytical operator — uses LINEAR_TRIANGLE model (sg_)
    Integral_dNT_rhsop_dN_dV<2U> analytical( sg_->Database(),
                                              "diffusivity",
                                              "fluid pressure" );

    // numerical operator — uses ISOPARAMETRIC_LINEAR_TRIANGLE model (sg_iso)
    // for linear simplex elements the numerical operator uses dN_AtBaryCenter
    // which gives the same result as the analytical formula
    NumIntegral_dNT_rhsop_dN_dV<2U> numerical( sg_iso.Database(),
                                                "diffusivity",
                                                "fluid pressure" );

    std::vector<double> rhs_anal( sg_->Region("Model").Nodes(),    0.0 );
    std::vector<double> rhs_num(  sg_iso.Region("Model").Nodes(),  0.0 );

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
      Region<2>& domain( sg_iso.Region("Model") );
      for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
        {
          numerical.GetOperands(        *(*it) );
          numerical.ComputeContribution(*(*it) );
          numerical.AssignToGlobal(     *(*it), rhs_num );
        }
    }

    if ( verbose_ )
      {
        cout << "\nmatchesNumericalTest: analytical vs numerical"
             << " on ISOPARAMETRIC_LINEAR_TRIANGLE mesh:\n";
        for ( size_t i{0U}; i < rhs_anal.size(); ++i )
          cout << "  analytical[" << i << "] = " << rhs_anal[i]
               << "  numerical["  << i << "] = " << rhs_num[i]  << "\n";
      }

    // for linear simplex elements both operators must agree exactly
    // any discrepancy indicates a bug in one of the two implementations
    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    assert( rhs_anal.size() == rhs_num.size() );
    for ( size_t i{0U}; i < rhs_anal.size(); ++i )
      _equal( rhs_anal[i], rhs_num[i], t );
 }


// =============================================================================
// Helper methods
// =============================================================================

void Integral_dNT_rhsop_dN_dV_Test::setNodeVariable( const vector<double>& val, const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Nodes() )
      throw csmp::Exception( ERROR,
          "Integral_dNT_rhsop_dN_dV_Test::setNodeVariable",
          var_name, "Variable size does not match number of nodes." );

    uint32_t i{0U};
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it )
      (*it)->Store( key, makeScalar( PLAIN, val[i++] ) );
 }


void Integral_dNT_rhsop_dN_dV_Test::setElementScalar(
    const vector<double>& val, const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR,
          "Integral_dNT_rhsop_dN_dV_Test::setElementScalar",
          var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }


void Integral_dNT_rhsop_dN_dV_Test::calculateRHS( std::vector<double>& rhs,
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

