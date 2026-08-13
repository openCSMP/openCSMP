//
//  NumIntegral_NT_rhsop_N_dV_Test.cpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 6/7/2026.
//  Copyright © 2026 Stephan Matthai. All rights reserved.
//

// NumIntegral_NT_rhsop_N_dV_Test.cpp
#include "NumIntegral_NT_rhsop_N_dV_Test.h"
#include "TRIANGLE_Interface.h"
#include "meshManagementUtilities.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

NumIntegral_NT_rhsop_N_dV_Test::NumIntegral_NT_rhsop_N_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nNumIntegral_NT_rhsop_N_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    if ( verbose_ )
      printNodeCoordinates<2>( sg_->Region("Model").NodesBegin(),
                               sg_->Region("Model").NodesEnd() );

    sg_->InputPropertyValue( "fluid pressure",            makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "diffusivity",               makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "fluid volume source",       makeScalar( PLAIN, 1.0 ) );
 }


NumIntegral_NT_rhsop_N_dV_Test::~NumIntegral_NT_rhsop_N_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void NumIntegral_NT_rhsop_N_dV_Test::run()
 {
    consistentTest();
    lumpedTest();
    lumpedEqualsConsistentRowSumTest();
 }


// =============================================================================
// consistentTest — verifies RHS vector against analytical row sums of mass matrix
// =============================================================================

void NumIntegral_NT_rhsop_N_dV_Test::consistentTest()
 {
    NumIntegral_NT_rhsop_N_dV<2U> integral( sg_->Database(), "diffusivity", "fluid pressure" );
    integral.LumpedFormulation( false );

    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    std::vector<double> rhs;
    calculateRHS( rhs, integral );

    if ( verbose_ )
      {
        cout << "\nconsistentTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // The consistent RHS vector is the row sum of the consistent mass matrix M:
    //   RHS[j] = sum_k M_jk
    //
    // From NumIntegral_NT_lhsop_N_dV_Test::consistentTest we know:
    //   M_jj^(e) = E * A/6
    //   M_jk^(e) = E * A/12   j != k
    //
    // Mesh: A=62500, A/6=10416.67, A/12=5208.33
    // Connectivity: e0={0,1,4} E=1, e1={1,2,4} E=2, e2={2,3,4} E=3, e3={3,0,4} E=4
    //
    // RHS[0] = M_00 + M_01 + M_03 + M_04
    //        = (A/6)*(1+4) + (A/12)*1 + (A/12)*4 + (A/12)*(1+4)
    //        = 5*A/6 + A/12 + 4*A/12 + 5*A/12
    //        = 10*A/12 + 10*A/12 = 20*A/12
    //
    // RHS[1] = M_11 + M_10 + M_12 + M_14
    //        = (A/6)*(1+2) + (A/12)*1 + (A/12)*2 + (A/12)*(1+2)
    //        = 6*A/12 + A/12 + 2*A/12 + 3*A/12 = 12*A/12
    //
    // RHS[2] = M_22 + M_21 + M_23 + M_24
    //        = (A/6)*(2+3) + (A/12)*2 + (A/12)*3 + (A/12)*(2+3)
    //        = 10*A/12 + 2*A/12 + 3*A/12 + 5*A/12 = 20*A/12
    //
    // RHS[3] = M_33 + M_30 + M_32 + M_34
    //        = (A/6)*(3+4) + (A/12)*4 + (A/12)*3 + (A/12)*(3+4)
    //        = 14*A/12 + 4*A/12 + 3*A/12 + 7*A/12 = 28*A/12
    //
    // RHS[4] = M_44 + M_40 + M_41 + M_42 + M_43
    //        = (A/6)*(1+2+3+4) + (A/12)*(1+4) + (A/12)*(1+2) + (A/12)*(2+3) + (A/12)*(3+4)
    //        = 20*A/12 + 5*A/12 + 3*A/12 + 5*A/12 + 7*A/12 = 40*A/12
    // -----------------------------------------------------------------------

    constexpr double A{   62500.0 };
    constexpr double A12{ A / 12.0 };
    const double t{ tol_ * A12 };

    _equal( rhs[0], 20.0*A12, t );
    _equal( rhs[1], 12.0*A12, t );
    _equal( rhs[2], 20.0*A12, t );
    _equal( rhs[3], 28.0*A12, t );
    _equal( rhs[4], 40.0*A12, t );
 }


// =============================================================================
// lumpedTest — verifies lumped RHS vector against analytical values
// =============================================================================

void NumIntegral_NT_rhsop_N_dV_Test::lumpedTest()
 {
    NumIntegral_NT_rhsop_N_dV<2U> lumped( sg_->Database(), "diffusivity", "fluid pressure" );
    lumped.LumpedFormulation( true );

    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    std::vector<double> rhs;
    calculateRHS( rhs, lumped );

    if ( verbose_ )
      {
        cout << "\nlumpedTest: RHS vector after assembly:\n";
        for ( size_t i{0U}; i < rhs.size(); ++i )
          cout << "  RHS[" << i << "] = " << rhs[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values for lumped formulation
    // -----------------------------------------------------------------------
    // Lumped: RHS[j] = sum_e E_e * (A/n_nodes) = E_e * A/3 per element
    //
    // Node 0: in e0(E=1) and e3(E=4)
    //   RHS[0] = (A/3)*(1+4) = 5*A/3 = 20*A/12
    //
    // Node 1: in e0(E=1) and e1(E=2)
    //   RHS[1] = (A/3)*(1+2) = 3*A/3 = 12*A/12  (wait: 3*A/3 = A = 12*A/12 ✓)
    //
    // Node 2: in e1(E=2) and e2(E=3)
    //   RHS[2] = (A/3)*(2+3) = 5*A/3 = 20*A/12
    //
    // Node 3: in e2(E=3) and e3(E=4)
    //   RHS[3] = (A/3)*(3+4) = 7*A/3 = 28*A/12
    //
    // Node 4: in all 4 elements
    //   RHS[4] = (A/3)*(1+2+3+4) = 10*A/3 = 40*A/12
    //
    // Note: lumped RHS equals consistent RHS for constant element operand —
    //       this is verified explicitly in lumpedEqualsConsistentRowSumTest
    // -----------------------------------------------------------------------

    constexpr double A{   62500.0 };
    constexpr double A12{ A / 12.0 };
    const double t{ tol_ * A12 };

    _equal( rhs[0], 20.0*A12, t );
    _equal( rhs[1], 12.0*A12, t );
    _equal( rhs[2], 20.0*A12, t );
    _equal( rhs[3], 28.0*A12, t );
    _equal( rhs[4], 40.0*A12, t );
 }


// =============================================================================
// lumpedEqualsConsistentRowSumTest — for constant element operand, lumped RHS
//                                    must equal consistent RHS entry by entry
// =============================================================================

void NumIntegral_NT_rhsop_N_dV_Test::lumpedEqualsConsistentRowSumTest()
 {
    NumIntegral_NT_rhsop_N_dV<2U> lumped( sg_->Database(), "diffusivity", "fluid pressure" );
    lumped.LumpedFormulation( true );

    NumIntegral_NT_rhsop_N_dV<2U> consistent( sg_->Database(), "diffusivity", "fluid pressure" );
    consistent.LumpedFormulation( false );

    // use non-uniform diffusivity to make the test more discriminating
    setElementVariable( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    std::vector<double> rhs_lumped, rhs_consistent;
    calculateRHS( rhs_lumped,     lumped     );
    calculateRHS( rhs_consistent, consistent );

    if ( verbose_ )
      {
        cout << "\nlumpedEqualsConsistentRowSumTest: RHS vectors:\n";
        for ( size_t i{0U}; i < rhs_lumped.size(); ++i )
          cout << "  lumped[" << i << "] = " << rhs_lumped[i]
               << "  consistent[" << i << "] = " << rhs_consistent[i] << "\n";
      }

    // for a piecewise-constant element operand the lumped and consistent
    // formulations must produce identical RHS vectors — this follows from
    // the fact that the consistent RHS is the row sum of the mass matrix
    // and the lumped RHS is the diagonal of the lumped mass matrix which
    // equals the row sum of the consistent mass matrix
    assert( rhs_lumped.size() == rhs_consistent.size() );
    for ( size_t i{0U}; i < rhs_lumped.size(); ++i )
      _equal( rhs_lumped[i], rhs_consistent[i], tol_ * rhs_consistent[i] );
 }


// =============================================================================
// Helper methods
// =============================================================================

void NumIntegral_NT_rhsop_N_dV_Test::setElementVariable( const vector<double>& var,
                                                         const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( var.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "NumIntegral_NT_rhsop_N_dV_Test::setElementVariable",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, var[i++] ) );
 }


void NumIntegral_NT_rhsop_N_dV_Test::calculateRHS( vector<double>& rhs,
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
