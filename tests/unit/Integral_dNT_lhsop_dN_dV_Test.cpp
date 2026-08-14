// Integral_dNT_lhsop_dN_dV_Test.cpp
#include "Integral_dNT_lhsop_dN_dV_Test.h"

#include "Integral_dNT_op_dN_dV.h"

#include "Model.h"
#include "SparseMatrix.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

Integral_dNT_lhsop_dN_dV_Test::Integral_dNT_lhsop_dN_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nIntegral_dNT_lhsop_dN_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
 }

Integral_dNT_lhsop_dN_dV_Test::~Integral_dNT_lhsop_dN_dV_Test()
 { delete sg_; }

void Integral_dNT_lhsop_dN_dV_Test::run()
 {
    scalarTest();
    symmetryTest();
    rowSumTest();
 }

void Integral_dNT_lhsop_dN_dV_Test::scalarTest()
 {
    Integral_dNT_op_dN_dV<2U> stiffness( sg_->Database(),
                                          "diffusivity",
                                          "fluid pressure",
                                          "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

    if ( verbose_ )
      {
        cout << "\nIntegral_dNT_op_dN_dV scalarTest:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Identical reference values to NumIntegral_dNT_lhsop_dN_dV_Test::scalarTest
    // s = 0.5, E={1,2,3,4}
    // -----------------------------------------------------------------------
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

    _equal( sm.At(4,0), -5.0*s, t );
    _equal( sm.At(4,1), -3.0*s, t );
    _equal( sm.At(4,2), -5.0*s, t );
    _equal( sm.At(4,3), -7.0*s, t );

    const std::vector<std::pair<size_t,size_t>> zero_pairs{
        {0U,1U},{0U,2U},{0U,3U},{1U,2U},{1U,3U},{2U,3U}
    };
    for ( const auto& [i,j] : zero_pairs )
      if ( sm.HasEntry(i,j) ) _equal( sm.At(i,j), 0.0, t );
 }

void Integral_dNT_lhsop_dN_dV_Test::symmetryTest()
 {
    Integral_dNT_op_dN_dV<2U> stiffness( sg_->Database(),
                                          "diffusivity",
                                          "fluid pressure",
                                          "fluid pressure" );
    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

    const std::vector<std::pair<size_t,size_t>> pairs{
        {0U,4U},{1U,4U},{2U,4U},{3U,4U}
    };
    for ( const auto& [i,j] : pairs )
      _equal( sm.At(i,j), sm.At(j,i), tol_ );
 }

void Integral_dNT_lhsop_dN_dV_Test::rowSumTest()
 {
    Integral_dNT_op_dN_dV<2U> stiffness( sg_->Database(),
                                          "diffusivity",
                                          "fluid pressure",
                                          "fluid pressure" );
    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

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
          if ( sm.HasEntry(i,j) ) row_sum += sm.At(i,j);
        _equal( row_sum, 0.0, tol_ * sm.At(i,i) );
      }
 }

void Integral_dNT_lhsop_dN_dV_Test::setElementScalar( const vector<double>& val,
                                                    const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "Integral_dNT_lhsop_dN_dV_Test::setElementScalar",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }

void Integral_dNT_lhsop_dN_dV_Test::calculateGlobalMatrix( SparseMatrix& sm,
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

