// Integral_NT_lhsop_N_dV_Test.cpp
#include "Integral_NT_lhsop_N_dV_Test.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

Integral_NT_lhsop_N_dV_Test::Integral_NT_lhsop_N_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nIntegral_NT_lhsop_N_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 0.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
 }

Integral_NT_lhsop_N_dV_Test::~Integral_NT_lhsop_N_dV_Test()
 { delete sg_; }
 

void Integral_NT_lhsop_N_dV_Test::run()
 {
    consistentTest();
    lumpedTest();
    rowSumTest();
    symmetryTest();
 }
 

void Integral_NT_lhsop_N_dV_Test::consistentTest()
 {
    Integral_NT_lhsop_N_dV<2U> mass( sg_->Database(),
                                      "diffusivity",
                                      "fluid pressure",
                                      "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, mass );

    if ( verbose_ )
      {
        cout << "\nIntegral_NT_lhsop_N_dV consistentTest:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Identical reference values to NumIntegral_NT_lhsop_N_dV_Test::consistentTest
    // A/6=10416.67, A/12=5208.33, E={1,2,3,4}
    // -----------------------------------------------------------------------
    constexpr double A{   62500.0 };
    constexpr double A6{  A / 6.0  };
    constexpr double A12{ A / 12.0 };
    const double t{ tol_ * A12 };

    _equal( sm.At(0,0),  5.0*A6,  t );
    _equal( sm.At(1,1),  3.0*A6,  t );
    _equal( sm.At(2,2),  5.0*A6,  t );
    _equal( sm.At(3,3),  7.0*A6,  t );
    _equal( sm.At(4,4), 10.0*A6,  t );

    _equal( sm.At(0,1),  1.0*A12, t );
    _equal( sm.At(0,3),  4.0*A12, t );
    _equal( sm.At(0,4),  5.0*A12, t );

    _equal( sm.At(1,0),  1.0*A12, t );
    _equal( sm.At(1,2),  2.0*A12, t );
    _equal( sm.At(1,4),  3.0*A12, t );

    _equal( sm.At(2,1),  2.0*A12, t );
    _equal( sm.At(2,3),  3.0*A12, t );
    _equal( sm.At(2,4),  5.0*A12, t );

    _equal( sm.At(3,0),  4.0*A12, t );
    _equal( sm.At(3,2),  3.0*A12, t );
    _equal( sm.At(3,4),  7.0*A12, t );

    _equal( sm.At(4,0),  5.0*A12, t );
    _equal( sm.At(4,1),  3.0*A12, t );
    _equal( sm.At(4,2),  5.0*A12, t );
    _equal( sm.At(4,3),  7.0*A12, t );
 }



void Integral_NT_lhsop_N_dV_Test::lumpedTest()
 {
    Integral_NT_lhsop_N_dV<2U> mass( sg_->Database(),
                                      "diffusivity",
                                      "fluid pressure",
                                      "fluid pressure" );
    mass.LumpedFormulation( true );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, mass );

    if ( verbose_ )
      {
        cout << "\nIntegral_NT_lhsop_N_dV lumpedTest:\n";
        sm.Out();
      }

    constexpr double A{   62500.0 };
    constexpr double A12{ A / 12.0 };
    const double t{ tol_ * A12 };

    _equal( sm.At(0,0), 20.0*A12, t );
    _equal( sm.At(1,1), 12.0*A12, t );
    _equal( sm.At(2,2), 20.0*A12, t );
    _equal( sm.At(3,3), 28.0*A12, t );
    _equal( sm.At(4,4), 40.0*A12, t );

    // lumped matrix must have exactly 5 non-zero entries
    _equal( static_cast<double>( sm.Entries() ), 5.0, 0.0 );
 }



void Integral_NT_lhsop_N_dV_Test::rowSumTest()
 {
    Integral_NT_lhsop_N_dV<2U> lumped( sg_->Database(),
                                        "diffusivity",
                                        "fluid pressure",
                                        "fluid pressure" );
    lumped.LumpedFormulation( true );

    Integral_NT_lhsop_N_dV<2U> consistent( sg_->Database(),
                                            "diffusivity",
                                            "fluid pressure",
                                            "fluid pressure" );
    consistent.LumpedFormulation( false );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm_lumped( n ), sm_consistent( n );
    calculateGlobalMatrix( sm_lumped,     lumped     );
    calculateGlobalMatrix( sm_consistent, consistent );

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
          if ( sm_consistent.HasEntry(i,j) ) row_sum += sm_consistent.At(i,j);
        _equal( sm_lumped.At(i,i), row_sum, tol_ * row_sum );
      }
 }



void Integral_NT_lhsop_N_dV_Test::symmetryTest()
 {
    Integral_NT_lhsop_N_dV<2U> mass( sg_->Database(),
                                      "diffusivity",
                                      "fluid pressure",
                                      "fluid pressure" );
    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, mass );

    const std::vector<std::pair<size_t,size_t>> pairs{
        {0U,1U},{0U,3U},{0U,4U},
        {1U,2U},{1U,4U},
        {2U,3U},{2U,4U},
        {3U,4U}
    };
    for ( const auto& [i,j] : pairs )
      _equal( sm.At(i,j), sm.At(j,i), tol_ );
 }



void Integral_NT_lhsop_N_dV_Test::setElementScalar( const vector<double>& val,
                                                     const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "Integral_NT_lhsop_N_dV_Test::setElementScalar",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }
 
 

void Integral_NT_lhsop_N_dV_Test::calculateGlobalMatrix( SparseMatrix& sm,
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

