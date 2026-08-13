// PropertyHandle_MathTest.cpp

#include "PropertyHandle_MathTest.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "Model.h"

using namespace std;

namespace csmp {

// ============================================================================
//  Constructor / destructor
// ============================================================================

PropertyHandle_MathTest::PropertyHandle_MathTest(
    double tolerance, bool verbose )
    : model_( new ANSYS_Model3D( "BoxHalfs3D",
                                 "PropertyHandle_Test-variables.txt" ) ),
      TOLERANCE( tolerance ),
      verbose_( verbose ),
      sc_elem( *model_, "math sc elem", SCALAR, ELEMENT ),
      vc_elem( *model_, "math vc elem", VECTOR, ELEMENT ),
      ts_elem( *model_, "math ts elem", TENSOR, ELEMENT ),
      ar_elem( *model_, "math ar elem", ARRAY,        ELEMENT, 4U ),
      fa_elem( *model_, "math fa elem", FLAGGEDARRAY, ELEMENT, 4U )
{
}

PropertyHandle_MathTest::~PropertyHandle_MathTest()
{
    delete model_;
}

// ============================================================================
//  run
// ============================================================================

void PropertyHandle_MathTest::run()
{
    if ( verbose_ )
    {
        cout << "\n=================================" << endl;
        cout << "Testing PropertyHandle math methods" << endl;
        cout << "================================="  << endl;
    }

    TestScalarMath();
    TestVectorMath();
    TestTensorMath();
    TestArrayMath();
    TestFlaggedArrayMath();
}

// ============================================================================
//  TestScalarMath
// ============================================================================

void PropertyHandle_MathTest::TestScalarMath()
{
    if ( verbose_ ) cout << "\n--- SCALAR math ---" << endl;

    ScalarVariable sc;
    const Index key = model_->Database().StorageKey( "math sc elem" );
    auto firstCell  = model_->Region("Model").CellsBegin();

    auto read = [&]() -> double
    {
        ( *firstCell )->Read( key, sc );
        return sc();
    };

    // Squared: 3^2 = 9
    sc_elem = 3.0;
    sc_elem.Squared();
    _equal( read(), 9.0, TOLERANCE );

    // Sqrt: sqrt(9) = 3
    sc_elem.Sqrt();
    _equal( read(), 3.0, TOLERANCE );

    // Ln: ln(3)
    sc_elem.Ln();
    _equal( read(), std::log(3.0), TOLERANCE );

    // Exp: exp(ln(3)) = 3
    sc_elem.Exp();
    _equal( read(), 3.0, TOLERANCE );

    // Log10: log10(3)
    sc_elem.Log10();
    _equal( read(), std::log10(3.0), TOLERANCE );

    // Pow: 2^4 = 16
    sc_elem = 2.0;
    sc_elem.Pow( 4.0 );
    _equal( read(), 16.0, TOLERANCE );

    // Abs: abs(-16) = 16
    sc_elem = -16.0;
    sc_elem.Abs();
    _equal( read(), 16.0, TOLERANCE );

    // ZapNAN: NaN -> 7
    sc_elem = numeric_limits<double>::quiet_NaN();
    sc_elem.ZapNAN( 7.0 );
    _equal( read(), 7.0, TOLERANCE );

    // Clamp: 7 -> [0,5] = 5
    sc_elem.Clamp( 0.0, 5.0 );
    _equal( read(), 5.0, TOLERANCE );

    // Clamp: -3 -> [0,5] = 0
    sc_elem = -3.0;
    sc_elem.Clamp( 0.0, 5.0 );
    _equal( read(), 0.0, TOLERANCE );

    // Sin: sin(30 deg) = 0.5
    sc_elem = 30.0;
    sc_elem.Sin();
    _equal( read(), 0.5, TOLERANCE );

    // Cos: cos(60 deg) = 0.5
    sc_elem = 60.0;
    sc_elem.Cos();
    _equal( read(), 0.5, TOLERANCE );

    // Tan: tan(45 deg) = 1
    sc_elem = 45.0;
    sc_elem.Tan();
    _equal( read(), 1.0, TOLERANCE );

    // Acos: acos(0.5) = 60 deg
    sc_elem = 0.5;
    sc_elem.Acos();
    _equal( read(), 60.0, TOLERANCE );

    // Asin: asin(0.5) = 30 deg
    sc_elem = 0.5;
    sc_elem.Asin();
    _equal( read(), 30.0, TOLERANCE );

    // Atan: atan(1) = 45 deg
    sc_elem = 1.0;
    sc_elem.Atan();
    _equal( read(), 45.0, TOLERANCE );

    if ( verbose_ ) cout << "  SCALAR math: passed" << endl;
}

// ============================================================================
//  TestVectorMath
// ============================================================================

void PropertyHandle_MathTest::TestVectorMath()
{
    if ( verbose_ ) cout << "\n--- VECTOR math ---" << endl;

    VectorVariable<3> vc;
    const Index key = model_->Database().StorageKey( "math vc elem" );
    auto firstCell  = model_->Region("Model").CellsBegin();

    // Helper: read all three components and verify they are equal.
    auto readAll = [&]( double expected )
    {
        ( *firstCell )->Read( key, vc );
        _equal( vc(0), expected, TOLERANCE );
        _equal( vc(1), expected, TOLERANCE );
        _equal( vc(2), expected, TOLERANCE );
    };

    // Squared: each component 3^2 = 9
    vc_elem = 3.0;
    vc_elem.Squared();
    readAll( 9.0 );

    // Sqrt: sqrt(9) = 3
    vc_elem.Sqrt();
    readAll( 3.0 );

    // Ln: ln(3)
    vc_elem.Ln();
    readAll( std::log(3.0) );

    // Exp: exp(ln(3)) = 3
    vc_elem.Exp();
    readAll( 3.0 );

    // Log10: log10(3)
    vc_elem.Log10();
    readAll( std::log10(3.0) );

    // Pow: 2^3 = 8 per component
    vc_elem = 2.0;
    vc_elem.Pow( 3.0 );
    readAll( 8.0 );

    // Abs: abs(-8) = 8 per component
    vc_elem = -8.0;
    vc_elem.Abs();
    readAll( 8.0 );

    // ZapNAN: NaN -> 5 per component
    vc_elem = numeric_limits<double>::quiet_NaN();
    vc_elem.ZapNAN( 5.0 );
    readAll( 5.0 );

    // Clamp: 5 -> [0,3] = 3 per component
    vc_elem.Clamp( 0.0, 3.0 );
    readAll( 3.0 );

    // Sin: sin(30 deg) = 0.5 per component
    vc_elem = 30.0;
    vc_elem.Sin();
    readAll( 0.5 );

    // Cos: cos(60 deg) = 0.5 per component
    vc_elem = 60.0;
    vc_elem.Cos();
    readAll( 0.5 );

    // Tan: tan(45 deg) = 1 per component
    vc_elem = 45.0;
    vc_elem.Tan();
    readAll( 1.0 );

    // Acos: acos(0.5) = 60 deg per component
    vc_elem = 0.5;
    vc_elem.Acos();
    readAll( 60.0 );

    // Asin: asin(0.5) = 30 deg per component
    vc_elem = 0.5;
    vc_elem.Asin();
    readAll( 30.0 );

    // Atan: atan(1) = 45 deg per component
    vc_elem = 1.0;
    vc_elem.Atan();
    readAll( 45.0 );

    if ( verbose_ ) cout << "  VECTOR math: passed" << endl;
}

// ============================================================================
//  TestTensorMath
// ============================================================================

void PropertyHandle_MathTest::TestTensorMath()
{
    if ( verbose_ ) cout << "\n--- TENSOR math ---" << endl;

    TensorVariable<3> ts;
    const Index key = model_->Database().StorageKey( "math ts elem" );
    auto firstCell  = model_->Region("Model").CellsBegin();

    auto readAll = [&]( double expected )
    {
        ( *firstCell )->Read( key, ts );
        for ( uint32_t i = 0; i < 3; ++i )
            for ( uint32_t j = 0; j < 3; ++j )
                _equal( ts(i,j), expected, TOLERANCE );
    };

    // ------------------------------------------------------------------
    //  Squared: component-wise T_ij^2.
    //  For all entries = 2: result = 4.
    // ------------------------------------------------------------------
    ts_elem = 2.0;
    ts_elem.Squared();
    readAll( 4.0 );

    // ------------------------------------------------------------------
    //  MatrixSquared: matrix product T*T.
    //  For uniform 3x3 with all entries = 2:
    //  result(i,j) = sum_k 2*2 = 3*4 = 12.
    // ------------------------------------------------------------------
    ts_elem = 2.0;
    ts_elem.MatrixSquared();
    readAll( 12.0 );

    // ------------------------------------------------------------------
    //  DoubleContraction: T:T = sum_ij T_ij^2.
    //  For uniform 3x3 with all entries = 2:
    //  T:T = 9 * 4 = 36.
    // ------------------------------------------------------------------
    {
        PropertyHandle<3> sc_result( *model_, "math ts dc result",
                                      SCALAR, ELEMENT );
        ts_elem = 2.0;
        sc_result.DoubleContraction( ts_elem );

        ScalarVariable sc;
        const Index skey =
            model_->Database().StorageKey( "math ts dc result" );
        ( *firstCell )->Read( skey, sc );
        _equal( sc(), 36.0, TOLERANCE );
    }

    // ------------------------------------------------------------------
    //  Sqrt: component-wise sqrt(4) = 2.
    // ------------------------------------------------------------------
    ts_elem = 4.0;
    ts_elem.Sqrt();
    readAll( 2.0 );

    // ------------------------------------------------------------------
    //  Ln: component-wise ln(2).
    // ------------------------------------------------------------------
    ts_elem.Ln();
    readAll( std::log(2.0) );

    // ------------------------------------------------------------------
    //  Exp: component-wise exp(ln(2)) = 2.
    // ------------------------------------------------------------------
    ts_elem.Exp();
    readAll( 2.0 );

    // ------------------------------------------------------------------
    //  Log10: component-wise log10(2).
    // ------------------------------------------------------------------
    ts_elem.Log10();
    readAll( std::log10(2.0) );

    // ------------------------------------------------------------------
    //  Pow: 2^3 = 8 component-wise.
    // ------------------------------------------------------------------
    ts_elem = 2.0;
    ts_elem.Pow( 3.0 );
    readAll( 8.0 );

    // ------------------------------------------------------------------
    //  Abs: abs(-8) = 8 component-wise.
    // ------------------------------------------------------------------
    ts_elem = -8.0;
    ts_elem.Abs();
    readAll( 8.0 );

    // ------------------------------------------------------------------
    //  ZapNAN: NaN -> 4 component-wise.
    // ------------------------------------------------------------------
    ts_elem = std::numeric_limits<double>::quiet_NaN();
    ts_elem.ZapNAN( 4.0 );
    readAll( 4.0 );

    // ------------------------------------------------------------------
    //  Clamp: 4 -> [0,2] = 2 component-wise.
    // ------------------------------------------------------------------
    ts_elem.Clamp( 0.0, 2.0 );
    readAll( 2.0 );

    // ------------------------------------------------------------------
    //  Sin: sin(30 deg) = 0.5 component-wise.
    // ------------------------------------------------------------------
    ts_elem = 30.0;
    ts_elem.Sin();
    readAll( 0.5 );

    // ------------------------------------------------------------------
    //  Cos: cos(60 deg) = 0.5 component-wise.
    // ------------------------------------------------------------------
    ts_elem = 60.0;
    ts_elem.Cos();
    readAll( 0.5 );

    // ------------------------------------------------------------------
    //  Tan: tan(45 deg) = 1 component-wise.
    // ------------------------------------------------------------------
    ts_elem = 45.0;
    ts_elem.Tan();
    readAll( 1.0 );

    // ------------------------------------------------------------------
    //  Acos: acos(0.5) = 60 deg component-wise.
    // ------------------------------------------------------------------
    ts_elem = 0.5;
    ts_elem.Acos();
    readAll( 60.0 );

    // ------------------------------------------------------------------
    //  Asin: asin(0.5) = 30 deg component-wise.
    // ------------------------------------------------------------------
    ts_elem = 0.5;
    ts_elem.Asin();
    readAll( 30.0 );

    // ------------------------------------------------------------------
    //  Atan: atan(1) = 45 deg component-wise.
    // ------------------------------------------------------------------
    ts_elem = 1.0;
    ts_elem.Atan();
    readAll( 45.0 );

    if ( verbose_ ) cout << "  TENSOR math: passed" << endl;
}

// ============================================================================
//  TestArrayMath
// ============================================================================

void PropertyHandle_MathTest::TestArrayMath()
{
    if ( verbose_ ) cout << "\n--- ARRAY math ---" << endl;

    ArrayVariable av;
    const Index key = model_->Database().StorageKey( "math ar elem" );
    auto firstCell  = model_->Region("Model").CellsBegin();

    // Helper: read all four elements and verify they are equal.
    auto readAll = [&]( double expected )
    {
        ( *firstCell )->Read( key, av );
        for ( uint32_t i = 0; i < 4; ++i )
            _equal( av[i], expected, TOLERANCE );
    };

    // Squared: 3^2 = 9 per element
    ar_elem = 3.0;
    ar_elem.Squared();
    readAll( 9.0 );

    // Sqrt: sqrt(9) = 3 per element
    ar_elem.Sqrt();
    readAll( 3.0 );

    // Ln: ln(3) per element
    ar_elem.Ln();
    readAll( std::log(3.0) );

    // Exp: exp(ln(3)) = 3 per element
    ar_elem.Exp();
    readAll( 3.0 );

    // Log10: log10(3) per element
    ar_elem.Log10();
    readAll( std::log10(3.0) );

    // Pow: 2^4 = 16 per element
    ar_elem = 2.0;
    ar_elem.Pow( 4.0 );
    readAll( 16.0 );

    // Abs: abs(-16) = 16 per element
    ar_elem = -16.0;
    ar_elem.Abs();
    readAll( 16.0 );

    // ZapNAN: NaN -> 7 per element
    ar_elem = numeric_limits<double>::quiet_NaN();
    ar_elem.ZapNAN( 7.0 );
    readAll( 7.0 );

    // Clamp: 7 -> [0,5] = 5 per element
    ar_elem.Clamp( 0.0, 5.0 );
    readAll( 5.0 );

    // Sin: sin(30 deg) = 0.5 per element
    ar_elem = 30.0;
    ar_elem.Sin();
    readAll( 0.5 );

    // Cos: cos(60 deg) = 0.5 per element
    ar_elem = 60.0;
    ar_elem.Cos();
    readAll( 0.5 );

    // Tan: tan(45 deg) = 1 per element
    ar_elem = 45.0;
    ar_elem.Tan();
    readAll( 1.0 );

    // Acos: acos(0.5) = 60 deg per element
    ar_elem = 0.5;
    ar_elem.Acos();
    readAll( 60.0 );

    // Asin: asin(0.5) = 30 deg per element
    ar_elem = 0.5;
    ar_elem.Asin();
    readAll( 30.0 );

    // Atan: atan(1) = 45 deg per element
    ar_elem = 1.0;
    ar_elem.Atan();
    readAll( 45.0 );

    if ( verbose_ ) cout << "  ARRAY math: passed" << endl;
}

// ============================================================================
//  TestFlaggedArrayMath
// ============================================================================

void PropertyHandle_MathTest::TestFlaggedArrayMath()
{
    if ( verbose_ ) cout << "\n--- FLAGGEDARRAY math ---" << endl;

    // ------------------------------------------------------------------
    //  Setup: 4-element FLAGGEDARRAY, element 0 = DIRICH (protected),
    //  elements 1-3 = ANY (modifiable).
    //  All elements initialised to the same value so we can verify
    //  that element 0 is unchanged after every operation.
    // ------------------------------------------------------------------
    FlaggedArrayVariable fav;
    const Index key = model_->Database().StorageKey( "math fa elem" );
    auto firstCell  = model_->Region("Model").CellsBegin();

    auto setDirich = [&]()
    {
        fa_elem.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
        {
            f.Flag( 0, DIRICH );
        });
    };

    // Helper: read element 0 (DIRICH, must be unchanged) and
    // element 1 (ANY, must reflect the operation).
    auto readCheck = [&]( double protected_val, double modified_val )
    {
        ( *firstCell )->Read( key, fav );
        _equal( fav[0], protected_val, TOLERANCE );  // DIRICH: unchanged
        _equal( fav[1], modified_val,  TOLERANCE );  // ANY: modified
        _equal( fav[2], modified_val,  TOLERANCE );  // ANY: modified
        _equal( fav[3], modified_val,  TOLERANCE );  // ANY: modified
    };

    // Helper: reset to known value and set element 0 to DIRICH.
    auto reset = [&]( double val )
    {
        // First reset ALL flags to ANY so operator=(double) sets all elements.
        fa_elem.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
        {
            for ( uint32_t i = 0; i < f.Size(); ++i )
                f.Flag( i, ANY );
        });
        fa_elem = val;
        // Now set element 0 to DIRICH.
        fa_elem.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
        {
            f.Flag( 0, DIRICH );
        });
    };

    // Squared: 3^2 = 9 for ANY; element 0 (DIRICH) = 3 unchanged.
    reset( 3.0 );
    fa_elem.Squared();
    readCheck( 3.0, 9.0 );

    // Sqrt: sqrt(4) = 2 for ANY; element 0 (DIRICH) = 4 unchanged.
    reset( 4.0 );
    fa_elem.Sqrt();
    readCheck( 4.0, 2.0 );

    // Ln: ln(e) = 1 for ANY; element 0 (DIRICH) = e unchanged.
    reset( std::exp(1.0) );
    fa_elem.Ln();
    readCheck( std::exp(1.0), 1.0 );

    // Exp: exp(1) = e for ANY; element 0 (DIRICH) = 1 unchanged.
    reset( 1.0 );
    fa_elem.Exp();
    readCheck( 1.0, std::exp(1.0) );

    // Log10: log10(10) = 1 for ANY; element 0 (DIRICH) = 10 unchanged.
    reset( 10.0 );
    fa_elem.Log10();
    readCheck( 10.0, 1.0 );

    // Pow: 2^3 = 8 for ANY; element 0 (DIRICH) = 2 unchanged.
    reset( 2.0 );
    fa_elem.Pow( 3.0 );
    readCheck( 2.0, 8.0 );

    // Abs: abs(-4) = 4 for ANY; element 0 (DIRICH) = -4 unchanged.
    reset( -4.0 );
    fa_elem.Abs();
    readCheck( -4.0, 4.0 );

    // ZapNAN: NaN -> 6 for ANY; element 0 (DIRICH) = NaN unchanged.
    fa_elem.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
    {
        for ( uint32_t i = 0; i < f.Size(); ++i )
            f.Flag( i, ANY );
    });
    fa_elem = numeric_limits<double>::quiet_NaN();
    fa_elem.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
    {
        f.Flag( 0, DIRICH );
    });
    fa_elem.ZapNAN( 6.0 );
    ( *firstCell )->Read( key, fav );
    _equal( std::isnan( fav[0] ), true, 0.5 );
    _equal( fav[1], 6.0, TOLERANCE );
    _equal( fav[2], 6.0, TOLERANCE );
    _equal( fav[3], 6.0, TOLERANCE );

    // Clamp: clamp(8,[0,5])=5 for ANY; element 0 (DIRICH) = 8 unchanged.
    reset( 8.0 );
    fa_elem.Clamp( 0.0, 5.0 );
    readCheck( 8.0, 5.0 );

    // Sin: sin(30 deg) = 0.5 for ANY; element 0 (DIRICH) = 30 unchanged.
    reset( 30.0 );
    fa_elem.Sin();
    readCheck( 30.0, 0.5 );

    // Cos: cos(60 deg) = 0.5 for ANY; element 0 (DIRICH) = 60 unchanged.
    reset( 60.0 );
    fa_elem.Cos();
    readCheck( 60.0, 0.5 );

    // Tan: tan(45 deg) = 1 for ANY; element 0 (DIRICH) = 45 unchanged.
    reset( 45.0 );
    fa_elem.Tan();
    readCheck( 45.0, 1.0 );

    // Acos: acos(0.5) = 60 deg for ANY; element 0 (DIRICH) = 0.5 unchanged.
    reset( 0.5 );
    fa_elem.Acos();
    readCheck( 0.5, 60.0 );

    // Asin: asin(0.5) = 30 deg for ANY; element 0 (DIRICH) = 0.5 unchanged.
    reset( 0.5 );
    fa_elem.Asin();
    readCheck( 0.5, 30.0 );

    // Atan: atan(1) = 45 deg for ANY; element 0 (DIRICH) = 1 unchanged.
    reset( 1.0 );
    fa_elem.Atan();
    readCheck( 1.0, 45.0 );

    if ( verbose_ ) cout << "  FLAGGEDARRAY math: passed" << endl;
}

} // namespace csmp

