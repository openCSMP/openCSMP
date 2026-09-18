#include "TensorVar_Test.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "CSMP_definitions.h"

using std::vector;
using std::cout;
using std::endl;

namespace csmp {

TensorVariable_Test::TensorVariable_Test()
{
}

TensorVariable_Test::~TensorVariable_Test()
{
}

void TensorVariable_Test::run()
{
    Assignment_Operator();
    Addition_Operator();
    Subtraction_Operator();
    Multiplication_Operator();
    Division_Operator();
    Addition_Assignment_Operator();
    Subtraction_Assignment_Operator();
    Multiplication_Assignment_Operator();
    Division_Assignment_Operator();
    Equality_Operator();
    LessThan_Operator();          // was missing from run()
    MinElement_Function();
    MaxElement_Function();
    IsWithinRange_Function();
    Adjoint_Function();
    Identity_Function();
    Transposed_Function();
    Determinant_Function();
    Inverse_Function();
    DiagonalValues_Function();
    EigenValues_Function();
    Trace_Function();
    AssignToRow_Function();
    AssignToColumn_Function();
    Row_Function();
    Column_Function();
    HadamardSquared_Function();   // new
    MatrixSquared_Function();     // new
    DoubleContraction_Function(); // new
    TestEigenMethods();
}

// ============================================================================

void TensorVariable_Test::Assignment_Operator()
{
    ScalarVariable     sc( PLAIN, 5.0 );
    VectorVariable<3U> v( PLAIN, ROBIN, ANY, 10.0, 10.0, 10.0 );
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0,
                           4.0, 5.0, 6.0,
                           7.0, 8.0, 9.0 );
    TensorVariable<3U> t2;

    // VectorVariable assignment — diagonal set, off-diagonal zeroed.
    t2 = v;
    _equal( t2(0,0), 10.0, fTolerance ); _equal( t2(0,1), 0.0, fTolerance );
    _equal( t2(0,2),  0.0, fTolerance ); _equal( t2(1,0), 0.0, fTolerance );
    _equal( t2(1,1), 10.0, fTolerance ); _equal( t2(1,2), 0.0, fTolerance );
    _equal( t2(2,0),  0.0, fTolerance ); _equal( t2(2,1), 0.0, fTolerance );
    _equal( t2(2,2), 10.0, fTolerance );
    _test( t2.Flag(0) == PLAIN );
    _test( t2.Flag(1) == ROBIN );
    _test( t2.Flag(2) == ANY   );

    // TensorVariable assignment — all elements and flags copied.
    t2 = t1;
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(0,2), 3.0, fTolerance ); _equal( t2(1,0), 4.0, fTolerance );
    _equal( t2(1,1), 5.0, fTolerance ); _equal( t2(1,2), 6.0, fTolerance );
    _equal( t2(2,0), 7.0, fTolerance ); _equal( t2(2,1), 8.0, fTolerance );
    _equal( t2(2,2), 9.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );
    _test( t2.Flag(2) == NEUMANN   );

    // double assignment — all elements set, flags unchanged.
    t2 = 10.0;
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( t2(i,j), 10.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );
    _test( t2.Flag(2) == NEUMANN   );

    // ScalarVariable assignment — all elements and all flags set.
    t2 = sc;
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( t2(i,j), 5.0, fTolerance );
    _test( t2.Flag(0) == PLAIN );
    _test( t2.Flag(1) == PLAIN );
    _test( t2.Flag(2) == PLAIN );
}

// ============================================================================

void TensorVariable_Test::Addition_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2( PLAIN, PLAIN, PLAIN,
                           9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 );
    TensorVariable<3U> t3;

    // Tensor + tensor — element-by-element, flags from lhs.
    t3 = t1 + t2;
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( t3(i,j), 10.0, fTolerance );
    _test( t3.Flag(0) == INIT_GUESS );
    _test( t3.Flag(1) == PLAIN     );
    _test( t3.Flag(2) == NEUMANN   );

    // Tensor + double.
    t3 = t1 + 10.0;
    _equal( t3(0,0), 11.0, fTolerance );
    _equal( t3(1,1), 15.0, fTolerance );
    _equal( t3(2,2), 19.0, fTolerance );
    _test( t3.Flag(0) == INIT_GUESS );
    _test( t3.Flag(1) == PLAIN     );
    _test( t3.Flag(2) == NEUMANN   );
}

// ============================================================================

void TensorVariable_Test::Subtraction_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2( PLAIN, PLAIN, PLAIN,
                           9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 );
    TensorVariable<3U> t3;

    t3 = t1 - t2;
    _equal( t3(0,0), -8.0, fTolerance ); _equal( t3(0,1), -6.0, fTolerance );
    _equal( t3(0,2), -4.0, fTolerance ); _equal( t3(1,0), -2.0, fTolerance );
    _equal( t3(1,1),  0.0, fTolerance ); _equal( t3(1,2),  2.0, fTolerance );
    _equal( t3(2,0),  4.0, fTolerance ); _equal( t3(2,1),  6.0, fTolerance );
    _equal( t3(2,2),  8.0, fTolerance );
    _test( t3.Flag(0) == INIT_GUESS );
    _test( t3.Flag(1) == PLAIN     );
    _test( t3.Flag(2) == NEUMANN   );

    // Tensor - double.
    t3 = t1 - 10.0;
    _equal( t3(0,0), -9.0, fTolerance );
    _equal( t3(1,1), -5.0, fTolerance );
    _equal( t3(2,2), -1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Multiplication_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2;
    TensorVariable<3U> t3( PLAIN, PLAIN, PLAIN,
                           9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 );
    VectorVariable<3U> v1( PLAIN, PERIODIC, ANY, 1.0, 1.0, 1.0 );
    VectorVariable<3U> v2;

    // Vector-matrix multiplication: v^T * A.
    v2 = v1 * t1;
    _equal( v2(0), 12.0, fTolerance );
    _equal( v2(1), 15.0, fTolerance );
    _equal( v2(2), 18.0, fTolerance );

    // Matrix-vector multiplication: A * v.
    v2 = t1 * v1;
    _equal( v2(0),  6.0, fTolerance );
    _equal( v2(1), 15.0, fTolerance );
    _equal( v2(2), 24.0, fTolerance );
    _test( v2.Flag(0) == PLAIN    );
    _test( v2.Flag(1) == PERIODIC );
    _test( v2.Flag(2) == ANY      );

    // Matrix-matrix multiplication.
    t2 = t1 * t3;
    _equal( t2(0,0),  30.0, fTolerance ); _equal( t2(0,1),  24.0, fTolerance );
    _equal( t2(0,2),  18.0, fTolerance ); _equal( t2(1,0),  84.0, fTolerance );
    _equal( t2(1,1),  69.0, fTolerance ); _equal( t2(1,2),  54.0, fTolerance );
    _equal( t2(2,0), 138.0, fTolerance ); _equal( t2(2,1), 114.0, fTolerance );
    _equal( t2(2,2),  90.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );
    _test( t2.Flag(2) == NEUMANN   );

    // Tensor * double.
    t2 = t1 * 10.0;
    _equal( t2(0,0), 10.0, fTolerance );
    _equal( t2(1,1), 50.0, fTolerance );
    _equal( t2(2,2), 90.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Division_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t3( PERIODIC, PLAIN, PLAIN,
                           10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
                           70.0, 80.0, 90.0 );
    TensorVariable<3U> t2;

    // Element-by-element division.
    t2 = t3 / t1;
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( t2(i,j), 10.0, fTolerance );
    _test( t2.Flag(0) == PERIODIC );
    _test( t2.Flag(1) == PLAIN    );
    _test( t2.Flag(2) == PLAIN    );

    // Tensor / double.
    t2 = t3 / 10.0;
    _equal( t2(0,0),  1.0, fTolerance );
    _equal( t2(1,1),  5.0, fTolerance );
    _equal( t2(2,2),  9.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Addition_Assignment_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2( PLAIN, PLAIN, PLAIN,
                           9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 );
    TensorVariable<3U> t3;
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor += tensor.
    t3 = t1; t3 += t2;
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( t3(i,j), 10.0, fTolerance );

    // Tensor += scalar.
    t3 = t1; t3 += sc;
    _equal( t3(0,0),  6.0, fTolerance ); _equal( t3(0,1),  7.0, fTolerance );
    _equal( t3(0,2),  8.0, fTolerance ); _equal( t3(1,0),  9.0, fTolerance );
    _equal( t3(1,1), 10.0, fTolerance ); _equal( t3(1,2), 11.0, fTolerance );
    _equal( t3(2,0), 12.0, fTolerance ); _equal( t3(2,1), 13.0, fTolerance );
    _equal( t3(2,2), 14.0, fTolerance );

    // Tensor += double.
    t3 = t1; t3 += 10.0;
    _equal( t3(0,0), 11.0, fTolerance );
    _equal( t3(1,1), 15.0, fTolerance );
    _equal( t3(2,2), 19.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Subtraction_Assignment_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2( PLAIN, PLAIN, PLAIN,
                           9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 );
    TensorVariable<3U> t3;
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor -= tensor.
    t3 = t1; t3 -= t2;
    _equal( t3(0,0), -8.0, fTolerance ); _equal( t3(1,1), 0.0, fTolerance );
    _equal( t3(2,2),  8.0, fTolerance );

    // Tensor -= scalar.
    t3 = t1; t3 -= sc;
    _equal( t3(0,0), -4.0, fTolerance ); _equal( t3(1,1), 0.0, fTolerance );
    _equal( t3(2,2),  4.0, fTolerance );

    // Tensor -= double.
    t3 = t1; t3 -= 10.0;
    _equal( t3(0,0), -9.0, fTolerance );
    _equal( t3(1,1), -5.0, fTolerance );
    _equal( t3(2,2), -1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Multiplication_Assignment_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2;
    TensorVariable<3U> t3( PLAIN, PLAIN, PLAIN,
                           9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 );
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor *= tensor (matrix product).
    t2 = t1; t2 *= t3;
    _equal( t2(0,0),  30.0, fTolerance ); _equal( t2(1,1),  69.0, fTolerance );
    _equal( t2(2,2),  90.0, fTolerance );

    // Tensor *= scalar.
    t2 = t1; t2 *= sc;
    _equal( t2(0,0),  5.0, fTolerance ); _equal( t2(1,1), 25.0, fTolerance );
    _equal( t2(2,2), 45.0, fTolerance );

    // Tensor *= double.
    t2 = t1; t2 *= 10.0;
    _equal( t2(0,0), 10.0, fTolerance );
    _equal( t2(1,1), 50.0, fTolerance );
    _equal( t2(2,2), 90.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Division_Assignment_Operator()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
    TensorVariable<3U> t2( PLAIN, PLAIN, PLAIN,
                           10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
                           70.0, 80.0, 90.0 );
    TensorVariable<3U> t3;
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor /= tensor (element-by-element).
    t3 = t1; t3 /= t2;
    _equal( t3(0,0), 0.1, fTolerance ); _equal( t3(1,1), 0.1, fTolerance );
    _equal( t3(2,2), 0.1, fTolerance );

    // Tensor /= scalar.
    t3 = t1; t3 /= sc;
    _equal( t3(0,0), 0.2, fTolerance ); _equal( t3(1,1), 1.0, fTolerance );
    _equal( t3(2,2), 1.8, fTolerance );

    // Tensor /= double.
    t3 = t1; t3 /= 10.0;
    _equal( t3(0,0), 0.1, fTolerance );
    _equal( t3(1,1), 0.5, fTolerance );
    _equal( t3(2,2), 0.9, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Equality_Operator()
{
    TensorVariable<3U> t1( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> t2( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> t3( PLAIN, 9.0,8.0,7.0, 6.0,5.0,4.0, 3.0,2.0,1.0 );
    TensorVariable<3U> t4( ROBIN, 9.0,8.0,7.0, 6.0,5.0,4.0, 3.0,2.0,1.0 );

    _test(  t1 == t2 );
    _test(  t1 != t3 );
    _test(  t4 != t3 );  // different flags
    _test( !( t1 != t2 ) );
    _test( !( t1 == t3 ) );
}

// ============================================================================

void TensorVariable_Test::LessThan_Operator()
{
    // operator< compares pointer addresses — satisfies strict weak ordering
    // for STL associative containers. Does not compare values.
    TensorVariable<3U> t1( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> t2( PLAIN, 9.0,8.0,7.0, 6.0,5.0,4.0, 3.0,2.0,1.0 );

    // Exactly one of (t1 < t2) or (t2 < t1) must be true.
    const bool lt = ( t1 < t2 );
    const bool gt = ( t2 < t1 );
    _test(  lt != gt );       // exactly one is true
    _test( !( t1 < t1 ) );   // irreflexive
    _test( !( t2 < t2 ) );   // irreflexive
}

// ============================================================================

void TensorVariable_Test::MinElement_Function()
{
    TensorVariable<3U> t( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    _equal( t.MinElement(), 1.0, fTolerance );

    // Negative values.
    TensorVariable<3U> t2( PLAIN, -5.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    _equal( t2.MinElement(), -5.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::MaxElement_Function()
{
    TensorVariable<3U> t( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    _equal( t.MaxElement(), 9.0, fTolerance );

    // Negative values.
    TensorVariable<3U> t2( PLAIN, -9.0,-8.0,-7.0, -6.0,-5.0,-4.0,
                                   -3.0,-2.0,-1.0 );
    _equal( t2.MaxElement(), -1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::IsWithinRange_Function()
{
    TensorVariable<3U> t( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );

    // IsWithinRange checks only the diagonal (1, 5, 9).
    _test( t.IsWithinRange( -5.0,  5.0 ) == false );  // 9 > 5
    _test( t.IsWithinRange(  5.0, 15.0 ) == false );  // 1 < 5
    _test( t.IsWithinRange(  1.0,  9.0 ) == true  );  // diagonal: 1,5,9
    _test( t.IsWithinRange(  4.0,  6.0 ) == false );  // 1 < 4
    _test( t.IsWithinRange( -1.0, 10.0 ) == true  );  // diagonal: 1,5,9

    // Verify off-diagonal elements are NOT checked.
    // Diagonal = (2,2,2), off-diagonal = 100 — should pass.
    TensorVariable<3U> t2( PLAIN, 2.0, 100.0, 100.0,
                                   100.0, 2.0, 100.0,
                                   100.0, 100.0, 2.0 );
    _test( t2.IsWithinRange( 1.0, 3.0 ) == true );
}

// ============================================================================

void TensorVariable_Test::Adjoint_Function()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           9.0, 1.0, 0.0, 4.0, 3.0, 8.0, 1.0, -1.0, 1.0 );
    TensorVariable<3U> expected( INIT_GUESS, PLAIN, NEUMANN,
                                 11.0,  4.0,  -7.0,
                                 -1.0,  9.0,  10.0,
                                  8.0, -72.0, 23.0 );

    TensorVariable<3U> result = t1.Adjoint();
    _test( result == expected );
}

// ============================================================================

void TensorVariable_Test::Identity_Function()
{
    TensorVariable<3U> t( INIT_GUESS, PLAIN, NEUMANN,
                          1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> expected( INIT_GUESS, PLAIN, NEUMANN,
                                 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 );
    t.Identity();
    _test( t == expected );
}

// ============================================================================

void TensorVariable_Test::Transposed_Function()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> expected( INIT_GUESS, PLAIN, NEUMANN,
                                 1.0,4.0,7.0, 2.0,5.0,8.0, 3.0,6.0,9.0 );

    TensorVariable<3U> result = t1.Transposed();
    _test( result == expected );

    // Transpose of transpose = original.
    _test( result.Transposed() == t1 );
}

// ============================================================================

void TensorVariable_Test::Determinant_Function()
{
    TensorVariable<3U> t1( PLAIN, 9.0,1.0,0.0, 4.0,3.0,8.0, 1.0,-1.0,1.0 );
    TensorVariable<3U> t2( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );

    _equal( t1.Determinant(), 103.0, fTolerance );
    _equal( t2.Determinant(),   0.0, fTolerance );  // singular

    // Identity has determinant 1.
    TensorVariable<3U> id( PLAIN, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 );
    _equal( id.Determinant(), 1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Inverse_Function()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           9.0,1.0,0.0, 4.0,3.0,8.0, 1.0,-1.0,1.0 );

    // Expected inverse = (1/103) * adjoint^T
    TensorVariable<3U> expected( INIT_GUESS, PLAIN, NEUMANN,
                                 11.0, -1.0,  8.0,
                                  4.0,  9.0, -72.0,
                                 -7.0, 10.0,  23.0 );
    expected /= 103.0;

    TensorVariable<3U> inv = t1.Inverse();
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( inv(i,j), expected(i,j), fTolerance );

    // A * A^-1 = I.
    TensorVariable<3U> product = t1 * inv;
    _equal( product(0,0), 1.0, fTolerance );
    _equal( product(1,1), 1.0, fTolerance );
    _equal( product(2,2), 1.0, fTolerance );
    _equal( product(0,1), 0.0, fTolerance );
    _equal( product(0,2), 0.0, fTolerance );
    _equal( product(1,0), 0.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::DiagonalValues_Function()
{
    TensorVariable<3U> t( INIT_GUESS, PLAIN, NEUMANN,
                          1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    VectorVariable<3U> v( PLAIN, PLAIN, PLAIN, 10.0, 10.0, 10.0 );
    vector<double>     vec = { 100.0, 100.0, 100.0 };

    // DiagonalValues( std::vector ) — only diagonal changed.
    t.DiagonalValues( vec );
    _equal( t(0,0), 100.0, fTolerance );
    _equal( t(1,1), 100.0, fTolerance );
    _equal( t(2,2), 100.0, fTolerance );
    _equal( t(0,1),   2.0, fTolerance );  // off-diagonal unchanged

    // DiagonalValues( VectorVariable ) — diagonal and flags changed.
    t.DiagonalValues( v );
    _equal( t(0,0), 10.0, fTolerance );
    _equal( t(1,1), 10.0, fTolerance );
    _equal( t(2,2), 10.0, fTolerance );
    _test( t.Flag(0) == PLAIN );
    _test( t.Flag(1) == PLAIN );
    _test( t.Flag(2) == PLAIN );

    // DiagonalValues( double, double, double ).
    t.DiagonalValues( 2.0, 3.0, 4.0 );
    _equal( t(0,0), 2.0, fTolerance );
    _equal( t(1,1), 3.0, fTolerance );
    _equal( t(2,2), 4.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::EigenValues_Function()
{
    // Test example verified in Maple (Karim Heinz Muxi, soil mechanics).
    // Tolerance widened to 1e-4 since the expected values are given to
    // 6 significant digits only.
    TensorVariable<3U> t( PLAIN,  40., -20.,  10.,
                                  -20., -80.,   5.,
                                   10.,   5.,  60. );

    VectorVariable<3U> expected( PLAIN, PLAIN, PLAIN,
                                 64.213203, 39.331419, -83.544623 );
    VectorVariable<3U> result;
    t.EigenValues( result );

    _equal( result(0), expected(0), 1.0e-4 );
    _equal( result(1), expected(1), 1.0e-4 );
    _equal( result(2), expected(2), 1.0e-4 );
}

// ============================================================================

void TensorVariable_Test::Trace_Function()
{
    TensorVariable<3U> t( PLAIN, 9.0,1.0,0.0, 4.0,3.0,8.0, 1.0,-1.0,1.0 );
    _equal( t.Trace(), 13.0, fTolerance );

    // Trace of identity = 3.
    TensorVariable<3U> id( PLAIN, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 );
    _equal( id.Trace(), 3.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::AssignToRow_Function()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> t2;
    VectorVariable<3U> zero( PLAIN, PLAIN, PLAIN, 0.0, 0.0, 0.0 );

    // Assign to row 0.
    t2 = t1; t2.AssignToRow( 0, zero );
    _equal( t2(0,0), 0.0, fTolerance ); _equal( t2(0,1), 0.0, fTolerance );
    _equal( t2(0,2), 0.0, fTolerance ); _equal( t2(1,0), 4.0, fTolerance );
    _equal( t2(1,1), 5.0, fTolerance ); _equal( t2(1,2), 6.0, fTolerance );
    _equal( t2(2,0), 7.0, fTolerance ); _equal( t2(2,1), 8.0, fTolerance );
    _equal( t2(2,2), 9.0, fTolerance );

    // Assign to row 1.
    t2 = t1; t2.AssignToRow( 1, zero );
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(0,2), 3.0, fTolerance ); _equal( t2(1,0), 0.0, fTolerance );
    _equal( t2(1,1), 0.0, fTolerance ); _equal( t2(1,2), 0.0, fTolerance );
    _equal( t2(2,0), 7.0, fTolerance ); _equal( t2(2,1), 8.0, fTolerance );
    _equal( t2(2,2), 9.0, fTolerance );

    // Assign to row 2.
    t2 = t1; t2.AssignToRow( 2, zero );
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(0,2), 3.0, fTolerance ); _equal( t2(1,0), 4.0, fTolerance );
    _equal( t2(1,1), 5.0, fTolerance ); _equal( t2(1,2), 6.0, fTolerance );
    _equal( t2(2,0), 0.0, fTolerance ); _equal( t2(2,1), 0.0, fTolerance );
    _equal( t2(2,2), 0.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::AssignToColumn_Function()
{
    TensorVariable<3U> t1( INIT_GUESS, PLAIN, NEUMANN,
                           1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    TensorVariable<3U> t2;
    VectorVariable<3U> zero( PLAIN, PLAIN, PLAIN, 0.0, 0.0, 0.0 );

    // Assign to column 0.
    t2 = t1; t2.AssignToColumn( 0, zero );
    _equal( t2(0,0), 0.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(0,2), 3.0, fTolerance ); _equal( t2(1,0), 0.0, fTolerance );
    _equal( t2(1,1), 5.0, fTolerance ); _equal( t2(1,2), 6.0, fTolerance );
    _equal( t2(2,0), 0.0, fTolerance ); _equal( t2(2,1), 8.0, fTolerance );
    _equal( t2(2,2), 9.0, fTolerance );

    // Assign to column 1.
    t2 = t1; t2.AssignToColumn( 1, zero );
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 0.0, fTolerance );
    _equal( t2(0,2), 3.0, fTolerance ); _equal( t2(1,0), 4.0, fTolerance );
    _equal( t2(1,1), 0.0, fTolerance ); _equal( t2(1,2), 6.0, fTolerance );
    _equal( t2(2,0), 7.0, fTolerance ); _equal( t2(2,1), 0.0, fTolerance );
    _equal( t2(2,2), 9.0, fTolerance );

    // Assign to column 2.
    t2 = t1; t2.AssignToColumn( 2, zero );
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(0,2), 0.0, fTolerance ); _equal( t2(1,0), 4.0, fTolerance );
    _equal( t2(1,1), 5.0, fTolerance ); _equal( t2(1,2), 0.0, fTolerance );
    _equal( t2(2,0), 7.0, fTolerance ); _equal( t2(2,1), 8.0, fTolerance );
    _equal( t2(2,2), 0.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Row_Function()
{
    TensorVariable<3U> t( INIT_GUESS, PLAIN, NEUMANN,
                          1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    VectorVariable<3U> v;

    v = t.Row(0);
    _equal( v(0), 1.0, fTolerance );
    _equal( v(1), 2.0, fTolerance );
    _equal( v(2), 3.0, fTolerance );

    v = t.Row(1);
    _equal( v(0), 4.0, fTolerance );
    _equal( v(1), 5.0, fTolerance );
    _equal( v(2), 6.0, fTolerance );

    v = t.Row(2);
    _equal( v(0), 7.0, fTolerance );
    _equal( v(1), 8.0, fTolerance );
    _equal( v(2), 9.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::Column_Function()
{
    TensorVariable<3U> t( INIT_GUESS, PLAIN, NEUMANN,
                          1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    VectorVariable<3U> v;

    v = t.Column(0);
    _equal( v(0), 1.0, fTolerance );
    _equal( v(1), 4.0, fTolerance );
    _equal( v(2), 7.0, fTolerance );

    v = t.Column(1);
    _equal( v(0), 2.0, fTolerance );
    _equal( v(1), 5.0, fTolerance );
    _equal( v(2), 8.0, fTolerance );

    v = t.Column(2);
    _equal( v(0), 3.0, fTolerance );
    _equal( v(1), 6.0, fTolerance );
    _equal( v(2), 9.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::HadamardSquared_Function()
{
    // Component-wise squaring: T_ij^2.
    TensorVariable<3U> t( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    t.HadamardSquared();
    _equal( t(0,0),  1.0, fTolerance );
    _equal( t(0,1),  4.0, fTolerance );
    _equal( t(0,2),  9.0, fTolerance );
    _equal( t(1,0), 16.0, fTolerance );
    _equal( t(1,1), 25.0, fTolerance );
    _equal( t(1,2), 36.0, fTolerance );
    _equal( t(2,0), 49.0, fTolerance );
    _equal( t(2,1), 64.0, fTolerance );
    _equal( t(2,2), 81.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::MatrixSquared_Function()
{
    // Matrix product T*T.
    // For uniform 3x3 with all entries = 2:
    // result(i,j) = sum_k 2*2 = 3*4 = 12.
    TensorVariable<3U> t( PLAIN, 2.0,2.0,2.0, 2.0,2.0,2.0, 2.0,2.0,2.0 );
    t.MatrixSquared();
    for ( uint32_t i = 0; i < 3; ++i )
        for ( uint32_t j = 0; j < 3; ++j )
            _equal( t(i,j), 12.0, fTolerance );

    // Identity squared = identity.
    TensorVariable<3U> id( PLAIN, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 );
    id.MatrixSquared();
    _equal( id(0,0), 1.0, fTolerance ); _equal( id(0,1), 0.0, fTolerance );
    _equal( id(1,1), 1.0, fTolerance ); _equal( id(2,2), 1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test::DoubleContraction_Function()
{
    // T:T = sum_ij T_ij^2.
    // For T = diag(1,2,3): T:T = 1+4+9 = 14.
    TensorVariable<3U> t( PLAIN, 1.0,0.0,0.0, 0.0,2.0,0.0, 0.0,0.0,3.0 );
    _equal( t.DoubleContraction(), 14.0, fTolerance );

    // For uniform tensor with all entries = 1: T:T = 9.
    TensorVariable<3U> t2( PLAIN, 1.0,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,1.0 );
    _equal( t2.DoubleContraction(), 9.0, fTolerance );

    // DoubleContraction = sum of HadamardSquared elements.
    TensorVariable<3U> t3( PLAIN, 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0 );
    const double expected = 1+4+9+16+25+36+49+64+81;  // = 285
    _equal( t3.DoubleContraction(), expected, fTolerance );
}

// ============================================================================
//  TestEigenMethods — unchanged from previous version except:
//  - boolean _equal calls use tolerance 0.5
//  - asymmetry tolerances are explicit
// ============================================================================

void TensorVariable_Test::TestEigenMethods()
{
    if ( verbose_ )
    {
        cout << "\n================================" << endl;
        cout << "Testing Eigen methods"             << endl;
        cout << "================================"  << endl;
    }

    auto residual2D = []( const TensorVariable<2U>& A,
                          const VectorVariable<2U>& v,
                          double lambda ) -> double
    {
        double r0 = A(0,0)*v(0) + A(0,1)*v(1) - lambda*v(0);
        double r1 = A(1,0)*v(0) + A(1,1)*v(1) - lambda*v(1);
        return std::sqrt( r0*r0 + r1*r1 );
    };

    auto residual3D = []( const TensorVariable<3U>& A,
                          const VectorVariable<3U>& v,
                          double lambda ) -> double
    {
        double r0 = A(0,0)*v(0)+A(0,1)*v(1)+A(0,2)*v(2) - lambda*v(0);
        double r1 = A(1,0)*v(0)+A(1,1)*v(1)+A(1,2)*v(2) - lambda*v(1);
        double r2 = A(2,0)*v(0)+A(2,1)*v(1)+A(2,2)*v(2) - lambda*v(2);
        return std::sqrt( r0*r0 + r1*r1 + r2*r2 );
    };

    auto col2D = []( const TensorVariable<2U>& M, uint32_t j ) -> VectorVariable<2U>
    {
        VectorVariable<2U> v; v(0)=M(0,j); v(1)=M(1,j); return v;
    };

    auto col3D = []( const TensorVariable<3U>& M, uint32_t j ) -> VectorVariable<3U>
    {
        VectorVariable<3U> v; v(0)=M(0,j); v(1)=M(1,j); v(2)=M(2,j); return v;
    };

    auto dot2D = []( const VectorVariable<2U>& a, const VectorVariable<2U>& b ) -> double
    { return a(0)*b(0) + a(1)*b(1); };

    auto dot3D = []( const VectorVariable<3U>& a, const VectorVariable<3U>& b ) -> double
    { return a(0)*b(0) + a(1)*b(1) + a(2)*b(2); };

    // ==================================================================
    //  2D TESTS
    // ==================================================================

    if ( verbose_ ) cout << "\n--- 2D EigenValues ---" << endl;

    // Test 1: Diagonal 2x2 — eigenvalues 5, 2.
    {
        TensorVariable<2U> A; A(0,0)=5.; A(0,1)=0.; A(1,0)=0.; A(1,1)=2.;
        VectorVariable<2U> ev;
        _equal( A.EigenValues(ev), true, 0.5 );
        _equal( ev(0), 5.0, fTolerance );
        _equal( ev(1), 2.0, fTolerance );
        if ( verbose_ ) cout << "  Test 1 (diagonal 2D): passed" << endl;
    }

    // Test 2: Symmetric off-diagonal — eigenvalues 4, 2.
    {
        TensorVariable<2U> A; A(0,0)=3.; A(0,1)=1.; A(1,0)=1.; A(1,1)=3.;
        VectorVariable<2U> ev;
        _equal( A.EigenValues(ev), true, 0.5 );
        _equal( ev(0), 4.0, fTolerance );
        _equal( ev(1), 2.0, fTolerance );
        if ( verbose_ ) cout << "  Test 2 (symmetric 2D): passed" << endl;
    }

    // Test 3: Complex roots — returns false, NaN.
    {
        TensorVariable<2U> A; A(0,0)=0.; A(0,1)=1.; A(1,0)=-1.; A(1,1)=0.;
        VectorVariable<2U> ev;
        _equal( A.EigenValues(ev), false, 0.5 );
        _equal( std::isnan(ev(0)), true, 0.5 );
        _equal( std::isnan(ev(1)), true, 0.5 );
        if ( verbose_ ) cout << "  Test 3 (complex roots 2D): passed" << endl;
    }

    if ( verbose_ ) cout << "\n--- 2D EigenSymmetric ---" << endl;

    // Test 4: Diagonal 2x2 eigenvectors.
    {
        TensorVariable<2U> A; A(0,0)=5.; A(0,1)=0.; A(1,0)=0.; A(1,1)=2.;
        VectorVariable<2U> ev; TensorVariable<2U> evecs;
        _equal( A.EigenSymmetric(ev, evecs, true), true, 0.5 );
        VectorVariable<2U> v0=col2D(evecs,0), v1=col2D(evecs,1);
        _equal( residual2D(A,v0,ev(0)), 0.0, fTolerance );
        _equal( residual2D(A,v1,ev(1)), 0.0, fTolerance );
        _equal( std::fabs(dot2D(v0,v1)), 0.0, fTolerance );
        _equal( dot2D(v0,v0), 1.0, fTolerance );
        _equal( dot2D(v1,v1), 1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 4 (diagonal 2D eigenvectors): passed" << endl;
    }

    // Test 5: Symmetric off-diagonal 2x2 eigenvectors.
    {
        TensorVariable<2U> A; A(0,0)=3.; A(0,1)=1.; A(1,0)=1.; A(1,1)=3.;
        VectorVariable<2U> ev; TensorVariable<2U> evecs;
        _equal( A.EigenSymmetric(ev, evecs, true), true, 0.5 );
        VectorVariable<2U> v0=col2D(evecs,0), v1=col2D(evecs,1);
        _equal( residual2D(A,v0,ev(0)), 0.0, fTolerance );
        _equal( residual2D(A,v1,ev(1)), 0.0, fTolerance );
        _equal( std::fabs(dot2D(v0,v1)), 0.0, fTolerance );
        _equal( dot2D(v0,v0), 1.0, fTolerance );
        _equal( dot2D(v1,v1), 1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 5 (symmetric off-diagonal 2D): passed" << endl;
    }

    if ( verbose_ ) cout << "\n--- 2D EigenWeaklyNonSymmetric ---" << endl;

    // Test 6: Weakly non-symmetric 2x2 — within tolerance 1e-3.
    {
        TensorVariable<2U> A;
        A(0,0)=3.; A(0,1)=1.0001; A(1,0)=0.9999; A(1,1)=3.;
        VectorVariable<2U> ev; TensorVariable<2U> evecs;
        _equal( A.EigenWeaklyNonSymmetric(ev, evecs, 1.0e-3), true, 0.5 );
        _equal( ev(0), 4.0, 1.0e-3 );
        _equal( ev(1), 2.0, 1.0e-3 );
        if ( verbose_ ) cout << "  Test 6 (weakly non-symmetric 2D): passed" << endl;
    }

    // Test 7: Strongly non-symmetric 2x2 — exceeds default tolerance.
    {
        TensorVariable<2U> A;
        A(0,0)=3.; A(0,1)=5.; A(1,0)=0.; A(1,1)=3.;
        VectorVariable<2U> ev; TensorVariable<2U> evecs;
        _equal( A.EigenWeaklyNonSymmetric(ev, evecs, 1.0e-6), false, 0.5 );
        if ( verbose_ ) cout << "  Test 7 (strongly non-symmetric 2D): passed" << endl;
    }

    // ==================================================================
    //  3D TESTS
    // ==================================================================

    if ( verbose_ ) cout << "\n--- 3D EigenValuesPositiveDefiniteSymmetric ---" << endl;

    // Test 8: Diagonal 3x3 — eigenvalues 6, 3, 1.
    {
        TensorVariable<3U> A;
        A(0,0)=6.; A(0,1)=0.; A(0,2)=0.;
        A(1,0)=0.; A(1,1)=3.; A(1,2)=0.;
        A(2,0)=0.; A(2,1)=0.; A(2,2)=1.;
        double e0,e1,e2;
        _equal( A.EigenValuesPositiveDefiniteSymmetricMatrix(e0,e1,e2), true, 0.5 );
        _equal( e0, 6.0, fTolerance );
        _equal( e1, 3.0, fTolerance );
        _equal( e2, 1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 8 (diagonal 3D): passed" << endl;
    }

    // Test 9: Symmetric tridiagonal — eigenvalues 2+sqrt(2), 2, 2-sqrt(2).
    {
        TensorVariable<3U> A;
        A(0,0)=2.; A(0,1)=1.; A(0,2)=0.;
        A(1,0)=1.; A(1,1)=2.; A(1,2)=1.;
        A(2,0)=0.; A(2,1)=1.; A(2,2)=2.;
        double e0,e1,e2;
        _equal( A.EigenValuesPositiveDefiniteSymmetricMatrix(e0,e1,e2), true, 0.5 );
        _equal( e0, 2.0+std::sqrt(2.0), fTolerance );
        _equal( e1, 2.0,                fTolerance );
        _equal( e2, 2.0-std::sqrt(2.0), fTolerance );
        if ( verbose_ ) cout << "  Test 9 (symmetric tridiagonal 3D): passed" << endl;
    }

    // Test 10: Two equal eigenvalues — diag(3,3,1).
    {
        TensorVariable<3U> A;
        A(0,0)=3.; A(0,1)=0.; A(0,2)=0.;
        A(1,0)=0.; A(1,1)=3.; A(1,2)=0.;
        A(2,0)=0.; A(2,1)=0.; A(2,2)=1.;
        double e0,e1,e2;
        _equal( A.EigenValuesPositiveDefiniteSymmetricMatrix(e0,e1,e2), true, 0.5 );
        _equal( e0, 3.0, fTolerance );
        _equal( e1, 3.0, fTolerance );
        _equal( e2, 1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 10 (two equal eigenvalues 3D): passed" << endl;
    }

    if ( verbose_ ) cout << "\n--- 3D EigenSymmetric ---" << endl;

    // Test 11: Diagonal 3x3 eigenvectors.
    {
        TensorVariable<3U> A;
        A(0,0)=6.; A(0,1)=0.; A(0,2)=0.;
        A(1,0)=0.; A(1,1)=3.; A(1,2)=0.;
        A(2,0)=0.; A(2,1)=0.; A(2,2)=1.;
        VectorVariable<3U> ev; TensorVariable<3U> evecs;
        _equal( A.EigenSymmetric(ev, evecs), true, 0.5 );
        VectorVariable<3U> v0=col3D(evecs,0), v1=col3D(evecs,1), v2=col3D(evecs,2);
        _equal( residual3D(A,v0,ev(0)), 0.0, fTolerance );
        _equal( residual3D(A,v1,ev(1)), 0.0, fTolerance );
        _equal( residual3D(A,v2,ev(2)), 0.0, fTolerance );
        _equal( std::fabs(dot3D(v0,v1)), 0.0, fTolerance );
        _equal( std::fabs(dot3D(v0,v2)), 0.0, fTolerance );
        _equal( std::fabs(dot3D(v1,v2)), 0.0, fTolerance );
        _equal( dot3D(v0,v0), 1.0, fTolerance );
        _equal( dot3D(v1,v1), 1.0, fTolerance );
        _equal( dot3D(v2,v2), 1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 11 (diagonal 3D eigenvectors): passed" << endl;
    }

    // Test 12: Symmetric tridiagonal 3x3 eigenvectors.
    {
        TensorVariable<3U> A;
        A(0,0)=2.; A(0,1)=1.; A(0,2)=0.;
        A(1,0)=1.; A(1,1)=2.; A(1,2)=1.;
        A(2,0)=0.; A(2,1)=1.; A(2,2)=2.;
        VectorVariable<3U> ev; TensorVariable<3U> evecs;
        _equal( A.EigenSymmetric(ev, evecs), true, 0.5 );
        VectorVariable<3U> v0=col3D(evecs,0), v1=col3D(evecs,1), v2=col3D(evecs,2);
        _equal( residual3D(A,v0,ev(0)), 0.0, fTolerance );
        _equal( residual3D(A,v1,ev(1)), 0.0, fTolerance );
        _equal( residual3D(A,v2,ev(2)), 0.0, fTolerance );
        _equal( std::fabs(dot3D(v0,v1)), 0.0, fTolerance );
        _equal( std::fabs(dot3D(v0,v2)), 0.0, fTolerance );
        _equal( std::fabs(dot3D(v1,v2)), 0.0, fTolerance );
        _equal( dot3D(v0,v0), 1.0, fTolerance );
        _equal( dot3D(v1,v1), 1.0, fTolerance );
        _equal( dot3D(v2,v2), 1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 12 (symmetric tridiagonal 3D): passed" << endl;
    }

    // Test 13: EigenSymmetric and SPD agree.
    {
        TensorVariable<3U> A;
        A(0,0)=4.; A(0,1)=2.; A(0,2)=0.;
        A(1,0)=2.; A(1,1)=3.; A(1,2)=1.;
        A(2,0)=0.; A(2,1)=1.; A(2,2)=2.;
        double e0s,e1s,e2s;
        A.EigenValuesPositiveDefiniteSymmetricMatrix(e0s,e1s,e2s);
        VectorVariable<3U> ev; TensorVariable<3U> evecs;
        A.EigenSymmetric(ev, evecs);
        _equal( ev(0), e0s, fTolerance );
        _equal( ev(1), e1s, fTolerance );
        _equal( ev(2), e2s, fTolerance );
        if ( verbose_ ) cout << "  Test 13 (EigenSymmetric vs SPD): passed" << endl;
    }

    if ( verbose_ ) cout << "\n--- 3D EigenWeaklyNonSymmetric ---" << endl;

    // Test 14: Weakly non-symmetric 3x3 — within tolerance 1e-4.
    {
        TensorVariable<3U> A;
        A(0,0)=2.;      A(0,1)=1.00001; A(0,2)=0.;
        A(1,0)=0.99999; A(1,1)=2.;      A(1,2)=1.00001;
        A(2,0)=0.;      A(2,1)=0.99999; A(2,2)=2.;
        VectorVariable<3U> ev; TensorVariable<3U> evecs;
        _equal( A.EigenWeaklyNonSymmetric(ev, evecs, 1.0e-4), true, 0.5 );
        _equal( ev(0), 2.0+std::sqrt(2.0), fTolerance );
        _equal( ev(1), 2.0,                fTolerance );
        _equal( ev(2), 2.0-std::sqrt(2.0), fTolerance );
        if ( verbose_ ) cout << "  Test 14 (weakly non-symmetric 3D): passed" << endl;
    }

    // Test 15: Strongly non-symmetric 3x3 — exceeds default tolerance.
    {
        TensorVariable<3U> A;
        A(0,0)=2.; A(0,1)=5.; A(0,2)=3.;
        A(1,0)=0.; A(1,1)=2.; A(1,2)=4.;
        A(2,0)=0.; A(2,1)=0.; A(2,2)=2.;
        VectorVariable<3U> ev; TensorVariable<3U> evecs;
        _equal( A.EigenWeaklyNonSymmetric(ev, evecs, 1.0e-6), false, 0.5 );
        if ( verbose_ ) cout << "  Test 15 (strongly non-symmetric 3D): passed" << endl;
    }

    // Test 16: Non-positive-definite — EigenSymmetric handles negative eigenvalues.
    {
        TensorVariable<3U> A;
        A(0,0)=-1.; A(0,1)=0.; A(0,2)=0.;
        A(1,0)=0.;  A(1,1)=2.; A(1,2)=0.;
        A(2,0)=0.;  A(2,1)=0.; A(2,2)=3.;
        VectorVariable<3U> ev; TensorVariable<3U> evecs;
        _equal( A.EigenSymmetric(ev, evecs), true, 0.5 );
        _equal( ev(0),  3.0, fTolerance );
        _equal( ev(1),  2.0, fTolerance );
        _equal( ev(2), -1.0, fTolerance );
        if ( verbose_ ) cout << "  Test 16 (non-positive-definite 3D): passed" << endl;
    }

    if ( verbose_ ) cout << "\nAll Eigen method tests completed." << endl;
}

} // namespace csmp

