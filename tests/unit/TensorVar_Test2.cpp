#include "TensorVar_Test2.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "CSMP_definitions.h"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace csmp {

TensorVariable_Test2::TensorVariable_Test2()
{
    fTolerance = 1.0e-6;
}

TensorVariable_Test2::~TensorVariable_Test2()
{
}

void TensorVariable_Test2::run()
{
    Assignment_Operator2();
    Addition_Operator2();
    Subtraction_Operator2();
    Multiplication_Operator2();
    Division_Operator2();
    Addition_Assignment_Operator2();
    Subtraction_Assignment_Operator2();
    Multiplication_Assignment_Operator2();
    Division_Assignment_Operator2();
    Equality_Operator2();
    LessThan_Operator2();           // was missing from run()
    MinElement_Function2();
    MaxElement_Function2();
    IsWithinRange_Function2();
    Adjoint_Function2();
    Identity_Function2();
    Transposed_Function2();
    Determinant_Function2();
    Inverse_Function2();
    DiagonalValues_Function2();
    EigenValues_Function2();
    Trace_Function2();
    AssignToRow_Function2();
    AssignToColumn_Function2();
    Row_Function2();
    Column_Function2();
    HadamardSquared_Function2();    // new
    MatrixSquared_Function2();      // new
    DoubleContraction_Function2();  // new
}

// ============================================================================

void TensorVariable_Test2::Assignment_Operator2()
{
    ScalarVariable     sc( PLAIN, 5.0 );
    VectorVariable<2U> v( PLAIN, ROBIN, 10.0, 10.0 );
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2;

    // VectorVariable assignment — diagonal set, off-diagonal zeroed.
    t2 = v;
    _equal( t2(0,0), 10.0, fTolerance ); _equal( t2(0,1), 0.0,  fTolerance );
    _equal( t2(1,0),  0.0, fTolerance ); _equal( t2(1,1), 10.0, fTolerance );
    _test( t2.Flag(0) == PLAIN );
    _test( t2.Flag(1) == ROBIN );

    // TensorVariable assignment — all elements and flags copied.
    t2 = t1;
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(1,0), 3.0, fTolerance ); _equal( t2(1,1), 4.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );

    // double assignment — all elements set, flags unchanged.
    t2 = 10.0;
    _equal( t2(0,0), 10.0, fTolerance ); _equal( t2(0,1), 10.0, fTolerance );
    _equal( t2(1,0), 10.0, fTolerance ); _equal( t2(1,1), 10.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );

    // ScalarVariable assignment — all elements and all flags set.
    t2 = sc;
    _equal( t2(0,0), 5.0, fTolerance ); _equal( t2(0,1), 5.0, fTolerance );
    _equal( t2(1,0), 5.0, fTolerance ); _equal( t2(1,1), 5.0, fTolerance );
    _test( t2.Flag(0) == PLAIN );
    _test( t2.Flag(1) == PLAIN );
}

// ============================================================================

void TensorVariable_Test2::Addition_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, PLAIN,      9.0, 8.0, 7.0, 6.0 );
    TensorVariable<2U> t3;

    // Tensor + tensor — flags from lhs.
    t3 = t1 + t2;
    _equal( t3(0,0), 10.0, fTolerance ); _equal( t3(0,1), 10.0, fTolerance );
    _equal( t3(1,0), 10.0, fTolerance ); _equal( t3(1,1), 10.0, fTolerance );
    _test( t3.Flag(0) == INIT_GUESS );
    _test( t3.Flag(1) == PLAIN     );

    // Tensor + double.
    t3 = t1 + 10.0;
    _equal( t3(0,0), 11.0, fTolerance ); _equal( t3(0,1), 12.0, fTolerance );
    _equal( t3(1,0), 13.0, fTolerance ); _equal( t3(1,1), 14.0, fTolerance );
    _test( t3.Flag(0) == INIT_GUESS );
    _test( t3.Flag(1) == PLAIN     );
}

// ============================================================================

void TensorVariable_Test2::Subtraction_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, PLAIN,      9.0, 8.0, 7.0, 6.0 );
    TensorVariable<2U> t3;

    // Tensor - tensor.
    t3 = t1 - t2;
    _equal( t3(0,0), -8.0, fTolerance ); _equal( t3(0,1), -6.0, fTolerance );
    _equal( t3(1,0), -4.0, fTolerance ); _equal( t3(1,1), -2.0, fTolerance );
    _test( t3.Flag(0) == INIT_GUESS );
    _test( t3.Flag(1) == PLAIN     );

    // Tensor - double.
    t3 = t1 - 10.0;
    _equal( t3(0,0), -9.0, fTolerance ); _equal( t3(0,1), -8.0, fTolerance );
    _equal( t3(1,0), -7.0, fTolerance ); _equal( t3(1,1), -6.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Multiplication_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2;
    TensorVariable<2U> t3( PLAIN, PLAIN, 9.0, 8.0, 7.0, 6.0 );
    VectorVariable<2U> v1( PLAIN, PERIODIC, 1.0, 1.0 );
    VectorVariable<2U> v2;

    // Tensor * double.
    t2 = t1 * 10.0;
    _equal( t2(0,0), 10.0, fTolerance ); _equal( t2(0,1), 20.0, fTolerance );
    _equal( t2(1,0), 30.0, fTolerance ); _equal( t2(1,1), 40.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );

    // Vector-matrix multiplication: v^T * A.
    v2 = v1 * t1;
    _equal( v2(0), 4.0, fTolerance );
    _equal( v2(1), 6.0, fTolerance );

    // Matrix-vector multiplication: A * v.
    v2 = t1 * v1;
    _equal( v2(0), 3.0, fTolerance );
    _equal( v2(1), 7.0, fTolerance );
    _test( v2.Flag(0) == PLAIN    );
    _test( v2.Flag(1) == PERIODIC );

    // Matrix-matrix multiplication.
    t2 = t1 * t3;
    _equal( t2(0,0), 23.0, fTolerance ); _equal( t2(0,1), 20.0, fTolerance );
    _equal( t2(1,0), 55.0, fTolerance ); _equal( t2(1,1), 48.0, fTolerance );
    _test( t2.Flag(0) == INIT_GUESS );
    _test( t2.Flag(1) == PLAIN     );
}

// ============================================================================

void TensorVariable_Test2::Division_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t3( PERIODIC, PLAIN,  10.0, 20.0, 30.0, 40.0 );
    TensorVariable<2U> t2;

    // Tensor / double.
    t2 = t1 / 10.0;
    _equal( t2(0,0), 0.1, fTolerance ); _equal( t2(0,1), 0.2, fTolerance );
    _equal( t2(1,0), 0.3, fTolerance ); _equal( t2(1,1), 0.4, fTolerance );

    // Element-by-element division — all four elements checked.
    t2 = t3 / t1;
    _equal( t2(0,0), 10.0, fTolerance ); _equal( t2(0,1), 10.0, fTolerance );
    _equal( t2(1,0), 10.0, fTolerance ); _equal( t2(1,1), 10.0, fTolerance );
    _test( t2.Flag(0) == PERIODIC );
    _test( t2.Flag(1) == PLAIN    );
}

// ============================================================================

void TensorVariable_Test2::Addition_Assignment_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, PLAIN,      9.0, 8.0, 7.0, 6.0 );
    TensorVariable<2U> t3;
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor += tensor.
    t3 = t1; t3 += t2;
    _equal( t3(0,0), 10.0, fTolerance ); _equal( t3(0,1), 10.0, fTolerance );
    _equal( t3(1,0), 10.0, fTolerance ); _equal( t3(1,1), 10.0, fTolerance );

    // Tensor += double.
    t3 = t1; t3 += 10.0;
    _equal( t3(0,0), 11.0, fTolerance ); _equal( t3(0,1), 12.0, fTolerance );
    _equal( t3(1,0), 13.0, fTolerance ); _equal( t3(1,1), 14.0, fTolerance );

    // Tensor += scalar.
    t3 = t1; t3 += sc;
    _equal( t3(0,0), 6.0, fTolerance ); _equal( t3(0,1), 7.0, fTolerance );
    _equal( t3(1,0), 8.0, fTolerance ); _equal( t3(1,1), 9.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Subtraction_Assignment_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, PLAIN,      9.0, 8.0, 7.0, 6.0 );
    TensorVariable<2U> t3;
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor -= tensor.
    t3 = t1; t3 -= t2;
    _equal( t3(0,0), -8.0, fTolerance ); _equal( t3(0,1), -6.0, fTolerance );
    _equal( t3(1,0), -4.0, fTolerance ); _equal( t3(1,1), -2.0, fTolerance );

    // Tensor -= double.
    t3 = t1; t3 -= 10.0;
    _equal( t3(0,0), -9.0, fTolerance ); _equal( t3(0,1), -8.0, fTolerance );
    _equal( t3(1,0), -7.0, fTolerance ); _equal( t3(1,1), -6.0, fTolerance );

    // Tensor -= scalar.
    t3 = t1; t3 -= sc;
    _equal( t3(0,0), -4.0, fTolerance ); _equal( t3(0,1), -3.0, fTolerance );
    _equal( t3(1,0), -2.0, fTolerance ); _equal( t3(1,1), -1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Multiplication_Assignment_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2;
    TensorVariable<2U> t3( PLAIN, PLAIN, 9.0, 8.0, 7.0, 6.0 );
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor *= double.
    t2 = t1; t2 *= 10.0;
    _equal( t2(0,0), 10.0, fTolerance ); _equal( t2(0,1), 20.0, fTolerance );
    _equal( t2(1,0), 30.0, fTolerance ); _equal( t2(1,1), 40.0, fTolerance );

    // Tensor *= tensor (matrix product).
    t2 = t1; t2 *= t3;
    _equal( t2(0,0), 23.0, fTolerance ); _equal( t2(0,1), 20.0, fTolerance );
    _equal( t2(1,0), 55.0, fTolerance ); _equal( t2(1,1), 48.0, fTolerance );

    // Tensor *= scalar.
    t2 = t1; t2 *= sc;
    _equal( t2(0,0),  5.0, fTolerance ); _equal( t2(0,1), 10.0, fTolerance );
    _equal( t2(1,0), 15.0, fTolerance ); _equal( t2(1,1), 20.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Division_Assignment_Operator2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, PLAIN,     10.0, 20.0, 30.0, 40.0 );
    TensorVariable<2U> t3;
    ScalarVariable     sc( PLAIN, 5.0 );

    // Tensor /= double.
    t3 = t1; t3 /= 10.0;
    _equal( t3(0,0), 0.1, fTolerance ); _equal( t3(0,1), 0.2, fTolerance );
    _equal( t3(1,0), 0.3, fTolerance ); _equal( t3(1,1), 0.4, fTolerance );

    // Tensor /= tensor (element-by-element).
    t3 = t1; t3 /= t2;
    _equal( t3(0,0), 0.1, fTolerance ); _equal( t3(0,1), 0.1, fTolerance );
    _equal( t3(1,0), 0.1, fTolerance ); _equal( t3(1,1), 0.1, fTolerance );

    // Tensor /= scalar.
    t3 = t1; t3 /= sc;
    _equal( t3(0,0), 0.2, fTolerance ); _equal( t3(0,1), 0.4, fTolerance );
    _equal( t3(1,0), 0.6, fTolerance ); _equal( t3(1,1), 0.8, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Equality_Operator2()
{
    TensorVariable<2U> t1( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, 1.0, 2.0, 3.0, 4.0 );  // same
    TensorVariable<2U> t3( PLAIN, 9.0, 8.0, 7.0, 6.0 );  // different values
    TensorVariable<2U> t4( ROBIN, 1.0, 2.0, 3.0, 4.0 );  // different flags

    _test(  t1 == t2 );
    _test(  t1 != t3 );
    _test(  t1 != t4 );  // different flags — was missing
    _test( !( t1 != t2 ) );
    _test( !( t1 == t3 ) );
}

// ============================================================================

void TensorVariable_Test2::LessThan_Operator2()
{
    // operator< compares pointer addresses — satisfies strict weak ordering
    // for STL associative containers. Does not compare values.
    TensorVariable<2U> t1( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2( PLAIN, 9.0, 8.0, 7.0, 6.0 );

    // Exactly one of (t1 < t2) or (t2 < t1) must be true.
    const bool lt = ( t1 < t2 );
    const bool gt = ( t2 < t1 );
    _test(  lt != gt );       // exactly one is true
    _test( !( t1 < t1 ) );   // irreflexive
    _test( !( t2 < t2 ) );   // irreflexive
}

// ============================================================================

void TensorVariable_Test2::MinElement_Function2()
{
    TensorVariable<2U> t( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    _equal( t.MinElement(), 1.0, fTolerance );

    // Negative values.
    TensorVariable<2U> t2( PLAIN, -5.0, 2.0, 3.0, 4.0 );
    _equal( t2.MinElement(), -5.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::MaxElement_Function2()
{
    TensorVariable<2U> t( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    _equal( t.MaxElement(), 4.0, fTolerance );

    // All negative.
    TensorVariable<2U> t2( PLAIN, -4.0, -3.0, -2.0, -1.0 );
    _equal( t2.MaxElement(), -1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::IsWithinRange_Function2()
{
    // 2D IsWithinRange checks ALL four elements (unlike 3D which checks
    // only the diagonal).
    TensorVariable<2U> t( PLAIN, 1.0, 2.0, 3.0, 4.0 );

    _test( t.IsWithinRange( -5.0,  2.0 ) == false );  // 3 and 4 > 2
    _test( t.IsWithinRange(  3.0, 15.0 ) == false );  // 1 and 2 < 3
    _test( t.IsWithinRange(  1.0,  4.0 ) == true  );  // all in [1,4]
    _test( t.IsWithinRange(  2.0,  3.0 ) == false );  // 1 < 2 and 4 > 3
    _test( t.IsWithinRange( -1.0, 10.0 ) == true  );  // all in range

    // Verify all four elements are checked (not just diagonal).
    // Diagonal = (1,4), off-diagonal = (100,100) — should fail.
    TensorVariable<2U> t2( PLAIN, 1.0, 100.0, 100.0, 4.0 );
    _test( t2.IsWithinRange( 0.0, 5.0 ) == false );
}

// ============================================================================

void TensorVariable_Test2::Adjoint_Function2()
{
    // Adjoint of [[1,2],[3,4]] = [[4,-3],[-2,1]]
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> expected( INIT_GUESS, PLAIN, 4.0, -3.0, -2.0, 1.0 );

    TensorVariable<2U> result = t1.Adjoint();
    _test( result == expected );
}

// ============================================================================

void TensorVariable_Test2::Identity_Function2()
{
    TensorVariable<2U> t( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> expected( INIT_GUESS, PLAIN, 1.0, 0.0, 0.0, 1.0 );

    t.Identity();
    _test( t == expected );
}

// ============================================================================

void TensorVariable_Test2::Transposed_Function2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> expected( INIT_GUESS, PLAIN, 1.0, 3.0, 2.0, 4.0 );

    TensorVariable<2U> result = t1.Transposed();
    _test( result == expected );

    // Transpose of transpose = original.
    _test( result.Transposed() == t1 );
}

// ============================================================================

void TensorVariable_Test2::Determinant_Function2()
{
    // det([[1,2],[3,4]]) = 4 - 6 = -2
    TensorVariable<2U> t1( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    _equal( t1.Determinant(), -2.0, fTolerance );

    // Singular matrix.
    TensorVariable<2U> t2( PLAIN, 1.0, 2.0, 2.0, 4.0 );
    _equal( t2.Determinant(), 0.0, fTolerance );

    // Identity has determinant 1.
    TensorVariable<2U> id( PLAIN, 1.0, 0.0, 0.0, 1.0 );
    _equal( id.Determinant(), 1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Inverse_Function2()
{
    // Inverse of [[1,2],[3,4]] = (1/-2) * [[4,-2],[-3,1]]
    //                          = [[-2,1],[1.5,-0.5]]
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> expected( INIT_GUESS, PLAIN,
                                 -2.0, 1.0, 1.5, -0.5 );

    TensorVariable<2U> inv = t1.Inverse();
    _test( inv == expected );

    // A * A^-1 = I.
    TensorVariable<2U> product = t1 * inv;
    _equal( product(0,0), 1.0, fTolerance ); _equal( product(0,1), 0.0, fTolerance );
    _equal( product(1,0), 0.0, fTolerance ); _equal( product(1,1), 1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::DiagonalValues_Function2()
{
    TensorVariable<2U> t( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    VectorVariable<2U> v( PLAIN, PLAIN, 10.0, 10.0 );
    vector<double>     vec = { 100.0, 100.0 };

    // DiagonalValues( std::vector ) — only diagonal changed.
    t.DiagonalValues( vec );
    _equal( t(0,0), 100.0, fTolerance ); _equal( t(0,1),   2.0, fTolerance );
    _equal( t(1,0),   3.0, fTolerance ); _equal( t(1,1), 100.0, fTolerance );

    // DiagonalValues( VectorVariable ) — diagonal and flags changed.
    t.DiagonalValues( v );
    _equal( t(0,0), 10.0, fTolerance ); _equal( t(0,1),  2.0, fTolerance );
    _equal( t(1,0),  3.0, fTolerance ); _equal( t(1,1), 10.0, fTolerance );
    _test( t.Flag(0) == PLAIN );
    _test( t.Flag(1) == PLAIN );

    // DiagonalValues( double, double ).
    t.DiagonalValues( 5.0, 5.0 );
    _equal( t(0,0), 5.0, fTolerance ); _equal( t(0,1), 2.0, fTolerance );
    _equal( t(1,0), 3.0, fTolerance ); _equal( t(1,1), 5.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::EigenValues_Function2()
{
    // t = [[1,2],[3,4]]
    // a1 = -(1+4) = -5,  a0 = 1*4 - 2*3 = -2
    // D = 25 - 4*(-2) = 33
    // lambda = (5 +/- sqrt(33)) / 2
    TensorVariable<2U> t( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );

    const double sqrtD    = std::sqrt( 33.0 );
    const double lambda0  = ( 5.0 + sqrtD ) / 2.0;
    const double lambda1  = ( 5.0 - sqrtD ) / 2.0;

    VectorVariable<2U> result;
    _equal( t.EigenValues( result ), true, 0.5 );
    _equal( result(0), lambda0, fTolerance );
    _equal( result(1), lambda1, fTolerance );

    // Eigenvalues must be in descending order.
    _test( result(0) > result(1) );
}

// ============================================================================

void TensorVariable_Test2::Trace_Function2()
{
    TensorVariable<2U> t( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 7.0 );
    _equal( t.Trace(), 8.0, fTolerance );

    // Trace of identity = 2.
    TensorVariable<2U> id( PLAIN, 1.0, 0.0, 0.0, 1.0 );
    _equal( id.Trace(), 2.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::AssignToRow_Function2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2;
    VectorVariable<2U> zero( PLAIN, PLAIN, 0.0, 0.0 );

    // Assign to row 0.
    t2 = t1; t2.AssignToRow( 0, zero );
    _equal( t2(0,0), 0.0, fTolerance ); _equal( t2(0,1), 0.0, fTolerance );
    _equal( t2(1,0), 3.0, fTolerance ); _equal( t2(1,1), 4.0, fTolerance );

    // Assign to row 1.
    t2 = t1; t2.AssignToRow( 1, zero );
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(1,0), 0.0, fTolerance ); _equal( t2(1,1), 0.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::AssignToColumn_Function2()
{
    TensorVariable<2U> t1( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    TensorVariable<2U> t2;
    VectorVariable<2U> zero( PLAIN, PLAIN, 0.0, 0.0 );

    // Assign to column 0.
    t2 = t1; t2.AssignToColumn( 0, zero );
    _equal( t2(0,0), 0.0, fTolerance ); _equal( t2(0,1), 2.0, fTolerance );
    _equal( t2(1,0), 0.0, fTolerance ); _equal( t2(1,1), 4.0, fTolerance );

    // Assign to column 1.
    t2 = t1; t2.AssignToColumn( 1, zero );
    _equal( t2(0,0), 1.0, fTolerance ); _equal( t2(0,1), 0.0, fTolerance );
    _equal( t2(1,0), 3.0, fTolerance ); _equal( t2(1,1), 0.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Row_Function2()
{
    TensorVariable<2U> t( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    VectorVariable<2U> v;

    v = t.Row(0);
    _equal( v(0), 1.0, fTolerance ); _equal( v(1), 2.0, fTolerance );

    v = t.Row(1);
    _equal( v(0), 3.0, fTolerance ); _equal( v(1), 4.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::Column_Function2()
{
    TensorVariable<2U> t( INIT_GUESS, PLAIN, 1.0, 2.0, 3.0, 4.0 );
    VectorVariable<2U> v;

    v = t.Column(0);
    _equal( v(0), 1.0, fTolerance ); _equal( v(1), 3.0, fTolerance );

    v = t.Column(1);
    _equal( v(0), 2.0, fTolerance ); _equal( v(1), 4.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::HadamardSquared_Function2()
{
    // Component-wise squaring: T_ij^2.
    TensorVariable<2U> t( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    t.HadamardSquared();
    _equal( t(0,0),  1.0, fTolerance ); _equal( t(0,1),  4.0, fTolerance );
    _equal( t(1,0),  9.0, fTolerance ); _equal( t(1,1), 16.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::MatrixSquared_Function2()
{
    // Matrix product T*T for [[1,2],[3,4]]:
    // result = [[1*1+2*3, 1*2+2*4],[3*1+4*3, 3*2+4*4]]
    //        = [[7, 10],[15, 22]]
    TensorVariable<2U> t( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    t.MatrixSquared();
    _equal( t(0,0),  7.0, fTolerance ); _equal( t(0,1), 10.0, fTolerance );
    _equal( t(1,0), 15.0, fTolerance ); _equal( t(1,1), 22.0, fTolerance );

    // Identity squared = identity.
    TensorVariable<2U> id( PLAIN, 1.0, 0.0, 0.0, 1.0 );
    id.MatrixSquared();
    _equal( id(0,0), 1.0, fTolerance ); _equal( id(0,1), 0.0, fTolerance );
    _equal( id(1,0), 0.0, fTolerance ); _equal( id(1,1), 1.0, fTolerance );
}

// ============================================================================

void TensorVariable_Test2::DoubleContraction_Function2()
{
    // T:T = sum_ij T_ij^2.
    // For T = [[1,2],[3,4]]: T:T = 1+4+9+16 = 30.
    TensorVariable<2U> t( PLAIN, 1.0, 2.0, 3.0, 4.0 );
    _equal( t.DoubleContraction(), 30.0, fTolerance );

    // For identity: T:T = 1+0+0+1 = 2.
    TensorVariable<2U> id( PLAIN, 1.0, 0.0, 0.0, 1.0 );
    _equal( id.DoubleContraction(), 2.0, fTolerance );

    // For diagonal tensor diag(2,3): T:T = 4+0+0+9 = 13.
    TensorVariable<2U> diag( PLAIN, 2.0, 0.0, 0.0, 3.0 );
    _equal( diag.DoubleContraction(), 13.0, fTolerance );
}

} // namespace csmp

