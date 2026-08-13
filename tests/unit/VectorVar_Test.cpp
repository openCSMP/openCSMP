#include "VectorVar_Test.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"

using std::vector;
using std::cout;
using std::endl;

namespace csmp {

VectorVariable_Test::VectorVariable_Test()
{
    fTolerance = 1.0e-6;
}

VectorVariable_Test::~VectorVariable_Test()
{
}

void VectorVariable_Test::run()
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
    LessThan_Operator();
    Pow_Function();
    Length_Function();
    EuclideanNormalize_Function();
    DotProduct_Function();
    CrossProduct_Function();
    IsWithinRange_Function();
    Flip_Function();
    AngleTo_Function();
    ProjectOnto_Function();
    Invert_Function();
}

// ============================================================================

void VectorVariable_Test::Assignment_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2;
    ScalarVariable     sc( DIRICH, 10.0 );

    // VectorVariable assignment — values and flags copied.
    v2 = v1;
    _equal( v2(0), 1.0, fTolerance );
    _equal( v2(1), 2.0, fTolerance );
    _equal( v2(2), 3.0, fTolerance );
    _test(  v2.Flag(0) == ANY   );
    _test(  v2.Flag(1) == PLAIN );
    _test(  v2.Flag(2) == ROBIN );

    // double assignment — values set, flags unchanged.
    v2 = 5.0;
    _equal( v2(0), 5.0, fTolerance );
    _equal( v2(1), 5.0, fTolerance );
    _equal( v2(2), 5.0, fTolerance );
    _test(  v2.Flag(0) == ANY   );
    _test(  v2.Flag(1) == PLAIN );
    _test(  v2.Flag(2) == ROBIN );

    // ScalarVariable assignment — values and flag copied.
    v2 = sc;
    _equal( v2(0), 10.0, fTolerance );
    _equal( v2(1), 10.0, fTolerance );
    _equal( v2(2), 10.0, fTolerance );
    _test(  v2.Flag(0) == DIRICH );
    _test(  v2.Flag(1) == DIRICH );
    _test(  v2.Flag(2) == DIRICH );

    // Point assignment — values set, flags unchanged.
    vector<double> pts = { 15.0, 15.0, 15.0 };
    Point<3U> p( pts );
    v2 = p;
    _equal( v2(0), 15.0, fTolerance );
    _equal( v2(1), 15.0, fTolerance );
    _equal( v2(2), 15.0, fTolerance );
    _test(  v2.Flag(0) == DIRICH );
    _test(  v2.Flag(1) == DIRICH );
    _test(  v2.Flag(2) == DIRICH );
}

// ============================================================================

void VectorVariable_Test::Addition_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;

    // VectorVariable + VectorVariable — flags from lhs.
    v3 = v1 + v2;
    _equal( v3(0), 5.0, fTolerance );
    _equal( v3(1), 7.0, fTolerance );
    _equal( v3(2), 9.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    // VectorVariable + double — flags from lhs.
    v3 = v1 + 10.0;
    _equal( v3(0), 11.0, fTolerance );
    _equal( v3(1), 12.0, fTolerance );
    _equal( v3(2), 13.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    // Point + VectorVariable.
    vector<double> pts = { 15.0, 15.0, 15.0 };
    Point<3U> p( pts );
    Point<3U> result = p + v1;
    _equal( result[0], 16.0, fTolerance );
    _equal( result[1], 17.0, fTolerance );
    _equal( result[2], 18.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::Subtraction_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;

    v3 = v1 - v2;
    _equal( v3(0), -3.0, fTolerance );
    _equal( v3(1), -3.0, fTolerance );
    _equal( v3(2), -3.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1 - 10.0;
    _equal( v3(0), -9.0, fTolerance );
    _equal( v3(1), -8.0, fTolerance );
    _equal( v3(2), -7.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    vector<double> pts = { 15.0, 15.0, 15.0 };
    Point<3U> p( pts );
    Point<3U> result = p - v1;
    _equal( result[0], 14.0, fTolerance );
    _equal( result[1], 13.0, fTolerance );
    _equal( result[2], 12.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::Multiplication_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;

    v3 = v1 * v2;
    _equal( v3(0),  4.0, fTolerance );
    _equal( v3(1), 10.0, fTolerance );
    _equal( v3(2), 18.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1 * 10.0;
    _equal( v3(0), 10.0, fTolerance );
    _equal( v3(1), 20.0, fTolerance );
    _equal( v3(2), 30.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    vector<double> pts = { 5.0, 5.0, 5.0 };
    Point<3U> p( pts );
    Point<3U> result = p * v1;
    _equal( result[0],  5.0, fTolerance );
    _equal( result[1], 10.0, fTolerance );
    _equal( result[2], 15.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::Division_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;

    v3 = v1 / v2;
    _equal( v3(0), 0.25, fTolerance );
    _equal( v3(1), 0.4,  fTolerance );
    _equal( v3(2), 0.5,  fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1 / 10.0;
    _equal( v3(0), 0.1, fTolerance );
    _equal( v3(1), 0.2, fTolerance );
    _equal( v3(2), 0.3, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    vector<double> pts = { 2.0, 2.0, 3.0 };
    Point<3U> p( pts );
    Point<3U> result = p / v1;
    _equal( result[0], 2.0, fTolerance );
    _equal( result[1], 1.0, fTolerance );
    _equal( result[2], 1.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::Addition_Assignment_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 += v2;
    _equal( v3(0), 5.0, fTolerance );
    _equal( v3(1), 7.0, fTolerance );
    _equal( v3(2), 9.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 += sc;
    _equal( v3(0), 6.0, fTolerance );
    _equal( v3(1), 7.0, fTolerance );
    _equal( v3(2), 8.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 += 10.0;
    _equal( v3(0), 11.0, fTolerance );
    _equal( v3(1), 12.0, fTolerance );
    _equal( v3(2), 13.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );
}

// ============================================================================

void VectorVariable_Test::Subtraction_Assignment_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 -= v2;
    _equal( v3(0), -3.0, fTolerance );
    _equal( v3(1), -3.0, fTolerance );
    _equal( v3(2), -3.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 -= sc;
    _equal( v3(0), -4.0, fTolerance );
    _equal( v3(1), -3.0, fTolerance );
    _equal( v3(2), -2.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 -= 10.0;
    _equal( v3(0), -9.0, fTolerance );
    _equal( v3(1), -8.0, fTolerance );
    _equal( v3(2), -7.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );
}

// ============================================================================

void VectorVariable_Test::Multiplication_Assignment_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 *= v2;
    _equal( v3(0),  4.0, fTolerance );
    _equal( v3(1), 10.0, fTolerance );
    _equal( v3(2), 18.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 *= sc;
    _equal( v3(0),  5.0, fTolerance );
    _equal( v3(1), 10.0, fTolerance );
    _equal( v3(2), 15.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 *= 10.0;
    _equal( v3(0), 10.0, fTolerance );
    _equal( v3(1), 20.0, fTolerance );
    _equal( v3(2), 30.0, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );
}

// ============================================================================

void VectorVariable_Test::Division_Assignment_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY,  4.0, 5.0, 6.0 );
    VectorVariable<3U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 /= v2;
    _equal( v3(0), 0.25, fTolerance );
    _equal( v3(1), 0.4,  fTolerance );
    _equal( v3(2), 0.5,  fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 /= sc;
    _equal( v3(0), 0.2, fTolerance );
    _equal( v3(1), 0.4, fTolerance );
    _equal( v3(2), 0.6, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );

    v3 = v1; v3 /= 10.0;
    _equal( v3(0), 0.1, fTolerance );
    _equal( v3(1), 0.2, fTolerance );
    _equal( v3(2), 0.3, fTolerance );
    _test(  v3.Flag(0) == ANY   );
    _test(  v3.Flag(1) == PLAIN );
    _test(  v3.Flag(2) == ROBIN );
}

// ============================================================================

void VectorVariable_Test::Equality_Operator()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );  // same
    VectorVariable<3U> v3( ANY, PLAIN, ROBIN, 4.0, 5.0, 6.0 );  // different values
    VectorVariable<3U> v4( ROBIN, PLAIN, ANY, 1.0, 2.0, 3.0 );  // different flags

    _test(  v1 == v2 );
    _test( !( v1 == v3 ) );
    _test( !( v1 != v2 ) );
    _test(  v1 != v3 );
    _test(  v1 != v4 );  // different flags — was missing
}

// ============================================================================

void VectorVariable_Test::LessThan_Operator()
{
    // operator< compares pointer addresses — satisfies strict weak ordering
    // for use in STL associative containers. It does not compare values.
    // We verify it is consistent: a < b implies !(b < a).
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> v2( ANY, PLAIN, ROBIN, 4.0, 5.0, 6.0 );

    // Exactly one of (v1 < v2) or (v2 < v1) must be true,
    // and (v1 < v1) must be false (irreflexivity).
    const bool lt = ( v1 < v2 );
    const bool gt = ( v2 < v1 );
    _test(  lt != gt );       // exactly one is true
    _test( !( v1 < v1 ) );   // irreflexive
    _test( !( v2 < v2 ) );   // irreflexive
}

// ============================================================================

void VectorVariable_Test::Pow_Function()
{
    // Replaces the missing test for Pow(double) — operator^ was removed.
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 2.0, 3.0, 4.0 );

    // Pow(2): [4, 9, 16]
    VectorVariable<3U> result = v1.Pow( 2.0 );
    _equal( result(0),  4.0, fTolerance );
    _equal( result(1),  9.0, fTolerance );
    _equal( result(2), 16.0, fTolerance );
    _test(  result.Flag(0) == ANY   );
    _test(  result.Flag(1) == PLAIN );
    _test(  result.Flag(2) == ROBIN );

    // Pow(0.5) = sqrt — round-trip.
    VectorVariable<3U> roundtrip = result.Pow( 0.5 );
    _equal( roundtrip(0), 2.0, fTolerance );
    _equal( roundtrip(1), 3.0, fTolerance );
    _equal( roundtrip(2), 4.0, fTolerance );

    // Pow(3): [8, 27, 64]
    VectorVariable<3U> cubed = v1.Pow( 3.0 );
    _equal( cubed(0),  8.0, fTolerance );
    _equal( cubed(1), 27.0, fTolerance );
    _equal( cubed(2), 64.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::Length_Function()
{
    // ||(1,2,3)|| = sqrt(14)
    VectorVariable<3U> v( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    _equal( v.Length(), std::sqrt(14.0), fTolerance );

    // Known Pythagorean triple in 3D: ||(1,2,2)|| = 3
    VectorVariable<3U> v2( ANY, ANY, ANY, 1.0, 2.0, 2.0 );
    _equal( v2.Length(), 3.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::EuclideanNormalize_Function()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    const double norm = std::sqrt(14.0);
    VectorVariable<3U> expected( ANY, PLAIN, ROBIN,
        1.0/norm, 2.0/norm, 3.0/norm );

    v1.EuclideanNormalize();
    _test( v1 == expected );
    _test( v1.Flag(0) == ANY   );
    _test( v1.Flag(1) == PLAIN );
    _test( v1.Flag(2) == ROBIN );

    // Normalised vector must have unit length.
    _equal( v1.Length(), 1.0, fTolerance );

    // Zero vector — no-op.
    VectorVariable<3U> zero( ANY, ANY, ANY, 0.0, 0.0, 0.0 );
    zero.EuclideanNormalize();
    _equal( zero(0), 0.0, fTolerance );
    _equal( zero(1), 0.0, fTolerance );
    _equal( zero(2), 0.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::DotProduct_Function()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );

    // DotProduct with Point: 1*2 + 2*4 + 3*6 = 28.
    vector<double> pts = { 2.0, 4.0, 6.0 };
    Point<3U> p( pts );
    _equal( v1.DotProduct( p ), 28.0, fTolerance );

    // DotProduct with VectorVariable: 1*1 + 2*2 + 3*3 = 14.
    VectorVariable<3U> v2( ANY, ANY, ANY, 1.0, 2.0, 3.0 );
    _equal( v1.DotProduct( v2 ), 14.0, fTolerance );

    // Orthogonal vectors: dot product = 0.
    VectorVariable<3U> vx( ANY, ANY, ANY, 1.0, 0.0, 0.0 );
    VectorVariable<3U> vy( ANY, ANY, ANY, 0.0, 1.0, 0.0 );
    _equal( vx.DotProduct( vy ), 0.0, fTolerance );

    // operator& is the same as DotProduct.
    _equal( v1 & v2, 14.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::CrossProduct_Function()
{
    // v1=(1,2,3), p=(3,1,5)
    // cross = (2*5-3*1, 3*3-1*5, 1*1-2*3) = (7, 4, -5)
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> expected( ANY, PLAIN, ROBIN, 7.0, 4.0, -5.0 );

    vector<double> pts = { 3.0, 1.0, 5.0 };
    Point<3U> p( pts );

    // CrossProduct with Point.
    VectorVariable<3U> result = v1.CrossProduct( p );
    _test( result == expected );
    _test( result.Flag(0) == ANY   );
    _test( result.Flag(1) == PLAIN );
    _test( result.Flag(2) == ROBIN );

    // CrossProduct with VectorVariable — same result.
    VectorVariable<3U> v2( ANY, ANY, ANY, 3.0, 1.0, 5.0 );
    VectorVariable<3U> result2 = v1.CrossProduct( v2 );
    _equal( result2(0),  7.0, fTolerance );
    _equal( result2(1),  4.0, fTolerance );
    _equal( result2(2), -5.0, fTolerance );

    // operator% is the same as CrossProduct.
    VectorVariable<3U> result3 = v1 % v2;
    _equal( result3(0),  7.0, fTolerance );
    _equal( result3(1),  4.0, fTolerance );
    _equal( result3(2), -5.0, fTolerance );

    // Cross product of parallel vectors is zero.
    VectorVariable<3U> vx( ANY, ANY, ANY, 1.0, 0.0, 0.0 );
    VectorVariable<3U> vx2( ANY, ANY, ANY, 2.0, 0.0, 0.0 );
    VectorVariable<3U> zero = vx.CrossProduct( vx2 );
    _equal( zero(0), 0.0, fTolerance );
    _equal( zero(1), 0.0, fTolerance );
    _equal( zero(2), 0.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test::IsWithinRange_Function()
{
    VectorVariable<3U> v( ANY, PLAIN, ROBIN, 1.0, 5.0, 9.0 );

    _test( v.IsWithinRange( -5.0,  5.0 ) == false );  // 9 > 5
    _test( v.IsWithinRange(  5.0, 15.0 ) == false );  // 1 < 5
    _test( v.IsWithinRange(  1.0,  9.0 ) == true  );  // exact bounds
    _test( v.IsWithinRange(  4.0,  6.0 ) == false );  // 1 and 9 out
    _test( v.IsWithinRange( -1.0, 10.0 ) == true  );  // all in
}

// ============================================================================

void VectorVariable_Test::Flip_Function()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN, 1.0, 2.0, 3.0 );
    VectorVariable<3U> expected( ROBIN, PLAIN, ANY, 3.0, 2.0, 1.0 );

    VectorVariable<3U> result = v1.Flip();
    _test( result == expected );
    _test( result.Flag(0) == ROBIN );
    _test( result.Flag(1) == PLAIN );
    _test( result.Flag(2) == ANY   );

    // Original must be unchanged.
    _equal( v1(0), 1.0, fTolerance );
    _equal( v1(1), 2.0, fTolerance );
    _equal( v1(2), 3.0, fTolerance );
    _test(  v1.Flag(0) == ANY   );
    _test(  v1.Flag(1) == PLAIN );
    _test(  v1.Flag(2) == ROBIN );
}

// ============================================================================

void VectorVariable_Test::AngleTo_Function()
{
    VectorVariable<3U> v1( ANY, PLAIN, ROBIN,  1.0,  1.0,  1.0 );
    VectorVariable<3U> v2( ROBIN, PLAIN, ANY, -1.0, -1.0, -1.0 );
    VectorVariable<3U> v3( ANY, PLAIN, ROBIN,  1.0,  1.0,  1.0 );
    VectorVariable<3U> v4( ROBIN, PLAIN, ANY,  1.0,  1.0,  1.0 );
    VectorVariable<3U> v5( ANY, PLAIN, ROBIN,  1.0,  1.0,  0.0 );
    VectorVariable<3U> v6( ROBIN, PLAIN, ANY,  0.0,  0.0,  1.0 );
    VectorVariable<3U> v7( ANY, PLAIN, ROBIN,  2.0,  1.0,  5.0 );
    VectorVariable<3U> v8( ROBIN, PLAIN, ANY, -1.0,  4.0, -2.0 );

    _equal( v1.AngleTo( v2 ), 180.0, fTolerance );
    _equal( v3.AngleTo( v4 ),   0.0, fTolerance );
    _equal( v5.AngleTo( v6 ),  90.0, fTolerance );

    // v7.v8 = -2+4-10 = -8
    // |v7| = sqrt(4+1+25) = sqrt(30)
    // |v8| = sqrt(1+16+4) = sqrt(21)
    // angle = acos(-8 / sqrt(630)) in degrees
    const double expected_angle =
        ( 180.0 / 3.14159265358979323846 ) *
        std::acos( -8.0 / std::sqrt( 30.0 * 21.0 ) );
    _equal( v7.AngleTo( v8 ), expected_angle, fTolerance );

    // Angle is symmetric.
    _equal( v5.AngleTo( v6 ), v6.AngleTo( v5 ), fTolerance );
}

// ============================================================================

void VectorVariable_Test::ProjectOnto_Function()
{
    // v = (2,1,5), u = (-1,4,-2)
    // ratio = v.u / u.u = (-2+4-10) / (1+16+4) = -8/21
    // result = (-8/21) * (-1,4,-2) = (8/21, -32/21, 16/21)
    VectorVariable<3U> v( DIRICH, PLAIN, ROBIN, 2.0, 1.0, 5.0 );
    VectorVariable<3U> u( ROBIN,  PLAIN, DIRICH, -1.0, 4.0, -2.0 );

    const double expected_x =  8.0 / 21.0;
    const double expected_y = -32.0 / 21.0;
    const double expected_z =  16.0 / 21.0;

    // ProjectOnto VectorVariable.
    VectorVariable<3U> result = v.ProjectOnto( u );
    _equal( result(0), expected_x, fTolerance );
    _equal( result(1), expected_y, fTolerance );
    _equal( result(2), expected_z, fTolerance );
    _test(  result.Flag(0) == DIRICH );
    _test(  result.Flag(1) == PLAIN  );
    _test(  result.Flag(2) == ROBIN  );

    // ProjectOnto std::vector — same result.
    vector<double> uv = { -1.0, 4.0, -2.0 };
    VectorVariable<3U> result2 = v.ProjectOnto( uv );
    _equal( result2(0), expected_x, fTolerance );
    _equal( result2(1), expected_y, fTolerance );
    _equal( result2(2), expected_z, fTolerance );
    _test(  result2.Flag(0) == DIRICH );
    _test(  result2.Flag(1) == PLAIN  );
    _test(  result2.Flag(2) == ROBIN  );
}

// ============================================================================

void VectorVariable_Test::Invert_Function()
{
    VectorVariable<3U> v( ANY, PLAIN, ROBIN, 1.0, -2.0, 3.0 );

    v.Invert();
    _equal( v(0), -1.0, fTolerance );
    _equal( v(1),  2.0, fTolerance );
    _equal( v(2), -3.0, fTolerance );

    // Flags must be unchanged.
    _test( v.Flag(0) == ANY   );
    _test( v.Flag(1) == PLAIN );
    _test( v.Flag(2) == ROBIN );

    // Double invert returns to original.
    v.Invert();
    _equal( v(0),  1.0, fTolerance );
    _equal( v(1), -2.0, fTolerance );
    _equal( v(2),  3.0, fTolerance );
}

} // namespace csmp

