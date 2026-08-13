#include "VectorVar_Test2.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace csmp {

VectorVariable_Test2::VectorVariable_Test2()
{
    fTolerance = 1.0e-6;
}

VectorVariable_Test2::~VectorVariable_Test2()
{
}

void VectorVariable_Test2::run()
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
    Pow_Function2();          // replaces Power_Operator2
    Length_Function2();
    EuclideanNormalize_Function2();
    DotProduct_Function2();
    CrossProduct_Function2();
    IsWithinRange_Function2();
    Flip_Function2();
    AngleTo_Function2();
    ProjectOnto_Function2();
    Invert_Function2();       // was missing
}

// ============================================================================

void VectorVariable_Test2::Assignment_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2;
    ScalarVariable     sc( DIRICH, 10.0 );

    // VectorVariable assignment — values and flags copied.
    v2 = v1;
    _equal( v2(0), 1.0, fTolerance );
    _equal( v2(1), 2.0, fTolerance );
    _test(  v2.Flag(0) == PLAIN );
    _test(  v2.Flag(1) == ROBIN );

    // double assignment — values set, flags unchanged.
    v2 = 5.0;
    _equal( v2(0), 5.0, fTolerance );
    _equal( v2(1), 5.0, fTolerance );
    _test(  v2.Flag(0) == PLAIN );
    _test(  v2.Flag(1) == ROBIN );

    // ScalarVariable assignment — values and flag copied.
    v2 = sc;
    _equal( v2(0), 10.0, fTolerance );
    _equal( v2(1), 10.0, fTolerance );
    _test(  v2.Flag(0) == DIRICH );
    _test(  v2.Flag(1) == DIRICH );

    // Point assignment — values set, flags unchanged.
    vector<double> pts = { 15.0, 15.0 };
    Point<2U> p( pts );
    v2 = p;
    _equal( v2(0), 15.0, fTolerance );
    _equal( v2(1), 15.0, fTolerance );
    _test(  v2.Flag(0) == DIRICH );
    _test(  v2.Flag(1) == DIRICH );
}

// ============================================================================

void VectorVariable_Test2::Addition_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 6.0 );
    VectorVariable<2U> v3;

    // VectorVariable + VectorVariable — flags from lhs.
    v3 = v1 + v2;
    _equal( v3(0), 6.0,  fTolerance );
    _equal( v3(1), 8.0,  fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    // VectorVariable + double — flags from lhs.
    v3 = v1 + 10.0;
    _equal( v3(0), 11.0, fTolerance );
    _equal( v3(1), 12.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    // Point + VectorVariable.
    vector<double> pts = { 15.0, 15.0 };
    Point<2U> p( pts );
    Point<2U> result = p + v1;
    _equal( result[0], 16.0, fTolerance );
    _equal( result[1], 17.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::Subtraction_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 6.0 );
    VectorVariable<2U> v3;

    v3 = v1 - v2;
    _equal( v3(0), -4.0, fTolerance );
    _equal( v3(1), -4.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1 - 10.0;
    _equal( v3(0), -9.0, fTolerance );
    _equal( v3(1), -8.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    vector<double> pts = { 15.0, 15.0 };
    Point<2U> p( pts );
    Point<2U> result = p - v1;
    _equal( result[0], 14.0, fTolerance );
    _equal( result[1], 13.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::Multiplication_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 6.0 );
    VectorVariable<2U> v3;

    v3 = v1 * v2;
    _equal( v3(0),  5.0, fTolerance );
    _equal( v3(1), 12.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1 * 10.0;
    _equal( v3(0), 10.0, fTolerance );
    _equal( v3(1), 20.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    vector<double> pts = { 5.0, 5.0 };
    Point<2U> p( pts );
    Point<2U> result = p * v1;
    _equal( result[0],  5.0, fTolerance );
    _equal( result[1], 10.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::Division_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 4.0 );
    VectorVariable<2U> v3;

    v3 = v1 / v2;
    _equal( v3(0), 0.2, fTolerance );
    _equal( v3(1), 0.5, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1 / 10.0;
    _equal( v3(0), 0.1, fTolerance );
    _equal( v3(1), 0.2, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    vector<double> pts = { 2.0, 2.0 };
    Point<2U> p( pts );
    Point<2U> result = p / v1;
    _equal( result[0], 2.0, fTolerance );
    _equal( result[1], 1.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::Addition_Assignment_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 6.0 );
    VectorVariable<2U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 += v2;
    _equal( v3(0),  6.0, fTolerance );
    _equal( v3(1),  8.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 += sc;
    _equal( v3(0),  6.0, fTolerance );
    _equal( v3(1),  7.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 += 10.0;
    _equal( v3(0), 11.0, fTolerance );
    _equal( v3(1), 12.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );
}

// ============================================================================

void VectorVariable_Test2::Subtraction_Assignment_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 6.0 );
    VectorVariable<2U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 -= v2;
    _equal( v3(0), -4.0, fTolerance );
    _equal( v3(1), -4.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 -= sc;
    _equal( v3(0), -4.0, fTolerance );
    _equal( v3(1), -3.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 -= 10.0;
    _equal( v3(0), -9.0, fTolerance );
    _equal( v3(1), -8.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );
}

// ============================================================================

void VectorVariable_Test2::Multiplication_Assignment_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 6.0 );
    VectorVariable<2U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 *= v2;
    _equal( v3(0),  5.0, fTolerance );
    _equal( v3(1), 12.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 *= sc;
    _equal( v3(0),  5.0, fTolerance );
    _equal( v3(1), 10.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 *= 10.0;
    _equal( v3(0), 10.0, fTolerance );
    _equal( v3(1), 20.0, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );
}

// ============================================================================

void VectorVariable_Test2::Division_Assignment_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 5.0, 4.0 );
    VectorVariable<2U> v3;
    ScalarVariable     sc( DIRICH, 5.0 );

    v3 = v1; v3 /= v2;
    _equal( v3(0), 0.2, fTolerance );
    _equal( v3(1), 0.5, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 /= sc;
    _equal( v3(0), 0.2, fTolerance );
    _equal( v3(1), 0.4, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );

    v3 = v1; v3 /= 10.0;
    _equal( v3(0), 0.1, fTolerance );
    _equal( v3(1), 0.2, fTolerance );
    _test(  v3.Flag(0) == PLAIN );
    _test(  v3.Flag(1) == ROBIN );
}

// ============================================================================

void VectorVariable_Test2::Equality_Operator2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> v2( ROBIN, PLAIN, 1.0, 2.0 );  // different flags
    VectorVariable<2U> v3( PLAIN, ROBIN, 1.0, 2.0 );  // same as v1
    VectorVariable<2U> v4( PLAIN, ROBIN, 5.0, 4.0 );  // different values

    _test( v1 != v2 );  // different flags
    _test( v1 != v4 );  // different values
    _test( v1 == v3 );  // identical
}

// ============================================================================

void VectorVariable_Test2::Pow_Function2()
{
    // Replaces Power_Operator2 which used the removed operator^.
    VectorVariable<2U> v1( PLAIN, ROBIN, 3.0, 2.0 );
    VectorVariable<2U> v2( PLAIN, ROBIN, 9.0, 4.0 );

    // Pow(2) should give [9, 4].
    VectorVariable<2U> result = v1.Pow( 2.0 );
    _test( result == v2 );
    _test( result.Flag(0) == PLAIN );
    _test( result.Flag(1) == ROBIN );

    // Pow(0.5) = sqrt: sqrt(9)=3, sqrt(4)=2 — round-trip.
    VectorVariable<2U> roundtrip = result.Pow( 0.5 );
    _equal( roundtrip(0), 3.0, fTolerance );
    _equal( roundtrip(1), 2.0, fTolerance );

    // Pow(3): 3^3=27, 2^3=8.
    VectorVariable<2U> cubed = v1.Pow( 3.0 );
    _equal( cubed(0), 27.0, fTolerance );
    _equal( cubed(1),  8.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::Length_Function2()
{
    // 3-4-5 right triangle.
    VectorVariable<2U> v( PLAIN, ROBIN, 3.0, 4.0 );
    _equal( v.Length(), 5.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::EuclideanNormalize_Function2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 3.0, 4.0 );
    VectorVariable<2U> v2( PLAIN, ROBIN, 0.6, 0.8 );

    v1.EuclideanNormalize();
    _test( v1 == v2 );
    _test( v1.Flag(0) == PLAIN );
    _test( v1.Flag(1) == ROBIN );

    // Normalised vector should have unit length.
    _equal( v1.Length(), 1.0, fTolerance );

    // Zero vector — no-op.
    VectorVariable<2U> zero( ANY, ANY, 0.0, 0.0 );
    zero.EuclideanNormalize();
    _equal( zero(0), 0.0, fTolerance );
    _equal( zero(1), 0.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::DotProduct_Function2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );

    // DotProduct with Point: 1*2 + 2*4 = 10.
    vector<double> pts = { 2.0, 4.0 };
    Point<2U> p( pts );
    _equal( v1.DotProduct( p ), 10.0, fTolerance );

    // DotProduct with VectorVariable: 1*3 + 2*5 = 13.
    VectorVariable<2U> v2( ANY, ANY, 3.0, 5.0 );
    _equal( v1.DotProduct( v2 ), 13.0, fTolerance );

    // Orthogonal vectors: dot product = 0.
    VectorVariable<2U> vx( ANY, ANY, 1.0, 0.0 );
    VectorVariable<2U> vy( ANY, ANY, 0.0, 1.0 );
    _equal( vx.DotProduct( vy ), 0.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::CrossProduct_Function2()
{
    // 2D cross product: result = (0, v.x*p.y - v.y*p.x)
    // v1=(1,2), p=(3,1): y = 1*1 - 2*3 = -5
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> expected( PLAIN, ROBIN, 0.0, -5.0 );

    vector<double> pts = { 3.0, 1.0 };
    Point<2U> p( pts );

    VectorVariable<2U> result = v1.CrossProduct( p );
    _test( result == expected );
    _test( result.Flag(0) == PLAIN );
    _test( result.Flag(1) == ROBIN );

    // CrossProduct with VectorVariable: same result.
    VectorVariable<2U> v2( ANY, ANY, 3.0, 1.0 );
    VectorVariable<2U> result2 = v1.CrossProduct( v2 );
    _equal( result2(0), 0.0,  fTolerance );
    _equal( result2(1), -5.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::IsWithinRange_Function2()
{
    VectorVariable<2U> v( PLAIN, ROBIN, 1.0, 9.0 );

    _test( v.IsWithinRange( -5.0,  5.0 ) == false );  // 9 > 5
    _test( v.IsWithinRange(  5.0, 15.0 ) == false );  // 1 < 5
    _test( v.IsWithinRange(  1.0,  9.0 ) == true  );  // exact bounds
    _test( v.IsWithinRange(  4.0,  6.0 ) == false );  // both out
    _test( v.IsWithinRange( -1.0, 10.0 ) == true  );  // both in
}

// ============================================================================

void VectorVariable_Test2::Flip_Function2()
{
    VectorVariable<2U> v1( PLAIN, ROBIN, 1.0, 2.0 );
    VectorVariable<2U> expected( ROBIN, PLAIN, 2.0, 1.0 );

    VectorVariable<2U> result = v1.Flip();
    _test( result == expected );
    _test( result.Flag(0) == ROBIN );
    _test( result.Flag(1) == PLAIN );

    // Original must be unchanged.
    _equal( v1(0), 1.0, fTolerance );
    _equal( v1(1), 2.0, fTolerance );
    _test(  v1.Flag(0) == PLAIN );
    _test(  v1.Flag(1) == ROBIN );
}

// ============================================================================

void VectorVariable_Test2::AngleTo_Function2()
{
    VectorVariable<2U> vx(  PLAIN, ROBIN,  1.0,  0.0 );
    VectorVariable<2U> v45( PLAIN, ROBIN,  1.0,  1.0 );
    VectorVariable<2U> vy(  PLAIN, ROBIN,  0.0,  1.0 );
    VectorVariable<2U> vnx( PLAIN, ROBIN, -1.0,  0.0 );
    VectorVariable<2U> vny( PLAIN, ROBIN,  0.0, -1.0 );

    _equal( vx.AngleTo( v45 ),  45.0, fTolerance );
    _equal( vx.AngleTo( vy  ),  90.0, fTolerance );
    _equal( vx.AngleTo( vnx ), 180.0, fTolerance );
    _equal( vx.AngleTo( vny ),  90.0, fTolerance );

    // Angle is symmetric.
    _equal( vy.AngleTo( vx ), 90.0, fTolerance );
}

// ============================================================================

void VectorVariable_Test2::ProjectOnto_Function2()
{
    // v = (2,1) projected onto u = (-1,4).
    // ratio = (v.u) / (u.u) = (-2+4) / (1+16) = 2/17
    // result = ratio * u = (2/17)*(-1,4) = (-2/17, 8/17)
    VectorVariable<2U> v( DIRICH, PLAIN, 2.0, 1.0 );
    VectorVariable<2U> u( ROBIN,  PLAIN, -1.0, 4.0 );

    const double expected_x = -2.0 / 17.0;
    const double expected_y =  8.0 / 17.0;

    // ProjectOnto VectorVariable.
    VectorVariable<2U> result = v.ProjectOnto( u );
    _equal( result(0), expected_x, fTolerance );
    _equal( result(1), expected_y, fTolerance );
    _test(  result.Flag(0) == DIRICH );
    _test(  result.Flag(1) == PLAIN  );

    // ProjectOnto std::vector — same result.
    vector<double> uv = { -1.0, 4.0 };
    VectorVariable<2U> result2 = v.ProjectOnto( uv );
    _equal( result2(0), expected_x, fTolerance );
    _equal( result2(1), expected_y, fTolerance );
    _test(  result2.Flag(0) == DIRICH );
    _test(  result2.Flag(1) == PLAIN  );
}

// ============================================================================

void VectorVariable_Test2::Invert_Function2()
{
    VectorVariable<2U> v( PLAIN, ROBIN, 3.0, -5.0 );

    v.Invert();
    _equal( v(0), -3.0, fTolerance );
    _equal( v(1),  5.0, fTolerance );

    // Flags must be unchanged.
    _test( v.Flag(0) == PLAIN );
    _test( v.Flag(1) == ROBIN );

    // Double invert returns to original.
    v.Invert();
    _equal( v(0),  3.0, fTolerance );
    _equal( v(1), -5.0, fTolerance );
}

} // namespace csmp

