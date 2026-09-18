#include "ScalarVar_Test.h"
#include "ScalarVariable.h"

#include <cmath>
#include <limits>
#include <iostream>

using std::cout;
using std::endl;

namespace csmp {

ScalarVariable_Test::ScalarVariable_Test()
{
    fTolerance = 1.0e-14;
}

ScalarVariable_Test::~ScalarVariable_Test()
{
}

void ScalarVariable_Test::run()
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
    Less_Than_Operator();
    Greater_Than_Operator();
    Less_Than_Or_Equal_To_Operator();
    Greater_Than_Or_Equal_To_Operator();
    Equality_Operator();
    IsWithinRange_Function();
    HasNaN_Function();
}

// ============================================================================

void ScalarVariable_Test::Assignment_Operator()
{
    // Default constructor — value should be NaN, flag ANY.
    ScalarVariable sc_default;
    _test( std::isnan( sc_default() ) );
    _test( sc_default.Flag() == ANY );

    ScalarVariable sc1( ANY, 5.0 );
    ScalarVariable sc2;

    // Scalar-to-scalar assignment — value and flag copied.
    sc2 = sc1;
    _equal( sc2(), sc1(), fTolerance );
    _test(  sc2.Flag() == sc1.Flag() );
    _test(  sc2.Flag() == ANY );

    // double assignment — value set, flag unchanged.
    sc2 = 99.0;
    _equal( sc2(), 99.0, fTolerance );
    _test(  sc2.Flag() == ANY );  // flag must not change

    // Assign different flag scalar — flag is also copied.
    ScalarVariable sc3( DIRICH, 7.0 );
    sc2 = sc3;
    _equal( sc2(), 7.0, fTolerance );
    _test(  sc2.Flag() == DIRICH );
}

// ============================================================================

void ScalarVariable_Test::Addition_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 15.0 );
    ScalarVariable sc3;

    // Scalar + double — flag from lhs.
    sc3 = sc1 + 10.0;
    _equal( sc3(), 15.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // double + scalar — flag from rhs.
    sc3 = 10.0 + sc1;
    _equal( sc3(), 15.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Scalar + scalar — flag from lhs.
    sc3 = sc1 + sc2;
    _equal( sc3(), 20.0, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Subtraction_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 15.0 );
    ScalarVariable sc3;

    // Scalar - double.
    sc3 = sc1 - 10.0;
    _equal( sc3(), -5.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // double - scalar.
    sc3 = 10.0 - sc1;
    _equal( sc3(), 5.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Scalar - scalar.
    sc3 = sc1 - sc2;
    _equal( sc3(), -10.0, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Multiplication_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 15.0 );
    ScalarVariable sc3;

    // Scalar * double.
    sc3 = sc1 * 10.0;
    _equal( sc3(), 50.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // double * scalar.
    sc3 = 10.0 * sc1;
    _equal( sc3(), 50.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Scalar * scalar.
    sc3 = sc1 * sc2;
    _equal( sc3(), 75.0, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Division_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3;

    // Scalar / double.
    sc3 = sc1 / 10.0;
    _equal( sc3(), 0.5, fTolerance );
    _test(  sc3.Flag() == ANY );

    // double / scalar.
    sc3 = 10.0 / sc1;
    _equal( sc3(), 2.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Scalar / scalar.
    sc3 = sc1 / sc2;
    _equal( sc3(), 0.25, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Addition_Assignment_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 15.0 );
    ScalarVariable sc3;

    // Compound += double via raw accessor (tests double& access).
    sc3 = sc1;
    sc3() += 10.0;
    _equal( sc3(), 15.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Compound += scalar.
    sc3 = sc1;
    sc3 += sc2;
    _equal( sc3(), 20.0, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Subtraction_Assignment_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 15.0 );
    ScalarVariable sc3;

    // Compound -= double via raw accessor.
    sc3 = sc1;
    sc3() -= 10.0;
    _equal( sc3(), -5.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Compound -= scalar.
    sc3 = sc1;
    sc3 -= sc2;
    _equal( sc3(), -10.0, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Multiplication_Assignment_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 15.0 );
    ScalarVariable sc3;

    // Compound *= double via raw accessor.
    sc3 = sc1;
    sc3() *= 10.0;
    _equal( sc3(), 50.0, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Compound *= scalar.
    sc3 = sc1;
    sc3 *= sc2;
    _equal( sc3(), 75.0, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Division_Assignment_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3;

    // Compound /= double via raw accessor.
    sc3 = sc1;
    sc3() /= 10.0;
    _equal( sc3(), 0.5, fTolerance );
    _test(  sc3.Flag() == ANY );

    // Compound /= scalar.
    sc3 = sc1;
    sc3 /= sc2;
    _equal( sc3(), 0.25, fTolerance );
    _test(  sc3.Flag() == ANY );
}

// ============================================================================

void ScalarVariable_Test::Less_Than_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3( PLAIN,  5.0 );
    ScalarVariable sc4( PLAIN,  2.0 );

    // Comparison with double (via raw accessor).
    _test(  sc1() < 10.0 );
    _test( !( sc1() < 5.0 ) );
    _test( !( sc1() < 2.0 ) );

    // Comparison with scalar.
    _test(  sc1 < sc2 );
    _test( !( sc1 < sc3 ) );
    _test( !( sc1 < sc4 ) );
}

// ============================================================================

void ScalarVariable_Test::Greater_Than_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3( PLAIN,  5.0 );
    ScalarVariable sc4( PLAIN,  2.0 );

    // Comparison with double.
    _test(  sc1() > 2.0 );
    _test( !( sc1() > 5.0 ) );
    _test( !( sc1() > 10.0 ) );

    // Comparison with scalar.
    _test(  sc2 > sc1 );
    _test( !( sc3 > sc1 ) );
    _test( !( sc4 > sc1 ) );
}

// ============================================================================

void ScalarVariable_Test::Less_Than_Or_Equal_To_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3( PLAIN,  5.0 );
    ScalarVariable sc4( PLAIN,  2.0 );

    // Comparison with double.
    _test(  sc1() <= 10.0 );
    _test(  sc1() <=  5.0 );
    _test( !( sc1() <= 2.0 ) );

    // Comparison with scalar.
    _test(  sc1 <= sc2 );
    _test(  sc1 <= sc3 );
    _test( !( sc1 <= sc4 ) );
}

// ============================================================================

void ScalarVariable_Test::Greater_Than_Or_Equal_To_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3( PLAIN,  5.0 );
    ScalarVariable sc4( PLAIN,  2.0 );

    // Comparison with double.
    _test(  sc1() >= 2.0 );
    _test(  sc1() >= 5.0 );
    _test( !( sc1() >= 10.0 ) );

    // Comparison with scalar.
    _test(  sc2 >= sc1 );
    _test(  sc3 >= sc1 );
    _test( !( sc4 >= sc1 ) );
}

// ============================================================================

void ScalarVariable_Test::Equality_Operator()
{
    ScalarVariable sc1( ANY,   5.0 );
    ScalarVariable sc2( PLAIN, 20.0 );
    ScalarVariable sc3( PLAIN,  5.0 );  // same value as sc1, different flag
    ScalarVariable sc4( PLAIN,  2.0 );
    ScalarVariable sc5( ANY,    5.0 );  // same value and flag as sc1

    // operator== compares both value and flag.
    _test( !( sc2 == sc1 ) );  // different value and flag
    _test( !( sc3 == sc1 ) );  // same value, different flag
    _test( !( sc4 == sc1 ) );  // different value and flag
    _test(    sc5 == sc1   );  // same value and flag

    // operator!=
    _test(  sc2 != sc1 );
    _test(  sc3 != sc1 );
    _test(  sc4 != sc1 );
    _test( !( sc5 != sc1 ) );
}

// ============================================================================

void ScalarVariable_Test::IsWithinRange_Function()
{
    ScalarVariable sc( ANY, 5.0 );

    _test( !sc.IsWithinRange( 2.0, 3.0 ) );  // 5 > 3
    _test(  sc.IsWithinRange( 2.0, 5.0 ) );  // 5 == upper bound
    _test(  sc.IsWithinRange( 2.0, 7.0 ) );  // 5 in (2,7)
    _test(  sc.IsWithinRange( 5.0, 7.0 ) );  // 5 == lower bound
    _test( !sc.IsWithinRange( 6.0, 9.0 ) );  // 5 < 6
}

// ============================================================================

void ScalarVariable_Test::HasNaN_Function()
{
    // Normal value — not NaN.
    ScalarVariable sc1( ANY, 5.0 );
    _test( !sc1.Has_NaN_Values() );

    // Default constructor — NaN.
    ScalarVariable sc2;
    _test( sc2.Has_NaN_Values() );

    // Explicitly set to NaN.
    ScalarVariable sc3( ANY, std::numeric_limits<double>::quiet_NaN() );
    _test( sc3.Has_NaN_Values() );

    // After ZapNAN equivalent — set to real value.
    sc3() = 42.0;
    _test( !sc3.Has_NaN_Values() );
}

} // namespace csmp

