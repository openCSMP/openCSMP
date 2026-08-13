// ArrayVariable_Test.cpp

#include "ArrayVariable_Test.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "ScalarVariable.h"
#include "PropertyDatabase.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>

using namespace std;

namespace csmp {

void ArrayVariable_Test::run()
{
    TestConstruction();
    TestElementAccess();
    TestFlag();
    TestAssignment();
    TestArithmetic();
    TestCompoundAssignment();
    TestComparison();
    TestQueries();
    TestRangeInterface();
    TestBinaryIO();
}

// ============================================================================
//  TestConstruction
// ============================================================================

void ArrayVariable_Test::TestConstruction()
{
    // Default constructor — zero elements.
    ArrayVariable av0;
    _test( av0.Size() == 0 );
    _test( av0.Flag() == ANY );

    // Size + value + flag constructor.
    ArrayVariable av1( 5, 3.0, DIRICH );
    _test( av1.Size() == 5 );
    _test( av1.Flag() == DIRICH );
    for ( uint32_t i = 0; i < 5; ++i )
        _equal( av1[i], 3.0, TOL );

    // Size + value (default flag = ANY).
    ArrayVariable av2( 3, 7.0 );
    _test( av2.Size() == 3 );
    _test( av2.Flag() == ANY );
    _equal( av2[0], 7.0, TOL );

    // Initializer list constructor.
    ArrayVariable av3( { 1.0, 2.0, 3.0, 4.0 }, PLAIN );
    _test( av3.Size() == 4 );
    _test( av3.Flag() == PLAIN );
    _equal( av3[0], 1.0, TOL );
    _equal( av3[1], 2.0, TOL );
    _equal( av3[2], 3.0, TOL );
    _equal( av3[3], 4.0, TOL );

    // Copy constructor.
    ArrayVariable av4( av3 );
    _test( av4 == av3 );
    _test( av4.Flag() == PLAIN );

    // Move constructor.
    ArrayVariable av5( ArrayVariable( 3, 9.0, ROBIN ) );
    _test( av5.Size() == 3 );
    _test( av5.Flag() == ROBIN );
    _equal( av5[0], 9.0, TOL );

    // PropertyDatabase constructor.
    PropertyDatabase<3> pdb( "PropertyDatabase_Test-variables.txt" );
    ArrayVariable av6( "element array 2", pdb, 1.3, DIRICH );
    _test( av6.Size() == 22 );
    _test( av6.Flag() == DIRICH );
    for ( uint32_t i = 0; i < av6.Size(); ++i )
        _equal( av6[i], 1.3, TOL );

    // Index constructor.
    const csmp::Index key = pdb.StorageKey( "element array 2" );
    ArrayVariable av7( key, 2.5, PLAIN );
    _test( av7.Size() == 22 );
    _test( av7.Flag() == PLAIN );
    _equal( av7[0], 2.5, TOL );
}

// ============================================================================
//  TestElementAccess
// ============================================================================

void ArrayVariable_Test::TestElementAccess()
{
    ArrayVariable av( 5, 0.0 );

    // operator() write, operator[] read.
    av(0) = 10.0;
    av(2) = 20.0;
    av(4) = 30.0;
    _equal( av[0], 10.0, TOL );
    _equal( av[1],  0.0, TOL );
    _equal( av[2], 20.0, TOL );
    _equal( av[3],  0.0, TOL );
    _equal( av[4], 30.0, TOL );

    // Component mutator and accessor.
    av.Component( 1, 99.0 );
    _equal( av.Component(1), 99.0, TOL );
    _equal( av[1], 99.0, TOL );

    // Size.
    _test( av.Size() == 5 );

    // Resize — grow: new elements get default NaN.
    av.Resize( 8 );
    _test( av.Size() == 8 );
    _equal( av[0], 10.0, TOL );  // existing values preserved
    _equal( av[2], 20.0, TOL );
    _test( std::isnan( av[5] ) );  // new elements are NaN

    // Resize — shrink: existing values preserved up to new size.
    av.Resize( 3 );
    _test( av.Size() == 3 );
    _equal( av[0], 10.0, TOL );
    _equal( av[1], 99.0, TOL );
    _equal( av[2], 20.0, TOL );

    // Resize with explicit fill value.
    av.Resize( 6, -1.0 );
    _test( av.Size() == 6 );
    _equal( av[3], -1.0, TOL );
    _equal( av[4], -1.0, TOL );
    _equal( av[5], -1.0, TOL );
}

// ============================================================================
//  TestFlag
// ============================================================================

void ArrayVariable_Test::TestFlag()
{
    ArrayVariable av( 4, 1.0, ANY );

    // Initial flag.
    _test( av.Flag() == ANY );

    // Set flag.
    av.Flag( DIRICH );
    _test( av.Flag() == DIRICH );

    // Non-const reference.
    av.Flag() = PLAIN;
    _test( av.Flag() == PLAIN );

    // Flag is a single value for the whole array — not per-element.
    // Changing it affects all reads.
    av.Flag( ROBIN );
    _test( av.Flag() == ROBIN );
    _test( av.Size() == 4 );  // size unaffected by flag change
}

// ============================================================================
//  TestAssignment
// ============================================================================

void ArrayVariable_Test::TestAssignment()
{
    ArrayVariable av( 4, 0.0, ANY );

    // operator=(double) — sets all elements, flag unchanged.
    av = 5.0;
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( av[i], 5.0, TOL );
    _test( av.Flag() == ANY );

    // operator=(ScalarVariable) — sets all elements to scalar value,
    // flag unchanged.
    ScalarVariable sc( DIRICH, 3.0 );
    av = sc;
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( av[i], 3.0, TOL );
    _test( av.Flag() == ANY );  // flag NOT copied from scalar

    // operator=(ArrayVariable) — values and flag copied.
    ArrayVariable av2( 4, 7.0, ROBIN );
    av = av2;
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( av[i], 7.0, TOL );
    _test( av.Flag() == ROBIN );

    // Move assignment.
    ArrayVariable av3( 4, 2.0, PLAIN );
    av = std::move( av3 );
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( av[i], 2.0, TOL );
    _test( av.Flag() == PLAIN );

    // CopyValuesOnly from FlaggedArrayVariable — values copied, flag unchanged.
    FlaggedArrayVariable fav( 4, 9.0, DIRICH );
    av.Flag( ROBIN );
    av.CopyValuesOnly( fav );
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( av[i], 9.0, TOL );
    _test( av.Flag() == ROBIN );  // flag NOT changed by CopyValuesOnly
}

// ============================================================================
//  TestArithmetic
// ============================================================================

void ArrayVariable_Test::TestArithmetic()
{
    ArrayVariable av( { 1.0, 2.0, 3.0, 4.0 }, ANY );

    // operator+(double) — returns new, original unchanged.
    ArrayVariable result = av + 10.0;
    _equal( result[0], 11.0, TOL ); _equal( result[1], 12.0, TOL );
    _equal( result[2], 13.0, TOL ); _equal( result[3], 14.0, TOL );
    _equal( av[0], 1.0, TOL );  // original unchanged

    // operator-(double).
    result = av - 1.0;
    _equal( result[0], 0.0, TOL ); _equal( result[1], 1.0, TOL );
    _equal( result[2], 2.0, TOL ); _equal( result[3], 3.0, TOL );

    // operator*(double).
    result = av * 3.0;
    _equal( result[0],  3.0, TOL ); _equal( result[1],  6.0, TOL );
    _equal( result[2],  9.0, TOL ); _equal( result[3], 12.0, TOL );

    // operator/(double).
    result = av / 2.0;
    _equal( result[0], 0.5, TOL ); _equal( result[1], 1.0, TOL );
    _equal( result[2], 1.5, TOL ); _equal( result[3], 2.0, TOL );

    // operator+(ArrayVariable) — element-by-element.
    ArrayVariable av2( { 4.0, 3.0, 2.0, 1.0 }, ANY );
    result = av + av2;
    _equal( result[0], 5.0, TOL ); _equal( result[1], 5.0, TOL );
    _equal( result[2], 5.0, TOL ); _equal( result[3], 5.0, TOL );

    // operator-(ArrayVariable).
    result = av - av2;
    _equal( result[0], -3.0, TOL ); _equal( result[1], -1.0, TOL );
    _equal( result[2],  1.0, TOL ); _equal( result[3],  3.0, TOL );

    // operator*(ArrayVariable) — element-by-element.
    result = av * av2;
    _equal( result[0],  4.0, TOL ); _equal( result[1],  6.0, TOL );
    _equal( result[2],  6.0, TOL ); _equal( result[3],  4.0, TOL );

    // operator/(ArrayVariable) — element-by-element.
    result = av / av2;
    _equal( result[0], 0.25, TOL ); _equal( result[1], 2.0/3.0, TOL );
    _equal( result[2], 1.5,  TOL ); _equal( result[3], 4.0,     TOL );

    // Pow(double).
    result = av.Pow( 2.0 );
    _equal( result[0],  1.0, TOL ); _equal( result[1],  4.0, TOL );
    _equal( result[2],  9.0, TOL ); _equal( result[3], 16.0, TOL );
    _equal( av[0], 1.0, TOL );  // original unchanged

    // Pow(0.5) = sqrt — round-trip.
    ArrayVariable sq = av.Pow( 2.0 );
    ArrayVariable rt = sq.Pow( 0.5 );
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( rt[i], av[i], TOL );

    // Flags are preserved in results.
    ArrayVariable avf( { 1.0, 2.0 }, DIRICH );
    result = avf + 1.0;
    _test( result.Flag() == DIRICH );
}

// ============================================================================
//  TestCompoundAssignment
// ============================================================================

void ArrayVariable_Test::TestCompoundAssignment()
{
    ArrayVariable av( { 2.0, 4.0, 6.0, 8.0 }, ANY );

    // += double.
    av += 1.0;
    _equal( av[0], 3.0, TOL ); _equal( av[1], 5.0, TOL );
    _equal( av[2], 7.0, TOL ); _equal( av[3], 9.0, TOL );

    // -= double.
    av -= 1.0;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // *= double.
    av *= 2.0;
    _equal( av[0],  4.0, TOL ); _equal( av[1],  8.0, TOL );
    _equal( av[2], 12.0, TOL ); _equal( av[3], 16.0, TOL );

    // /= double.
    av /= 2.0;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // += ScalarVariable.
    ScalarVariable sc( ANY, 10.0 );
    av += sc;
    _equal( av[0], 12.0, TOL ); _equal( av[1], 14.0, TOL );
    _equal( av[2], 16.0, TOL ); _equal( av[3], 18.0, TOL );

    // -= ScalarVariable.
    av -= sc;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // *= ScalarVariable.
    ScalarVariable sc2( ANY, 3.0 );
    av *= sc2;
    _equal( av[0],  6.0, TOL ); _equal( av[1], 12.0, TOL );
    _equal( av[2], 18.0, TOL ); _equal( av[3], 24.0, TOL );

    // /= ScalarVariable.
    av /= sc2;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // += ArrayVariable.
    ArrayVariable av2( { 1.0, 1.0, 1.0, 1.0 }, ANY );
    av += av2;
    _equal( av[0], 3.0, TOL ); _equal( av[1], 5.0, TOL );
    _equal( av[2], 7.0, TOL ); _equal( av[3], 9.0, TOL );

    // -= ArrayVariable.
    av -= av2;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // *= ArrayVariable (element-by-element).
    av *= av2;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // /= ArrayVariable.
    av /= av2;
    _equal( av[0], 2.0, TOL ); _equal( av[1], 4.0, TOL );
    _equal( av[2], 6.0, TOL ); _equal( av[3], 8.0, TOL );

    // Flags unchanged by compound assignment.
    ArrayVariable avf( { 1.0, 2.0 }, DIRICH );
    avf += 1.0;
    _test( avf.Flag() == DIRICH );
}

// ============================================================================
//  TestComparison
// ============================================================================

void ArrayVariable_Test::TestComparison()
{
    ArrayVariable av1( { 1.0, 2.0, 3.0 }, ANY );
    ArrayVariable av2( { 1.0, 2.0, 3.0 }, ANY );   // same values and flag
    ArrayVariable av3( { 1.0, 2.0, 4.0 }, ANY );   // different values
    ArrayVariable av4( { 1.0, 2.0, 3.0 }, DIRICH ); // different flag

    // operator==: compares values and flag.
    _test(  av1 == av2 );
    _test( !( av1 == av3 ) );
    _test( !( av1 == av4 ) );  // different flag

    // operator!=.
    _test( !( av1 != av2 ) );
    _test(  av1 != av3 );
    _test(  av1 != av4 );

    // operator<: lexicographic on values.
    ArrayVariable avA( { 1.0, 2.0, 3.0 }, ANY );
    ArrayVariable avB( { 1.0, 2.0, 4.0 }, ANY );
    _test(  avA < avB );
    _test( !( avB < avA ) );
    _test( !( avA < avA ) );  // irreflexive
}

// ============================================================================
//  TestQueries
// ============================================================================

void ArrayVariable_Test::TestQueries()
{
    // IsWithinRange.
    ArrayVariable av( { 3.0, 5.0, 7.0 }, ANY );
    _test(  av.IsWithinRange( 2.0, 8.0 ) );   // all in range
    _test(  av.IsWithinRange( 3.0, 7.0 ) );   // exact bounds
    _test( !av.IsWithinRange( 4.0, 8.0 ) );   // 3 < 4
    _test( !av.IsWithinRange( 2.0, 6.0 ) );   // 7 > 6
    _test( !av.IsWithinRange( 4.0, 6.0 ) );   // 3 < 4 and 7 > 6

    // Has_NaN_Values.
    ArrayVariable avNoNaN( { 1.0, 2.0, 3.0 }, ANY );
    _test( !avNoNaN.Has_NaN_Values() );

    ArrayVariable avNaN( 3, std::numeric_limits<double>::quiet_NaN() );
    _test( avNaN.Has_NaN_Values() );

    // Mixed: one NaN.
    ArrayVariable avMixed( { 1.0, std::numeric_limits<double>::quiet_NaN(), 3.0 }, ANY );
    _test( avMixed.Has_NaN_Values() );

    // MinMax.
    ArrayVariable avMM( { 3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0 }, ANY );
    double mn, mx;
    avMM.MinMax( mn, mx );
    _equal( mn, 1.0, TOL );
    _equal( mx, 9.0, TOL );

    // MinMax on single element.
    ArrayVariable avSingle( 1, 42.0 );
    avSingle.MinMax( mn, mx );
    _equal( mn, 42.0, TOL );
    _equal( mx, 42.0, TOL );

    // Sort.
    ArrayVariable avSort( { 5.0, 2.0, 8.0, 1.0, 3.0 }, ANY );
    avSort.Sort();
    _equal( avSort[0], 1.0, TOL );
    _equal( avSort[1], 2.0, TOL );
    _equal( avSort[2], 3.0, TOL );
    _equal( avSort[3], 5.0, TOL );
    _equal( avSort[4], 8.0, TOL );

    // NextLargestEntry.
    ArrayVariable avNL( { 1.0, 3.0, 5.0, 7.0, 9.0 }, ANY );
    _equal( avNL.NextLargestEntry( 4.0 ), 5.0, TOL );
    _equal( avNL.NextLargestEntry( 0.0 ), 1.0, TOL );
    _equal( avNL.NextLargestEntry( 8.0 ), 9.0, TOL );
    // No larger entry — returns max double.
    _equal( avNL.NextLargestEntry( 9.0 ),
            std::numeric_limits<double>::max(), TOL );

    // HasLargerEntry.
    _test(  avNL.HasLargerEntry( 4.0 ) );
    _test(  avNL.HasLargerEntry( 0.0 ) );
    _test( !avNL.HasLargerEntry( 9.0 ) );
    _test( !avNL.HasLargerEntry( 10.0 ) );
}

// ============================================================================
//  TestRangeInterface
// ============================================================================

void ArrayVariable_Test::TestRangeInterface()
{
    ArrayVariable av( { 1.0, 2.0, 3.0, 4.0, 5.0 }, ANY );

    // begin/end — range-based for loop.
    double sum = 0.0;
    for ( double v : av ) sum += v;
    _equal( sum, 15.0, TOL );

    // std::accumulate via begin/end.
    double acc = std::accumulate( av.begin(), av.end(), 0.0 );
    _equal( acc, 15.0, TOL );

    // cbegin/cend — const iteration.
    const ArrayVariable& cav = av;
    double csum = 0.0;
    for ( auto it = cav.cbegin(); it != cav.cend(); ++it )
        csum += *it;
    _equal( csum, 15.0, TOL );

    // Legacy uppercase Begin/End.
    double lsum = 0.0;
    for ( auto it = av.Begin(); it != av.End(); ++it )
        lsum += *it;
    _equal( lsum, 15.0, TOL );

    // std::sort via begin/end.
    ArrayVariable avSort( { 5.0, 3.0, 1.0, 4.0, 2.0 }, ANY );
    std::sort( avSort.begin(), avSort.end() );
    _equal( avSort[0], 1.0, TOL );
    _equal( avSort[4], 5.0, TOL );
}

// ============================================================================
//  TestBinaryIO
// ============================================================================

void ArrayVariable_Test::TestBinaryIO()
{
    // Write.
    ArrayVariable avOut( 500, 999.0, DIRICH );
    {
        fstream fp( "ArrayVariableBinaryIO_Test",
                    ios::out | ios::binary );
        if ( !fp.is_open() )
        {
            cerr << "\nArrayVariable_Test: binary file could not be created.\n";
            _test( false );
            return;
        }
        avOut.Out( fp );
    }

    // Read back.
    ArrayVariable avIn;
    {
        fstream fp( "ArrayVariableBinaryIO_Test",
                    ios::in | ios::binary );
        if ( !fp.is_open() )
        {
            cerr << "\nArrayVariable_Test: binary file could not be opened.\n";
            _test( false );
            return;
        }
        avIn.In( fp );
    }

    // Verify round-trip.
    _test( avIn == avOut );
    _test( avIn.Size() == 500 );
    _test( avIn.Flag() == DIRICH );
    _equal( avIn[0],   999.0, TOL );
    _equal( avIn[499], 999.0, TOL );

    // Text file output — verify file is created and non-empty.
    ArrayVariable avText( { 1.0, 2.0, 3.0 }, ANY );
    _test( avText.Out( "ArrayVariableTextIO_Test.txt" ) );

    // Screen output — just verify it does not crash.
    avText.Out( 3 );
}

} // namespace csmp

