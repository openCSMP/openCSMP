// FlaggedArrayVariable_Test.cpp

#include "FlaggedArrayVariable_Test.h"
#include "FlaggedArrayVariable.h"
#include "ArrayVariable.h"
#include "ScalarVariable.h"
#include "PropertyDatabase.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>

using namespace std;

namespace csmp {

void FlaggedArrayVariable_Test::run()
{
    TestConstruction();
    TestElementAccess();
    TestPerElementFlags();
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

void FlaggedArrayVariable_Test::TestConstruction()
{
    // Default constructor — zero elements.
    FlaggedArrayVariable fav0;
    _test( fav0.Size() == 0 );

    // Size + value + flag constructor.
    FlaggedArrayVariable fav1( 5, 3.0, DIRICH );
    _test( fav1.Size() == 5 );
    for ( uint32_t i = 0; i < 5; ++i )
    {
        _equal( fav1[i],    3.0,   TOL );
        _test(  fav1.Flag(i) == DIRICH );
    }

    // Size + value (default flag = ANY).
    FlaggedArrayVariable fav2( 3, 7.0 );
    _test( fav2.Size() == 3 );
    for ( uint32_t i = 0; i < 3; ++i )
    {
        _equal( fav2[i],    7.0, TOL );
        _test(  fav2.Flag(i) == ANY );
    }

    // Copy constructor — values and flags copied.
    FlaggedArrayVariable fav3( fav1 );
    _test( fav3 == fav1 );
    _test( fav3.Size() == 5 );
    for ( uint32_t i = 0; i < 5; ++i )
        _test( fav3.Flag(i) == DIRICH );

    // Move constructor.
    FlaggedArrayVariable fav4( FlaggedArrayVariable( 3, 9.0, ROBIN ) );
    _test( fav4.Size() == 3 );
    for ( uint32_t i = 0; i < 3; ++i )
    {
        _equal( fav4[i],    9.0,  TOL );
        _test(  fav4.Flag(i) == ROBIN );
    }

    // Index constructor.
    PropertyDatabase<3> pdb( "PropertyDatabase_Test-variables.txt" );
    const csmp::Index key = pdb.StorageKey( "element array 2" );
    FlaggedArrayVariable fav5( key, 2.5, PLAIN );
    _test( fav5.Size() == 22 );
    for ( uint32_t i = 0; i < fav5.Size(); ++i )
    {
        _equal( fav5[i],    2.5,   TOL );
        _test(  fav5.Flag(i) == PLAIN );
    }

    // PropertyDatabase constructor.
    FlaggedArrayVariable fav6( "element array 2", pdb, 1.3, DIRICH );
    _test( fav6.Size() == 22 );
    for ( uint32_t i = 0; i < fav6.Size(); ++i )
    {
        _equal( fav6[i],    1.3,   TOL );
        _test(  fav6.Flag(i) == DIRICH );
    }
}

// ============================================================================
//  TestElementAccess
// ============================================================================

void FlaggedArrayVariable_Test::TestElementAccess()
{
    FlaggedArrayVariable fav( 5, 0.0, ANY );

    // operator() write, operator[] read.
    fav(0) = 10.0;
    fav(2) = 20.0;
    fav(4) = 30.0;
    _equal( fav[0], 10.0, TOL );
    _equal( fav[1],  0.0, TOL );
    _equal( fav[2], 20.0, TOL );
    _equal( fav[3],  0.0, TOL );
    _equal( fav[4], 30.0, TOL );

    // Component mutator and accessor.
    fav.Component( 1, 99.0 );
    _equal( fav.Component(1), 99.0, TOL );
    _equal( fav[1], 99.0, TOL );

    // Size.
    _test( fav.Size() == 5 );

    // Resize — grow: new elements get NaN, new flags get ANY.
    fav.Resize( 8 );
    _test( fav.Size() == 8 );
    _equal( fav[0], 10.0, TOL );  // existing values preserved
    _equal( fav[2], 20.0, TOL );
    _test( std::isnan( fav[5] ) );  // new elements are NaN
    _test( fav.Flag(5) == ANY );    // new flags are ANY

    // Resize — shrink: existing values preserved up to new size.
    fav.Resize( 3 );
    _test( fav.Size() == 3 );
    _equal( fav[0], 10.0, TOL );
    _equal( fav[1], 99.0, TOL );
    _equal( fav[2], 20.0, TOL );
}

// ============================================================================
//  TestPerElementFlags
// ============================================================================

void FlaggedArrayVariable_Test::TestPerElementFlags()
{
    FlaggedArrayVariable fav( 5, 1.0, ANY );

    // All flags initially ANY.
    for ( uint32_t i = 0; i < 5; ++i )
        _test( fav.Flag(i) == ANY );

    // Set individual flags independently.
    fav.Flag( 0, DIRICH );
    fav.Flag( 2, PLAIN  );
    fav.Flag( 4, ROBIN  );

    _test( fav.Flag(0) == DIRICH );
    _test( fav.Flag(1) == ANY    );  // unchanged
    _test( fav.Flag(2) == PLAIN  );
    _test( fav.Flag(3) == ANY    );  // unchanged
    _test( fav.Flag(4) == ROBIN  );

    // Non-const reference — modify via reference.
    fav.Flag(1) = NEUMANN;
    _test( fav.Flag(1) == NEUMANN );

    // Setting one flag does not affect others.
    fav.Flag( 3, CONSTANT_FLUX );
    _test( fav.Flag(0) == DIRICH       );
    _test( fav.Flag(1) == NEUMANN      );
    _test( fav.Flag(2) == PLAIN        );
    _test( fav.Flag(3) == CONSTANT_FLUX );
    _test( fav.Flag(4) == ROBIN        );

    // Flags survive Resize (existing flags preserved, new flags = ANY).
    fav.Resize( 7 );
    _test( fav.Flag(0) == DIRICH       );
    _test( fav.Flag(4) == ROBIN        );
    _test( fav.Flag(5) == ANY          );  // new
    _test( fav.Flag(6) == ANY          );  // new

    // Flags survive shrink.
    fav.Resize( 3 );
    _test( fav.Flag(0) == DIRICH );
    _test( fav.Flag(1) == NEUMANN );
    _test( fav.Flag(2) == PLAIN  );
}

// ============================================================================
//  TestAssignment
// ============================================================================

void FlaggedArrayVariable_Test::TestAssignment()
{
    FlaggedArrayVariable fav( 4, 0.0, ANY );

    // operator=(double) — sets all values, flags unchanged.
    fav.Flag( 0, DIRICH );
    fav = 5.0;
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( fav[i], 5.0, TOL );
    _test( fav.Flag(0) == DIRICH );  // flag NOT changed
    _test( fav.Flag(1) == ANY    );

    // operator=(ScalarVariable) — sets all values, flags unchanged.
    ScalarVariable sc( ROBIN, 3.0 );
    fav = sc;
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( fav[i], 3.0, TOL );
    _test( fav.Flag(0) == DIRICH );  // flags NOT changed by scalar assignment
    _test( fav.Flag(1) == ANY    );

    // operator=(FlaggedArrayVariable) — values and flags copied.
    FlaggedArrayVariable fav2( 4, 7.0, ROBIN );
    fav2.Flag( 1, PLAIN );
    fav = fav2;
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( fav[i], 7.0, TOL );
    _test( fav.Flag(0) == ROBIN );
    _test( fav.Flag(1) == PLAIN );
    _test( fav.Flag(2) == ROBIN );

    // Move assignment.
    FlaggedArrayVariable fav3( 4, 2.0, PLAIN );
    fav3.Flag( 2, DIRICH );
    fav = std::move( fav3 );
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( fav[i], 2.0, TOL );
    _test( fav.Flag(0) == PLAIN );
    _test( fav.Flag(2) == DIRICH );

    // CopyValuesOnly from ArrayVariable — values copied, flags unchanged.
    ArrayVariable av( { 9.0, 8.0, 7.0, 6.0 }, ANY );
    fav.Flag( 0, ROBIN );
    fav.CopyValuesOnly( av );
    _equal( fav[0], 9.0, TOL );
    _equal( fav[1], 8.0, TOL );
    _equal( fav[2], 7.0, TOL );
    _equal( fav[3], 6.0, TOL );
    _test( fav.Flag(0) == ROBIN  );  // flags NOT changed by CopyValuesOnly
    _test( fav.Flag(2) == DIRICH );
}

// ============================================================================
//  TestArithmetic
// ============================================================================

void FlaggedArrayVariable_Test::TestArithmetic()
{
    // Arithmetic operators return new objects with flags from *this.
    FlaggedArrayVariable fav( 4, 0.0, ANY );
    fav(0) = 1.0; fav(1) = 2.0; fav(2) = 3.0; fav(3) = 4.0;
    fav.Flag( 0, DIRICH );
    fav.Flag( 2, PLAIN  );

    // operator+(double) — returns new, original unchanged.
    FlaggedArrayVariable result = fav + 10.0;
    _equal( result[0], 11.0, TOL ); _equal( result[1], 12.0, TOL );
    _equal( result[2], 13.0, TOL ); _equal( result[3], 14.0, TOL );
    _test(  result.Flag(0) == DIRICH );  // flags from *this
    _test(  result.Flag(2) == PLAIN  );
    _equal( fav[0], 1.0, TOL );  // original unchanged

    // operator-(double).
    result = fav - 1.0;
    _equal( result[0], 0.0, TOL ); _equal( result[1], 1.0, TOL );
    _equal( result[2], 2.0, TOL ); _equal( result[3], 3.0, TOL );

    // operator*(double).
    result = fav * 3.0;
    _equal( result[0],  3.0, TOL ); _equal( result[1],  6.0, TOL );
    _equal( result[2],  9.0, TOL ); _equal( result[3], 12.0, TOL );

    // operator/(double).
    result = fav / 2.0;
    _equal( result[0], 0.5, TOL ); _equal( result[1], 1.0, TOL );
    _equal( result[2], 1.5, TOL ); _equal( result[3], 2.0, TOL );

    // operator+(FlaggedArrayVariable) — element-by-element, flags from *this.
    FlaggedArrayVariable fav2( 4, 0.0, ROBIN );
    fav2(0)=4.0; fav2(1)=3.0; fav2(2)=2.0; fav2(3)=1.0;
    result = fav + fav2;
    _equal( result[0], 5.0, TOL ); _equal( result[1], 5.0, TOL );
    _equal( result[2], 5.0, TOL ); _equal( result[3], 5.0, TOL );
    _test(  result.Flag(0) == DIRICH );  // flags from lhs (*this)
    _test(  result.Flag(1) == ANY    );

    // operator-(FlaggedArrayVariable).
    result = fav - fav2;
    _equal( result[0], -3.0, TOL ); _equal( result[1], -1.0, TOL );
    _equal( result[2],  1.0, TOL ); _equal( result[3],  3.0, TOL );

    // operator*(FlaggedArrayVariable).
    result = fav * fav2;
    _equal( result[0],  4.0, TOL ); _equal( result[1],  6.0, TOL );
    _equal( result[2],  6.0, TOL ); _equal( result[3],  4.0, TOL );

    // operator/(FlaggedArrayVariable).
    result = fav / fav2;
    _equal( result[0], 0.25,     TOL );
    _equal( result[1], 2.0/3.0,  TOL );
    _equal( result[2], 1.5,      TOL );
    _equal( result[3], 4.0,      TOL );

    // Pow(double) — returns new, flags from *this.
    result = fav.Pow( 2.0 );
    _equal( result[0],  1.0, TOL ); _equal( result[1],  4.0, TOL );
    _equal( result[2],  9.0, TOL ); _equal( result[3], 16.0, TOL );
    _test(  result.Flag(0) == DIRICH );
    _equal( fav[0], 1.0, TOL );  // original unchanged

    // Pow round-trip.
    FlaggedArrayVariable sq  = fav.Pow( 2.0 );
    FlaggedArrayVariable rt  = sq.Pow( 0.5 );
    for ( uint32_t i = 0; i < 4; ++i )
        _equal( rt[i], fav[i], TOL );
}

// ============================================================================
//  TestCompoundAssignment
// ============================================================================

void FlaggedArrayVariable_Test::TestCompoundAssignment()
{
    // Compound assignment modifies values only — flags are NOT modified.
    FlaggedArrayVariable fav( 4, 2.0, ANY );
    fav.Flag( 0, DIRICH );
    fav.Flag( 2, PLAIN  );

    // += double.
    fav += 1.0;
    _equal( fav[0], 3.0, TOL ); _equal( fav[1], 3.0, TOL );
    _equal( fav[2], 3.0, TOL ); _equal( fav[3], 3.0, TOL );
    _test( fav.Flag(0) == DIRICH );  // flags unchanged
    _test( fav.Flag(1) == ANY    );

    // -= double.
    fav -= 1.0;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );

    // *= double.
    fav *= 3.0;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 6.0, TOL );

    // /= double.
    fav /= 3.0;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );

    // += ScalarVariable.
    ScalarVariable sc( ANY, 10.0 );
    fav += sc;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 12.0, TOL );
    _test( fav.Flag(0) == DIRICH );  // flags unchanged

    // -= ScalarVariable.
    fav -= sc;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );

    // *= ScalarVariable.
    ScalarVariable sc2( ANY, 4.0 );
    fav *= sc2;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 8.0, TOL );

    // /= ScalarVariable.
    fav /= sc2;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );

    // += FlaggedArrayVariable.
    FlaggedArrayVariable fav2( 4, 1.0, ROBIN );
    fav += fav2;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 3.0, TOL );
    _test( fav.Flag(0) == DIRICH );  // flags of *this unchanged

    // -= FlaggedArrayVariable.
    fav -= fav2;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );

    // *= FlaggedArrayVariable.
    fav *= fav2;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );

    // /= FlaggedArrayVariable.
    fav /= fav2;
    for ( uint32_t i = 0; i < 4; ++i ) _equal( fav[i], 2.0, TOL );
}

// ============================================================================
//  TestComparison
// ============================================================================

void FlaggedArrayVariable_Test::TestComparison()
{
    // Comparison is value-only — flags are NOT compared.
    FlaggedArrayVariable fav1( 3, 0.0, ANY );
    fav1(0)=1.0; fav1(1)=2.0; fav1(2)=3.0;

    FlaggedArrayVariable fav2( 3, 0.0, DIRICH );  // different flags
    fav2(0)=1.0; fav2(1)=2.0; fav2(2)=3.0;

    FlaggedArrayVariable fav3( 3, 0.0, ANY );
    fav3(0)=1.0; fav3(1)=2.0; fav3(2)=4.0;  // different values

    // operator==: value-only (flags ignored).
    _test(  fav1 == fav2 );  // same values, different flags — equal
    _test( !( fav1 == fav3 ) );

    // operator!=.
    _test( !( fav1 != fav2 ) );
    _test(  fav1 != fav3 );

    // operator<: lexicographic on values.
    _test(  fav1 < fav3 );
    _test( !( fav3 < fav1 ) );
    _test( !( fav1 < fav1 ) );  // irreflexive

    // operator<=.
    _test(  fav1 <= fav2 );  // equal
    _test(  fav1 <= fav3 );  // less
    _test( !( fav3 <= fav1 ) );

    // operator>.
    _test(  fav3 > fav1 );
    _test( !( fav1 > fav3 ) );
    _test( !( fav1 > fav1 ) );

    // operator>=.
    _test(  fav2 >= fav1 );  // equal
    _test(  fav3 >= fav1 );  // greater
    _test( !( fav1 >= fav3 ) );
}

// ============================================================================
//  TestQueries
// ============================================================================

void FlaggedArrayVariable_Test::TestQueries()
{
    // IsWithinRange — checks all values regardless of flags.
    FlaggedArrayVariable fav( 3, 0.0, ANY );
    fav(0)=3.0; fav(1)=5.0; fav(2)=7.0;
    fav.Flag( 0, DIRICH );  // flag does not affect range check

    _test(  fav.IsWithinRange( 2.0, 8.0 ) );
    _test(  fav.IsWithinRange( 3.0, 7.0 ) );
    _test( !fav.IsWithinRange( 4.0, 8.0 ) );  // 3 < 4
    _test( !fav.IsWithinRange( 2.0, 6.0 ) );  // 7 > 6

    // Has_NaN_Values — checks all values regardless of flags.
    FlaggedArrayVariable favNoNaN( 3, 1.0, ANY );
    _test( !favNoNaN.Has_NaN_Values() );

    FlaggedArrayVariable favNaN( 3, std::numeric_limits<double>::quiet_NaN(), ANY );
    _test( favNaN.Has_NaN_Values() );

    // Mixed: one NaN, one DIRICH-flagged.
    FlaggedArrayVariable favMixed( 3, 1.0, ANY );
    favMixed(1) = std::numeric_limits<double>::quiet_NaN();
    favMixed.Flag( 1, DIRICH );
    _test( favMixed.Has_NaN_Values() );  // NaN found regardless of flag

    // MinMax — checks all values regardless of flags.
    FlaggedArrayVariable favMM( 5, 0.0, ANY );
    favMM(0)=3.0; favMM(1)=1.0; favMM(2)=4.0; favMM(3)=1.0; favMM(4)=9.0;
    favMM.Flag( 4, DIRICH );  // largest value is DIRICH — still checked
    double mn, mx;
    favMM.MinMax( mn, mx );
    _equal( mn, 1.0, TOL );
    _equal( mx, 9.0, TOL );

    // Sort — data and flags sorted together by data value.
    FlaggedArrayVariable favSort( 4, 0.0, ANY );
    favSort(0)=5.0; favSort(1)=2.0; favSort(2)=8.0; favSort(3)=1.0;
    favSort.Flag( 0, DIRICH );  // flag 0 is on value 5.0
    favSort.Flag( 2, PLAIN  );  // flag 2 is on value 8.0
    favSort.Sort();

    // After sort: values [1,2,5,8], flags follow their values.
    _equal( favSort[0], 1.0, TOL ); _test( favSort.Flag(0) == ANY   );
    _equal( favSort[1], 2.0, TOL ); _test( favSort.Flag(1) == ANY   );
    _equal( favSort[2], 5.0, TOL ); _test( favSort.Flag(2) == DIRICH );
    _equal( favSort[3], 8.0, TOL ); _test( favSort.Flag(3) == PLAIN  );

    // NextLargestEntry.
    FlaggedArrayVariable favNL( 5, 0.0, ANY );
    favNL(0)=1.0; favNL(1)=3.0; favNL(2)=5.0; favNL(3)=7.0; favNL(4)=9.0;
    _equal( favNL.NextLargestEntry( 4.0 ), 5.0, TOL );
    _equal( favNL.NextLargestEntry( 0.0 ), 1.0, TOL );
    _equal( favNL.NextLargestEntry( 8.0 ), 9.0, TOL );
    _equal( favNL.NextLargestEntry( 9.0 ),
            std::numeric_limits<double>::max(), TOL );

    // HasLargerEntry.
    _test(  favNL.HasLargerEntry( 4.0 ) );
    _test(  favNL.HasLargerEntry( 0.0 ) );
    _test( !favNL.HasLargerEntry( 9.0 ) );
    _test( !favNL.HasLargerEntry( 10.0 ) );
}

// ============================================================================
//  TestRangeInterface
// ============================================================================

void FlaggedArrayVariable_Test::TestRangeInterface()
{
    FlaggedArrayVariable fav( 5, 0.0, ANY );
    fav(0)=1.0; fav(1)=2.0; fav(2)=3.0; fav(3)=4.0; fav(4)=5.0;

    // begin/end — range-based for loop (values only).
    double sum = 0.0;
    for ( double v : fav ) sum += v;
    _equal( sum, 15.0, TOL );

    // std::accumulate via begin/end.
    double acc = std::accumulate( fav.begin(), fav.end(), 0.0 );
    _equal( acc, 15.0, TOL );

    // cbegin/cend — const iteration.
    const FlaggedArrayVariable& cfav = fav;
    double csum = 0.0;
    for ( auto it = cfav.cbegin(); it != cfav.cend(); ++it )
        csum += *it;
    _equal( csum, 15.0, TOL );

    // Legacy uppercase Begin/End.
    double lsum = 0.0;
    for ( auto it = fav.Begin(); it != fav.End(); ++it )
        lsum += *it;
    _equal( lsum, 15.0, TOL );

    // std::sort via begin/end — note: this sorts data only, not flags.
    // Use Sort() member to sort data and flags together.
    FlaggedArrayVariable favSort( 5, 0.0, ANY );
    favSort(0)=5.0; favSort(1)=3.0; favSort(2)=1.0; favSort(3)=4.0; favSort(4)=2.0;
    std::sort( favSort.begin(), favSort.end() );
    _equal( favSort[0], 1.0, TOL );
    _equal( favSort[4], 5.0, TOL );
}

// ============================================================================
//  TestBinaryIO
// ============================================================================

void FlaggedArrayVariable_Test::TestBinaryIO()
{
    // Construct with mixed flags to verify they survive round-trip.
    FlaggedArrayVariable favOut( 6, 0.0, ANY );
    favOut(0)=1.0; favOut(1)=2.0; favOut(2)=3.0;
    favOut(3)=4.0; favOut(4)=5.0; favOut(5)=6.0;
    favOut.Flag( 0, DIRICH );
    favOut.Flag( 2, PLAIN  );
    favOut.Flag( 4, ROBIN  );

    // Write.
    {
        fstream fp( "FlaggedArrayVariableBinaryIO_Test",
                    ios::out | ios::binary );
        if ( !fp.is_open() )
        {
            cerr << "\nFlaggedArrayVariable_Test: "
                    "binary file could not be created.\n";
            _test( false );
            return;
        }
        favOut.Out( fp );
    }

    // Read back.
    FlaggedArrayVariable favIn;
    {
        fstream fp( "FlaggedArrayVariableBinaryIO_Test",
                    ios::in | ios::binary );
        if ( !fp.is_open() )
        {
            cerr << "\nFlaggedArrayVariable_Test: "
                    "binary file could not be opened.\n";
            _test( false );
            return;
        }
        favIn.In( fp );
    }

    // Verify round-trip — values.
    _test( favIn.Size() == 6 );
    for ( uint32_t i = 0; i < 6; ++i )
        _equal( favIn[i], favOut[i], TOL );

    // Verify round-trip — flags.
    _test( favIn.Flag(0) == DIRICH );
    _test( favIn.Flag(1) == ANY    );
    _test( favIn.Flag(2) == PLAIN  );
    _test( favIn.Flag(3) == ANY    );
    _test( favIn.Flag(4) == ROBIN  );
    _test( favIn.Flag(5) == ANY    );

    // Text file output — verify file is created.
    FlaggedArrayVariable favText( 3, 1.0, ANY );
    _test( favText.Out( "FlaggedArrayVariableTextIO_Test.txt" ) );

    // Screen output — verify it does not crash.
    favText.Out( 3 );
}

} // namespace csmp

