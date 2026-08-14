// FlaggedArrayVariable_Test.h

#ifndef FLAGGED_ARRAY_VARIABLE_TEST_H
#define FLAGGED_ARRAY_VARIABLE_TEST_H

#include "Test.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/**
@brief Comprehensive unit test for FlaggedArrayVariable.

@section coverage Coverage

  Construction:
    Default, size+value+flag, Index, PropertyDatabase, copy, move.

  Element access:
    operator[], operator(), Component, Size, Resize.

  Per-element flags:
    Flag(i) get/set/reference — one flag per element.
    Flags are independent: setting one does not affect others.

  Assignment:
    operator=(double), operator=(ScalarVariable),
    operator=(FlaggedArrayVariable), CopyValuesOnly(ArrayVariable),
    move assignment.

  Arithmetic (returning new object, flags from *this):
    operator+(double), operator-(double), operator*(double), operator/(double),
    operator+(FlaggedArrayVariable), operator-(FlaggedArrayVariable),
    operator*(FlaggedArrayVariable), operator/(FlaggedArrayVariable),
    Pow(double).

  Compound assignment (modifies values, flags not modified):
    operator+=(double), operator-=(double), operator*=(double), operator/=(double),
    operator+=(ScalarVariable), operator-=(ScalarVariable),
    operator*=(ScalarVariable), operator/=(ScalarVariable),
    operator+=(FlaggedArrayVariable), operator-=(FlaggedArrayVariable),
    operator*=(FlaggedArrayVariable), operator/=(FlaggedArrayVariable).

  Comparison (value-only; flags not compared):
    operator==, operator!=, operator<, operator<=, operator>, operator>=.

  Queries:
    IsWithinRange, Has_NaN_Values, MinMax, Sort (data and flags together),
    NextLargestEntry, HasLargerEntry.

  Range interface:
    begin/end/cbegin/cend, Begin/End.

  I/O:
    Out(screen), Out(filename), Out(fstream), In(fstream).
*/
class FlaggedArrayVariable_Test : public Test
{
public:
    virtual void run();

private:
    void TestConstruction();
    void TestElementAccess();
    void TestPerElementFlags();
    void TestAssignment();
    void TestArithmetic();
    void TestCompoundAssignment();
    void TestComparison();
    void TestQueries();
    void TestRangeInterface();
    void TestBinaryIO();

    static constexpr double TOL = 1.0e-14;
};

} // namespace csmp

#endif // FLAGGED_ARRAY_VARIABLE_TEST_H

