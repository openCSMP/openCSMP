// ArrayVariable_Test.h

#ifndef ARRAY_VARIABLE_TEST_H
#define ARRAY_VARIABLE_TEST_H

#include "Test.h"
#include "ArrayVariable.h"

namespace csmp {

/**
@brief Comprehensive unit test for ArrayVariable.

@section coverage Coverage

  Construction:
    Default, size+value+flag, Index, PropertyDatabase, initializer_list,
    copy, move.

  Element access:
    operator[], operator(), Component, Size, Resize.

  Flag:
    Flag() get/set, single flag for whole array.

  Assignment:
    operator=(double), operator=(ScalarVariable), operator=(ArrayVariable),
    CopyValuesOnly(FlaggedArrayVariable), move assignment.

  Arithmetic (returning new object):
    operator+(double), operator-(double), operator*(double), operator/(double),
    operator+(ArrayVariable), operator-(ArrayVariable),
    operator*(ArrayVariable), operator/(ArrayVariable),
    Pow(double).

  Compound assignment:
    operator+=(double), operator-=(double), operator*=(double), operator/=(double),
    operator+=(ScalarVariable), operator-=(ScalarVariable),
    operator*=(ScalarVariable), operator/=(ScalarVariable),
    operator+=(ArrayVariable), operator-=(ArrayVariable),
    operator*=(ArrayVariable), operator/=(ArrayVariable).

  Comparison:
    operator==, operator!=, operator<.

  Queries:
    IsWithinRange, Has_NaN_Values, MinMax, Sort,
    NextLargestEntry, HasLargerEntry.

  Range interface:
    begin/end/cbegin/cend, Begin/End.

  I/O:
    Out(screen), Out(filename), Out(fstream), In(fstream).

  PropertyDatabase construction:
    Verified against known size from variable file.
*/
class ArrayVariable_Test : public Test
{
public:
    virtual void run();

private:
    void TestConstruction();
    void TestElementAccess();
    void TestFlag();
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

#endif // ARRAY_VARIABLE_TEST_H
