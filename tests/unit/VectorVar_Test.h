#ifndef CSMP_VECTOR_VARIABLE_TEST_H
#define CSMP_VECTOR_VARIABLE_TEST_H

#include "CSMP_definitions.h"
#include "VectorVariable.h"
#include "Test.h"

namespace csmp
{

class VectorVariable_Test : public Test
{
public:
    VectorVariable_Test();
    ~VectorVariable_Test();

    /** Runs all tests. */
    void run();

    // --- Arithmetic operators ---
    void Assignment_Operator();
    void Addition_Operator();
    void Subtraction_Operator();
    void Multiplication_Operator();
    void Division_Operator();

    // --- Compound assignment operators ---
    void Addition_Assignment_Operator();
    void Subtraction_Assignment_Operator();
    void Multiplication_Assignment_Operator();
    void Division_Assignment_Operator();

    // --- Comparison operators ---
    void Equality_Operator();

    /**
    Tests operator< which compares pointer addresses to satisfy
    strict weak ordering for STL associative containers.
    Does not compare values.
    */
    void LessThan_Operator();

    // --- Math functions ---

    /**
    Tests Pow(double exponent) — raises each element to the given power.
    Replaces the removed operator^ which had ambiguous C++ precedence.
    */
    void Pow_Function();

    void Length_Function();
    void EuclideanNormalize_Function();
    void DotProduct_Function();
    void CrossProduct_Function();
    void IsWithinRange_Function();

    // --- Geometric operations ---
    void Flip_Function();
    void AngleTo_Function();
    void ProjectOnto_Function();

    /**
    Tests Invert() which multiplies all elements by -1.
    Was missing from the original test.
    */
    void Invert_Function();

private:
    double fTolerance;
};

} // namespace csmp

#endif // CSMP_VECTOR_VARIABLE_TEST_H

