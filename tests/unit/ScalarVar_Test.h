#ifndef CSMP_SCALAR_VARIABLE_TEST_H
#define CSMP_SCALAR_VARIABLE_TEST_H

#include "CSMP_definitions.h"
#include "ScalarVariable.h"
#include "Test.h"

namespace csmp
{

class ScalarVariable_Test : public Test
{
public:
    ScalarVariable_Test();
    ~ScalarVariable_Test();

    virtual void run();

    // --- Operators ---
    void Assignment_Operator();           // fixed spelling
    void Addition_Operator();
    void Subtraction_Operator();
    void Multiplication_Operator();
    void Division_Operator();
    void Addition_Assignment_Operator();  // fixed spelling
    void Subtraction_Assignment_Operator();
    void Multiplication_Assignment_Operator();
    void Division_Assignment_Operator();
    void Less_Than_Operator();
    void Greater_Than_Operator();
    void Less_Than_Or_Equal_To_Operator();
    void Greater_Than_Or_Equal_To_Operator();
    void Equality_Operator();

    // --- Functions ---
    void IsWithinRange_Function();

    /** Tests Has_NaN_Values() — was missing from original test. */
    void HasNaN_Function();

private:
    double fTolerance;
};

} // namespace csmp

#endif // CSMP_SCALAR_VARIABLE_TEST_H

