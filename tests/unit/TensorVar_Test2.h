#ifndef CSMP_TENSOR_VARIABLE_TEST2_H
#define CSMP_TENSOR_VARIABLE_TEST2_H

#include "TensorVariable.h"
#include "Test.h"

namespace csmp
{

class TensorVariable_Test2 : public Test
{
public:
    TensorVariable_Test2();
    ~TensorVariable_Test2();

    void run();

    void Assignment_Operator2();
    void Addition_Operator2();
    void Subtraction_Operator2();
    void Multiplication_Operator2();
    void Division_Operator2();
    void Addition_Assignment_Operator2();
    void Subtraction_Assignment_Operator2();
    void Multiplication_Assignment_Operator2();
    void Division_Assignment_Operator2();
    void Equality_Operator2();

    /**
    Tests operator< which compares pointer addresses to satisfy
    strict weak ordering for STL associative containers.
    Was implemented but never called in run().
    */
    void LessThan_Operator2();

    void MinElement_Function2();
    void MaxElement_Function2();
    void IsWithinRange_Function2();
    void Adjoint_Function2();
    void Identity_Function2();
    void Transposed_Function2();
    void Determinant_Function2();
    void Inverse_Function2();
    void DiagonalValues_Function2();
    void EigenValues_Function2();
    void Trace_Function2();
    void AssignToRow_Function2();
    void AssignToColumn_Function2();
    void Row_Function2();
    void Column_Function2();

    /** Tests HadamardSquared() — component-wise T_ij^2. */
    void HadamardSquared_Function2();

    /** Tests MatrixSquared() — matrix product T*T. */
    void MatrixSquared_Function2();

    /** Tests DoubleContraction() — scalar T:T = sum_ij T_ij^2. */
    void DoubleContraction_Function2();

private:
    double fTolerance;
};

} // namespace csmp

#endif // CSMP_TENSOR_VARIABLE_TEST2_H

