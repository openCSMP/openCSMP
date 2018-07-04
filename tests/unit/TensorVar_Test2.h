#ifndef CSMP_TENSORVARTEST2_H
#define CSMP_TENSORVARTEST2_H

#include "TensorVariable.h"
#include "Test.h"

namespace csmp 
{

/// @todo (1-T) Put flag tests back in
class TensorVariable_Test2 : public Test 
{

  public:
	TensorVariable_Test2();
	~TensorVariable_Test2();
	void run(); // runs all the tests for the class (register other methods)
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
	void Zero_Function2();
	void Average_Function2();
	void MinElement_Function2();
	void MaxElement_Function2();
	void IsWithinRange_Function2();
	void Fabs_Function2();
	void Ln_Function2();
	void Log10_Function2();
	void Sqrt_Function2();
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
	void LessThan_Operator2();
	void Row_Function2();
	void Column_Function2();
  
  private:
  double64 fTolerance;
  
}; //end class

} //end csmp

#endif