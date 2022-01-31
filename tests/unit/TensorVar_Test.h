#ifndef CSMP_TENSORVARTEST_H
#define CSMP_TENSORVARTEST_H

#include "TensorVariable.h"
#include "Test.h"

namespace csmp 
{

/// @todo (1-T) Put flag tests back in
class TensorVariable_Test : public Test
{

  public:
  TensorVariable_Test();
  ~TensorVariable_Test();
	void run(); // runs all the tests for the class (register other methods)
	void Assignment_Operator();
	void Addition_Operator();
	void Subtraction_Operator();
	void Multiplication_Operator();
	void Division_Operator();
	void Addition_Assignment_Operator();
	void Subtraction_Assignment_Operator();
	void Multiplication_Assignment_Operator();
	void Division_Assignment_Operator();
	void Equality_Operator();
	void Zero_Function();
	void Average_Function();
	void MinElement_Function();
	void MaxElement_Function();
	void IsWithinRange_Function();
	void Fabs_Function();
	void Ln_Function();
	void Log10_Function();
	void Sqrt_Function();
	void Adjoint_Function();
	void Identity_Function();
	void Transposed_Function();
	void Determinant_Function();
	void Inverse_Function();
	void DiagonalValues_Function();
	void EigenValues_Function();
	void Trace_Function();
	void AssignToRow_Function();
	void AssignToColumn_Function();
	void LessThan_Operator();
	void Row_Function();
	void Column_Function();
  
  private:
  double fTolerance;
  
}; //end class

} //end csmp

#endif
