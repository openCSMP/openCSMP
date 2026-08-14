#ifndef CSMP_TENSORVARTEST1_H
#define CSMP_TENSORVARTEST1_H

#include "TensorVariable.h"
#include "Test.h"

namespace csmp 
{

class TensorVariable_Test1 : public Test 
{

  public:
	TensorVariable_Test1();
	~TensorVariable_Test1();
	void run(); // runs all the tests for the class (register other methods)
	void Assignment_Operator1();
	void Addition_Operator1();
	void Subtraction_Operator1();
	void Multiplication_Operator1();
	void Division_Operator1();
	void Addition_Assignment_Operator1();
	void Subtraction_Assignment_Operator1();
	void Multiplication_Assignment_Operator1();
	void Division_Assignment_Operator1();
	void Equality_Operator1();
	void Zero_Function1();
	void Average_Function1();
	void MinElement_Function1();
	void MaxElement_Function1();
	void IsWithinRange_Function1();
	void Fabs_Function1();
	void Ln_Function1();
	void Log10_Function1();
	void Sqrt_Function1();
	void Adjoint_Function1();
	void Identity_Function1();
	void Transposed_Function1();
	void Determinant_Function1();
	void Inverse_Function1();
	void AssignToRow_Function1();
	void AssignToColumn_Function1();
	void LessThan_Operator1();
	void Row_Function1();
	void Column_Function1();
  
  private:
  double fTolerance;
  
}; //end class

} //end csmp

#endif
