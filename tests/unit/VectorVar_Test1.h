#ifndef CSMP_VECTORVARTEST1_H
#define CSMP_VECTORVARTEST1_H

#include "CSMP_definitions.h"
#include "VectorVariable.h"
#include "Test.h"

namespace csmp 
{

class VectorVariable_Test1 : public Test 
{

  public:
	VectorVariable_Test1();
	~VectorVariable_Test1();
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
	void Ampersand_Operator1();
	void LessThan_Operator1();
	void Power_Operator1();
	void CrossProduct_Operator1();
	void Zero_Function1();
	void Average_Function1();
	void Length_Function1();
	void DotProduct_Function1();
	void CrossProduct_Function1();
	void IsWithinRange_Function1();
	void Fabs_Function1();
	void Ln_Function1();
	void Log10_Function1();
	void Sqrt_Function1();
	void Flip_Function1();
	void ProjectOnto_Function1();

  private:
  double64 fTolerance;
  
}; //end class

} //end csmp

#endif