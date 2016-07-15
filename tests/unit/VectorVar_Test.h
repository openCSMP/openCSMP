#ifndef CSMP_VECTORVARTEST_H
#define CSMP_VECTORVARTEST_H

#include "VectorVariable.h"
#include "Test.h"

namespace csmp 
{

class VectorVariable_Test : public Test 
{

  public:
	VectorVariable_Test();
	~VectorVariable_Test();
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
	void Ampersand_Operator();
	void LessThan_Operator();
	void Power_Operator();
	void CrossProduct_Operator();
	void Zero_Function();
	void Average_Function();
	void Length_Function();
	void EuclideanNormalize_Function();
	void DotProduct_Function();
	void CrossProduct_Function();
	void IsWithinRange_Function();
	void Fabs_Function();
	void Ln_Function();
	void Log10_Function();
	void Sqrt_Function();
	void Flip_Function();
	void AngleTo_Function();
	void ProjectOnto_Function();

  private:
  double64 fTolerance;
  
}; //end class

} //end csmp

#endif