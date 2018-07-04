#ifndef CSMP_VECTORVARTEST2_H
#define CSMP_VECTORVARTEST2_H

#include "CSMP_definitions.h"
#include "VectorVariable.h"
#include "Test.h"

namespace csmp 
{

class VectorVariable_Test2 : public Test 
{

  public:
	VectorVariable_Test2();
	~VectorVariable_Test2();
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
	void Ampersand_Operator2();
	void LessThan_Operator2();
	void Power_Operator2();
	void CrossProduct_Operator2();
	void Zero_Function2();
	void Average_Function2();
	void Length_Function2();
	void EuclideanNormalize_Function2();
	void DotProduct_Function2();
	void CrossProduct_Function2();
	void IsWithinRange_Function2();
	void Fabs_Function2();
	void Ln_Function2();
	void Log10_Function2();
	void Sqrt_Function2();
	void Flip_Function2();
	void AngleTo_Function2();
	void ProjectOnto_Function2();

  private:
  double64 fTolerance;
  
}; //end class

} //end csmp

#endif