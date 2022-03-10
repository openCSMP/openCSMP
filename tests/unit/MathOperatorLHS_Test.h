#ifndef CSMP_MATH_OPERATOR_LHS_TEST_H
#define CSMP_MATH_OPERATOR_LHS_TEST_H

#include "Test.h"
#include "Model1D.h"
#include "PropertyDatabase.h"

namespace csmp 
{

class MathOperatorLHS_Test : public Test 
{

public:
	MathOperatorLHS_Test();
	~MathOperatorLHS_Test();
  
	void run(); // runs all the tests for the class (register other methods below)
  
  void MathOperatorLHS_Ctor();
  void MathOperatorLHS_CopyCtor();
  void MathOperatorLHS_Equal();
  void MathOperatorLHS_Name();
  void MathOperatorLHS_OperandName();
  void MathOperatorLHS_BasicOperandName();
  void MathOperatorLHS_TestFunctionName();
  void MathOperatorLHS_Add();
  void MathOperatorLHS_Multiply();
  void MathOperatorLHS_AddLater();
  void MathOperatorLHS_MultiplyWithTimeIncrement();
  void MathOperatorLHS_MultiplyBy();
  void MathOperatorLHS_LumpedFormulation();
  void MathOperatorLHS_BasicOffset();
  void MathOperatorLHS_TestOffset();
  void MathOperatorLHS_ApplicationCycles();
  void MathOperatorLHS_ApplicationCycle();
  void MathOperatorLHS_MaterialOperand();
  void MathOperatorLHS_BasicOperand();
  void MathOperatorLHS_TestFunctionOperand();
    
private:
  double               fTolerance_;
  Model1D<1U>*         model_ = nullptr;
  PropertyDatabase<1>& database_;
  
}; //end class

} //end csmp

#endif
