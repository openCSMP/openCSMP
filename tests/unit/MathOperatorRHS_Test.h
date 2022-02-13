#ifndef CSMP_MATH_OPERATOR_RHS_TEST_H
#define CSMP_MATH_OPERATOR_RHS_TEST_H

#include "Test.h"
#include "Model1D.h"
#include "PropertyDatabase.h"

namespace csmp 
{

class MathOperatorRHS_Test : public Test 
{

public:
	MathOperatorRHS_Test();
	~MathOperatorRHS_Test();
	void run(); // runs all the tests for the class (register other methods below)
  void MathOperatorRHS_Ctor();
  void MathOperatorRHS_CopyCtor();
  void MathOperatorRHS_Equal();
  void MathOperatorRHS_Name();
  void MathOperatorRHS_OperandName();
  void MathOperatorRHS_BasicOperandName();
  void MathOperatorRHS_TestFunctionName();
  void MathOperatorRHS_Subtract();
  void MathOperatorRHS_Add();
  void MathOperatorRHS_Multiply();
  void MathOperatorRHS_AddLater();
  void MathOperatorRHS_SubtractLater();
  void MathOperatorRHS_MultiplyWithTimeIncrement();
  void MathOperatorRHS_MultiplyBy();
  void MathOperatorRHS_LumpedFormulation();
  void MathOperatorRHS_ApplicationCycles();
  void MathOperatorRHS_ApplicationCycle();
  void MathOperatorRHS_Offset();
  void MathOperatorRHS_MaterialOperand();
  void MathOperatorRHS_BasicOperand();
  void MathOperatorRHS_TestFunctionOperand();
    
private:
  double fTolerance;
  Model1D<1U>* model_;
  PropertyDatabase<1>& database_;
  
}; //end class

} //end csmp

#endif
