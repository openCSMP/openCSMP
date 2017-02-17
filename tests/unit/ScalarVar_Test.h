#ifndef CSMP_SCALARVARTEST_H
#define CSMP_SCALARVARTEST_H

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
  
    virtual void run(); // runs all the tests for the class (register other methods)
  
    void Assigment_Operator();
    void Addition_Operator();
    void Subtraction_Operator();
    void Multiplication_Operator();
    void Division_Operator();
    void Addition_Assigment_Operator();
    void Subtraction_Assigment_Operator();
    void Multiplication_Assigment_Operator();
    void Division_Assigment_Operator();
    void Less_Than_Operator();
    void Greater_Than_Operator();
    void Less_Than_Or_Equal_To_Operator();
    void Greater_Than_Or_Equal_To_Operator();
    void Equality_Operator();
    void IsWithinRange_Function();
	
  private:
    double64 fTolerance;
  
}; //end class

} //end csmp

#endif
