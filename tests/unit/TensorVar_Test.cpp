#include "TensorVar_Test.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "CSMP_definitions.h"

#include <vector>
using std::vector;
using namespace std;
namespace csmp
{

TensorVariable_Test::TensorVariable_Test()
{
  fTolerance = 1.e-5;
}
	
TensorVariable_Test::~TensorVariable_Test()
{
}
	
void TensorVariable_Test::run()
{
	Assignment_Operator();
	Addition_Operator();
	Subtraction_Operator();
	Multiplication_Operator();
	Division_Operator();
	Addition_Assignment_Operator();
	Subtraction_Assignment_Operator();
	Multiplication_Assignment_Operator();
	Division_Assignment_Operator();
	Equality_Operator();
  MinElement_Function();
	MaxElement_Function();
	IsWithinRange_Function();
	Adjoint_Function();
	Identity_Function();
	Transposed_Function();
	Determinant_Function();
	Inverse_Function();
	DiagonalValues_Function();
	EigenValues_Function();
	Trace_Function();
	AssignToRow_Function();
	AssignToColumn_Function();
	Row_Function();
	Column_Function();
}

void TensorVariable_Test::Assignment_Operator()
{
   ScalarVariable scalar1( PLAIN, 5 );
   VectorVariable<3U> vector1( PLAIN, ROBIN, ANY,
                                          10.0, 10.0, 10.0 );
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2;
   
   // Test for Vector assignment
   tensor2 = vector1;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 0.0, fTolerance);
   _equal(tensor2( 0, 2 ), 0.0, fTolerance);
   _equal(tensor2( 1, 0 ), 0.0, fTolerance);
   _equal(tensor2( 1, 1 ), 10.0, fTolerance);
   _equal(tensor2( 1, 2 ), 0.0, fTolerance);
   _equal(tensor2( 2, 0 ), 0.0, fTolerance);
   _equal(tensor2( 2, 1 ), 0.0, fTolerance);
   _equal(tensor2( 2, 2 ), 10.0, fTolerance); 
   _test( tensor2.Flag(0) == PLAIN );
   _test( tensor2.Flag(1) == ROBIN );
   _test( tensor2.Flag(2) == ANY );


   
   // Test for Tensor assignment
   tensor2 = tensor1;
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 0, 2 ), 3.0, fTolerance);
   _equal(tensor2( 1, 0 ), 4.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 2 ), 6.0, fTolerance);
   _equal(tensor2( 2, 0 ), 7.0, fTolerance);
   _equal(tensor2( 2, 1 ), 8.0, fTolerance);
   _equal(tensor2( 2, 2 ), 9.0, fTolerance);
   _test( tensor2.Flag(0) == INIT_GUESS );
   _test( tensor2.Flag(1) == PLAIN );
   _test( tensor2.Flag(2) == NEUMANN );
   
   //Test for constant_value assignment  
   tensor2 = 10;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 10.0, fTolerance);
   _equal(tensor2( 0, 2 ), 10.0, fTolerance);
   _equal(tensor2( 1, 0 ), 10.0, fTolerance);
   _equal(tensor2( 1, 1 ), 10.0, fTolerance);
   _equal(tensor2( 1, 2 ), 10.0, fTolerance);
   _equal(tensor2( 2, 0 ), 10.0, fTolerance);
   _equal(tensor2( 2, 1 ), 10.0, fTolerance);
   _equal(tensor2( 2, 2 ), 10.0, fTolerance);  
   _test( tensor2.Flag(0) == INIT_GUESS );
   _test( tensor2.Flag(1) == PLAIN );
   _test( tensor2.Flag(2) == NEUMANN );

   
   //Test for scalar assignment
   tensor2 = scalar1;
   _equal(tensor2( 0, 0 ), 5.0, fTolerance);
   _equal(tensor2( 0, 1 ), 5.0, fTolerance);
   _equal(tensor2( 0, 2 ), 5.0, fTolerance);
   _equal(tensor2( 1, 0 ), 5.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 2 ), 5.0, fTolerance);
   _equal(tensor2( 2, 0 ), 5.0, fTolerance);
   _equal(tensor2( 2, 1 ), 5.0, fTolerance);
   _equal(tensor2( 2, 2 ), 5.0, fTolerance); 
   _test( tensor2.Flag(0) == PLAIN );
   _test( tensor2.Flag(1) == PLAIN );
   _test( tensor2.Flag(2) == PLAIN );

}

void TensorVariable_Test::Addition_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2 ( PLAIN, PLAIN, PLAIN,
                                           9.0, 8.0, 7.0,
                                           6.0, 5.0, 4.0,
                                           3.0, 2.0, 1.0);
   TensorVariable<3U> tensor3;
   
   //Test for tensor-tensor addition
   tensor3 = tensor1 + tensor2;
   _equal(tensor3( 0, 0 ), 10.0, fTolerance);
   _equal(tensor3( 0, 1 ), 10.0, fTolerance);
   _equal(tensor3( 0, 2 ), 10.0, fTolerance);
   _equal(tensor3( 1, 0 ), 10.0, fTolerance);
   _equal(tensor3( 1, 1 ), 10.0, fTolerance);
   _equal(tensor3( 1, 2 ), 10.0, fTolerance);
   _equal(tensor3( 2, 0 ), 10.0, fTolerance);
   _equal(tensor3( 2, 1 ), 10.0, fTolerance);
   _equal(tensor3( 2, 2 ), 10.0, fTolerance);
   _test( tensor3.Flag(0) == INIT_GUESS );
   _test( tensor3.Flag(1) == PLAIN );
   _test( tensor3.Flag(2) == NEUMANN );
/*
   //Test for tensor-constant_value addition
   tensor3 = tensor1 + 10;
   _equal(tensor3( 0, 0 ), 11.0, fTolerance);
   _equal(tensor3( 0, 1 ), 12.0, fTolerance);
   _equal(tensor3( 0, 2 ), 13.0, fTolerance);
   _equal(tensor3( 1, 0 ), 14.0, fTolerance);
   _equal(tensor3( 1, 1 ), 15.0, fTolerance);
   _equal(tensor3( 1, 2 ), 16.0, fTolerance);
   _equal(tensor3( 2, 0 ), 17.0, fTolerance);
   _equal(tensor3( 2, 1 ), 18.0, fTolerance);
   _equal(tensor3( 2, 2 ), 19.0, fTolerance); 
   _test( tensor3.Flag(0) == INIT_GUESS );
   _test( tensor3.Flag(1) == PLAIN );
   _test( tensor3.Flag(2) == NEUMANN );
*/
}


void TensorVariable_Test::Subtraction_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2 ( PLAIN, PLAIN, PLAIN,
                                           9.0, 8.0, 7.0,
                                           6.0, 5.0, 4.0,
                                           3.0, 2.0, 1.0);
   TensorVariable<3U> tensor3;
   
   //Test for tensor-tensor subtraction
   tensor3 = tensor1 - tensor2;     
   _equal(tensor3( 0, 0 ), -8.0, fTolerance);
   _equal(tensor3( 0, 1 ), -6.0, fTolerance);
   _equal(tensor3( 0, 2 ), -4.0, fTolerance);
   _equal(tensor3( 1, 0 ), -2.0, fTolerance);
   _equal(tensor3( 1, 1 ), 0.0, fTolerance);
   _equal(tensor3( 1, 2 ), 2.0, fTolerance);
   _equal(tensor3( 2, 0 ), 4.0, fTolerance);
   _equal(tensor3( 2, 1 ), 6.0, fTolerance);
   _equal(tensor3( 2, 2 ), 8.0, fTolerance);
   _test( tensor3.Flag(0) == INIT_GUESS );
   _test( tensor3.Flag(1) == PLAIN );
   _test( tensor3.Flag(2) == NEUMANN );   
/*
   //Test for tensor-constant_value subtraction
   tensor3 = tensor1 - 10;
   _equal(tensor3( 0, 0 ), -9.0, fTolerance);
   _equal(tensor3( 0, 1 ), -8.0, fTolerance);
   _equal(tensor3( 0, 2 ), -7.0, fTolerance);
   _equal(tensor3( 1, 0 ), -6.0, fTolerance);
   _equal(tensor3( 1, 1 ), -5.0, fTolerance);
   _equal(tensor3( 1, 2 ), -4.0, fTolerance);
   _equal(tensor3( 2, 0 ), -3.0, fTolerance);
   _equal(tensor3( 2, 1 ), -2.0, fTolerance);
   _equal(tensor3( 2, 2 ), -1.0, fTolerance);
   _test( tensor3.Flag(0) == INIT_GUESS );
   _test( tensor3.Flag(1) == PLAIN );
   _test( tensor3.Flag(2) == NEUMANN );
*/
}


void TensorVariable_Test::Multiplication_Operator()
{

   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2;
   TensorVariable<3U> tensor3 ( PLAIN, PLAIN, PLAIN,
                                           9.0, 8.0, 7.0,
                                           6.0, 5.0, 4.0,
                                           3.0, 2.0, 1.0);
   VectorVariable<3U> vector1( PLAIN, PERIODIC, ANY,
                                          1.0, 1.0, 1.0 );
   VectorVariable<3U> vector2;
      
   //Testing constant_value-tensor multiplication
   /*
   tensor2 = tensor1 * 10;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 20.0, fTolerance);
   _equal(tensor2( 0, 2 ), 30.0, fTolerance);
   _equal(tensor2( 1, 0 ), 40.0, fTolerance);
   _equal(tensor2( 1, 1 ), 50.0, fTolerance);
   _equal(tensor2( 1, 2 ), 60.0, fTolerance);
   _equal(tensor2( 2, 0 ), 70.0, fTolerance);
   _equal(tensor2( 2, 1 ), 80.0, fTolerance);
   _equal(tensor2( 2, 2 ), 90.0, fTolerance);
   _test( tensor2.Flag(0) == INIT_GUESS );
   _test( tensor2.Flag(1) == PLAIN );
   _test( tensor2.Flag(2) == NEUMANN );
   */
   //Testing vector-tensor multiplication
   vector2 = vector1 * tensor1;  
   _equal(vector2( 0 ), 12.0, fTolerance);
   _equal(vector2( 1 ), 15.0, fTolerance);
   _equal(vector2( 2 ), 18.0, fTolerance);
      
   
   //Testing tensor-vector multiplication
   vector2 = tensor1 * vector1;  
   _equal(vector2( 0 ), 6.0, fTolerance);
   _equal(vector2( 1 ), 15.0, fTolerance);
   _equal(vector2( 2 ), 24.0, fTolerance);
   _test( vector2.Flag(0) == PLAIN );
   _test( vector2.Flag(1) == PERIODIC );
   _test( vector2.Flag(2) == ANY );

   //Testing tensor-tensor multiplication
   tensor2 = tensor1 * tensor3;
   _equal(tensor2( 0, 0 ), 30.0, fTolerance);
   _equal(tensor2( 0, 1 ), 24.0, fTolerance);
   _equal(tensor2( 0, 2 ), 18.0, fTolerance);
   _equal(tensor2( 1, 0 ), 84.0, fTolerance);
   _equal(tensor2( 1, 1 ), 69.0, fTolerance);
   _equal(tensor2( 1, 2 ), 54.0, fTolerance);
   _equal(tensor2( 2, 0 ), 138.0, fTolerance);
   _equal(tensor2( 2, 1 ), 114.0, fTolerance);
   _equal(tensor2( 2, 2 ), 90.0, fTolerance);
   _test( tensor2.Flag(0) == INIT_GUESS );
   _test( tensor2.Flag(1) == PLAIN );
   _test( tensor2.Flag(2) == NEUMANN );
}


void TensorVariable_Test::Division_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor3( PERIODIC, PLAIN, PLAIN,
                                          10.0, 20.0, 30.0,
                                          40.0, 50.0, 60.0,
                                          70.0, 80.0, 90.0 );
   TensorVariable<3U> tensor2;
   
   //Testing element by element division
   tensor2 = tensor3 / tensor1;
   _equal(tensor2( 0, 0 ), 10, fTolerance);
   _equal(tensor2( 0, 1 ), 10, fTolerance);
   _equal(tensor2( 0, 2 ), 10, fTolerance);
   _equal(tensor2( 1, 0 ), 10, fTolerance);
   _equal(tensor2( 1, 1 ), 10, fTolerance);
   _equal(tensor2( 1, 2 ), 10, fTolerance);
   _equal(tensor2( 2, 0 ), 10, fTolerance);
   _equal(tensor2( 2, 1 ), 10, fTolerance);
   _equal(tensor2( 2, 2 ), 10, fTolerance);   
   _test( tensor2.Flag(0) == PERIODIC );
   _test( tensor2.Flag(1) == PLAIN );
   _test( tensor2.Flag(2) == PLAIN );
}


void TensorVariable_Test::Addition_Assignment_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2 ( PLAIN, PLAIN, PLAIN,
                                           9.0, 8.0, 7.0,
                                           6.0, 5.0, 4.0,
                                           3.0, 2.0, 1.0);
   TensorVariable<3U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Test for tensor-tensor addition
   tensor3 = tensor1;
   tensor3 += tensor2;
   _equal(tensor3( 0, 0 ), 10.0, fTolerance);
   _equal(tensor3( 0, 1 ), 10.0, fTolerance);
   _equal(tensor3( 0, 2 ), 10.0, fTolerance);
   _equal(tensor3( 1, 0 ), 10.0, fTolerance);
   _equal(tensor3( 1, 1 ), 10.0, fTolerance);
   _equal(tensor3( 1, 2 ), 10.0, fTolerance);
   _equal(tensor3( 2, 0 ), 10.0, fTolerance);
   _equal(tensor3( 2, 1 ), 10.0, fTolerance);
   _equal(tensor3( 2, 2 ), 10.0, fTolerance);
      
   
   //Test for tensor-constant_value addition
   tensor3 = tensor1;
   /*
   tensor3 += 10;
   _equal(tensor3( 0, 0 ), 11.0, fTolerance);
   _equal(tensor3( 0, 1 ), 12.0, fTolerance);
   _equal(tensor3( 0, 2 ), 13.0, fTolerance);
   _equal(tensor3( 1, 0 ), 14.0, fTolerance);
   _equal(tensor3( 1, 1 ), 15.0, fTolerance);
   _equal(tensor3( 1, 2 ), 16.0, fTolerance);
   _equal(tensor3( 2, 0 ), 17.0, fTolerance);
   _equal(tensor3( 2, 1 ), 18.0, fTolerance);
   _equal(tensor3( 2, 2 ), 19.0, fTolerance); 
   */

   //Test for tensor-scalar addition
   tensor3 = tensor1;
   tensor3 += scalar1;
   _equal(tensor3( 0, 0 ), 6.0, fTolerance);
   _equal(tensor3( 0, 1 ), 7.0, fTolerance);
   _equal(tensor3( 0, 2 ), 8.0, fTolerance);
   _equal(tensor3( 1, 0 ), 9.0, fTolerance);
   _equal(tensor3( 1, 1 ), 10.0, fTolerance);
   _equal(tensor3( 1, 2 ), 11.0, fTolerance);
   _equal(tensor3( 2, 0 ), 12.0, fTolerance);
   _equal(tensor3( 2, 1 ), 13.0, fTolerance);
   _equal(tensor3( 2, 2 ), 14.0, fTolerance);
      
}

void TensorVariable_Test::Subtraction_Assignment_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2 ( PLAIN, PLAIN, PLAIN,
                                           9.0, 8.0, 7.0,
                                           6.0, 5.0, 4.0,
                                           3.0, 2.0, 1.0);
   TensorVariable<3U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Test for tensor-tensor subtraction
   tensor3 = tensor1;
   tensor3 -= tensor2;
   _equal(tensor3( 0, 0 ), -8.0, fTolerance);
   _equal(tensor3( 0, 1 ), -6.0, fTolerance);
   _equal(tensor3( 0, 2 ), -4.0, fTolerance);
   _equal(tensor3( 1, 0 ), -2.0, fTolerance);
   _equal(tensor3( 1, 1 ), 0.0, fTolerance);
   _equal(tensor3( 1, 2 ), 2.0, fTolerance);
   _equal(tensor3( 2, 0 ), 4.0, fTolerance);
   _equal(tensor3( 2, 1 ), 6.0, fTolerance);
   _equal(tensor3( 2, 2 ), 8.0, fTolerance);
      
   
   //Test for tensor-scalar addition
   tensor3 = tensor1;
   tensor3 -= scalar1;
   _equal(tensor3( 0, 0 ), -4.0, fTolerance);
   _equal(tensor3( 0, 1 ), -3.0, fTolerance);
   _equal(tensor3( 0, 2 ), -2.0, fTolerance);
   _equal(tensor3( 1, 0 ), -1.0, fTolerance);
   _equal(tensor3( 1, 1 ), 0.0, fTolerance);
   _equal(tensor3( 1, 2 ), 1.0, fTolerance);
   _equal(tensor3( 2, 0 ), 2.0, fTolerance);
   _equal(tensor3( 2, 1 ), 3.0, fTolerance);
   _equal(tensor3( 2, 2 ), 4.0, fTolerance);
      
}

void TensorVariable_Test::Multiplication_Assignment_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2;
   TensorVariable<3U> tensor3( PLAIN, PLAIN, PLAIN,
                                          9.0, 8.0, 7.0,
                                          6.0, 5.0, 4.0,
                                          3.0, 2.0, 1.0);
   ScalarVariable scalar1( PLAIN, 5 );

   //Testing tensor-tensor multiplication
   tensor2 = tensor1;
   tensor2 *= tensor3;
   _equal(tensor2( 0, 0 ), 30.0, fTolerance);
   _equal(tensor2( 0, 1 ), 24.0, fTolerance);
   _equal(tensor2( 0, 2 ), 18.0, fTolerance);
   _equal(tensor2( 1, 0 ), 84.0, fTolerance);
   _equal(tensor2( 1, 1 ), 69.0, fTolerance);
   _equal(tensor2( 1, 2 ), 54.0, fTolerance);
   _equal(tensor2( 2, 0 ), 138.0, fTolerance);
   _equal(tensor2( 2, 1 ), 114.0, fTolerance);
   _equal(tensor2( 2, 2 ), 90.0, fTolerance);
      
   
   //Test for tensor-scalar addition
   tensor2 = tensor1;
   tensor2 *= scalar1;
   _equal(tensor2( 0, 0 ), 5.0, fTolerance);
   _equal(tensor2( 0, 1 ), 10.0, fTolerance);
   _equal(tensor2( 0, 2 ), 15.0, fTolerance);
   _equal(tensor2( 1, 0 ), 20.0, fTolerance);
   _equal(tensor2( 1, 1 ), 25.0, fTolerance);
   _equal(tensor2( 1, 2 ), 30.0, fTolerance);
   _equal(tensor2( 2, 0 ), 35.0, fTolerance);
   _equal(tensor2( 2, 1 ), 40.0, fTolerance);
   _equal(tensor2( 2, 2 ), 45.0, fTolerance);
   
}

void TensorVariable_Test::Division_Assignment_Operator()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2( PLAIN, PLAIN, PLAIN,
                                          10.0, 20.0, 30.0,
                                          40.0, 50.0, 60.0,
                                          70.0, 80.0, 90.0 );
   TensorVariable<3U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Testing element by element division
   tensor3 = tensor1;
   tensor3 /= tensor2;
   _equal(tensor3( 0, 0 ), 0.1, fTolerance);
   _equal(tensor3( 0, 1 ), 0.1, fTolerance);
   _equal(tensor3( 0, 2 ), 0.1, fTolerance);
   _equal(tensor3( 1, 0 ), 0.1, fTolerance);
   _equal(tensor3( 1, 1 ), 0.1, fTolerance);
   _equal(tensor3( 1, 2 ), 0.1, fTolerance);
   _equal(tensor3( 2, 0 ), 0.1, fTolerance);
   _equal(tensor3( 2, 1 ), 0.1, fTolerance);
   _equal(tensor3( 2, 2 ), 0.1, fTolerance);
      
   
   //Testing tensor-scalar division
   tensor3 = tensor1;
   tensor3 /= scalar1;
   _equal(tensor3( 0, 0 ), 0.2, fTolerance);
   _equal(tensor3( 0, 1 ), 0.4, fTolerance);
   _equal(tensor3( 0, 2 ), 0.6, fTolerance);
   _equal(tensor3( 1, 0 ), 0.8, fTolerance);
   _equal(tensor3( 1, 1 ), 1.0, fTolerance);
   _equal(tensor3( 1, 2 ), 1.2, fTolerance);
   _equal(tensor3( 2, 0 ), 1.4, fTolerance);
   _equal(tensor3( 2, 1 ), 1.6, fTolerance);
   _equal(tensor3( 2, 2 ), 1.8, fTolerance);
   
}

void TensorVariable_Test::Equality_Operator()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor3 ( PLAIN,
                                          9.0, 8.0, 7.0,
                                          6.0, 5.0, 4.0,
                                          3.0, 2.0, 1.0);
   TensorVariable<3U> tensor4 ( ROBIN,
                                          9.0, 8.0, 7.0,
                                          6.0, 5.0, 4.0,
                                          3.0, 2.0, 1.0);
                                                                                                                             
   _test( tensor1 == tensor2 );
   _test( tensor1 != tensor3 );
   _test( tensor4 != tensor3 );
}




void TensorVariable_Test::MinElement_Function()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   _equal(tensor1.MinElement(), 1.0, fTolerance);
   
}


void TensorVariable_Test::MaxElement_Function()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   _equal(tensor1.MaxElement(), 9.0, fTolerance);
   
}


void TensorVariable_Test::IsWithinRange_Function()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   _test( tensor1.IsWithinRange( -5.0, 5.0) == false );
   _test( tensor1.IsWithinRange( 5.0, 15.0) == false );
   _test( tensor1.IsWithinRange( 1.0, 9.0) == true );
   _test( tensor1.IsWithinRange( 4.0, 6.0) == false );
   _test( tensor1.IsWithinRange( -1.0, 10.0) == true );
   
}	







void TensorVariable_Test::Adjoint_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          9.0, 1.0, 0.0,
                                          4.0, 3.0, 8.0,
                                          1.0, -1.0, 1.0 );
   TensorVariable<3U> tensor2( INIT_GUESS, PLAIN, NEUMANN,
                                          11.0, 4.0, -7.0,
                                          -1.0, 9.0, 10.0,
                                          8.0, -72.0, 23.0 );
   TensorVariable<3U> tensor3;
   
   tensor3 = tensor1.Adjoint();
   _test( tensor3 == tensor2 );
      
//   cout << tensor3 << endl;
//   cout << tensor2 << endl;
}


void TensorVariable_Test::Identity_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0 );
   
   tensor1.Identity();
   _test( tensor1 == tensor2 );
}


void TensorVariable_Test::Transposed_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 4.0, 7.0,
                                          2.0, 5.0, 8.0,
                                          3.0, 6.0, 9.0 );
   TensorVariable<3U> tensor3;
   
   tensor3 = tensor1.Transposed();
   _test( tensor3 == tensor2 );
      
}


void TensorVariable_Test::Determinant_Function()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          9.0, 1.0, 0.0,
                                          4.0, 3.0, 8.0,
                                          1.0, -1.0, 1.0 );
   TensorVariable<3U> tensor2( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   
   _test( tensor1.Determinant() == 103 );
   _test( tensor2.Determinant() == 0 );
   
}

void TensorVariable_Test::Inverse_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          9.0, 1.0, 0.0,
                                          4.0, 3.0, 8.0,
                                          1.0, -1.0, 1.0 );
  
   TensorVariable<3U> tensor2( INIT_GUESS, PLAIN, NEUMANN,
                                          11.0, -1.0, 8.0,
                                          4.0, 9.0, -72.0,
                                          -7.0, 10.0, 23.0 );
   tensor2 /= 103.;
   TensorVariable<3U> tensor3 = tensor1.Inverse();
   _equal(tensor3( 0, 0 ), tensor2( 0, 0 ), fTolerance);
   _equal(tensor3( 0, 1 ), tensor2( 0, 1 ), fTolerance);
   _equal(tensor3( 0, 2 ), tensor2( 0, 2 ), fTolerance);
   _equal(tensor3( 1, 0 ), tensor2( 1, 0 ), fTolerance);
   _equal(tensor3( 1, 1 ), tensor2( 1, 1 ), fTolerance);
   _equal(tensor3( 1, 2 ), tensor2( 1, 2 ), fTolerance);
   _equal(tensor3( 2, 0 ), tensor2( 2, 0 ), fTolerance);
   _equal(tensor3( 2, 1 ), tensor2( 2, 1 ), fTolerance);
   _equal(tensor3( 2, 2 ), tensor2( 2, 2 ), fTolerance);
   
}


void TensorVariable_Test::DiagonalValues_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   VectorVariable<3U> vectorvariable1( PLAIN, PLAIN, PLAIN,
                                          10.0, 10.0, 10.0 );
   vector<double64> vector1;
   
   //Testing for vector
   vector1.push_back(100);
   vector1.push_back(100);
   vector1.push_back(100);
   tensor1.DiagonalValues( vector1 );
   _equal(tensor1( 0, 0 ), 100.0, fTolerance);
   _equal(tensor1( 1, 1 ), 100.0, fTolerance);
   _equal(tensor1( 2, 2 ), 100.0, fTolerance);
   

   //Testing for vectorvariable
   tensor1.DiagonalValues( vectorvariable1 );
   _equal(tensor1( 0, 0 ), 10.0, fTolerance);
   _equal(tensor1( 1, 1 ), 10.0, fTolerance);
   _equal(tensor1( 2, 2 ), 10.0, fTolerance);
      

   tensor1.DiagonalValues( 2.0, 2.0, 2.0 );
   _equal(tensor1( 0, 0 ), 2.0, fTolerance);
   _equal(tensor1( 1, 1 ), 2.0, fTolerance);
   _equal(tensor1( 2, 2 ), 2.0, fTolerance); 

}

void TensorVariable_Test::EigenValues_Function()
{
   // test example verified in Maple (Karim Heinz Muxi, soil mechanics)
   TensorVariable<3U> tensor1( PLAIN,  40., -20.,  10.,
                                      -20., -80.,   5.,
                                       10.,   5.,  60. );
   // Eigen values
   VectorVariable<3U> vector1( PLAIN, PLAIN, PLAIN, 64.213203, 39.331419, -83.544623 );
  
   VectorVariable<3U> vector2;
   
   tensor1.EigenValues( vector2 );
   
   _equal(vector1( 0 ), vector2( 0 ), fTolerance);
   _equal(vector1( 1 ), vector2( 1 ), fTolerance);
   _equal(vector1( 2 ), vector2( 2 ), fTolerance);
}


void TensorVariable_Test::Trace_Function()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          9.0, 1.0, 0.0,
                                          4.0, 3.0, 8.0,
                                          1.0, -1.0, 1.0 );
   
   _test( tensor1.Trace() == 13 );
}


void TensorVariable_Test::AssignToRow_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2;
   VectorVariable<3U> vector1( PLAIN, PLAIN, PLAIN,
                                          0.0, 0.0, 0.0 );
   
   tensor2 = tensor1;
   tensor2.AssignToRow( 0, vector1);
   _equal(tensor2( 0, 0 ), 0.0, fTolerance);
   _equal(tensor2( 0, 1 ), 0.0, fTolerance);
   _equal(tensor2( 0, 2 ), 0.0, fTolerance);
   _equal(tensor2( 1, 0 ), 4.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 2 ), 6.0, fTolerance);
   _equal(tensor2( 2, 0 ), 7.0, fTolerance);
   _equal(tensor2( 2, 1 ), 8.0, fTolerance);
   _equal(tensor2( 2, 2 ), 9.0, fTolerance);
      
   tensor2 = tensor1;
   tensor2.AssignToRow( 1, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 0, 2 ), 3.0, fTolerance);
   _equal(tensor2( 1, 0 ), 0.0, fTolerance);
   _equal(tensor2( 1, 1 ), 0.0, fTolerance);
   _equal(tensor2( 1, 2 ), 0.0, fTolerance);
   _equal(tensor2( 2, 0 ), 7.0, fTolerance);
   _equal(tensor2( 2, 1 ), 8.0, fTolerance);
   _equal(tensor2( 2, 2 ), 9.0, fTolerance);
      
   tensor2 = tensor1;
   tensor2.AssignToRow( 2, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 0, 2 ), 3.0, fTolerance);
   _equal(tensor2( 1, 0 ), 4.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 2 ), 6.0, fTolerance);
   _equal(tensor2( 2, 0 ), 0.0, fTolerance);
   _equal(tensor2( 2, 1 ), 0.0, fTolerance);
   _equal(tensor2( 2, 2 ), 0.0, fTolerance);
            
}


void TensorVariable_Test::AssignToColumn_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   TensorVariable<3U> tensor2;
   VectorVariable<3U> vector1( PLAIN, PLAIN, PLAIN,
                                          0.0, 0.0, 0.0 );
   
   tensor2 = tensor1;
   tensor2.AssignToColumn( 0, vector1);
   _equal(tensor2( 0, 0 ), 0.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 0, 2 ), 3.0, fTolerance);
   _equal(tensor2( 1, 0 ), 0.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 2 ), 6.0, fTolerance);
   _equal(tensor2( 2, 0 ), 0.0, fTolerance);
   _equal(tensor2( 2, 1 ), 8.0, fTolerance);
   _equal(tensor2( 2, 2 ), 9.0, fTolerance);
      
   tensor2 = tensor1;
   tensor2.AssignToColumn( 1, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 0.0, fTolerance);
   _equal(tensor2( 0, 2 ), 3.0, fTolerance);
   _equal(tensor2( 1, 0 ), 4.0, fTolerance);
   _equal(tensor2( 1, 1 ), 0.0, fTolerance);
   _equal(tensor2( 1, 2 ), 6.0, fTolerance);
   _equal(tensor2( 2, 0 ), 7.0, fTolerance);
   _equal(tensor2( 2, 1 ), 0.0, fTolerance);
   _equal(tensor2( 2, 2 ), 9.0, fTolerance);
      
   tensor2 = tensor1;
   tensor2.AssignToColumn( 2, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 0, 2 ), 0.0, fTolerance);
   _equal(tensor2( 1, 0 ), 4.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 2 ), 0.0, fTolerance);
   _equal(tensor2( 2, 0 ), 7.0, fTolerance);
   _equal(tensor2( 2, 1 ), 8.0, fTolerance);
   _equal(tensor2( 2, 2 ), 0.0, fTolerance);
   
}


void TensorVariable_Test::LessThan_Operator()
{
   TensorVariable<3U> tensor1( PLAIN,
                                          9.0, 1.0, 0.0,
                                          4.0, 3.0, 8.0,
                                          1.0, -1.0, 1.0 );
   TensorVariable<3U> tensor2( PLAIN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   
   _test( tensor2 < tensor1 );

}


void TensorVariable_Test::Row_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   VectorVariable<3U> vector1;
   
   vector1 = tensor1.Row( 0 );
   _equal(vector1( 0 ), 1.0, fTolerance);
   _equal(vector1( 1 ), 2.0, fTolerance);
   _equal(vector1( 2 ), 3.0, fTolerance);
   
   vector1 = tensor1.Row( 1 );
   _equal(vector1( 0 ), 4.0, fTolerance);
   _equal(vector1( 1 ), 5.0, fTolerance);
   _equal(vector1( 2 ), 6.0, fTolerance);
   
   vector1 = tensor1.Row( 2 );
   _equal(vector1( 0 ), 7.0, fTolerance);
   _equal(vector1( 1 ), 8.0, fTolerance);
   _equal(vector1( 2 ), 9.0, fTolerance);
            
}


void TensorVariable_Test::Column_Function()
{
   TensorVariable<3U> tensor1( INIT_GUESS, PLAIN, NEUMANN,
                                          1.0, 2.0, 3.0,
                                          4.0, 5.0, 6.0,
                                          7.0, 8.0, 9.0 );
   VectorVariable<3U> vector1;
   
   vector1 = tensor1.Column( 0 );
   _equal(vector1( 0 ), 1.0, fTolerance);
   _equal(vector1( 1 ), 4.0, fTolerance);
   _equal(vector1( 2 ), 7.0, fTolerance);
   
   vector1 = tensor1.Column( 1 );
   _equal(vector1( 0 ), 2.0, fTolerance);
   _equal(vector1( 1 ), 5.0, fTolerance);
   _equal(vector1( 2 ), 8.0, fTolerance);
   
   vector1 = tensor1.Column( 2 );
   _equal(vector1( 0 ), 3.0, fTolerance);
   _equal(vector1( 1 ), 6.0, fTolerance);
   _equal(vector1( 2 ), 9.0, fTolerance);
   
}
	

} //end namespace csmp
