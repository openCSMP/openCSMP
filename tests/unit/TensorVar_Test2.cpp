#include "TensorVar_Test2.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "CSMP_definitions.h"

#include<vector>

using std::log;
using std::log10;
using std::sqrt;
using std::vector;
using std::cout;
using std::endl;

namespace csmp
{

TensorVariable_Test2::TensorVariable_Test2()
{
  fTolerance = 1.e-6;
}
	
TensorVariable_Test2::~TensorVariable_Test2()
{
}
	
void TensorVariable_Test2::run()
{
	Assignment_Operator2();
	Addition_Operator2();
	Subtraction_Operator2();
	Multiplication_Operator2();
	Division_Operator2();
	Addition_Assignment_Operator2();
	Subtraction_Assignment_Operator2();
	Multiplication_Assignment_Operator2();
	Division_Assignment_Operator2();
	Equality_Operator2();
  MinElement_Function2();
	MaxElement_Function2();
	IsWithinRange_Function2();
	Adjoint_Function2();
	Identity_Function2();
	Transposed_Function2();
	Determinant_Function2();
	Inverse_Function2();
	DiagonalValues_Function2();
	EigenValues_Function2();
	Trace_Function2();
	AssignToRow_Function2();
	AssignToColumn_Function2();
	Row_Function2();
	Column_Function2();
}


void TensorVariable_Test2::Assignment_Operator2()
{
   ScalarVariable scalar1( PLAIN, 5 );
   VectorVariable<2U> vector1( PLAIN, ROBIN,
                               10.0, 10.0 );
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                               1.0, 2.0,
                               3.0, 4.0 );
   TensorVariable<2U> tensor2;
   
   // Test for Vector assignment
   tensor2 = vector1;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 0.0, fTolerance);
   _equal(tensor2( 1, 0 ), 0.0, fTolerance);
   _equal(tensor2( 1, 1 ), 10.0, fTolerance);


   
   // Test for Tensor assignment
   tensor2 = tensor1;
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 1, 0 ), 3.0, fTolerance);
   _equal(tensor2( 1, 1 ), 4.0, fTolerance);

   
   //Test for constant_value assignment  
   tensor2 = 10;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 10.0, fTolerance);
   _equal(tensor2( 1, 0 ), 10.0, fTolerance);
   _equal(tensor2( 1, 1 ), 10.0, fTolerance);

   
   //Test for scalar assignment
   tensor2 = scalar1;
   _equal(tensor2( 0, 0 ), 5.0, fTolerance);
   _equal(tensor2( 0, 1 ), 5.0, fTolerance);
   _equal(tensor2( 1, 0 ), 5.0, fTolerance);
   _equal(tensor2( 1, 1 ), 5.0, fTolerance);      
}

void TensorVariable_Test2::Addition_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2 ( PLAIN, PLAIN,
                                           9.0, 8.0,
                                           7.0, 6.0 );
   TensorVariable<2U> tensor3;
   
   //Test for tensor-tensor addition
   tensor3 = tensor1 + tensor2;
   _equal(tensor3( 0, 0 ), 10.0, fTolerance);
   _equal(tensor3( 0, 1 ), 10.0, fTolerance);
   _equal(tensor3( 1, 0 ), 10.0, fTolerance);
   _equal(tensor3( 1, 1 ), 10.0, fTolerance);
      
   
   //Test for tensor-constant_value addition
   tensor3 = tensor1 + 10;
   _equal(tensor3( 0, 0 ), 11.0, fTolerance);
   _equal(tensor3( 0, 1 ), 12.0, fTolerance);
   _equal(tensor3( 1, 0 ), 13.0, fTolerance);
   _equal(tensor3( 1, 1 ), 14.0, fTolerance); 
   
}

void TensorVariable_Test2::Subtraction_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2 ( PLAIN, PLAIN,
                                           9.0, 8.0,
                                           7.0, 6.0 );
   TensorVariable<2U> tensor3;
   
   //Test for tensor-tensor subtraction
   tensor3 = tensor1 - tensor2;     
   _equal(tensor3( 0, 0 ), -8.0, fTolerance);
   _equal(tensor3( 0, 1 ), -6.0, fTolerance);
   _equal(tensor3( 1, 0 ), -4.0, fTolerance);
   _equal(tensor3( 1, 1 ), -2.0, fTolerance);
      
   
   //Test for tensor-constant_value subtraction
   tensor3 = tensor1 - 10;
   _equal(tensor3( 0, 0 ), -9.0, fTolerance);
   _equal(tensor3( 0, 1 ), -8.0, fTolerance);
   _equal(tensor3( 1, 0 ), -7.0, fTolerance);
   _equal(tensor3( 1, 1 ), -6.0, fTolerance);  
 	
}

void TensorVariable_Test2::Multiplication_Operator2()
{

   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2;
   TensorVariable<2U> tensor3 ( PLAIN, PLAIN,
                                           9.0, 8.0,
                                           7.0, 6.0 );
   VectorVariable<2U> vector1( PLAIN, PERIODIC,
                                          1.0, 1.0 );
   VectorVariable<2U> vector2;
      
   //Testing constant_value-tensor multiplication
   tensor2 = tensor1 * 10;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 20.0, fTolerance);
   _equal(tensor2( 1, 0 ), 30.0, fTolerance);
   _equal(tensor2( 1, 1 ), 40.0, fTolerance);
      
   
   //Testing vector-tensor multiplication
   vector2 = vector1 * tensor1;  
   _equal(vector2( 0 ), 4.0, fTolerance);
   _equal(vector2( 1 ), 6.0, fTolerance);
      
   
   //Testing tensor-vector multiplication
   vector2 = tensor1 * vector1;  
   _equal(vector2( 0 ), 3.0, fTolerance);
   _equal(vector2( 1 ), 7.0, fTolerance);

   //Testing tensor-tensor multiplication
   tensor2 = tensor1 * tensor3;
   _equal(tensor2( 0, 0 ), 23.0, fTolerance);
   _equal(tensor2( 0, 1 ), 20.0, fTolerance);
   _equal(tensor2( 1, 0 ), 55.0, fTolerance);
   _equal(tensor2( 1, 1 ), 48.0, fTolerance);

}

void TensorVariable_Test2::Division_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor3( PERIODIC, PLAIN,
                                          10.0, 20.0,
                                          30.0, 40.0 );
   TensorVariable<2U> tensor2;
   
   //Testing tensor-constant_value division
   tensor2 = tensor1 / 10;
   _equal(tensor2( 0, 0 ), 0.1, fTolerance);
   _equal(tensor2( 0, 1 ), 0.2, fTolerance);
   _equal(tensor2( 1, 0 ), 0.3, fTolerance);
   _equal(tensor2( 1, 1 ), 0.4, fTolerance);
      
   
   //Testing element by element division
   tensor2 = tensor3 / tensor1;
   _equal(tensor2( 0, 0 ), 10, fTolerance);
   _equal(tensor2( 0, 1 ), 10, fTolerance);
   _equal(tensor2( 1, 0 ), 10, fTolerance);
   
}


void TensorVariable_Test2::Addition_Assignment_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2 ( PLAIN, PLAIN,
                                           9.0, 8.0,
                                           7.0, 6.0 );
   TensorVariable<2U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Test for tensor-tensor addition
   tensor3 = tensor1;
   tensor3 += tensor2;
   _equal(tensor3( 0, 0 ), 10.0, fTolerance);
   _equal(tensor3( 0, 1 ), 10.0, fTolerance);
   _equal(tensor3( 1, 0 ), 10.0, fTolerance);
   _equal(tensor3( 1, 1 ), 10.0, fTolerance);
      
   
   //Test for tensor-constant_value addition
   tensor3 = tensor1;
   tensor3 += 10;
   _equal(tensor3( 0, 0 ), 11.0, fTolerance);
   _equal(tensor3( 0, 1 ), 12.0, fTolerance);
   _equal(tensor3( 1, 0 ), 13.0, fTolerance);
   _equal(tensor3( 1, 1 ), 14.0, fTolerance);
   

   //Test for tensor-scalar addition
   tensor3 = tensor1;
   tensor3 += scalar1;
   _equal(tensor3( 0, 0 ), 6.0, fTolerance);
   _equal(tensor3( 0, 1 ), 7.0, fTolerance);
   _equal(tensor3( 1, 0 ), 8.0, fTolerance);
   _equal(tensor3( 1, 1 ), 9.0, fTolerance);
      
}

void TensorVariable_Test2::Subtraction_Assignment_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2 ( PLAIN, PLAIN,
                                           9.0, 8.0,
                                           7.0, 6.0 );
   TensorVariable<2U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Test for tensor-tensor subtraction
   tensor3 = tensor1;
   tensor3 -= tensor2;
   _equal(tensor3( 0, 0 ), -8.0, fTolerance);
   _equal(tensor3( 0, 1 ), -6.0, fTolerance);
   _equal(tensor3( 1, 0 ), -4.0, fTolerance);
   _equal(tensor3( 1, 1 ), -2.0, fTolerance);
      
   
   //Test for tensor-constant_value subtraction
   tensor3 = tensor1;
   tensor3 -= 10;
   _equal(tensor3( 0, 0 ), -9.0, fTolerance);
   _equal(tensor3( 0, 1 ), -8.0, fTolerance);
   _equal(tensor3( 1, 0 ), -7.0, fTolerance);
   _equal(tensor3( 1, 1 ), -6.0, fTolerance);
   

   //Test for tensor-scalar addition
   tensor3 = tensor1;
   tensor3 -= scalar1;
   _equal(tensor3( 0, 0 ), -4.0, fTolerance);
   _equal(tensor3( 0, 1 ), -3.0, fTolerance);
   _equal(tensor3( 1, 0 ), -2.0, fTolerance);
   _equal(tensor3( 1, 1 ), -1.0, fTolerance);
      
}

void TensorVariable_Test2::Multiplication_Assignment_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2;
   TensorVariable<2U> tensor3( PLAIN, PLAIN,
                                          9.0, 8.0,
                                          7.0, 6.0 );
   ScalarVariable scalar1( PLAIN, 5 );
      
   //Testing constant_value-tensor multiplication
   tensor2 = tensor1;
   tensor2 *= 10;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2( 0, 1 ), 20.0, fTolerance);
   _equal(tensor2( 1, 0 ), 30.0, fTolerance);
   _equal(tensor2( 1, 1 ), 40.0, fTolerance);
   

   //Testing tensor-tensor multiplication
   tensor2 = tensor1;
   tensor2 *= tensor3;
   _equal(tensor2( 0, 0 ), 23.0, fTolerance);
   _equal(tensor2( 0, 1 ), 20.0, fTolerance);
   _equal(tensor2( 1, 0 ), 55.0, fTolerance);
   _equal(tensor2( 1, 1 ), 48.0, fTolerance);
      
   
   //Test for tensor-scalar addition
   tensor2 = tensor1;
   tensor2 *= scalar1;
   _equal(tensor2( 0, 0 ), 5.0, fTolerance);
   _equal(tensor2( 0, 1 ), 10.0, fTolerance);
   _equal(tensor2( 1, 0 ), 15.0, fTolerance);
   _equal(tensor2( 1, 1 ), 20.0, fTolerance);
   
}

void TensorVariable_Test2::Division_Assignment_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( PLAIN, PLAIN,
                                          10.0, 20.0,
                                          30.0, 40.0 );
   TensorVariable<2U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Testing tensor-constant_value division
   tensor3 = tensor1;
   tensor3 /= 10;
   _equal(tensor3( 0, 0 ), 0.1, fTolerance);
   _equal(tensor3( 0, 1 ), 0.2, fTolerance);
   _equal(tensor3( 1, 0 ), 0.3, fTolerance);
   _equal(tensor3( 1, 1 ), 0.4, fTolerance);
      
   
   //Testing element by element division
   tensor3 = tensor1;
   tensor3 /= tensor2;
   _equal(tensor3( 0, 0 ), 0.1, fTolerance);
   _equal(tensor3( 0, 1 ), 0.1, fTolerance);
   _equal(tensor3( 1, 0 ), 0.1, fTolerance);
   _equal(tensor3( 1, 1 ), 0.1, fTolerance);
   
   //Testing tensor-scalar division
   tensor3 = tensor1;
   tensor3 /= scalar1;
   _equal(tensor3( 0, 0 ), 0.2, fTolerance);
   _equal(tensor3( 0, 1 ), 0.4, fTolerance);
   _equal(tensor3( 1, 0 ), 0.6, fTolerance);
   _equal(tensor3( 1, 1 ), 0.8, fTolerance);
   
}

void TensorVariable_Test2::Equality_Operator2()
{
   TensorVariable<2U> tensor1( PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor3 ( PLAIN,
                                          9.0, 8.0,
                                          7.0, 6.0 );
                                                                                    
   _test( tensor1 == tensor2 );
   _test( tensor1 != tensor3 );
}


void TensorVariable_Test2::MinElement_Function2()
{
   TensorVariable<2U> tensor1( PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
                                          
   _equal(tensor1.MinElement(), 1.0, fTolerance);
   
}

void TensorVariable_Test2::MaxElement_Function2()
{
   TensorVariable<2U> tensor1( PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
                                          
   _equal(tensor1.MaxElement(), 4.0, fTolerance);
   
}

void TensorVariable_Test2::IsWithinRange_Function2()
{
   TensorVariable<2U> tensor1( PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
                                          
   _test( tensor1.IsWithinRange( -5.0, 2.0) == false );
   _test( tensor1.IsWithinRange( 3.0, 15.0) == false );
   _test( tensor1.IsWithinRange( 1.0, 4.0) == true );
   _test( tensor1.IsWithinRange( 2.0, 3.0) == false );
   _test( tensor1.IsWithinRange( -1.0, 10.0) == true );
   
}	


void TensorVariable_Test2::Adjoint_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( INIT_GUESS, PLAIN,
                                          4.0, -3.0,
                                          -2.0, 1.0 );
   TensorVariable<2U> tensor3;
   
   tensor3 = tensor1.Adjoint();
   _test( tensor3 == tensor2 );
      
}

void TensorVariable_Test2::Identity_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( INIT_GUESS, PLAIN,
                                          1.0, 0.0,
                                          0.0, 1.0 );
   
   tensor1.Identity();
   _test( tensor1 == tensor2 ); 
}

void TensorVariable_Test2::Transposed_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( INIT_GUESS, PLAIN,
                                          1.0, 3.0,
                                          2.0, 4.0 );
   TensorVariable<2U> tensor3;
   
   tensor3 = tensor1.Transposed();
   _test( tensor3 == tensor2 );
      
}

void TensorVariable_Test2::Determinant_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   
   _test( tensor1.Determinant() == -2.0 );
   
}

void TensorVariable_Test2::Inverse_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( INIT_GUESS, PLAIN,
                                          -2.0, 1.0,
                                          1.5, -0.5 );
   TensorVariable<2U> tensor3;
    
   tensor3 = tensor1.Inverse();
   _test( tensor3 == tensor2 );
   
}

void TensorVariable_Test2::DiagonalValues_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   VectorVariable<2U> vectorvariable1( PLAIN, PLAIN,
                                                  10.0, 10.0 );
   vector<double64> vector1;
   
   //Testing for vector
   vector1.push_back(100);
   vector1.push_back(100);
   tensor1.DiagonalValues( vector1 );
   _equal(tensor1( 0, 0 ), 100.0, fTolerance);
   _equal(tensor1( 0, 1 ), 2.0, fTolerance);
   _equal(tensor1( 1, 0 ), 3.0, fTolerance);
   _equal(tensor1( 1, 1 ), 100.0, fTolerance);
   

   //Testing for vectorvariable
   tensor1.DiagonalValues( vectorvariable1 );
   _equal(tensor1( 0, 0 ), 10.0, fTolerance);
   _equal(tensor1( 0, 1 ), 2.0, fTolerance);
   _equal(tensor1( 1, 0 ), 3.0, fTolerance);
   _equal(tensor1( 1, 1 ), 10.0, fTolerance);
      

   tensor1.DiagonalValues( 5.0, 5.0 );
   _equal(tensor1( 0, 0 ), 5.0, fTolerance);
   _equal(tensor1( 0, 1 ), 2.0, fTolerance);
   _equal(tensor1( 1, 0 ), 3.0, fTolerance);
   _equal(tensor1( 1, 1 ), 5.0, fTolerance);

}

void TensorVariable_Test2::EigenValues_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   VectorVariable<2U> vector1( PLAIN, PLAIN,
                                          5.372281324, -.372281324 );
   VectorVariable<2U> vector2;
   
   tensor1.EigenValues( vector2 );
   _equal(vector1( 0 ), vector2( 0 ), fTolerance);
   _equal(vector1( 1 ), vector2( 1 ), fTolerance);
 
   
}

void TensorVariable_Test2::Trace_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 7.0 );
   
   _test( tensor1.Trace() == 8 );
}

void TensorVariable_Test2::AssignToRow_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2;
   VectorVariable<2U> vector1( PLAIN, PLAIN,
                                          0.0, 0.0 );
   
   tensor2 = tensor1;
   tensor2.AssignToRow( 0, vector1);
   _equal(tensor2( 0, 0 ), 0.0, fTolerance);
   _equal(tensor2( 0, 1 ), 0.0, fTolerance);
   _equal(tensor2( 1, 0 ), 3.0, fTolerance);
   _equal(tensor2( 1, 1 ), 4.0, fTolerance);
      
   tensor2 = tensor1;
   tensor2.AssignToRow( 1, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 1, 0 ), 0.0, fTolerance);
   _equal(tensor2( 1, 1 ), 0.0, fTolerance);
            
}

void TensorVariable_Test2::AssignToColumn_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2;
   VectorVariable<2U> vector1( PLAIN, PLAIN,
                                          0.0, 0.0 );
   
   tensor2 = tensor1;
   tensor2.AssignToColumn( 0, vector1);
   _equal(tensor2( 0, 0 ), 0.0, fTolerance);
   _equal(tensor2( 0, 1 ), 2.0, fTolerance);
   _equal(tensor2( 1, 0 ), 0.0, fTolerance);
   _equal(tensor2( 1, 1 ), 4.0, fTolerance);
      
   tensor2 = tensor1;
   tensor2.AssignToColumn( 1, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2( 0, 1 ), 0.0, fTolerance);
   _equal(tensor2( 1, 0 ), 3.0, fTolerance);
   _equal(tensor2( 1, 1 ), 0.0, fTolerance);
   
}

void TensorVariable_Test2::LessThan_Operator2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   TensorVariable<2U> tensor2( INIT_GUESS, PLAIN,
                                          6.0, 1.0,
                                          3.0, 4.0 );
   
   _test( tensor1 < tensor2 );

}

void TensorVariable_Test2::Row_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   VectorVariable<2U> vector1;
   
   vector1 = tensor1.Row( 0 );
   _equal(vector1( 0 ), 1.0, fTolerance);
   _equal(vector1( 1 ), 2.0, fTolerance);
   
   vector1 = tensor1.Row( 1 );
   _equal(vector1( 0 ), 3.0, fTolerance);
   _equal(vector1( 1 ), 4.0, fTolerance);
            
}

void TensorVariable_Test2::Column_Function2()
{
   TensorVariable<2U> tensor1( INIT_GUESS, PLAIN,
                                          1.0, 2.0,
                                          3.0, 4.0 );
   VectorVariable<2U> vector1;
   
   vector1 = tensor1.Column( 0 );
   _equal(vector1( 0 ), 1.0, fTolerance);
   _equal(vector1( 1 ), 3.0, fTolerance);
   
   vector1 = tensor1.Column( 1 );
   _equal(vector1( 0 ), 2.0, fTolerance);
   _equal(vector1( 1 ), 4.0, fTolerance);
   
}	

} //end namespace csmp
