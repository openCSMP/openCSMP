#include"TensorVar_Test1.h"
#include"ScalarVariable.h"
#include"VectorVariable.h"

#include "CSMP_definitions.h"
#include<vector>


using std::log;
using std::log10;
using std::sqrt;
using std::vector;
using std::cout;
using std::endl;

using namespace std;
namespace csmp
{

TensorVariable_Test1::TensorVariable_Test1()
{
  fTolerance = 1.e-16;
}
	
TensorVariable_Test1::~TensorVariable_Test1()
{
}
	
void TensorVariable_Test1::run()
{
	Assignment_Operator1();
	Addition_Operator1();
	Subtraction_Operator1();
	Multiplication_Operator1();
	Division_Operator1();
	Addition_Assignment_Operator1();
	Subtraction_Assignment_Operator1();
	Multiplication_Assignment_Operator1();
	Division_Assignment_Operator1();
	Equality_Operator1();
	Zero_Function1();
	Average_Function1();
  MinElement_Function1();
	MaxElement_Function1();
	IsWithinRange_Function1();
	Fabs_Function1();
	Ln_Function1();
	Log10_Function1();
	Sqrt_Function1();
	Adjoint_Function1();
	Identity_Function1();
	Transposed_Function1();
	Determinant_Function1();
	Inverse_Function1();
	AssignToRow_Function1();
	AssignToColumn_Function1();
	LessThan_Operator1();
	Row_Function1();
	Column_Function1();
}

void TensorVariable_Test1::Assignment_Operator1()
{
   ScalarVariable scalar1( PLAIN, 5 );
   VectorVariable<1U> vector1( ROBIN, 10.0 );
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2;
   
   // Test for Vector assignment
   tensor2 = vector1;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2.Flag(), ROBIN, fTolerance);

   
   // Test for Tensor assignment
   tensor2 = tensor1;
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);

   
   //Test for constant_value assignment  
   tensor2 = 10;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);

   
   //Test for scalar assignment
   tensor2 = scalar1;
   _equal(tensor2( 0, 0 ), 5.0, fTolerance);
   _equal(tensor2.Flag(), PLAIN, fTolerance);
      
}

void TensorVariable_Test1::Addition_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2( PLAIN, 9.0 );
   TensorVariable<1U> tensor3;
   
   //Test for tensor-tensor addition
   tensor3 = tensor1 + tensor2;
   _equal(tensor3( 0, 0 ), 10.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      
   
   //Test for tensor-constant_value addition
   tensor3 = tensor1 + 10;
   _equal(tensor3( 0, 0 ), 11.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
   
}

void TensorVariable_Test1::Subtraction_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2( PLAIN, 9.0 );
   TensorVariable<1U> tensor3;
   
   //Test for tensor-tensor subtraction
   tensor3 = tensor1 - tensor2;
   _equal(tensor3( 0, 0 ), -8.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance); 
      
   
   //Test for tensor-constant_value subtraction
   tensor3 = tensor1 - 10;
   _equal(tensor3( 0, 0 ), -9.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
 	
}

void TensorVariable_Test1::Multiplication_Operator1()
{

   TensorVariable<1U> tensor1( DIRICH, 3.0 );
   TensorVariable<1U> tensor2;
   TensorVariable<1U> tensor3( PLAIN, 9.0 );
   VectorVariable<1U> vector1( ROBIN, 2.0 );
   VectorVariable<1U> vector2;
      
   //Testing constant_value-tensor multiplication
   tensor2 = tensor1 * 10;
   _equal(tensor2( 0, 0 ), 30.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);
      
   
   //Testing vector-tensor multiplication
   vector2 = vector1 * tensor1;  
   _equal(vector2( 0 ), 6.0, fTolerance);
   _equal(vector2.Flag(), ROBIN, fTolerance);
      
   
   //Testing tensor-vector multiplication
   vector2 = tensor1 * vector1;  
   _equal(vector2( 0 ), 6.0, fTolerance);
   _equal(vector2.Flag(), DIRICH, fTolerance);
   

   //Testing tensor-tensor multiplication
   tensor2 = tensor1 * tensor3;
   _equal(tensor2( 0, 0 ), 27.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);

}

void TensorVariable_Test1::Division_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 3.0 );
   TensorVariable<1U> tensor3( ROBIN, 6.0 );
   TensorVariable<1U> tensor2;
   
   //Testing tensor-constant_value division
   tensor2 = tensor1 / 10;
   _equal(tensor2( 0, 0 ), 0.3, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);
      
   
   //Testing element by element division
   tensor2 = tensor3 / tensor1;
   _equal(tensor2( 0, 0 ), 2, fTolerance);
   _equal(tensor2.Flag(), ROBIN, fTolerance);
   
}

void TensorVariable_Test1::Zero_Function1()
{
   TensorVariable<1U> tensor1( DIRICH, 3.0 );
   
   tensor1.Zero();
   _equal(tensor1( 0, 0 ), 0, fTolerance);
   _equal(tensor1.Flag(), DIRICH, fTolerance);
   
}

void TensorVariable_Test1::Addition_Assignment_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2( PLAIN, 9.0 );
   TensorVariable<1U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Test for tensor-tensor addition
   tensor3 = tensor1;
   tensor3 += tensor2;
   _equal(tensor3( 0, 0 ), 10.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      
   
   //Test for tensor-constant_value addition
   tensor3 = tensor1;
   tensor3 += 10;
   _equal(tensor3( 0, 0 ), 11.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      

   //Test for tensor-scalar addition
   tensor3 = tensor1;
   tensor3 += scalar1;
   _equal(tensor3( 0, 0 ), 6.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
         
}

void TensorVariable_Test1::Subtraction_Assignment_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2( PLAIN, 9.0 );
   TensorVariable<1U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Test for tensor-tensor subtraction
   tensor3 = tensor1;
   tensor3 -= tensor2;
   _equal(tensor3( 0, 0 ), -8.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      
   
   //Test for tensor-constant_value subtraction
   tensor3 = tensor1;
   tensor3 -= 10;
   _equal(tensor3( 0, 0 ), -9.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
   

   //Test for tensor-scalar addition
   tensor3 = tensor1;
   tensor3 -= scalar1;
   _equal(tensor3( 0, 0 ), -4.0, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      
}

void TensorVariable_Test1::Multiplication_Assignment_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 2.0 );
   TensorVariable<1U> tensor2;
   TensorVariable<1U> tensor3( PLAIN, 9.0 );
   ScalarVariable scalar1( PLAIN, 5 );
      
   //Testing constant_value-tensor multiplication
   tensor2 = tensor1;
   tensor2 *= 10;
   _equal(tensor2( 0, 0 ), 20.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);
   

   //Testing tensor-tensor multiplication
   tensor2 = tensor1;
   tensor2 *= tensor3;
   _equal(tensor2( 0, 0 ), 18.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);
      
   
   //Test for tensor-scalar addition
   tensor2 = tensor1;
   tensor2 *= scalar1;
   _equal(tensor2( 0, 0 ), 10.0, fTolerance);
   _equal(tensor2.Flag(), DIRICH, fTolerance);
   
}

void TensorVariable_Test1::Division_Assignment_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2( PLAIN, 20.0 );
   TensorVariable<1U> tensor3;
   ScalarVariable scalar1( PLAIN, 5 );
   
   //Testing tensor-constant_value division
   tensor3 = tensor1;
   tensor3 /= 10;
   _equal(tensor3( 0, 0 ), 0.1, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      
   
   //Testing element by element division
   tensor3 = tensor1;
   tensor3 /= tensor2;
   _equal(tensor3( 0, 0 ), 0.05, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
      
   
   //Testing tensor-scalar division
   tensor3 = tensor1;
   tensor3 /= scalar1;
   _equal(tensor3( 0, 0 ), 0.2, fTolerance);
   _equal(tensor3.Flag(), DIRICH, fTolerance);
   
}

void TensorVariable_Test1::Equality_Operator1()
{
   TensorVariable<1U> tensor1( DIRICH, 1.0 );
   TensorVariable<1U> tensor2( DIRICH, 1.0 );
   TensorVariable<1U> tensor3( DIRICH, 2.0 );
                                                                                    
   _test( tensor1 == tensor2 );
   _test( tensor1 != tensor3 );
}

void TensorVariable_Test1::Average_Function1()
{
   TensorVariable<1U> tensor1( PLAIN, 5.0 );
   
   _equal(tensor1.Average(), 5.0, fTolerance);

}

void TensorVariable_Test1::MinElement_Function1()
{
   TensorVariable<1U> tensor1( PLAIN, 5.0 );
   
   _equal(tensor1.MinElement(), 5.0, fTolerance);
   
}

void TensorVariable_Test1::MaxElement_Function1()
{
   TensorVariable<1U> tensor1( PLAIN, 5.0 );
   
   _equal(tensor1.MaxElement(), 5.0, fTolerance);
   
}

void TensorVariable_Test1::IsWithinRange_Function1()
{
   TensorVariable<1U> tensor1( PLAIN, 5.0 );
   
   _test( tensor1.IsWithinRange( -5.0, 2.0) == false );
   _test( tensor1.IsWithinRange( 3.0, 5.0) == true );
   _test( tensor1.IsWithinRange( 5.0, 10.0) == true );
   _test( tensor1.IsWithinRange( 7.0, 9.0) == false );
   _test( tensor1.IsWithinRange( -1.0, 10.0) == true );
   
}	

void TensorVariable_Test1::Fabs_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, -4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor3;
   
   tensor3 = tensor1;
   tensor3.Fabs();
   _test( tensor3 == tensor2 );
   _equal(tensor3.Flag(), INIT_GUESS, fTolerance);
      
   tensor3 = tensor2;
   tensor3.Fabs();
   _test( tensor3 == tensor2 );
   _equal(tensor3.Flag(), INIT_GUESS, fTolerance);
      
}

void TensorVariable_Test1::Ln_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, -4.0 );
   TensorVariable<1U> tensor3( INIT_GUESS, log( 4.0 ) );
   TensorVariable<1U> tensor5;
   
   tensor5 = tensor1;
   tensor5.Ln(true);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
   
   tensor5 = tensor1;
   tensor5.Ln(false);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
      
   tensor5 = tensor2;
   tensor5.Ln(true);
   _test( tensor5 == tensor3 );      
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
      
   tensor5 = tensor2;
   tensor5.Ln(false);
   _test( isnan(tensor5( 0, 0 )) == true );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
     
}

void TensorVariable_Test1::Log10_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, -4.0 );
   TensorVariable<1U> tensor3( INIT_GUESS, log10( 4.0 ) );
   TensorVariable<1U> tensor5;
   
   tensor5 = tensor1;
   tensor5.Log10(true);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);

   tensor5 = tensor1;
   tensor5.Log10(false);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
      
   tensor5 = tensor2;
   tensor5.Log10(true);
   _test( tensor5 == tensor3 );      
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
      
   tensor5 = tensor2;
   tensor5.Log10(false);
   _test( isnan(tensor5( 0, 0 )) == true );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
   
}

void TensorVariable_Test1::Sqrt_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, -4.0 );
   TensorVariable<1U> tensor3( INIT_GUESS, sqrt( 4.0 ) );
   TensorVariable<1U> tensor5;
   
   tensor5 = tensor1;
   tensor5.Sqrt(true);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
   
   tensor5 = tensor1;
   tensor5.Sqrt(false);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);

   tensor5 = tensor2;
   tensor5.Sqrt(true);
   _test( tensor5 == tensor3 );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
      
   tensor5 = tensor2;
   tensor5.Sqrt(false);
   _test( isnan(tensor5( 0, 0 )) == true );
   _equal(tensor5.Flag(), INIT_GUESS, fTolerance);
   
}

void TensorVariable_Test1::Adjoint_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, 1.0 );
   TensorVariable<1U> tensor3;
   
   tensor3 = tensor1.Adjoint();
   _test( tensor3 == tensor2 );
   _equal(tensor3.Flag(), INIT_GUESS, fTolerance);
      
}

void TensorVariable_Test1::Identity_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, 1.0 );
   
   tensor1.Identity();
   _test( tensor1 == tensor2 );
   _equal(tensor1.Flag(), INIT_GUESS, fTolerance);  
}

void TensorVariable_Test1::Transposed_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor3;
   
   tensor3 = tensor1.Transposed();
   _test( tensor3 == tensor2 );
   _equal(tensor3.Flag(), INIT_GUESS, fTolerance);
      
}

void TensorVariable_Test1::Determinant_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   
   _test( tensor1.Determinant() == 4.0 );
   
}

void TensorVariable_Test1::Inverse_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 5.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, 0.2 );
   TensorVariable<1U> tensor3;
    
   tensor3 = tensor1.Inverse();
   _test( tensor3 == tensor2 );
   _equal(tensor3.Flag(), INIT_GUESS, fTolerance);
   
}


void TensorVariable_Test1::AssignToRow_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2;
   VectorVariable<1U> vector1( ROBIN, 1.0 );
   
   tensor2 = tensor1;
   tensor2.AssignToRow( 0, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2.Flag(), ROBIN, fTolerance);
            
}

void TensorVariable_Test1::AssignToColumn_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2;
   VectorVariable<1U> vector1( ROBIN, 1.0 );
   
   tensor2 = tensor1;
   tensor2.AssignToColumn( 0, vector1);
   _equal(tensor2( 0, 0 ), 1.0, fTolerance);
   _equal(tensor2.Flag(), ROBIN, fTolerance);
   
}

void TensorVariable_Test1::LessThan_Operator1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   TensorVariable<1U> tensor2( INIT_GUESS, 7.0 );
   
   _test( tensor1 < tensor2 );

}

void TensorVariable_Test1::Row_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   VectorVariable<1U> vector1;
   
   vector1 = tensor1.Row( 0 );
   _equal(vector1( 0 ), 4.0, fTolerance);
   _equal(vector1.Flag(), INIT_GUESS, fTolerance);
            
}

void TensorVariable_Test1::Column_Function1()
{
   TensorVariable<1U> tensor1( INIT_GUESS, 4.0 );
   VectorVariable<1U> vector1;
   
   vector1 = tensor1.Column( 0 );
   _equal(vector1( 0 ), 4.0, fTolerance);
   _equal(vector1.Flag(), INIT_GUESS, fTolerance);

}	

} //end namespace csmp
