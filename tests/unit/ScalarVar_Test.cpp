#include"ScalarVar_Test.h"
#include"ScalarVariable.h"
#include<cmath>
#include<vector>
#include<iostream>
using std::pow;
using std::log;
using std::log10;
using std::sqrt;
using std::vector;
using std::cout;
using std::endl;

namespace csmp
{

ScalarVariable_Test::ScalarVariable_Test()
{
  fTolerance = 1.0e-16;
}
	
ScalarVariable_Test::~ScalarVariable_Test()
{
}
	
void ScalarVariable_Test::run()
{
	Assigment_Operator();
	Addition_Operator();
	Subtraction_Operator();
	Multiplication_Operator();
	Division_Operator();
	Addition_Assigment_Operator();
	Subtraction_Assigment_Operator();
	Multiplication_Assigment_Operator();
	Division_Assigment_Operator();
	Less_Than_Operator();
	Greater_Than_Operator();
	Less_Than_Or_Equal_To_Operator();
	Greater_Than_Or_Equal_To_Operator();
	Equality_Operator();
	IsWithinRange_Function();
}

void ScalarVariable_Test::Assigment_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2;
   
   //Testing for scalar assignment
   scalar2 = scalar1;
   _equal( scalar2(), scalar1(), fTolerance);
   _equal( scalar2.Flag(), scalar1.Flag(), fTolerance);

}

void ScalarVariable_Test::Addition_Operator()
{
   ScalarVariable scalar1( ANY, 5. );
   ScalarVariable scalar2( PLAIN, 15. );
   ScalarVariable scalar3;
   
   //Testing for constant_value addition
   scalar3 = scalar1 + 10.0;
   _equal( scalar3(), 15., fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   
   scalar3 = 10.0 + scalar1;
   _equal( scalar3(), 15., fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);


   //Testing for scalar addition
   scalar3 = scalar1 + scalar2;
   _equal( scalar3(), 20, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);
   
}

void ScalarVariable_Test::Subtraction_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 15 );
   ScalarVariable scalar3;
   
   //Testing for constant_value subtraction
   scalar3 = scalar1 - 10.0;
   _equal( scalar3(), -5, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   
   scalar3 = 10.0 - scalar1;
   _equal( scalar3(), 5, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);


   //Testing for scalar subtraction
   scalar3 = scalar1 - scalar2;
   _equal( scalar3(), -10, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Multiplication_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 15 );
   ScalarVariable scalar3;
   
   //Testing for constant_value multiplication
   scalar3 = scalar1 * 10.0;
   _equal( scalar3(), 50, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   
   scalar3 = 10.0 * scalar1;
   _equal( scalar3(), 50, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);


   //Testing for scalar multiplication
   scalar3 = scalar1 * scalar2;
   _equal( scalar3(), 75, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Division_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3;
   
   //Testing for constant_value division
   scalar3 = scalar1 / 10.0;
   _equal( scalar3(), 0.5, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   
   scalar3 = 10.0 / scalar1;
   _equal( scalar3(), 2, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);


   //Testing for scalar division
   scalar3 = scalar1 / scalar2;
   _equal( scalar3(), 0.25, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Addition_Assigment_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 15 );
   ScalarVariable scalar3;
   
   //Testing for constant_value addition
   scalar3 = scalar1;
   scalar3() += 10;
   _equal( scalar3(), 15, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   


   //Testing for scalar addition
   scalar3 = scalar1;
   scalar3 += scalar2;
   _equal( scalar3(), 20, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Subtraction_Assigment_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 15 );
   ScalarVariable scalar3;
   
   //Testing for constant_value subtraction
   scalar3 = scalar1;
   scalar3() -= 10;
   _equal( scalar3(), -5, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   


   //Testing for scalar subtraction
   scalar3 = scalar1;
   scalar3 -= scalar2;
   _equal( scalar3(), -10, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Multiplication_Assigment_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 15 );
   ScalarVariable scalar3;
   
   //Testing for constant_value multiplication
   scalar3 = scalar1;
   scalar3() *= 10;
   _equal( scalar3(), 50, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   


   //Testing for scalar multiplication
   scalar3 = scalar1;
   scalar3 *= scalar2;
   _equal( scalar3(), 75, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Division_Assigment_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3;
   
   //Testing for constant_value division
   scalar3    = scalar1;
   scalar3() /= 10;
   _equal( scalar3(), 0.5, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);   


   //Testing for scalar division
   scalar3 = scalar1;
   scalar3 /= scalar2;
   _equal( scalar3(), 0.25, fTolerance);
   _equal( scalar3.Flag(), ANY, fTolerance);

}

void ScalarVariable_Test::Less_Than_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3( PLAIN, 5 );
   ScalarVariable scalar4( PLAIN, 2 );
   
   //Testing for constant_value
   _test( scalar1() < 10 );
   _test( !(scalar1() < 5) );
   _test( !(scalar1() < 2) );
   
   //Testing for scalar
   _test( scalar1 < scalar2 );
   _test( !(scalar1 < scalar3) );
   _test( !(scalar1 < scalar4) );

}

void ScalarVariable_Test::Greater_Than_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3( PLAIN, 5 );
   ScalarVariable scalar4( PLAIN, 2 );
   
   //Testing for constant_value
   _test( scalar1() > 2 );
   _test( !(scalar1() > 5) );
   _test( !(scalar1() > 10) );
   
   //Testing for scalar
   _test( scalar2 > scalar1 );
   _test( !(scalar3 > scalar1) );
   _test( !(scalar4 > scalar1) );

}

void ScalarVariable_Test::Less_Than_Or_Equal_To_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3( PLAIN, 5 );
   ScalarVariable scalar4( PLAIN, 2 );
   
   //Testing for constant_value
   _test( scalar1() <= 10 );
   _test( scalar1() <= 5 );
   _test( !(scalar1() <= 2) );
   
   //Testing for scalar
   _test( scalar1 <= scalar2 );
   _test( scalar1 <= scalar3 );
   _test( !(scalar1 <= scalar4) );

}

void ScalarVariable_Test::Greater_Than_Or_Equal_To_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3( PLAIN, 5 );
   ScalarVariable scalar4( PLAIN, 2 );
   
   //Testing for constant_value
   _test( scalar1() >= 2 );
   _test( scalar1() >= 5 );
   _test( !(scalar1() >= 10) );
   
   //Testing for scalar
   _test( scalar2 >= scalar1 );
   _test( scalar3 >= scalar1 );
   _test( !(scalar4 >= scalar1) );

}

void ScalarVariable_Test::Equality_Operator()
{
   ScalarVariable scalar1( ANY, 5 );
   ScalarVariable scalar2( PLAIN, 20 );
   ScalarVariable scalar3( PLAIN, 5 );
   ScalarVariable scalar4( PLAIN, 2 );
   ScalarVariable scalar5( ANY, 5 );
         
   //Testing for scalar
   _test( !(scalar2 == scalar1) );
   _test( !(scalar3 == scalar1) );
   _test( !(scalar4 == scalar1) );
   _test( scalar5 == scalar1 );
   _test( scalar2 != scalar1 );
   _test( scalar3 != scalar1 );
   _test( scalar4 != scalar1 );
   _test( !(scalar5 != scalar1) );
   
}


void ScalarVariable_Test::IsWithinRange_Function()
{
   ScalarVariable scalar1( ANY, 5. );
   
   _test( !scalar1.IsWithinRange( 2, 3 ) );
   _test( scalar1.IsWithinRange( 2, 5 ) );
   _test( scalar1.IsWithinRange( 2, 7 ) );
   _test( scalar1.IsWithinRange( 5, 7 ) );
   _test( !scalar1.IsWithinRange( 6, 9 ) );

}

	
} //end namespace csmp
