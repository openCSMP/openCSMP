#include"VectorVar_Test2.h"
#include"ScalarVariable.h"
#include"VectorVariable.h"
#include<cmath>
#include<vector>
#include<iostream>
using std::log;
using std::log10;
using std::sqrt;
using std::vector;
using std::cout;
using std::endl;

namespace csmp
{

VectorVariable_Test2::VectorVariable_Test2()
{
  fTolerance = 1.e-6;
}
	
VectorVariable_Test2::~VectorVariable_Test2()
{
}
	
void VectorVariable_Test2::run()
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
	Ampersand_Operator2();
	LessThan_Operator2();
	Power_Operator2();
	CrossProduct_Operator2();
	Zero_Function2();
	Average_Function2();
	Length_Function2();
	EuclideanNormalize_Function2();
	DotProduct_Function2();
	CrossProduct_Function2();
	IsWithinRange_Function2();
	Fabs_Function2();
	Ln_Function2();
	Log10_Function2();
	Sqrt_Function2();
	Flip_Function2();
	AngleTo_Function2();
	ProjectOnto_Function2();

}


void VectorVariable_Test2::Assignment_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2;
   ScalarVariable scalar1 ( DIRICH, 10 );

   //Testing for csmpvector assignment
   csmpvector2 = csmpvector1;
   _equal(csmpvector2( 0 ), 1.0, fTolerance);
   _equal(csmpvector2( 1 ), 2.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 1 ), ROBIN, fTolerance);
   
   
   //Testing for const_value assignment
   csmpvector2 = 5.0;
   _equal(csmpvector2( 0 ), 5.0, fTolerance);
   _equal(csmpvector2( 1 ), 5.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 1 ), ROBIN, fTolerance);
   
   
   //Testing for scalar assignment
   csmpvector2 = scalar1;
   _equal(csmpvector2( 0 ), 10.0, fTolerance);
   _equal(csmpvector2( 1 ), 10.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector2.Flag( 1 ), DIRICH, fTolerance);
   
   
   //Testing for point assignment
   vector<double64> vector1;
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   Point<2U> point1( vector1 );
   csmpvector2 = point1;
   _equal(csmpvector2( 0 ), 15.0, fTolerance);
   _equal(csmpvector2( 1 ), 15.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector2.Flag( 1 ), DIRICH, fTolerance);
   
}

void VectorVariable_Test2::Addition_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 6.0 );
   VectorVariable<2U> csmpvector3;
   
   //Testing for csmpvector addition
   csmpvector3 = csmpvector1 + csmpvector2;
   _equal(csmpvector3( 0 ), 6.0, fTolerance);
   _equal(csmpvector3( 1 ), 8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   

   //Testing for constant_value addition
   csmpvector3 = csmpvector1 + 10;
   _equal(csmpvector3( 0 ), 11.0, fTolerance);
   _equal(csmpvector3( 1 ), 12.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   
   
   //Testing for point addition
   vector<double64> vector1;
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   Point<2U> point1( vector1 );
   Point<2U> point2;
   point2 = point1 + csmpvector1;
   _equal(point2[ 0 ], 16.0, fTolerance);
   _equal(point2[ 1 ], 17.0, fTolerance);
      
}

void VectorVariable_Test2::Subtraction_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 6.0 );
   VectorVariable<2U> csmpvector3;
   
   //Testing for csmpvector subtraction
   csmpvector3 = csmpvector1 - csmpvector2;
   _equal(csmpvector3( 0 ), -4.0, fTolerance);
   _equal(csmpvector3( 1 ), -4.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   

   //Testing for constant_value subtraction
   csmpvector3 = csmpvector1 - 10;
   _equal(csmpvector3( 0 ), -9.0, fTolerance);
   _equal(csmpvector3( 1 ), -8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   
   
   //Testing for point subtraction
   vector<double64> vector1;
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   Point<2U> point1( vector1 );
   Point<2U> point2;
   point2 = point1 - csmpvector1;
   _equal(point2[ 0 ], 14.0, fTolerance);
   _equal(point2[ 1 ], 13.0, fTolerance);
      
}

void VectorVariable_Test2::Multiplication_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 6.0 );
   VectorVariable<2U> csmpvector3;
   
   //Testing for csmpvector multiplication
   csmpvector3 = csmpvector1 * csmpvector2;
   _equal(csmpvector3( 0 ), 5.0, fTolerance);
   _equal(csmpvector3( 1 ), 12.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   

   //Testing for constant_value multiplication
   csmpvector3 = csmpvector1 * 10;
   _equal(csmpvector3( 0 ), 10.0, fTolerance);
   _equal(csmpvector3( 1 ), 20.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   
   
   //Testing for point multiplication
   vector<double64> vector1;
   vector1.push_back( 5 );
   vector1.push_back( 5 );
   Point<2U> point1( vector1 );
   Point<2U> point2;
   point2 = point1 * csmpvector1;
   _equal(point2[ 0 ], 5.0, fTolerance);
   _equal(point2[ 1 ], 10.0, fTolerance);
      
}

void VectorVariable_Test2::Division_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 4.0 );
   VectorVariable<2U> csmpvector3;
   
   //Testing for csmpvector division
   csmpvector3 = csmpvector1 / csmpvector2;
   _equal(csmpvector3( 0 ), 0.2, fTolerance);
   _equal(csmpvector3( 1 ), 0.5, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   

   //Testing for constant_value division
   csmpvector3 = csmpvector1 / 10;
   _equal(csmpvector3( 0 ), 0.1, fTolerance);
   _equal(csmpvector3( 1 ), 0.2, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   
   
   //Testing for point division
   vector<double64> vector1;
   vector1.push_back( 2 );
   vector1.push_back( 2 );
   Point<2U> point1( vector1 );
   Point<2U> point2;
   point2 = point1 / csmpvector1;
   _equal(point2[ 0 ], 2.0, fTolerance);
   _equal(point2[ 1 ], 1.0, fTolerance);
      
}

void VectorVariable_Test2::Addition_Assignment_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 6.0 );
   VectorVariable<2U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector addition
   csmpvector3 = csmpvector1;
   csmpvector3 += csmpvector2;
   _equal(csmpvector3( 0 ), 6.0, fTolerance);
   _equal(csmpvector3( 1 ), 8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for scalar addition
   csmpvector3 = csmpvector1;
   csmpvector3 += scalar1;
   _equal(csmpvector3( 0 ), 6.0, fTolerance);
   _equal(csmpvector3( 1 ), 7.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for constant_value addition
   csmpvector3 = csmpvector1;
   csmpvector3 += 10;
   _equal(csmpvector3( 0 ), 11.0, fTolerance);
   _equal(csmpvector3( 1 ), 12.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

}

void VectorVariable_Test2::Subtraction_Assignment_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 6.0 );
   VectorVariable<2U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= csmpvector2;
   _equal(csmpvector3( 0 ), -4.0, fTolerance);
   _equal(csmpvector3( 1 ), -4.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for scalar subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= scalar1;
   _equal(csmpvector3( 0 ), -4.0, fTolerance);
   _equal(csmpvector3( 1 ), -3.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for constant_value subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= 10;
   _equal(csmpvector3( 0 ), -9.0, fTolerance);
   _equal(csmpvector3( 1 ), -8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

}

void VectorVariable_Test2::Multiplication_Assignment_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 6.0 );
   VectorVariable<2U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= csmpvector2;
   _equal(csmpvector3( 0 ), 5.0, fTolerance);
   _equal(csmpvector3( 1 ), 12.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for scalar multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= scalar1;
   _equal(csmpvector3( 0 ), 5.0, fTolerance);
   _equal(csmpvector3( 1 ), 10.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for constant_value multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= 10;
   _equal(csmpvector3( 0 ), 10.0, fTolerance);
   _equal(csmpvector3( 1 ), 20.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

}

void VectorVariable_Test2::Division_Assignment_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 4.0 );
   VectorVariable<2U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector division
   csmpvector3 = csmpvector1;
   csmpvector3 /= csmpvector2;
   _equal(csmpvector3( 0 ), 0.2, fTolerance);
   _equal(csmpvector3( 1 ), 0.5, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for scalar division
   csmpvector3 = csmpvector1;
   csmpvector3 /= scalar1;
   _equal(csmpvector3( 0 ), 0.2, fTolerance);
   _equal(csmpvector3( 1 ), 0.4, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

   //Testing for constant_value division
   csmpvector3 = csmpvector1;
   csmpvector3 /= 10;
   _equal(csmpvector3( 0 ), 0.1, fTolerance);
   _equal(csmpvector3( 1 ), 0.2, fTolerance);
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);

}

void VectorVariable_Test2::Equality_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector3( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector4( PLAIN, ROBIN,
                                             5.0, 4.0 );
   
   _test( csmpvector1 != csmpvector2 );
   _test( csmpvector1 != csmpvector4 );
   _test( csmpvector1 == csmpvector3 );
   
}

void VectorVariable_Test2::Ampersand_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 4.0 );
   
   _test( (csmpvector1 & csmpvector2) == 13 );
}

void VectorVariable_Test2::LessThan_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             5.0, 4.0 );

   _test( csmpvector1 < csmpvector2 );
   _test( !(csmpvector1 < csmpvector1) );
   _test( !(csmpvector2 < csmpvector1) );
   
}

void VectorVariable_Test2::Power_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             3.0, 2.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             9.0, 4.0 );
   
   _test( (csmpvector1 ^ 2) == csmpvector2 );
   
}

void VectorVariable_Test2::CrossProduct_Operator2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             3.0, 1.0 );
   VectorVariable<2U> csmpvector3( PLAIN, ROBIN,
                                             0.0, -5.0 );
   VectorVariable<2U> csmpvector4;
   
   csmpvector4 = csmpvector1 % csmpvector2;
   _test( csmpvector4 == csmpvector3 );
   _equal(csmpvector4.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector4.Flag( 1 ), ROBIN, fTolerance);
      
}

void VectorVariable_Test2::Zero_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   
   csmpvector1.Zero();
   _equal(csmpvector1( 0 ), 0.0, fTolerance);
   _equal(csmpvector1( 1 ), 0.0, fTolerance);
   _equal(csmpvector1.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector1.Flag( 1 ), ROBIN, fTolerance);
   
}

void VectorVariable_Test2::Average_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   
   _test( csmpvector1.Average() == 1.5 );
   
}

void VectorVariable_Test2::Length_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             3.0, 4.0 );
   
   _test( csmpvector1.Length() == 5.0 );
}

void VectorVariable_Test2::EuclideanNormalize_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             3.0, 4.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             0.6, 0.8 );
   
   csmpvector1.EuclideanNormalize();
   _test( csmpvector1 == csmpvector2 );
   _equal(csmpvector1.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector1.Flag( 1 ), ROBIN, fTolerance);  
}

void VectorVariable_Test2::DotProduct_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   vector<double64> vector1;
   
   vector1.push_back( 2 );
   vector1.push_back( 4 );
   
   Point<2U> point1( vector1 );

   _test( csmpvector1.DotProduct( point1 ) == 10.0 );   
   
}

void VectorVariable_Test2::CrossProduct_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             0.0, -5.0 );
   VectorVariable<2U> csmpvector3;
   vector<double64> vector1;
   
   vector1.push_back( 3 );
   vector1.push_back( 1 );
   Point<2U> point1( vector1 );
   csmpvector3 = csmpvector1.CrossProduct( point1 );
   _test( csmpvector3 == csmpvector2 );
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
}

void VectorVariable_Test2::IsWithinRange_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 9.0 );
   
   _test( csmpvector1.IsWithinRange( -5.0, 5.0) == false );
   _test( csmpvector1.IsWithinRange( 5.0, 15.0) == false );
   _test( csmpvector1.IsWithinRange( 1.0, 9.0) == true );
   _test( csmpvector1.IsWithinRange( 4.0, 6.0) == false );
   _test( csmpvector1.IsWithinRange( -1.0, 10.0) == true );
   
}

void VectorVariable_Test2::Fabs_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             -1.0, -2.0 );
   VectorVariable<2U> csmpvector3;
   
   csmpvector3 = csmpvector1;
   csmpvector3.Fabs();
   _test( csmpvector3 == csmpvector1 );
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
   
   csmpvector3 = csmpvector2;
   csmpvector3.Fabs();
   _test( csmpvector3 == csmpvector1 );
   _equal(csmpvector3.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), ROBIN, fTolerance);
      
}

void VectorVariable_Test2::Ln_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             -1.0, -2.0 );
   VectorVariable<2U> csmpvector3( PLAIN, ROBIN,
                                             log(1.0), log(2.0) );
   
   csmpvector1.Ln();
   _test( csmpvector1 == csmpvector3 );
   _equal(csmpvector1.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector1.Flag( 1 ), ROBIN, fTolerance);
      
   csmpvector2.Ln();
   _test( isnan(csmpvector2(0)));
   _test( isnan(csmpvector2(1)));
   _test( isnan(csmpvector2(2)));
   _equal(csmpvector2.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 1 ), ROBIN, fTolerance);
      
}

void VectorVariable_Test2::Log10_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             -1.0, -2.0 );
   VectorVariable<2U> csmpvector3( PLAIN, ROBIN,
                                             log10(1.0), log10(2.0) );
   
   csmpvector1.Log10();
   _test( csmpvector1 == csmpvector3 );
   _equal(csmpvector1.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector1.Flag( 1 ), ROBIN, fTolerance);
      
   csmpvector2.Log10();
   _test( isnan(csmpvector2(0)));
   _test( isnan(csmpvector2(1)));
   _test( isnan(csmpvector2(2)));
   _equal(csmpvector2.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 1 ), ROBIN, fTolerance);
   
}

void VectorVariable_Test2::Sqrt_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             -1.0, -2.0 );
   VectorVariable<2U> csmpvector3( PLAIN, ROBIN,
                                             sqrt(1.0), sqrt(2.0) );
   
   csmpvector1.Sqrt();
   _test( csmpvector1 == csmpvector3 );
   _equal(csmpvector1.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector1.Flag( 1 ), ROBIN, fTolerance);
      
   csmpvector2.Sqrt();
   _test( isnan(csmpvector2(0)));
   _test( isnan(csmpvector2(1)));
   _test( isnan(csmpvector2(2)));
   _equal(csmpvector2.Flag( 0 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 1 ), ROBIN, fTolerance);
   
}

void VectorVariable_Test2::Flip_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 2.0 );
   VectorVariable<2U> csmpvector2( ROBIN, PLAIN,
                                             2.0, 1.0 );
   VectorVariable<2U> csmpvector3;
   
   csmpvector3 = csmpvector1.Flip();
   _test( csmpvector3 == csmpvector2 );
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);   
}

void VectorVariable_Test2::AngleTo_Function2()
{
   VectorVariable<2U> csmpvector1( PLAIN, ROBIN,
                                             1.0, 0.0 );
   VectorVariable<2U> csmpvector2( PLAIN, ROBIN,
                                             1.0, 1.0 );
   VectorVariable<2U> csmpvector3( PLAIN, ROBIN,
                                             0.0, 1.0 );
   VectorVariable<2U> csmpvector4( PLAIN, ROBIN,
                                             -1.0, 0.0 );
   VectorVariable<2U> csmpvector5( PLAIN, ROBIN,
                                             0.0, -1.0 );
   
   _equal( csmpvector1.AngleTo( csmpvector2 ), 45, fTolerance );
   _equal( csmpvector1.AngleTo( csmpvector3 ), 90, fTolerance );
   _equal( csmpvector1.AngleTo( csmpvector4 ), 180, fTolerance );
   _equal( csmpvector1.AngleTo( csmpvector5 ), 90, fTolerance );

}

void VectorVariable_Test2::ProjectOnto_Function2()
{
   VectorVariable<2U> csmpvectorv( DIRICH, PLAIN,
                                             2.0, 1.0 );
   VectorVariable<2U> csmpvectoru( ROBIN, PLAIN,
                                             -1.0, 4.0 );
   VectorVariable<2U> csmpvector;
   vector<double64> vectoru;
   
   vectoru.push_back(-1.0);
   vectoru.push_back(4.0);
   
   csmpvector = csmpvectorv.ProjectOnto(csmpvectoru);
   _equal(csmpvector( 0 ), -0.11764705879, fTolerance);
   _equal(csmpvector( 1 ), 0.47058823519, fTolerance);
   _equal(csmpvector.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector.Flag( 1 ), PLAIN, fTolerance);

   csmpvector = csmpvectorv.ProjectOnto(vectoru);
   _equal(csmpvector( 0 ), -0.11764705879, fTolerance);
   _equal(csmpvector( 1 ), 0.47058823519, fTolerance);
   _equal(csmpvector.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector.Flag( 1 ), PLAIN, fTolerance);

}

} //end namespace csmp
