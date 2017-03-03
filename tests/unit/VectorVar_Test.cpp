#include "VectorVar_Test.h"
#include "vectorOperations.h"

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

VectorVariable_Test::VectorVariable_Test()
{
  fTolerance = 1.e-6;
}
	
VectorVariable_Test::~VectorVariable_Test()
{
}
	
void VectorVariable_Test::run()
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
	Length_Function();
	EuclideanNormalize_Function();
	DotProduct_Function();
	CrossProduct_Function();
	IsWithinRange_Function();
	Flip_Function();
	AngleTo_Function();
	ProjectOnto_Function();

}


void VectorVariable_Test::Assignment_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2;
   ScalarVariable scalar1 ( DIRICH, 10 );

   //Testing for csmpvector assignment
   csmpvector2 = csmpvector1;
   _equal(csmpvector2( 0 ), 1.0, fTolerance);
   _equal(csmpvector2( 1 ), 2.0, fTolerance);
   _equal(csmpvector2( 2 ), 3.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector2.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 2 ), ROBIN, fTolerance);
   
   
   //Testing for const_value assignment
   csmpvector2 = 5.0;
   _equal(csmpvector2( 0 ), 5.0, fTolerance);
   _equal(csmpvector2( 1 ), 5.0, fTolerance);
   _equal(csmpvector2( 2 ), 5.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector2.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector2.Flag( 2 ), ROBIN, fTolerance);
   
   
   //Testing for scalar assignment
   csmpvector2 = scalar1;
   _equal(csmpvector2( 0 ), 10.0, fTolerance);
   _equal(csmpvector2( 1 ), 10.0, fTolerance);
   _equal(csmpvector2( 2 ), 10.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector2.Flag( 1 ), DIRICH, fTolerance);
   _equal(csmpvector2.Flag( 2 ), DIRICH, fTolerance); 
   
   
   //Testing for point assignment
   vector<double64> vector1;
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   Point<3U> point1( vector1 );
   csmpvector2 = point1;
   _equal(csmpvector2( 0 ), 15.0, fTolerance);
   _equal(csmpvector2( 1 ), 15.0, fTolerance);
   _equal(csmpvector2( 2 ), 15.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector2.Flag( 1 ), DIRICH, fTolerance);
   _equal(csmpvector2.Flag( 2 ), DIRICH, fTolerance);
   
}

void VectorVariable_Test::Addition_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   
   //Testing for csmpvector addition
   csmpvector3 = csmpvector1 + csmpvector2;
   _equal(csmpvector3( 0 ), 5.0, fTolerance);
   _equal(csmpvector3( 1 ), 7.0, fTolerance);
   _equal(csmpvector3( 2 ), 9.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   

   //Testing for constant_value addition
   csmpvector3 = csmpvector1 + 10;
   _equal(csmpvector3( 0 ), 11.0, fTolerance);
   _equal(csmpvector3( 1 ), 12.0, fTolerance);
   _equal(csmpvector3( 2 ), 13.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   
   
   //Testing for point addition
   vector<double64> vector1;
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   Point<3U> point1( vector1 );
   Point<3U> point2;
   point2 = point1 + csmpvector1;
   _equal(point2[ 0 ], 16.0, fTolerance);
   _equal(point2[ 1 ], 17.0, fTolerance);
   _equal(point2[ 2 ], 18.0, fTolerance);
      
}

void VectorVariable_Test::Subtraction_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   
   //Testing for csmpvector subtraction
   csmpvector3 = csmpvector1 - csmpvector2;
   _equal(csmpvector3( 0 ), -3.0, fTolerance);
   _equal(csmpvector3( 1 ), -3.0, fTolerance);
   _equal(csmpvector3( 2 ), -3.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   

   //Testing for constant_value subtraction
   csmpvector3 = csmpvector1 - 10;
   _equal(csmpvector3( 0 ), -9.0, fTolerance);
   _equal(csmpvector3( 1 ), -8.0, fTolerance);
   _equal(csmpvector3( 2 ), -7.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   
   
   //Testing for point subtraction
   vector<double64> vector1;
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   vector1.push_back( 15 );
   Point<3U> point1( vector1 );
   Point<3U> point2;
   point2 = point1 - csmpvector1;
   _equal(point2[ 0 ], 14.0, fTolerance);
   _equal(point2[ 1 ], 13.0, fTolerance);
   _equal(point2[ 2 ], 12.0, fTolerance);
      
}

void VectorVariable_Test::Multiplication_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   
   //Testing for csmpvector multiplication
   csmpvector3 = csmpvector1 * csmpvector2;
   _equal(csmpvector3( 0 ), 4.0, fTolerance);
   _equal(csmpvector3( 1 ), 10.0, fTolerance);
   _equal(csmpvector3( 2 ), 18.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   

   //Testing for constant_value multiplication
   csmpvector3 = csmpvector1 * 10;
   _equal(csmpvector3( 0 ), 10.0, fTolerance);
   _equal(csmpvector3( 1 ), 20.0, fTolerance);
   _equal(csmpvector3( 2 ), 30.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   
   
   //Testing for point multiplication
   vector<double64> vector1;
   vector1.push_back( 5 );
   vector1.push_back( 5 );
   vector1.push_back( 5 );
   Point<3U> point1( vector1 );
   Point<3U> point2;
   point2 = point1 * csmpvector1;
   _equal(point2[ 0 ], 5.0, fTolerance);
   _equal(point2[ 1 ], 10.0, fTolerance);
   _equal(point2[ 2 ], 15.0, fTolerance);
      
}

void VectorVariable_Test::Division_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   
   //Testing for csmpvector division
   csmpvector3 = csmpvector1 / csmpvector2;
   _equal(csmpvector3( 0 ), 0.25, fTolerance);
   _equal(csmpvector3( 1 ), 0.4, fTolerance);
   _equal(csmpvector3( 2 ), 0.5, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   

   //Testing for constant_value division
   csmpvector3 = csmpvector1 / 10;
   _equal(csmpvector3( 0 ), 0.1, fTolerance);
   _equal(csmpvector3( 1 ), 0.2, fTolerance);
   _equal(csmpvector3( 2 ), 0.3, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
   
   
   //Testing for point division
   vector<double64> vector1;
   vector1.push_back( 2 );
   vector1.push_back( 2 );
   vector1.push_back( 3 );
   Point<3U> point1( vector1 );
   Point<3U> point2;
   point2 = point1 / csmpvector1;
   _equal(point2[ 0 ], 2.0, fTolerance);
   _equal(point2[ 1 ], 1.0, fTolerance);
   _equal(point2[ 2 ], 1.0, fTolerance);
      
}

void VectorVariable_Test::Addition_Assignment_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector addition
   csmpvector3 = csmpvector1;
   csmpvector3 += csmpvector2;
   _equal(csmpvector3( 0 ), 5.0, fTolerance);
   _equal(csmpvector3( 1 ), 7.0, fTolerance);
   _equal(csmpvector3( 2 ), 9.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for scalar addition
   csmpvector3 = csmpvector1;
   csmpvector3 += scalar1;
   _equal(csmpvector3( 0 ), 6.0, fTolerance);
   _equal(csmpvector3( 1 ), 7.0, fTolerance);
   _equal(csmpvector3( 2 ), 8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for constant_value addition
   csmpvector3 = csmpvector1;
   csmpvector3 += 10;
   _equal(csmpvector3( 0 ), 11.0, fTolerance);
   _equal(csmpvector3( 1 ), 12.0, fTolerance);
   _equal(csmpvector3( 2 ), 13.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

}

void VectorVariable_Test::Subtraction_Assignment_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= csmpvector2;
   _equal(csmpvector3( 0 ), -3.0, fTolerance);
   _equal(csmpvector3( 1 ), -3.0, fTolerance);
   _equal(csmpvector3( 2 ), -3.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for scalar subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= scalar1;
   _equal(csmpvector3( 0 ), -4.0, fTolerance);
   _equal(csmpvector3( 1 ), -3.0, fTolerance);
   _equal(csmpvector3( 2 ), -2.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for constant_value subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= 10;
   _equal(csmpvector3( 0 ), -9.0, fTolerance);
   _equal(csmpvector3( 1 ), -8.0, fTolerance);
   _equal(csmpvector3( 2 ), -7.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

}

void VectorVariable_Test::Multiplication_Assignment_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= csmpvector2;
   _equal(csmpvector3( 0 ), 4.0, fTolerance);
   _equal(csmpvector3( 1 ), 10.0, fTolerance);
   _equal(csmpvector3( 2 ), 18.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for scalar multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= scalar1;
   _equal(csmpvector3( 0 ), 5.0, fTolerance);
   _equal(csmpvector3( 1 ), 10.0, fTolerance);
   _equal(csmpvector3( 2 ), 15.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for constant_value multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= 10;
   _equal(csmpvector3( 0 ), 10.0, fTolerance);
   _equal(csmpvector3( 1 ), 20.0, fTolerance);
   _equal(csmpvector3( 2 ), 30.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

}

void VectorVariable_Test::Division_Assignment_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             4.0, 5.0, 6.0 );
   VectorVariable<3U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector division
   csmpvector3 = csmpvector1;
   csmpvector3 /= csmpvector2;
   _equal(csmpvector3( 0 ), 0.25, fTolerance);
   _equal(csmpvector3( 1 ), 0.4, fTolerance);
   _equal(csmpvector3( 2 ), 0.5, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for scalar division
   csmpvector3 = csmpvector1;
   csmpvector3 /= scalar1;
   _equal(csmpvector3( 0 ), 0.2, fTolerance);
   _equal(csmpvector3( 1 ), 0.4, fTolerance);
   _equal(csmpvector3( 2 ), 0.6, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

   //Testing for constant_value division
   csmpvector3 = csmpvector1;
   csmpvector3 /= 10;
   _equal(csmpvector3( 0 ), 0.1, fTolerance);
   _equal(csmpvector3( 1 ), 0.2, fTolerance);
   _equal(csmpvector3( 2 ), 0.3, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);

}

void VectorVariable_Test::Equality_Operator()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector4( ANY, PLAIN, ROBIN,
                                             4.0, 5.0, 6.0 );
   
   _test( csmpvector1 == csmpvector2 );
   _test( !(csmpvector1 == csmpvector4) );
   _test( !(csmpvector1 != csmpvector2) );
   _test( csmpvector1 != csmpvector4 );
   
}


void VectorVariable_Test::Length_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   
   _test( csmpvector1.Length() == sqrt( 14. ) );
}

void VectorVariable_Test::EuclideanNormalize_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ANY, PLAIN, ROBIN,
                                             1.0 / sqrt( 14. ), 2.0 / sqrt( 14. ), 3.0 / sqrt( 14. ) );
   
   csmpvector1.EuclideanNormalize();
   _test( csmpvector1 == csmpvector2 );
   _equal(csmpvector1.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector1.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector1.Flag( 2 ), ROBIN, fTolerance);  
}

void VectorVariable_Test::DotProduct_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   vector<double64> vector1;
   
   vector1.push_back( 2 );
   vector1.push_back( 4 );
   vector1.push_back( 6 );
   
   Point<3U> point1( vector1 );

   _test( csmpvector1.DotProduct( point1 ) == 28.0 );   
   
}

void VectorVariable_Test::CrossProduct_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ANY, PLAIN, ROBIN,
                                             7.0, 4.0, -5.0 );
   VectorVariable<3U> csmpvector3;
   vector<double64> vector1;
   
   vector1.push_back( 3 );
   vector1.push_back( 1 );
   vector1.push_back( 5 );
   
   Point<3U> point1( vector1 );
   
   csmpvector3 = csmpvector1.CrossProduct( point1 );
   
   _test( csmpvector3 == csmpvector2 );
   _equal(csmpvector3.Flag( 0 ), ANY, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ROBIN, fTolerance);
}


void VectorVariable_Test::IsWithinRange_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 5.0, 9.0 );
   
   _test( csmpvector1.IsWithinRange( -5.0, 5.0) == false );
   _test( csmpvector1.IsWithinRange( 5.0, 15.0) == false );
   _test( csmpvector1.IsWithinRange( 1.0, 9.0) == true );
   _test( csmpvector1.IsWithinRange( 4.0, 6.0) == false );
   _test( csmpvector1.IsWithinRange( -1.0, 10.0) == true );
   
}


void VectorVariable_Test::Flip_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 2.0, 3.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             3.0, 2.0, 1.0 );
   VectorVariable<3U> csmpvector3;
   
   csmpvector3 = csmpvector1.Flip();
   _test( csmpvector3 == csmpvector2 );
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   _equal(csmpvector3.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector3.Flag( 2 ), ANY, fTolerance);   
}

void VectorVariable_Test::AngleTo_Function()
{
   VectorVariable<3U> csmpvector1( ANY, PLAIN, ROBIN,
                                             1.0, 1.0, 1.0 );
   VectorVariable<3U> csmpvector2( ROBIN, PLAIN, ANY,
                                             -1.0, -1.0, -1.0 );
   VectorVariable<3U> csmpvector3( ANY, PLAIN, ROBIN,
                                             1.0, 1.0, 1.0 );
   VectorVariable<3U> csmpvector4( ROBIN, PLAIN, ANY,
                                             1.0, 1.0, 1.0 );
   VectorVariable<3U> csmpvector5( ANY, PLAIN, ROBIN,
                                             1.0, 1.0, 0.0 );
   VectorVariable<3U> csmpvector6( ROBIN, PLAIN, ANY,
                                             0.0, 0.0, 1.0 );
   VectorVariable<3U> csmpvector7( ANY, PLAIN, ROBIN,
                                             2.0, 1.0, 5.0 );
   VectorVariable<3U> csmpvector8( ROBIN, PLAIN, ANY,
                                             -1.0, 4.0, -2.0 );
   
   _equal( csmpvector1.AngleTo( csmpvector2 ), 180, fTolerance );
   _equal( csmpvector3.AngleTo( csmpvector4 ), 0, fTolerance );
   _equal( csmpvector5.AngleTo( csmpvector6 ), 90, fTolerance );
   _equal( csmpvector7.AngleTo( csmpvector8 ), 108.5859947, fTolerance );

}

void VectorVariable_Test::ProjectOnto_Function()
{
   VectorVariable<3U> csmpvectorv( DIRICH, PLAIN, ROBIN,
                                             2.0, 1.0, 5.0 );
   VectorVariable<3U> csmpvectoru( ROBIN, PLAIN, DIRICH,
                                             -1.0, 4.0, -2.0 );
   VectorVariable<3U> csmpvector;
   vector<double64> vectoru;
   
   vectoru.push_back(-1.0);
   vectoru.push_back(4.0);
   vectoru.push_back(-2.0);
   
   csmpvector = csmpvectorv.ProjectOnto(csmpvectoru);
   _equal(csmpvector( 0 ), 0.38095238, fTolerance);
   _equal(csmpvector( 1 ), -1.52380952, fTolerance);
   _equal(csmpvector( 2 ), 0.761904762, fTolerance);
   _equal(csmpvector.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector.Flag( 2 ), ROBIN, fTolerance);

   csmpvector = csmpvectorv.ProjectOnto(vectoru);
   _equal(csmpvector( 0 ), 0.38095238, fTolerance);
   _equal(csmpvector( 1 ), -1.52380952, fTolerance);
   _equal(csmpvector( 2 ), 0.761904762, fTolerance);
   _equal(csmpvector.Flag( 0 ), DIRICH, fTolerance);
   _equal(csmpvector.Flag( 1 ), PLAIN, fTolerance);
   _equal(csmpvector.Flag( 2 ), ROBIN, fTolerance);      
}

} //end namespace csmp
