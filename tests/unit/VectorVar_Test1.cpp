#include"VectorVar_Test1.h"
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

VectorVariable_Test1::VectorVariable_Test1()
{
  fTolerance = 1.e-6;
}
	
VectorVariable_Test1::~VectorVariable_Test1()
{
}
	
void VectorVariable_Test1::run()
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
	//Power_Operator1();
	Length_Function1();
	DotProduct_Function1();
	CrossProduct_Function1();
	IsWithinRange_Function1();
	Flip_Function1();
	ProjectOnto_Function1();
}


void VectorVariable_Test1::Assignment_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2;
   ScalarVariable scalar1 ( DIRICH, 10 );

   //Testing for csmpvector assignment
   csmpvector2 = csmpvector1;
   _equal(csmpvector2( 0 ), 2.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), ROBIN, fTolerance);
   
   
   //Testing for const_value assignment
   csmpvector2 = 5.0;
   _equal(csmpvector2( 0 ), 5.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), ROBIN, fTolerance);
   
   
   //Testing for scalar assignment
   csmpvector2 = scalar1;
   _equal(csmpvector2( 0 ), 10.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), DIRICH, fTolerance);
   
   
   //Testing for point assignment
   vector<double64> vector1;
   vector1.push_back( 15 );
   Point<1U> point1( vector1 );
   csmpvector2 = point1;
   _equal(csmpvector2( 0 ), 15.0, fTolerance);
   _equal(csmpvector2.Flag( 0 ), DIRICH, fTolerance);
   
}

void VectorVariable_Test1::Addition_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   
   //Testing for csmpvector addition
   csmpvector3 = csmpvector1 + csmpvector2;
   _equal(csmpvector3( 0 ), 6.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   

   //Testing for constant_value addition
   csmpvector3 = csmpvector1 + 10;
   _equal(csmpvector3( 0 ), 12.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   
   
   //Testing for point addition
   vector<double64> vector1;
   vector1.push_back( 15 );
   Point<1U> point1( vector1 );
   Point<1U> point2;
   point2 = point1 + csmpvector1;
   _equal(point2[ 0 ], 17.0, fTolerance);
      
}

void VectorVariable_Test1::Subtraction_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   
   //Testing for csmpvector subtraction
   csmpvector3 = csmpvector1 - csmpvector2;
   _equal(csmpvector3( 0 ), -2.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   

   //Testing for constant_value subtraction
   csmpvector3 = csmpvector1 - 10;
   _equal(csmpvector3( 0 ), -8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   
   
   //Testing for point subtraction
   vector<double64> vector1;
   vector1.push_back( 15 );
   Point<1U> point1( vector1 );
   Point<1U> point2;
   point2 = point1 - csmpvector1;
   _equal(point2[ 0 ], 13.0, fTolerance);
      
}

void VectorVariable_Test1::Multiplication_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   
   //Testing for csmpvector multiplication
   csmpvector3 = csmpvector1 * csmpvector2;
   _equal(csmpvector3( 0 ), 8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   

   //Testing for constant_value multiplication
   csmpvector3 = csmpvector1 * 10;
   _equal(csmpvector3( 0 ), 20.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   
   
   //Testing for point multiplication
   vector<double64> vector1;
   vector1.push_back( 5 );
   Point<1U> point1( vector1 );
   Point<1U> point2;
   point2 = point1 * csmpvector1;
   _equal(point2[ 0 ], 10.0, fTolerance);
      
}

void VectorVariable_Test1::Division_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   
   //Testing for csmpvector division
   csmpvector3 = csmpvector1 / csmpvector2;
   _equal(csmpvector3( 0 ), 0.5, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   

   //Testing for constant_value division
   csmpvector3 = csmpvector1 / 10;
   _equal(csmpvector3( 0 ), 0.2, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   
   
   //Testing for point division
   vector<double64> vector1;
   vector1.push_back( 2 );
   Point<1U> point1( vector1 );
   Point<1U> point2;
   point2 = point1 / csmpvector1;
   _equal(point2[ 0 ], 1.0, fTolerance);
      
}

void VectorVariable_Test1::Addition_Assignment_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector addition
   csmpvector3 = csmpvector1;
   csmpvector3 += csmpvector2;
   _equal(csmpvector3( 0 ), 6.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for scalar addition
   csmpvector3 = csmpvector1;
   csmpvector3 += scalar1;
   _equal(csmpvector3( 0 ), 7.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for constant_value addition
   csmpvector3 = csmpvector1;
   csmpvector3 += 10;
   _equal(csmpvector3( 0 ), 12.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

}

void VectorVariable_Test1::Subtraction_Assignment_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= csmpvector2;
   _equal(csmpvector3( 0 ), -2.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for scalar subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= scalar1;
   _equal(csmpvector3( 0 ), -3.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for constant_value subtraction
   csmpvector3 = csmpvector1;
   csmpvector3 -= 10;
   _equal(csmpvector3( 0 ), -8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

}

void VectorVariable_Test1::Multiplication_Assignment_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= csmpvector2;
   _equal(csmpvector3( 0 ), 8.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for scalar multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= scalar1;
   _equal(csmpvector3( 0 ), 10.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for constant_value multiplication
   csmpvector3 = csmpvector1;
   csmpvector3 *= 10;
   _equal(csmpvector3( 0 ), 20.0, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

}

void VectorVariable_Test1::Division_Assignment_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( DIRICH, 4.0 );
   VectorVariable<1U> csmpvector3;
   ScalarVariable scalar1 ( DIRICH, 5 );

   //Testing for csmpvector division
   csmpvector3 = csmpvector1;
   csmpvector3 /= csmpvector2;
   _equal(csmpvector3( 0 ), 0.5, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for scalar division
   csmpvector3 = csmpvector1;
   csmpvector3 /= scalar1;
   _equal(csmpvector3( 0 ), 0.4, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

   //Testing for constant_value division
   csmpvector3 = csmpvector1;
   csmpvector3 /= 10;
   _equal(csmpvector3( 0 ), 0.2, fTolerance);
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);

}

void VectorVariable_Test1::Equality_Operator1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector3( DIRICH, 2.0 );
   VectorVariable<1U> csmpvector4( DIRICH, 4.0 );
   
   _test( csmpvector1 == csmpvector2 );
   _test( !(csmpvector1 == csmpvector4) );
   _test( !(csmpvector1 != csmpvector2) );
   _test( (csmpvector2 != csmpvector3) );
   _test( csmpvector1 != csmpvector4 );
   
}




void VectorVariable_Test1::Length_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, -2.0 );
   
   _test( csmpvector1.Length() == 2.0 );
}

void VectorVariable_Test1::DotProduct_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   vector<double64> vector1;
   
   vector1.push_back( 4 );
   
   Point<1U> point1( vector1 );

   _test( csmpvector1.DotProduct( point1 ) == 8.0 );   
   
}

void VectorVariable_Test1::CrossProduct_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, 0.0 );
   VectorVariable<1U> csmpvector3;
   vector<double64> vector1;
   
   vector1.push_back( 3 );
   Point<1U> point1( vector1 );
   csmpvector3 = csmpvector1.CrossProduct( point1 );
   _test( csmpvector3 == csmpvector2 );
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
}

void VectorVariable_Test1::IsWithinRange_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   
   _test( csmpvector1.IsWithinRange( -5.0, 5.0) == true );
   _test( csmpvector1.IsWithinRange( 2.0, 10.0) == true );
   _test( csmpvector1.IsWithinRange( 1.0, 2.0) == true );
   _test( csmpvector1.IsWithinRange( 4.0, 6.0) == false );
   _test( csmpvector1.IsWithinRange( -1.0, 1.0) == false );
   
}

<<<<<<< HEAD
=======
void VectorVariable_Test1::Fabs_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, -2.0 );
   VectorVariable<1U> csmpvector3;
   
   csmpvector3 = csmpvector1;
   csmpvector3.Fabs();
   _test( csmpvector3 == csmpvector1 );
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
   
   csmpvector3 = csmpvector2;
   csmpvector3.Fabs();
   _test( csmpvector3 == csmpvector1 );
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);
      
}

void VectorVariable_Test1::Ln_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, -2.0 );
   VectorVariable<1U> csmpvector3( ROBIN, log(2.0) );
   
   csmpvector1.Ln();
   _test( csmpvector1 == csmpvector3 );
   _equal(csmpvector1.Flag( 0 ), ROBIN, fTolerance);
      
   csmpvector2.Ln();
#if 0
   // XXX FIXME Behaviour not guaranteed
   _test( isnan(csmpvector2(0)));
   _test( isnan(csmpvector2(1)));
   _test( isnan(csmpvector2(2)));
#endif
   _equal(csmpvector2.Flag( 0 ), ROBIN, fTolerance);
      
}

void VectorVariable_Test1::Log10_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, -2.0 );
   VectorVariable<1U> csmpvector3( ROBIN, log10(2.0) );
   
   csmpvector1.Log10();
   _test( csmpvector1 == csmpvector3 );
   _equal(csmpvector1.Flag( 0 ), ROBIN, fTolerance);
      
   csmpvector2.Log10();
#if 0
   // XXX FIXME Behaviour not guaranteed
   _test( isnan(csmpvector2(0)));
   _test( isnan(csmpvector2(1)));
   _test( isnan(csmpvector2(2)));
#endif
   _equal(csmpvector2.Flag( 0 ), ROBIN, fTolerance);
   
}

void VectorVariable_Test1::Sqrt_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, -2.0 );
   VectorVariable<1U> csmpvector3( ROBIN, sqrt(2.0) );
   
   csmpvector1.Sqrt();
   _test( csmpvector1 == csmpvector3 );
   _equal(csmpvector1.Flag( 0 ), ROBIN, fTolerance);
      
   csmpvector2.Sqrt();
   _test( isnan(csmpvector2(0)));
   _test( isnan(csmpvector2(1)));
   _test( isnan(csmpvector2(2)));
   _equal(csmpvector2.Flag( 0 ), ROBIN, fTolerance);
   
}
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3

void VectorVariable_Test1::Flip_Function1()
{
   VectorVariable<1U> csmpvector1( ROBIN, 2.0 );
   VectorVariable<1U> csmpvector2( ROBIN, -2.0 );
   VectorVariable<1U> csmpvector3;
   
   csmpvector3 = csmpvector1.Flip();
   _test( csmpvector3 == csmpvector2 );
   _equal(csmpvector3.Flag( 0 ), ROBIN, fTolerance);   
}

void VectorVariable_Test1::ProjectOnto_Function1()
{
   VectorVariable<1U> csmpvectorv( ROBIN, 2.0 );
   VectorVariable<1U> csmpvectoru( DIRICH, -1.0 );
   VectorVariable<1U> csmpvector;
   vector<double64> vectoru;
   
   vectoru.push_back(-1.0);
   
   csmpvector = csmpvectorv.ProjectOnto(csmpvectoru);
   _equal(csmpvector( 0 ), -1.0, fTolerance);
   _equal(csmpvector.Flag( 0 ), ROBIN, fTolerance);

   csmpvector = csmpvectorv.ProjectOnto(vectoru);
   _equal(csmpvector( 0 ), -1.0, fTolerance);
   _equal(csmpvector.Flag( 0 ), ROBIN, fTolerance);

}

} //end namespace csmp
