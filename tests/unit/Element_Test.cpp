#include "Element_Test.h"


namespace csmp
{

Element_Test::Element_Test()
{
  fTolerance = 1.e-3;
}
	
Element_Test::~Element_Test()
{
}
	
void Element_Test::run()
{
  lengthMeasurement2D_Test();
  //lengthMeasurement3D_Test();
}
	
void Element_Test::lengthMeasurement2D_Test()
{
  FiniteElement * fptr = new IsoparametricLinearQuadrilateral(2U) ;
  csmp::Element<2U>  e( fptr );
 
  csmp::Node<2U> n1, n2, n3, n4;
  n1.x(8.55096); n1.y(8.59504);
  n2.x(12.3416); n2.y(0.396694);
  n3.x(5.37741); n3.y(-4.71625);
  n4.x(5.90634); n4.y(4.27548);
    
  
  e.Idx( 1 );
  e.Assign( 0, &n1 );
  e.Assign( 1, &n2 );
  e.Assign( 2, &n3 );
  e.Assign( 3, &n4 );

  VectorVariable<2U> direction(PLAIN, PLAIN);
  double64 fLength(0.);
  
  //test 1
  //direction
  direction(0)=4.31956-12.3416;
  direction(1)=8.8595+4.45179;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 9.763, fTolerance);
  
  //test 2
  //direction
  direction(0)=12.2643-8.36609;
  direction(1)=7.48967+5.60591;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 13.663, fTolerance);
  
  //test 3a
  //direction
  direction(0)=-1.85124-19.0413;
  direction(1)=7.09642-7.09642;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);
  
  //test 3b
  //direction
  direction(0)=4-10;
  direction(1)=-3+3;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);

  //test 4
  //direction
  direction(0)=0-0;
  direction(1)=3-0;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 13.311, fTolerance);

  delete fptr;

}
  
  /*
void Element_Test::lengthMeasurement3D_Test()
{
  csmp::FiniteElement = new IsoparametricLinearQuadrilateral(3U);
  csmp::Element<3U> e( fptr );

  Node<3U> n1, n2, n3, n4;
  n1.x(5.90634); n1.y(5       ); n1.z(7.);
  n2.x(10     ); n2.y(8.59504 ); n2.z(3.);
  n3.x(12.3416); n3.y(0.396694); n3.z(0.);
  n4.x(5.37741); n4.y(-4.71625); n4.z(0.);
 
  e.Idx( 1 );
  e.Assign( 0, &n1 );
  e.Assign( 1, &n2 );
  e.Assign( 2, &n3 );
  e.Assign( 3, &n4 );
  
  VectorVariable<3U> direction(PLAIN, PLAIN);
  double64 fLength(0.);
  
  //test 1
  //direction
  direction(0)=12+12;
  direction(1)=-18-18;
  direction(2)=-7-7;
  //compare with measured distance
  fLength = e.LengthOfElementInDirection(direction);
  std:: cout << "\n3d Length: " << fLength << std::endl;
  _equal(fLength, 9.568, fTolerance);
  
   //test 1
  //direction
  direction(0)=20;
  direction(1)=0;
  direction(2)=0;
  //compare with measured distance
  fLength = e.LengthOfElementInDirection(direction);
  std:: cout << "\n3d Length: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);
  delete fptr;
  
}
  */
  
} //end namespace csp
