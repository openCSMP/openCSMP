#include "Element_Test.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp {

/// TODO: @todo check why the tolerance must be so high for the test to pass
Element_Test::Element_Test( bool verbose )
 : fTolerance(1.0e-3), verbose_(verbose)
{
}
	
  
  
Element_Test::~Element_Test()
{
}
	
  
  
/**   
    Testing:
 
      1. LengthInDirection() - quadrilateral 2D
      2. LengthInDirection() - quadrilateral 3D
 
      3. UnitNormalToFace() - actually testing FEPolicy functionality
 
*/
void Element_Test::run()
{
  ElementLengthTest2D();
  ElementLengthTest3D();
  
  UnitNormalTest();
  
  // A. Bromage test that does not rely of FE functionality (takes very long)
  //PointInVolumeElementTest();
}
	
  
  
void Element_Test::ElementLengthTest2D()
{
  IsoparametricLinearQuadrilateral fe(2U);
  csmp::Element<2U>  e( &fe );
 
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
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 9.763, fTolerance);
  
  //test 2
  //direction
  direction(0)=12.2643-8.36609;
  direction(1)=7.48967+5.60591;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 13.663, fTolerance);
  
  //test 3a
  //direction
  direction(0)=-1.85124-19.0413;
  direction(1)=7.09642-7.09642;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);
  
  //test 3b
  //direction
  direction(0)=4-10;
  direction(1)=-3+3;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);

  //test 4
  //direction
  direction(0)=0-0;
  direction(1)=3-0;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\nLength: " << fLength << std::endl;
  _equal(fLength, 13.311, fTolerance);
}
  
  
  
  
void Element_Test::ElementLengthTest3D()
{
  IsoparametricLinearQuadrilateral fe(3U);
  csmp::Element<3U> e( &fe );

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
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\n3d Length: " << fLength << std::endl;
  _equal(fLength, 9.568, fTolerance);
  
   //test 1
  //direction
  direction(0)=20;
  direction(1)=0;
  direction(2)=0;
  //compare with measured distance
  fLength = e.LengthInDirection(direction);
  if ( verbose_ ) std:: cout << "\n3d Length: " << fLength << std::endl;
  _equal(fLength, 6.964, fTolerance);
}
  
  
  
/**
    Tests FiniteVolumePolicy:
       - UnitNormalToFace( size_t face, std::vector<double64>& );
 
    Tests prism_test model because it contains elements of all
    types.
 
    @todo 2D model has to be tested as well
*/
void Element_Test::UnitNormalTest()
 {
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
      string  model_name("prism_test");
      // TODO: UnitNormalTest does not require any variables; remove property file
      ANSYS_Model3D  model( model_name.c_str(), "example25.txt");

     // ------------------------------------------------------------
     // 2. looping over all highest-dimensional elements
     //    testing whether normals are aligned with vectors
     //    between barycenter and face barycenters
     // ------------------------------------------------------------
     std::vector<double64> unrml;
     const Region<3U>& model_domain(model.Region("Model"));
     if ( verbose_ ) cout <<"\nElement_Test::UnitNormalTest: testing normal directions...\n";
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
       {
          Point<3U> bctr((*it)->BaryCenter());
          // for all the faces of the element
          for ( size_t face=0U; face<(*it)->Faces(); ++face ) {
               // constructing a vector from element to face barycenter
               Point<3U> fbctr((*it)->FaceBaryCenter( face ));
               Point<3U> outward_vec(fbctr - bctr);
               // testing that the face unit normal is aligned with the outward
               // pointing vector
               (*it)->UnitNormalToFace( face, unrml );
               Point<3U> unitnormal(unrml);
               // the normals are aligned if dotproduct is positive
               double64 dotproduct = dotProduct( outward_vec, unitnormal );
               _test( dotproduct > 0. );
               if ( verbose_ and dotproduct < 0. ) {
                    cerr <<"\nunit normal to face "<< face <<" is inward pointing:";
                    (*it)->Out();
                 }
            }
       }
   
 } // end UnitNormalTest
  
  /**
   Tests pointInVolumeElement:
   
   Tests prism_test model because it contains elements of all
   types.
   */
  void Element_Test::PointInVolumeElementTest()
  {
    ANSYS_Model3D model( "prism_test", "CSMP-variables.txt", true, true, true, true );
    
    Point<3u> query(2434.0f, -1510.0f, 5400.0f);
    
    auto& gref = model.Region("Model");
    
    auto eend = gref.ElementsEnd();
    for (auto eit = gref.ElementsBegin(); eit != eend; ++eit) {
      if (!(*eit)->IsVolumeElement()) {
        continue;
      }
      const Point<3u> bctr = (*eit)->BaryCenter();
      
      Element<3u>* e = pointInVolumeElement(gref, bctr);
      _test(e == *eit);
    }
  }
  

  
  
  
} //end namespace csmp





