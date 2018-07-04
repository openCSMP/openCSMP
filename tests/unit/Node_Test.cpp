#include "Node_Test.h"
#include "Node.h"
#include "Element.h"

using namespace std;

namespace csmp{

void Node_Test::run()
{

  Element<3U> e1,e2,e3;
  e1.Idx( 999 );
  e2.Idx( 9999 );
  e3.Idx( 99999 );

  // .) CONSTRUCTORS
  Node<3U> n1;

  // .) PARENTS
  n1.ResizeParentStorage( 3 );
  n1.Assign( 0, &e1 );
  n1.Assign( 1, &e2 );
  n1.Assign( 1, &e3 );
  _test( n1.Parents() == 3 );
  _test( n1.Parent( 0 ) ==  &e1 );
  _test( n1.Parent( 1 ) ==  &e2 );
  _test( n1.Parent( 2 ) ==  &e3 );
  _test( n1.ParentNodeNumber( 0 ) == 0 );
  _test( n1.ParentNodeNumber( 1 ) == 1 );
  _test( n1.ParentNodeNumber( 2 ) == 1 );

  // .) CONSTRUCTORS
  Node<3U> n2( n1 );
  _test( n1.Parents() == 3 );
  _test( n1.Parent( 0 ) ==  &e1 );
  _test( n1.Parent( 1 ) ==  &e2 );
  _test( n1.Parent( 2 ) ==  &e3 );
  _test( n1.ParentNodeNumber( 0 ) == 0 );
  _test( n1.ParentNodeNumber( 1 ) == 1 );
  _test( n1.ParentNodeNumber( 2 ) == 1 );

  // .) ERASE PARENTS
  n2.EraseParents();
  _test( n2.Parents() == 0 );

  // .) IDX
  n1.Idx( 999 );
  _test( n1.Idx() == 999 );

  // .) BOUNDRAY
  n1.AtBoundary( LEFT );
  _test( n1.AtBoundary() == LEFT );

  // .) COORDINATE / POINT
  vector<double64> p1Coordinates( 3, 9.9 );
  p1Coordinates.at( 1 ) = 99.9;
  p1Coordinates.at( 2 ) = 999.9;
  Point<3> p1( p1Coordinates );
  n1.Coordinate( p1 );
  _test( n1.Coordinate() == p1 );
  _test( n1.x() == 9.9 );
  _test( n1.y() == 99.9 );
  _test( n1.z() == 999.9 );
  n1.x( 1.1 );
  n1.y( 11.1 );
  n1.z( 111.1 );
  _test( n1.x() == 1.1 );
  _test( n1.y() == 11.1 );
  _test( n1.z() == 111.1 );
  n1[0] = 2.2;
  n1[1] = 22.2;
  n1[2] = 222.2;
  double64 n1x( n1[0] ), n1y( n1[1] ), n1z( n1[2] );
  _test( n1x == n1.x() );
  _test( n1y == n1.y() );
  _test( n1z == n1.z() );

}

} // csmp
