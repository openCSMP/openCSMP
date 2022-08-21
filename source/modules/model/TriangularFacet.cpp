#include "TriangularFacet.h"
#include "Point.h"

namespace csmp {

double triangleArea( const Point<3U>& pt0, 
                       const Point<3U>& pt1, 
                       const Point<3U>& pt2 )
{
   return crossProduct(pt1-pt0, pt2-pt0).Length()/2.;
}


double triangleArea( const Point<2U>& pt0, 
                       const Point<2U>& pt1, 
                       const Point<2U>& pt2 )
{
   return std::fabs( 0.5 * ( pt1[0]*pt2[1] + pt0[0]*pt1[1] + 
                             pt0[1]*pt2[0] - pt2[1]*pt0[0] -
                             pt2[0]*pt1[1] - pt1[0]*pt0[1] ) );
}



/**
    creates the unit normal to the triangle defined by the points in 
    counter-clockwise order.
*/
Point<3U>  normalOfTriangle( const Point<3U>& pt0, 
                             const Point<3U>& pt1, 
                             const Point<3U>& pt2 )
{ 
  Point<3U> vecNormal( crossProduct(pt1-pt0,pt2-pt0) );

  vecNormal.NormalizeLengthTo(1.); 

  return vecNormal;
  
}

// TODO: clarify / fix this (suggestion 1(z-direction) if counter-clockwise nodes, else -1
/// normal is zero as it points into the coordinate direction that does not exist
Point<2U>  normalOfTriangle( const Point<2U>&,
                             const Point<2U>&,
                             const Point<2U>& )
 {
    return Point<2U>();
 }



} // end csmp
