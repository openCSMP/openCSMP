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
    Creates  unit normal to the triangle defined by the points in
    counter-clockwise order.
    
    For example, if the counter-clockwise points are located in the XY plane, the normal will
    be point to +Z, given the right-hand coordinate system that CSMP uses.
*/
Point<3U>  normalOfTriangle( const Point<3U>& pt0, 
                             const Point<3U>& pt1, 
                             const Point<3U>& pt2 )
{ 
  Point<3U> vecNormal( crossProduct(pt1-pt0,pt2-pt0) );

  vecNormal.NormalizeLengthTo(1.); 

  return vecNormal;
  
}


/// normal is NAN as it points into the coordinate direction that does not exist
Point<2U>  normalOfTriangle( const Point<2U>&,
                             const Point<2U>&,
                             const Point<2U>& )
 {
    return Point<2U>( std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::quiet_NaN() );
 }



} // end csmp
