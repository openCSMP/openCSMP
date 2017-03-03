#include "MJL_Point.h"
#include "MJL_Edge.h"

using namespace std;

namespace mjl {


LOCATION  Point::Classify( const Edge& e ) const
 {
    return Classify( e.org_, e.dest_ );
 }


// Laszlo, p. 76
LOCATION  Point::Classify( const Point& p0, const Point& p1 ) const
  {
      Point p2 = *this;
      Point pa = p1 - p0;
      Point pb = p2 - p0;
      double sa = pa.x_ * pb.y_;
      sa          -= pb.x_ * pa.y_;
      if ( sa > 0. ) return LEFT;
      if ( sa < 0. ) return RIGHT;
      if ( ( pa.x_ * pb.x_ < 0.) || ( pa.y_ * pb.y_ < 0.) ) return BEHIND;
      if ( pa.Length() < pb.Length() ) return BEYOND;
      if ( p0 == p2 ) return ORIGIN;
      if ( p1 == p2 ) return DESTINATION;
      return BETWEEN;
  }




double  Point::PolarAngle() const
  {
     if ( (x_ == 0.) && (y_ == 0.) ) return -1.;
     if (  x_ == 0. ) return ( (y_ > 0.) ? 90. : 270. );
     double theta = std::atan( y_/x_ );         // in radians
     theta *= 360. / (2. * 3.14159265358979324); // convert to degrees
     if ( x_ > 0. ) // quadrants 1 and 4
       return (( y_ >= 0. ) ? theta : 360. + theta);
     else // quadrants 2 and 3
       return( 180. + theta );   
  }



/**

@return 1 if positive oriented (counterclockwise)
@return -1 if clockwise
@return 0 if collinear
*/
int   Point::Orientation( Point& p0, Point& p1, Point& p2 )
 {
    Point a = p1 - p0;
    Point b = p2 - p0;
    double sa = a.x_ * b.y_; 
    sa          -= b.x_ * a.y_;
    if ( sa > 0. ) return  1;
    if ( sa < 0. ) return -1;
    return 0;
 }



 void  Point::Out() const
  {
     cout <<"\nPoint::Out: "<< x_ <<"  "<< y_ << endl;
  }
  
  
  
/**

compares points a, b with regard to their radial ordering
around p. polarCmp:

@return -1 if a > b
@return 1 if a < b
@return 0 if a = b
 */
//                   originPt          
int  polarCmp( const mjl::Point& p, const mjl::Point& a, const mjl::Point& b )
  {
     mjl::Point vp = a - p;
     mjl::Point vq = b - p;
     double pPolar = vp.PolarAngle();
     double qPolar = vq.PolarAngle();
     if ( pPolar < qPolar ) return -1;
     if ( pPolar > qPolar ) return  1;
     if ( vp.Length() < vq.Length() ) return -1;
     if ( vp.Length() > vq.Length() ) return  1;
     return 0;  
  }

} // end namespace mjl
