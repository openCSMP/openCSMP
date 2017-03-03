#include "MJL_Plane3D.h"
#include "MJL_Edge3D.h"
#include "MJL_Point3D.h"

using namespace std;

namespace mjl {



/**

@param e = a segment
@param pl = a plane {Point p; Vector vc;}
@return the intersect point (so it exists)
@return 0 = disjoint (no intersection)
@return 1 = intersection in the unique point *I0
@return  2 = the segment lies in the plane
 */
INTERSECTION  crossingPoint( double tolerance, const mjl::Edge3D& e, const mjl::Plane3D& pl, mjl::Point3D& p  )
 {
    
    mjl::Point3D u = e.Destination() - e.Origin();
    mjl::Point3D w =  e.Origin() - pl.PointP();
    mjl::Point3D e1(pl.NVector());
    double  d = u.DotProduct( e1 );
    double  nn = -w.DotProduct( e1 ); // e1=pl.NVector()

    if (fabs(d) < tolerance) {     // segment is parallel to plane
    if (fabs(nn)  < tolerance)     // segment lies in plane
        return COLLINEAR;
     else
        return PARALLEL;           // no intersection
    } 
    // If they are not parallel
    // compute intersect param
        double sI = nn / d;
    if (sI < 0. || sI > 1.)
        return SKEW_NO_CROSS;  // no intersection

    p = e.Origin() + sI * u;  // compute segment intersect point
    return SKEW_CROSS; //?
 } 


 void Plane3D::Out(std::ostream& os) const
 { 
  os <<"\nMJL_Point3D::Out: "<< p1_.x_ <<"  "<< p1_.y_ <<"  "<< p1_.z_ <<std:: endl;
  os <<"\n vector::Out: "<< vc_.x_ <<"  "<< vc_.y_ <<"  "<< vc_.z_ << std::endl;
 } 
 

} // end namespace mjl
