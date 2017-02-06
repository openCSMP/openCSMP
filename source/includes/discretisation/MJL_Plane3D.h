#ifndef MJL_PLANE_3D_H
#define MJL_PLANE_3D_H

#include "MJL_geometry.h"
#include "MJL_Point3D.h"
#include "MJL_Edge3D.h"

namespace mjl {

class Plane3D {
  public:
    
    Plane3D();
    Plane3D( const mjl::Point3D&, const mjl::Point3D& ); //vector and one point
    Plane3D( const mjl::Point3D&, const mjl::Point3D&, const mjl::Point3D& ); //three points
    Plane3D( const mjl::Plane3D& );
    ~Plane3D() {};
    Plane3D&  operator=( const mjl::Plane3D& m );
    mjl::Point3D PointP() const;
    mjl::Point3D NVector() const;
    mjl::Point3D p1_, p2_, p3_;
    mjl::Point3D vc_;
    void Out(std::ostream& os) ;
};

INTERSECTION  crossingPoint( double tolerance, const Edge3D&, const Plane3D&, Point3D& ); 
    

// inlined methods
inline mjl::Point3D  Plane3D::PointP() const { return p1_; }


inline mjl::Point3D  Plane3D::NVector() const { return vc_; }


inline Plane3D::Plane3D( const Point3D& vc, const mjl::Point3D& p1 )
   :  vc_(vc), p1_(p1) 
 {
 }


inline Plane3D::Plane3D( const Plane3D& m )
   : p1_(m.p1_), vc_(m.vc_)
 {
 }


inline Plane3D& Plane3D::operator=( const Plane3D& m )
 {
    if ( &m == this ) return *this;
    p1_ = m.p1_;
    vc_ = m.vc_;
    
    return *this;
 }


}

#endif






















