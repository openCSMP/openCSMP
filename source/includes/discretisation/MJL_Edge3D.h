#ifndef MJL_EDGE_3D_H
#define MJL_EDGE_3D_H

#include "MJL_Point3D.h"

namespace mjl {

class Edge3D {
  public:
    Edge3D();
    Edge3D( const Edge3D& );
    Edge3D( const Point3D&, const Point3D& );
    ~Edge3D();
    Edge3D&      operator=( const Edge3D& );
    bool         operator<( const Edge3D& ) const;
    void         Set( const Point3D&, const Point3D& );
    INTERSECTION Intersect( const Triangle3D&, double& ) const;
    Point3D      Origin() const;
    Point3D      Destination() const;
    Point3D      Point( double ) const;
    void         MidPoint( Point3D& ) const;
    double       AngleTo( const Edge3D& ) const;
    double       Length() const;
    void         Normalize();

    // 90o clockwise rotation and normalization if so requested
    Edge3D&     Rot( const Point3D& rot_axis, bool normalize=false );
    
    void Out() const;

    Point3D org_;
    Point3D dest_; 
};



inline Edge3D::Edge3D() {}

inline Edge3D::~Edge3D() {}


inline Edge3D::Edge3D( const Point3D& org, const Point3D& dest )
     : org_(org), dest_(dest) 
 {}




inline Edge3D&  Edge3D::operator=( const Edge3D& e )
 {
    if ( &e == this ) return *this;
    org_  = e.org_;
    dest_ = e.dest_;
    return *this;
 }


inline  bool Edge3D::operator<( const Edge3D& m ) const
 {
    return ((org_<m.org_) || (org_==m.org_ && (dest_<m.dest_)) );
 }


inline Edge3D::Edge3D( const Edge3D& e )
 {
    *this = e;
 }

inline mjl::Point3D  Edge3D::Origin() const { return org_; }

inline mjl::Point3D  Edge3D::Destination() const { return dest_; }


inline Point3D Edge3D::Point( double t ) const
 {
    return org_ + t * (dest_ - org_);
 }



inline void Edge3D::Set( const Point3D& p1, const Point3D& p2 )
 {
    org_  = p1;
    dest_ = p2;
 }


inline double Edge3D::Length() const
 {
    Point3D p = dest_ - org_;
    return p.Length();
 }


inline void Edge3D::Normalize()
 {
    dest_  = dest_ - org_;
    org_.x_ = org_.y_ = org_.z_ = 0.0;
    double len = dest_.Length();
    dest_ /= len;
 }




inline void Edge3D::MidPoint( Point3D& mp ) const
 {
    mp.x_ = (org_.x_ + dest_.x_) / 2.;
    mp.y_ = (org_.y_ + dest_.y_) / 2.;
    mp.z_ = (org_.z_ + dest_.z_) / 2.;
 }


}

#endif















