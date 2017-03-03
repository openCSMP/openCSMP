#include "MJL_Edge3D.h"
#include "MJL_Triangle3D.h"

using namespace std;

namespace mjl {

INTERSECTION  Edge3D::Intersect( const Triangle3D& p, double& t ) const
 {
     Point3D n = p.n();
     double denom = n.DotProduct(dest_ - org_);
     if ( denom == 0. )
       {
          int aclass = p.Classify(org_);
          if (aclass!=ON) return PARALLEL;
          else            return COLLINEAR;
       }
     double num = n.DotProduct(org_ - p[0]);
     t = -num / denom;
     return SKEW; 
 }



/** 90o clockwise rotation of edge

computes the edge normal as the cross product between the edge and the edge
defined by the rot-axis. The use can center the rotated unit normal=edge
on the edge midpoint by using the center_edge boolean option

@warning method works only for a normalized edge
*/
Edge3D&  Edge3D::Rot( const Point3D& ax, bool normalize )
 {
    // remembering the midpoint of the edge
    Point3D  midpoint( (org_.x_+dest_.x_)/2.,
                       (org_.y_+dest_.y_)/2.,
                       (org_.z_+dest_.z_)/2. );
                           
    double elength = Length();                  
                      
    Normalize();
    // computing the cross product
    Point3D rot( dest_.y_ * ax.z_ - ax.y_ * dest_.z_,
               -(dest_.x_ * ax.z_ - ax.x_ * dest_.z_),
                 dest_.x_ * ax.y_ - ax.x_ * dest_.y_ );

    if ( normalize ) {
         dest_ = rot;
         return *this;
      }

    // the rotated edge is centered on the old midpoint
    rot /= rot.Length();
    rot *= elength / 2.;
    org_  = midpoint + rot;
    dest_ = midpoint - rot;
    
    return *this;
 }





/// see "Mathematische Formeln" p. 170, Schnittwinkel...
double Edge3D::AngleTo( const Edge3D& v ) const 
 {
    Point3D a = dest_   - org_;
    Point3D b = v.dest_ - v.org_;
    
    // a b
    // ---
    double ab = a.x_*b.x_ + a.y_*b.y_ + a.z_*b.z_;
    
    // |a| . |b|
    // ---------
    double a_dot_b = std::sqrt( (a.x_*a.x_+a.y_*a.y_+a.z_*a.z_) *
                                   (b.x_*b.x_+b.y_*b.y_+b.z_*b.z_) );
                                   
    double cos_angle = ab/a_dot_b;

    // if zero intercept
    if ( cos_angle == 0. ) return 90.;
    // if outside of range of 'acos' function
    if ( cos_angle >  1. ) return   0.;
    if ( cos_angle < -1. ) return 180.;
        
    return static_cast<double>(180./3.1415926535897932384) * std::acos(cos_angle);
    
 } // end AngleTo




void Edge3D::Out(std::ostream& os) const
 {
    os <<"\nEdge3D::Out:"<< endl;
    os <<"Origin: ";
    org_.Out(os);
    os <<"\nDestination";
    dest_.Out(os);
 }

} // end mjl

