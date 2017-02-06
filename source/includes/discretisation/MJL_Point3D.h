#ifndef MJL_POINT_3D_H
#define MJL_POINT_3D_H

#include "MJL_geometry.h"

namespace mjl {

class Triangle3D;

class Point3D {
  public:
    Point3D( double, double, double );
    Point3D();
    Point3D( const Point3D& p ); 
    Point3D& operator=( const Point3D& p );
    bool operator<( const Point3D& p ) const;
    double   operator[](int i) const;
    double&  operator()(int i);

    double X() const;
    double Y() const;
    double Z() const;
    void   Set( double, double, double );
    
    Point3D  operator+( const Point3D& p ) const;
    Point3D  operator-( const Point3D& p ) const;
    Point3D  operator+=( const Point3D& p );
    Point3D  operator-=( const Point3D& p );
    Point3D  operator*=( const Point3D& p );
    Point3D  operator/=( const Point3D& p );
    Point3D  operator+=( double p );
    Point3D  operator-=( double p );
    Point3D  operator*=( double p );
    Point3D  operator/=( double p );
    Point3D         operator*( double ) const;
    friend Point3D  operator*( double, const Point3D& ); 
     
    bool     operator==( const Point3D& ) const;
    bool     operator!=( const Point3D& ) const;
     
    double  DotProduct( const Point3D& ) const;
    double  Length() const;
    bool    CoincidesWithWithinTolerance( const Point3D&, double tolerance=1.0e-5 ) const;

    void Out(std::ostream& os) const;

    double x_, y_, z_;
 };


Point3D  crossProduct( const Point3D&, const Point3D& );


// member functions

inline Point3D::Point3D( double x, double y, double z )
     : x_(x), y_(y), z_(z) {}


inline Point3D& Point3D::operator=( const Point3D& p )
 {
    if ( &p != this ) {
         x_=p.x_; y_=p.y_; z_=p.z_;
      }
    return *this;
 }
    
inline Point3D::Point3D( const Point3D& p ) 
 {
    *this = p;
 }
     
inline Point3D::Point3D() {}



inline Point3D  Point3D::operator+( const Point3D& p ) const
  { return Point3D( x_+p.x_, y_+p.y_, z_+p.z_ ); }
 
  

inline Point3D  Point3D::operator-( const Point3D& p ) const
     { return Point3D( x_-p.x_, y_-p.y_, z_-p.z_ ); }
 
     
inline Point3D  Point3D::operator+=( const Point3D& p )
     {
       x_ += p.x_;
       y_ += p.y_;
       z_ += p.z_;
       return *this;
     }
   
             
inline Point3D  Point3D::operator-=( const Point3D& p )
 {
   x_ -= p.x_;
   y_ -= p.y_;
   z_ -= p.z_;
   return *this;
 }
     
             
inline Point3D  Point3D::operator*=( const Point3D& p )
 {
   x_ *= p.x_;
   y_ *= p.y_;
   z_ *= p.z_;
   return *this;
 }
     
             
inline Point3D  Point3D::operator/=( const Point3D& p )
 {
   x_ /= p.x_;
   y_ /= p.y_;
   z_ /= p.z_;
   return *this;
 }
     
            
inline Point3D  Point3D::operator+=( double p )
 {
   x_ += p;
   y_ += p;
   z_ += p;
   return *this;
 }
     
          
inline Point3D  Point3D::operator-=( double p )
 {
   x_ -= p;
   y_ -= p;
   z_ -= p;
   return *this;
 }
     
          
inline Point3D  Point3D::operator*=( double p )
 {
   x_ *= p;
   y_ *= p;
   z_ *= p;
   return *this;
 }
         
          
inline Point3D  Point3D::operator/=( double p )
 {
   x_ /= p;
   y_ /= p;
   z_ /= p;
   return *this;
 }     
     

inline Point3D  Point3D::operator*( double s ) const
 { return Point3D(s*x_, s*y_, s*z_); }
 
     
inline bool Point3D::operator==( const Point3D& p ) const
 { return ((x_==p.x_) && (y_==p.y_) && (z_==p.z_)); } 
 
inline bool Point3D::operator!=( const Point3D& p ) const
 { return !(*this == p); } 


 
inline double Point3D::operator[](int i) const
 { return((i==0) ? x_ : ((i==1) ? y_ : z_)); }


 
inline double&  Point3D::operator()(int i)
 { return((i==0) ? x_ : ((i==1) ? y_ : z_)); }



inline double  Point3D::DotProduct( const Point3D& p ) const
 { return (x_*p.x_ + y_*p.y_ + z_*p.z_); }
 
 
 
inline double  Point3D::Length() const { return std::sqrt(x_*x_ + y_*y_ + z_*z_); }
  

inline void Point3D::Set( double ix, double iy, double iz )
  { x_=ix; y_=iy; z_=iz; }
     
     
inline bool Point3D::operator<( const Point3D& p ) const
 { 
    return (x_<p.x_ || (x_==p.x_ && y_<p.y_) || (x_==p.x_ && y_==p.y_ && z_<p.z_) ); 
 }
     

inline double Point3D::X() const { return x_; }


inline double Point3D::Y() const { return y_; }


inline double Point3D::Z() const { return z_; }
    



// non member functions

inline Point3D  crossProduct( const Point3D& a, const Point3D& b )
 {
    return Point3D( a.y_*b.z_ - a.z_*b.y_,
                    a.z_*b.x_ - a.x_*b.z_,
                    a.x_*b.y_ - a.y_*b.x_ );
 } 


inline Point3D  operator*( double s, const Point3D& p ) 
 { return Point3D(s*p.x_, s*p.y_, s*p.z_); }
    

   
} // end namespace mjl

#endif

























