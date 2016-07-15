#ifndef MJL_POINT_H
#define MJL_POINT_H

#include "MJL_geometry.h"

namespace mjl {

class Edge;

// Lazlo p. 73 ff.
class Point {
  public:
    Point();
    Point( double, double );
    Point( const Point& );
    ~Point() {};
    Point& operator=( const Point& );    
    Point  operator+( const Point& ) const;    
    Point  operator-( const Point& ) const;    
    Point& operator+=( const Point& );    
    Point& operator-=( const Point& );    
    Point& operator*=( const Point& );    
    Point& operator/=( const Point& );    
    Point& operator*=( double );    
    Point& operator/=( double );    
    friend Point operator*( double, const Point& );
    double& operator()(int);
    double  operator[](int) const;
    double  X() const;
    double  Y() const;
    bool       operator==( const Point& ) const;
    bool       operator!=( const Point& ) const;
    bool       operator<( const Point& )  const;
    bool       operator>( const Point& )  const;
    int        Orientation( Point&, Point&, Point& );
    double  PolarAngle() const;
    double  Length() const;
    LOCATION   Classify( const Point&, const Point& ) const; 
    LOCATION   Classify( const Edge& ) const; 
    void       Set( double, double );
    void       Move( double dx, double dy );
    void       Scale( double xfac, double yfac );
    void       Out() const;

  private:
    double x_, y_;
};



int  polarCmp( const Point&, const Point&, const Point& ); // MJL modified, p. 111





inline  Point::Point()
  : x_(0.), y_(0.) {}


inline  Point::Point( double px, double py )
  : x_(px), y_(py) {}
  

inline  Point&  Point::operator=( const Point& p )
 {
    if ( &p == this ) return *this;
    x_ = p.x_;
    y_ = p.y_;
    return *this;
 }    

inline  Point&  Point::operator+=( const Point& p )
 {
    x_ += p.x_;
    y_ += p.y_;
    return *this;
 }    

inline  Point&  Point::operator-=( const Point& p )
 {
    x_ -= p.x_;
    y_ -= p.y_;
    return *this;
 }    

inline  Point&  Point::operator/=( double p )
 {
    x_ /= p;
    y_ /= p;
    return *this;
 }    

inline  Point&  Point::operator*=( double p )
 {
    x_ *= p;
    y_ *= p;
    return *this;
 }    

inline  Point&  Point::operator*=( const Point& p )
 {
    x_ *= p.x_;
    y_ *= p.y_;
    return *this;
 }    

inline  Point&  Point::operator/=( const Point& p )
 {
    x_ /= p.x_;
    y_ /= p.y_;
    return *this;
 }    

inline  Point::Point( const Point& p )  
 { 
    *this = p; 
 }


inline  Point  Point::operator+( const Point& p ) const
 {
    return Point( x_ + p.x_, y_ + p.y_ );
 }   


inline  Point  Point::operator-( const Point& p ) const
 {
    return Point( x_ - p.x_, y_ - p.y_ );
 }   


inline Point  operator*( double s, const Point& p )
 {
    return Point( s * p.x_, s * p.y_ );
 }



inline  double&  Point::operator()(int i) 
 {
    return ( i == 0 ) ? x_ : y_;
 }


inline  double Point::operator[](int i) const 
 {
    return ( i == 0 ) ? x_ : y_;
 }


inline  bool  Point::operator==( const Point& p ) const
 {
    return ( x_ == p.x_ ) && ( y_ == p.y_ );
 }


inline  bool  Point::operator!=( const Point& p ) const
 {
    return !( *this == p );
 }


inline  bool  Point::operator<( const Point& p ) const
 {
    return (( x_ < p.x_ ) || (( x_ == p.x_ ) && ( y_ < p.y_ )));
 }


inline  bool  Point::operator>( const Point& p ) const
 {
    return (( x_ > p.x_ ) || (( x_ == p.x_ ) && ( y_ > p.y_ )));
 }


inline double Point::X() const { return x_; }
inline double Point::Y() const { return y_; }

inline  double  Point::Length() const
  {
      return std::sqrt( x_*x_ + y_*y_ );
  }


inline  void  Point::Set( double cx, double cy )
  {
      x_ = cx;
      y_ = cy;
  }


  
inline  void  Point::Move( double dx, double dy )
 {
     x_ += dx;
     y_ += dy;
 }

inline  void  Point::Scale( double xfac, double yfac )
 {
     x_ *= xfac;
     y_ *= yfac;
 }
  
  
}

#endif















