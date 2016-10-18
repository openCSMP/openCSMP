#ifndef CSMP_POINT_H
#define CSMP_POINT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <cstdlib>
#include <limits>
#include "CSMP_number_types.h"

namespace csmp {

/**
    @brief Generic class that is used in CSMP to stored Node coordinate values
    and to manipulate them efficiently.

    @author Stephan K. Matthai and Adriana Paluszny
    @date 2007
*/
template<size_t dim>
class Point {
  public:
    explicit Point( double64 = 0. ); // guarantees initialization to zero
    ~Point();
    /// construct point from an STL vector of coordinate values
    explicit Point( const std::vector<double64>& );
    Point( const Point& );
    Point( Point&& );
    Point&    operator=( const Point& );    
    Point&    operator=( Point&& );
    Point&    operator=( double64 );
    Point     operator+( const Point& ) const;    
    Point     operator-( const Point& ) const;    
    Point     operator*( const Point& ) const;    
    Point     operator/( const Point& ) const;    
    Point     operator+( double64 ) const;    
    Point     operator-( double64 ) const;    
    Point     operator*( double64 ) const;    
    Point     operator/( double64 ) const;    
    Point&    operator+=( const Point& );    
    Point&    operator-=( const Point& );    
    Point&    operator*=( const Point& );    
    Point&    operator/=( const Point& );    
    Point&    operator+=( double64 );    
    Point&    operator-=( double64 );    
    Point&    operator*=( double64 );    
    Point&    operator/=( double64 );
    /// accessor and mutator of point (0=x coordinate, 1=y...)
    double64& operator[](size_t);
    /// accessor of point (0=x coordinate, 1=y...)
    double64  operator[](size_t) const;
    /// compares points using epsilon from numeric_limits
    bool      operator==( const Point& ) const;
    bool      operator!=( const Point& ) const;
    bool      operator<( const Point& ) const;
    bool      operator>( const Point& ) const;
    /// change the coordinates of an existing point to those stored in the supplied STL vector
    void      Set( const std::vector<double64>& );
    /// returns the offset of th point from the origin of the coordinate system
    double64  Length() const;
    /// enforce offset of point from coordinate origin (when point is used to store a vector)
    void      NormalizeLengthTo( double64 len=1. );
    /// return distance between current and other point
    double64  DistanceTo( const Point& ) const;
    /// checks whether points coincide within the giving tolerance
    bool      CoincidesWithWithinTolerance( const Point&, double64 tolerance=1.0e-5 ) const;
    /// checks whether point lies on a straight line between the supplied to points
    bool      IsBetween( const Point& pt1, const Point& pt2 );
    /// returns point coordinates into an STL vector
    std::vector<double64> Coordinates() const;
    /// prints point cooordinates to screen
    void  Out() const;

  protected:
    double64 xyz_[dim];
};


/// subtracts coordinates of point (3nd arg) from double (1st arg)
template<size_t dim>
Point<dim>  operator-( double64, const Point<dim>& );

/// adds point coordinates
template<size_t dim>
Point<dim>  operator+( double64, const Point<dim>& );

/// multiplies the coordinates of the 2 points
template<size_t dim>
Point<dim>  operator*( double64, const Point<dim>& );

/// writes the point coordinates to an output stream
template<size_t dim>
std::ostream&  operator<<( std::ostream&, const Point<dim>& );

/// returns the midpoint of the 2 points
template<size_t dim>
Point<dim>  midPoint( const Point<dim>&, const Point<dim>& );

/// treating the points as vectors originating in the origin, computes their scalar product
template<size_t dim>
double64  dotProduct( const Point<dim>&, const Point<dim>& );

/// treating the points as vectors originating in the origin, computes their cross product vector
template<size_t dim>
Point<dim>  crossProduct( const Point<dim>&, const Point<dim>& );




// partial specializations
template<>
class Point<1U> {
  public:
    Point( double64 = 0. ); // explicit is not required because conversion is desired
    ~Point();
    explicit Point( const std::vector<double64>& );
    Point( const Point<1U>& );
    Point( Point<1U>&& );
    Point<1U>& operator=( const Point<1U>& );
    Point<1U>& operator=( Point<1U>&& );
    Point<1U>& operator=( double64 );
    Point<1U>  operator+( const Point<1U>& ) const;    
    Point<1U>  operator-( const Point<1U>& ) const;    
    Point<1U>  operator*( const Point<1U>& ) const;    
    Point<1U>  operator/( const Point<1U>& ) const;    
    Point<1U>  operator+( double64 ) const;    
    Point<1U>  operator-( double64 ) const;    
    Point<1U>  operator*( double64 ) const;    
    Point<1U>  operator/( double64 ) const;    
    Point<1U>& operator+=( const Point<1U>& );    
    Point<1U>& operator-=( const Point<1U>& );    
    Point<1U>& operator*=( const Point<1U>& );    
    Point<1U>& operator/=( const Point<1U>& );    
    Point<1U>& operator+=( double64 );    
    Point<1U>& operator-=( double64 );    
    Point<1U>& operator*=( double64 );    
    Point<1U>& operator/=( double64 );    
    double64&    operator[](size_t);
    double64     operator[](size_t) const;
    bool   operator==( const Point<1U>& ) const;
    bool   operator!=( const Point<1U>& ) const;
    bool   operator<( const Point<1U>& ) const;
    bool   operator>( const Point<1U>& ) const;
    void   Set( const std::vector<double64>& );
    double64  Length() const;
    void      NormalizeLengthTo( double64 len=1. );
    double64  DistanceTo( const Point<1U>& ) const;
    bool      CoincidesWithWithinTolerance( const Point<1U>&, double64 tolerance=1.0e-5 ) const;
    bool      IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double64> Coordinates() const;
    void                  Out() const;

    friend Point<1U> operator-( double64, const Point<1U>& );
    friend Point<1U> operator+( double64, const Point<1U>& );
    friend Point<1U> operator*( double64, const Point<1U>& );

  protected:
    double64 x_;
};


template<>
class Point<2U> {
  public:
    explicit Point( double64 = 0. );
    ~Point();
    Point( double64, double64 );
    explicit Point( const std::vector<double64>& );
    Point( const Point& );
    Point( Point&& );
    Point& operator=( const Point& );
    Point& operator=( Point&& );
    Point& operator=( double64 );
    Point  operator+( const Point& ) const;    
    Point  operator-( const Point& ) const;    
    Point  operator*( const Point& ) const;    
    Point  operator/( const Point& ) const;    
    Point  operator+( double64 ) const;    
    Point  operator-( double64 ) const;    
    Point  operator*( double64 ) const;    
    Point  operator/( double64 ) const;    
    Point& operator+=( const Point& );    
    Point& operator-=( const Point& );    
    Point& operator*=( const Point& );    
    Point& operator/=( const Point& );    
    Point& operator+=( double64 );    
    Point& operator-=( double64 );    
    Point& operator*=( double64 );    
    Point& operator/=( double64 );    
    double64&    operator[](size_t);
    double64     operator[](size_t) const;
    bool   operator==( const Point& ) const;
    bool   operator!=( const Point& ) const;
    bool   operator<( const Point& ) const;
    bool   operator>( const Point& ) const;
    void   Set( const std::vector<double64>& );
    void   Set( double64, double64 );
    double64     Length() const;
    void     NormalizeLengthTo( double64 len=1. );
    double64 DistanceTo( const Point& ) const;
    bool     CoincidesWithWithinTolerance( const Point&, double64 tolerance=1.0e-5 ) const;
    bool     IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double64> Coordinates() const;
    void            Out() const;
  
    friend Point<2U> operator-( double64, const Point<2U>& );
    friend Point<2U> operator+( double64, const Point<2U>& );
    friend Point<2U> operator*( double64, const Point<2U>& );

  protected:
    double64 x_, y_;
};



template<>
class Point<3U> {
  public:
    explicit Point( double64 = 0. );
    ~Point();
    Point( double64, double64, double64 );
    explicit Point( const std::vector<double64>& );
    Point( const Point& );
    Point( Point&& );
    Point& operator=( const Point& );
    Point& operator=( Point&& );
    Point& operator=( double64 );
    Point  operator+( const Point& ) const;    
    Point  operator-( const Point& ) const;    
    Point  operator*( const Point& ) const;    
    Point  operator/( const Point& ) const;    
    Point  operator+( double64 ) const;    
    Point  operator-( double64 ) const;    
    Point  operator*( double64 ) const;    
    Point  operator/( double64 ) const;    
    Point& operator+=( const Point& );    
    Point& operator-=( const Point& );    
    Point& operator*=( const Point& );    
    Point& operator/=( const Point& );    
    Point& operator+=( double64 );    
    Point& operator-=( double64 );    
    Point& operator*=( double64 );    
    Point& operator/=( double64 );    
    double64&    operator[](size_t);
    double64     operator[](size_t) const;
    bool   operator==( const Point& ) const;
    bool   operator!=( const Point& ) const;
    bool   operator<( const Point& ) const;
    bool   operator>( const Point& ) const;
    void   Set( const std::vector<double64>& );
    void   Set( double64, double64, double64 );
    double64 Length() const;
    void     NormalizeLengthTo( double64 len=1. );
    double64 DistanceTo( const Point& ) const;
    bool     CoincidesWithWithinTolerance( const Point&, double64 tolerance=1.0e-5 ) const;
    bool     IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double64> Coordinates() const;
    void                  Out() const;
 
    friend Point<3U> operator-( double64, const Point<3U>& );
    friend Point<3U> operator+( double64, const Point<3U>& );
    friend Point<3U> operator*( double64, const Point<3U>& );

  protected:
    double64 x_, y_, z_;
};

// ---------------------------------------------------------------------------

// MEMBER DEFINITIONS

// ---------------------------------------------------------------------------
template<size_t dim>
inline Point<dim>::~Point()
 {
 }


template<size_t dim>
inline double64& Point<dim>::operator[](size_t i)
 {
    return xyz_[i];
 }
 
 
template<size_t dim>
inline double64 Point<dim>::operator[](size_t i) const
 {
    return xyz_[i];
 }




// ----------------------------------------------------------------------------
//
//  1D SPECIALIZATION
//
// ----------------------------------------------------------------------------

inline Point<1U>::Point( double64 val ) : x_(val)
 {
 }

inline Point<1U>::Point( const Point<1U>& pt ) : x_(pt.x_)
 {
 }

inline Point<1U>::Point( Point<1U>&& pt ) : x_{pt.x_}
 {
 }

inline Point<1U>::Point( const std::vector<double64>& v )
 : x_(v[0])
 {
 }

inline Point<1U>::~Point()
 {
 }


inline double64& Point<1U>::operator[](size_t)
 {
    return x_;
 }

inline double64 Point<1U>::operator[](size_t) const
 {
    return x_;
 }


inline void Point<1U>::Set( const std::vector<double64>& v )
 {
    x_ = v[0];
 }


inline Point<1U>& Point<1U>::operator=( const Point<1U>& pt )
 {
    if ( &pt != this ) x_ = pt.x_;
    return *this;
 }

inline Point<1U>& Point<1U>::operator=( Point<1U>&& pt )
 {
    if ( &pt != this ) x_ = {pt.x_};
    return *this;
 }

inline Point<1U>& Point<1U>::operator=( double64 val )
 {
    x_ = val;
    return *this;
 }
    

inline Point<1U>  Point<1U>::operator+( const Point<1U>& pt ) const
 {
    return std::move(Point<1U>(x_ + pt.x_));
 }

inline Point<1U>  Point<1U>::operator-( const Point<1U>& pt ) const
 {
    return std::move(Point<1U>(x_ - pt.x_));
 }

inline Point<1U>  Point<1U>::operator*( const Point<1U>& pt ) const
 {
    return std::move(Point<1U>(x_ * pt.x_));
 }

inline Point<1U>  Point<1U>::operator/( const Point<1U>& pt ) const
 {
    return std::move(Point<1U>(x_ / pt.x_));
 }
    

inline Point<1U>  Point<1U>::operator+( double64 val ) const
 {
    return std::move(Point<1U>(x_ + val));
 }

inline Point<1U>  Point<1U>::operator-( double64 val ) const
 {
    return std::move(Point<1U>(x_ - val));
 }

inline Point<1U>  Point<1U>::operator*( double64 val ) const
 {
    return std::move(Point<1U>(x_ * val));
 }

inline Point<1U>  Point<1U>::operator/( double64 val ) const
 {
    return std::move(Point<1U>(x_ / val));
 }



inline Point<1U>& Point<1U>::operator+=( const Point<1U>& pt )
 {
    x_ += pt.x_;
    return *this;
 }

inline Point<1U>& Point<1U>::operator-=( const Point<1U>& pt )
 {
    x_ -= pt.x_;
    return *this;
 }

inline Point<1U>& Point<1U>::operator*=( const Point<1U>& pt )
 {
    x_ *= pt.x_;
    return *this;
 }

inline Point<1U>& Point<1U>::operator/=( const Point<1U>& pt )
 {
    x_ /= pt.x_;
    return *this;
 }
 

inline Point<1U>& Point<1U>::operator+=( double64 val )
 {
    x_ += val;
    return *this;
 }

inline Point<1U>& Point<1U>::operator-=( double64 val )
 {
    x_ -= val;
    return *this;
 }

inline Point<1U>& Point<1U>::operator*=( double64 val )
 {
    x_ *= val;
    return *this;
 }

inline Point<1U>& Point<1U>::operator/=( double64 val )
 {
    x_ /= val;
    return *this;
 }
 

inline bool Point<1U>::operator==( const Point<1U>& pt ) const
 {
    if ( x_ != pt.x_ ) return false;
    return true;
 }

inline bool Point<1U>::operator!=( const Point<1U>& pt ) const
 {
    if ( x_ != pt.x_ ) return true;
    return false;
 }

inline bool Point<1U>::operator<( const Point<1U>& pt ) const
 {
   if ( x_ < pt.x_ ) return true;
   return false;
 }

inline bool Point<1U>::operator>( const Point<1U>& pt ) const
 {
   if ( x_ > pt.x_ ) return true;
   return false;
 }


inline double64  Point<1U>::Length() const
 {
    // return std::sqrt( x_ * x_ );
    return x_; 
 }


inline void Point<1U>::NormalizeLengthTo( double64 len ) 
 {
    x_ = len; 
 }


inline double64  Point<1U>::DistanceTo( const Point<1U>& pt ) const
 {
    return std::fabs(pt.x_ - x_); 
 }


inline bool Point<1U>::CoincidesWithWithinTolerance( const Point<1U>& pt, 
                                                     double64 tolerance ) const
 {
    if ( std::fabs(pt.x_ - x_) > tolerance ) return false;
    return true;
 }

inline bool Point<1U>::IsBetween( const Point<1U>& pt1, const Point<1U>& pt2)
 {
    if ( (x_ > pt1.x_)  && (x_ < pt2.x_) ) return true;
    return false;
 }


inline std::vector<double64> Point<1U>::Coordinates() const
 {
    return std::vector<double64>(1U,x_);
 }


inline void Point<1U>::Out() const
 {
    std::cout <<"\nPoint<" << 1U;
    std::cout <<">::Out(): coordinates: "<< x_ << std::endl;
    std::cout.flush();
 }




// ----------------------------------------------------------------------------
//
//  2D SPECIALIZATION
//
// ----------------------------------------------------------------------------

inline Point<2U>::Point(double64 val)
 : x_(val), y_(val)
 {
 }

inline Point<2U>::Point( double64 px, double64 py ) : x_(px), y_(py)
 {
 }

inline Point<2U>::Point( const Point<2U>& pt ) : x_(pt.x_), y_(pt.y_)
 {
 }

inline Point<2U>::Point( Point<2U>&& pt ) : x_{pt.x_}, y_{pt.y_}
 {
 }

inline Point<2U>::Point( const std::vector<double64>& v )
 : x_(v[0]), y_(v[1])
 {
 }

inline Point<2U>::~Point()
 {
 }


inline double64& Point<2U>::operator[]( size_t i )
 {
    return (i == 0U) ? x_ : y_;
 }

inline double64 Point<2U>::operator[]( size_t i ) const
 {
    return (i == 0U) ? x_ : y_;
 }


inline void Point<2U>::Set( const std::vector<double64>& v )
 {
    x_ = v[0];
    y_ = v[1];
 }


inline void Point<2U>::Set( double64 px, double64 py )
 {
    x_ = px;
    y_ = py;
 }


inline Point<2U>& Point<2U>::operator=( const Point<2U>& pt )
 {
    if ( &pt != this ) {
         x_ = pt.x_;
         y_ = pt.y_;
      }
    return *this;
 }

inline Point<2U>& Point<2U>::operator=( Point<2U>&& pt )
 {
    if ( &pt != this ) {
         x_ = {pt.x_};
         y_ = {pt.y_};
      }
    return *this;
 }

inline Point<2U>& Point<2U>::operator=( double64 val )
 {
    x_ = val;
    y_ = val;
    return *this;
 }
    

inline Point<2U>  Point<2U>::operator+( const Point<2U>& pt ) const
 {
    return std::move(Point<2U>(x_ + pt.x_,y_ + pt.y_));
 }

inline Point<2U>  Point<2U>::operator-( const Point<2U>& pt ) const
 {
    return std::move(Point<2U>(x_ - pt.x_,y_ - pt.y_));
 }

inline Point<2U>  Point<2U>::operator*( const Point<2U>& pt ) const
 {
    return std::move(Point<2U>(x_ * pt.x_,y_ * pt.y_));
 }

inline Point<2U>  Point<2U>::operator/( const Point<2U>& pt ) const
 {
    return std::move(Point<2U>(x_ / pt.x_,y_ / pt.y_));
 }
    

inline Point<2U>  Point<2U>::operator+( double64 val ) const
 {
    return std::move(Point<2U>(x_ + val, y_ + val));
 }

inline Point<2U>  Point<2U>::operator-( double64 val ) const
 {
    return std::move(Point<2U>(x_ - val, y_ - val));
 }

inline Point<2U>  Point<2U>::operator*( double64 val ) const
 {
    return std::move(Point<2U>(x_ * val, y_ * val));
 }

inline Point<2U>  Point<2U>::operator/( double64 val ) const
 {
    return std::move(Point<2U>(x_ / val, y_ / val));
 }



inline Point<2U>& Point<2U>::operator+=( const Point<2U>& pt )
 {
    x_ += pt.x_;
    y_ += pt.y_;
    return *this;
 }

inline Point<2U>& Point<2U>::operator-=( const Point<2U>& pt )
 {
    x_ -= pt.x_;
    y_ -= pt.y_;
    return *this;
 }

inline Point<2U>& Point<2U>::operator*=( const Point<2U>& pt )
 {
    x_ *= pt.x_;
    y_ *= pt.y_;
    return *this;
 }

inline Point<2U>& Point<2U>::operator/=( const Point<2U>& pt )
 {
    x_ /= pt.x_;
    y_ /= pt.y_;
    return *this;
 }
 

inline Point<2U>& Point<2U>::operator+=( double64 val )
 {
    x_ += val;
    y_ += val;
    return *this;
 }

inline Point<2U>& Point<2U>::operator-=( double64 val )
 {
    x_ -= val;
    y_ -= val;
    return *this;
 }

inline Point<2U>& Point<2U>::operator*=( double64 val )
 {
    x_ *= val;
    y_ *= val;
    return *this;
 }

inline Point<2U>& Point<2U>::operator/=( double64 val )
 {
    x_ /= val;
    y_ /= val;
    return *this;
 }
 

inline bool Point<2U>::operator==( const Point<2U>& pt ) const
 {
    return (x_ == pt.x_) && (y_ == pt.y_);
 }

inline bool Point<2U>::operator!=( const Point<2U>& pt ) const
 {
    return !(*this == pt);
 }

inline bool Point<2U>::operator<( const Point<2U>& pt )  const
 {
    return (( x_ < pt.x_ ) || (( x_ == pt.x_ ) && ( y_ < pt.y_ )));
 }

inline bool Point<2U>::operator>( const Point<2U>& pt )  const
 {
    return (( x_ > pt.x_ ) || (( x_ == pt.x_ ) && ( y_ > pt.y_ )));
 }


inline double64  Point<2U>::Length() const
 {
    return std::sqrt( x_ * x_ + y_ * y_ );
 }


inline void Point<2U>::NormalizeLengthTo( double64 len )
 {
    *this /= Length();
    *this *= len; 
 }


inline double64  Point<2U>::DistanceTo( const Point<2U>& pt ) const
 {
    return Point<2U>(pt - *this).Length(); 
 }


inline bool Point<2U>::CoincidesWithWithinTolerance( const Point<2U>& pt, 
                                                     double64 tolerance ) const
 {
    if ( Point<2U>(pt - *this).Length() > tolerance ) return false;
    return true;
 }

inline bool Point<2U>::IsBetween( const Point<2U>& pt1, const Point<2U>& pt2 )
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (x_ < pt2.x_) && (y_ < pt2.y_) ) return true;
    return false;
}


inline std::vector<double64> Point<2U>::Coordinates() const
 {
    std::vector<double64> temp(2U,x_);
    temp[1U] = y_;
    return temp;
 }



inline void Point<2U>::Out() const
 {
    std::cout <<"\nPoint<"<< 2U;
    std::cout <<">::Out(): coordinates: "<< x_ <<","<< y_ << std::endl;
    std::cout.flush();
 }






// ----------------------------------------------------------------------------
//
//  3D SPECIALIZATION
//
// ----------------------------------------------------------------------------

inline Point<3U>::Point( double64 val )
 : x_(val), y_(val), z_(val)
 {
 }

inline Point<3U>::Point( double64 px, double64 py, double64 pz ) : x_(px), y_(py), z_(pz)
 {
 }

inline Point<3U>::Point( const Point<3U>& pt ) : x_(pt.x_), y_(pt.y_), z_(pt.z_)
 {
 }

inline Point<3U>::Point( Point<3U>&& pt ) : x_{pt.x_}, y_{pt.y_}, z_{pt.z_}
 {
//std::cerr <<"\nPoint<3>: called move constructor.";
 }

inline Point<3U>::Point( const std::vector<double64>& v )
 : x_(v[0]), y_(v[1]), z_(v[2])
 {
 }

inline Point<3U>::~Point()
 {
 }


inline double64& Point<3U>::operator[]( size_t i )
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }

inline double64 Point<3U>::operator[]( size_t i ) const
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }


inline void Point<3U>::Set( const std::vector<double64>& v )
 {
    x_ = v[0];
    y_ = v[1];
    z_ = v[2];
 }


inline void Point<3U>::Set( double64 px, double64 py, double64 pz )
 {
    x_ = px;
    y_ = py;
    z_ = pz;
 }


inline Point<3U>& Point<3U>::operator=( const Point<3U>& pt )
 {
    if ( &pt != this ) {
         x_ = pt.x_;
         y_ = pt.y_;
         z_ = pt.z_;
      }
    return *this;
 }

inline Point<3U>& Point<3U>::operator=( Point<3U>&& pt )
 {
    if ( &pt != this ) {
         x_ = {pt.x_};
         y_ = {pt.y_};
         z_ = {pt.z_};
      }
//std::cerr <<"\nPoint<3>::operator= called move assignment.";
    return *this;
 }

inline Point<3U>& Point<3U>::operator=( double64 val )
 {
    x_ = val;
    y_ = val;
    z_ = val;
    return *this;
 }
    

inline Point<3U>  Point<3U>::operator+( const Point<3U>& pt ) const
 {
    return std::move(Point<3U>(x_ + pt.x_, y_ + pt.y_, z_ + pt.z_));
 }

inline Point<3U>  Point<3U>::operator-( const Point<3U>& pt ) const
 {
    return std::move(Point<3U>(x_ - pt.x_,y_ - pt.y_, z_ - pt.z_));
 }

inline Point<3U>  Point<3U>::operator*( const Point<3U>& pt ) const
 {
    return std::move(Point<3U>(x_ * pt.x_,y_ * pt.y_, z_ * pt.z_));
 }

inline Point<3U>  Point<3U>::operator/( const Point<3U>& pt ) const
 {
    return std::move(Point<3U>(x_ / pt.x_,y_ / pt.y_, z_ / pt.z_));
 }
    

inline Point<3U>  Point<3U>::operator+( double64 val ) const
 {
    return std::move(Point<3U>(x_ + val, y_ + val, z_ + val));
 }

inline Point<3U>  Point<3U>::operator-( double64 val ) const
 {
    return std::move(Point<3U>(x_ - val, y_ - val, z_ - val));
 }

inline Point<3U>  Point<3U>::operator*( double64 val ) const
 {
    return std::move(Point<3U>(x_ * val, y_ * val, z_ * val));
 }

inline Point<3U>  Point<3U>::operator/( double64 val ) const
 {
    return std::move(Point<3U>(x_ / val, y_ / val, z_ / val));
 }


inline Point<3U>& Point<3U>::operator+=( const Point<3U>& pt )
 {
    x_ += pt.x_;
    y_ += pt.y_;
    z_ += pt.z_;
    return *this;
 }

inline Point<3U>& Point<3U>::operator-=( const Point<3U>& pt )
 {
    x_ -= pt.x_;
    y_ -= pt.y_;
    z_ -= pt.z_;
    return *this;
 }

inline Point<3U>& Point<3U>::operator*=( const Point<3U>& pt )
 {
    x_ *= pt.x_;
    y_ *= pt.y_;
    z_ *= pt.z_;
    return *this;
 }

inline Point<3U>& Point<3U>::operator/=( const Point<3U>& pt )
 {
    x_ /= pt.x_;
    y_ /= pt.y_;
    z_ /= pt.z_;
    return *this;
 }
 

inline Point<3U>& Point<3U>::operator+=( double64 val )
 {
    x_ += val;
    y_ += val;
    z_ += val;
    return *this;
 }

inline Point<3U>& Point<3U>::operator-=( double64 val )
 {
    x_ -= val;
    y_ -= val;
    z_ -= val;
    return *this;
 }

inline Point<3U>& Point<3U>::operator*=( double64 val )
 {
    x_ *= val;
    y_ *= val;
    z_ *= val;
    return *this;
 }

inline Point<3U>& Point<3U>::operator/=( double64 val )
 {
    x_ /= val;
    y_ /= val;
    z_ /= val;
    return *this;
 }
 

inline bool Point<3U>::operator==( const Point<3U>& pt ) const
 {
    return (x_ == pt.x_) && (y_ == pt.y_) && (z_ == pt.z_);
 }

inline bool Point<3U>::operator!=( const Point<3U>& pt ) const
 {
    return !(*this == pt);
 }

inline bool Point<3U>::operator<( const Point<3U>& p )  const
 {
    return (x_<p.x_ || (x_==p.x_ && y_<p.y_) || (x_==p.x_ &&  y_==p.y_ && z_<p.z_) );
 }

inline bool Point<3U>::operator>( const Point<3U>& p )  const
 {
    return (x_>p.x_ || (x_==p.x_ && y_>p.y_) || (x_==p.x_ && y_==p.y_ && z_>p.z_) );
 }


inline double64  Point<3U>::Length() const
 {
    return std::sqrt( x_ * x_ + y_ * y_ + z_ * z_ );
 }


inline void Point<3U>::NormalizeLengthTo( double64 len )
 {
    *this /= Length();
    *this *= len; 
 }


inline double64  Point<3U>::DistanceTo( const Point<3U>& pt ) const
 {
    return Point<3U>(pt - *this).Length(); 
 }


inline bool Point<3U>::CoincidesWithWithinTolerance( const Point<3U>& pt, 
                                                     double64 tolerance ) const
 {
    if ( Point<3U>(pt - *this).Length() > tolerance ) return false;
    return true;
 }


inline bool Point<3U>::IsBetween( const Point<3U>& pt1, const Point<3U>& pt2 )
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (z_ > pt1.z_) && (x_ < pt2.x_) && (y_ < pt2.y_) && (z_ < pt2.z_) ) return true;
    return false;
}


inline std::vector<double64> Point<3U>::Coordinates() const
 {
    std::vector<double64> temp(3U,x_);
    temp[1U] = y_;
    temp[2U] = z_;
    return temp;
 }



inline void Point<3U>::Out() const
 {
    std::cout <<"\nPoint<"<< 3U;
    std::cout <<">::Out(): coordinates: ";
    std::cout << x_ <<","<< y_ <<","<< z_ << std::endl;
    std::cout.flush();
 }





// -------------------------------------------------------------------------------
//
//    GENERIC INLINE FUNCTIONS
//
// -------------------------------------------------------------------------------
template<size_t dim>
inline double64  Point<dim>::DistanceTo( const Point& pt ) const
 {
    return Point<dim>( pt - *this ).Length();
 }

// general case
template<size_t dim>
inline Point<dim> operator-( double64 val, const Point<dim>& pt )
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) temp[i] = val - pt[i];  
    return std::move(temp);
 } 

template<size_t dim>
inline Point<dim> operator+( double64 val, const Point<dim>& pt )
 {
    return std::move(pt + val);
 } 
  
template<size_t dim>
inline Point<dim> operator*( double64 val, const Point<dim>& pt )
 {
    return std::move(pt * val);
 }   

// 1D

inline Point<1U> operator-( double64 val, const Point<1U>& pt )
 {
    return Point<1U>( val - pt.x_ );
 }   

inline Point<1U> operator+( double64 val, const Point<1U>& pt )
 {
    return Point<1U>( pt.x_ + val );
 }   

inline Point<1U> operator*( double64 val, const Point<1U>& pt )
 {
    return Point<1U>( pt.x_ * val );
 }   

// 2D

inline Point<2U> operator-( double64 val, const Point<2U>& pt )
 {
    return Point<2U>( val - pt.x_, val - pt.y_ );
 }   

inline Point<2U> operator+( double64 val, const Point<2U>& pt )
 {
    return Point<2U>( pt.x_ + val, pt.y_ + val );
 }   

inline Point<2U> operator*( double64 val, const Point<2U>& pt )
 {
    return Point<2U>( pt.x_ * val, pt.y_ * val );
 }   
 
// 3D

inline Point<3U> operator-( double64 val, const Point<3U>& pt )
 {
    return std::move(Point<3U>( val - pt.x_, val - pt.y_, val - pt.z_ ));
 }   

inline Point<3U> operator+( double64 val, const Point<3U>& pt )
 {
    return std::move(Point<3U>( pt.x_ + val, pt.y_ + val, pt.z_ + val ));
 }   

inline Point<3U> operator*( double64 val, const Point<3U>& pt )
 {
    return std::move(Point<3U>( pt.x_ * val, pt.y_ * val, pt.z_ * val ));
 }   




// generic copy constructor
template<size_t dim>
inline Point<dim>::Point( const Point<dim>& pt )
 {
    *this = pt;
 }

template<size_t dim>
inline Point<dim>  midPoint( const Point<dim>& p1, const Point<dim>& p2 )
 { 
    return std::move(Point<dim>( p1 + p2 ) / 2.);
 }

// dot = scalar product
template<size_t dim>
inline double64 dotProduct( const Point<dim>& p1, const Point<dim>& p2 )
  {
    return inner_product( p1.xyz_[0], p1.xyz_[0]+dim, p2.xyz_[0], p2.xyz_[0]+dim, 0. );
  }

inline double64 dotProduct( const Point<1U>& p1, const Point<1U>& p2 )
  {
    return p1[0U] * p2[0U];
  }

inline double64 dotProduct( const Point<2U>& p1, const Point<2U>& p2 )
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U];
  }

inline double64 dotProduct( const Point<3U>& p1, const Point<3U>& p2 )
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U] + p1[2U] * p2[2U];
  }
 

inline Point<1U> crossProduct( const Point<1U>& p1, const Point<1U>& p2 )
  {
     return std::move(Point<1U>( p1[0U] * p2[0U] ));
  }



/** crossProduct(2D)

Watch out in 2D the cross product is not uniquely defined:

CrossProductAnalog1(U,V)=(U.x*V.y-U.y*V.x)  is a scalar!
CrossProductAnalog2(U)=(U.y, -U.x)
(yes, second analog takes only one argument, and return orthogonal vector :)
First analog makes some physical and geometrical sense, 
second analog comes from "determinant rule", for determinant of 2x2 matrix,

|A B|
|C D| = AD-BC 

To avoid ambiguities, NAN is returned in the second component of Point.

*/
inline Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 )
  {
     return std::move(Point<2U>( p1[0U] * p2[1U] - p2[0U] * p1[1U], std::numeric_limits<double64>::quiet_NaN() ));
  }

 // tested: SKM O.K.
inline Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 )
  {
     return std::move(Point<3U>( p1[1U] * p2[2U] - p1[2U] * p2[1U],
                                 p1[2U] * p2[0U] - p1[0U] * p2[2U],
                                 p1[0U] * p2[1U] - p1[1U] * p2[0U] ));
  }

std::ostream&   operator<<( std::ostream& stream, const Point<1U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<2U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<3U>& pt );

} // end namespace csmp

#endif

