// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_POINT_H
#define CSMP_POINT_H

#include "CSMP_definitions.h"

namespace csmp {

/**
    @brief Generic class that is used in CSMP to stored Node coordinate values
    and to manipulate them efficiently.

    @author Stephan K. Matthai and Adriana Paluszny
    @date 2007
    @note refactored by SKM 20/12/2016
*/
template<uint32_t dim>
class Point {
  public:
    /// initialises point to default position of zero
    explicit Point( double = 0. ) noexcept;
  
    /// construct point from an STL vector of coordinate values
    explicit Point( const std::vector<double>& );
    explicit Point( const std::array<double,dim>& ) noexcept;
    
    // note: move semantics constructed automatically by compiler (see static_assert in .cpp)
    
    Point& operator=( double );
    
    Point operator+( const Point& ) const noexcept;
    Point operator-( const Point& ) const noexcept;
    Point operator*( const Point& ) const noexcept;
    Point operator/( const Point& ) const noexcept;
    
    Point operator+( double ) const noexcept;
    Point operator-( double ) const noexcept;
    Point operator*( double ) const noexcept;
    Point operator/( double ) const noexcept;
    
    Point& operator+=( const Point& ) noexcept;
    Point& operator-=( const Point& ) noexcept;
    Point& operator*=( const Point& ) noexcept;
    Point& operator/=( const Point& ) noexcept;
    
    Point& operator+=( double ) noexcept;
    Point& operator-=( double ) noexcept;
    Point& operator*=( double ) noexcept;
    Point& operator/=( double ) noexcept;
  
    /// accessor and mutator of point (0=x coordinate, 1=y...)
    double& operator[]( uint32_t );
  
    /// accessor of point (0=x coordinate, 1=y...)
    const double& operator[]( uint32_t ) const;
  
    /// compares points using epsilon from numeric_limits
    bool operator==( const Point& ) const noexcept;
    bool operator!=( const Point& ) const noexcept;
  
    /// comparison of position in space; ascertains spatial proximity ordering in associative containers
    bool operator<( const Point& ) const noexcept;

    /// comparison of position in space; ascertains spatial proximity ordering in associative containers
    bool operator>( const Point& ) const noexcept;
  
    /// change the coordinates of an existing point to those stored in the supplied STL vector
    void Set( const std::vector<double>& ) noexcept;
    void Set( const std::array<double,dim>& ) noexcept;
  
    /// returns the offset of th point from the origin of the coordinate system
    double Length() const noexcept;
  
    /// returns the square of the distance of the point from the origin of the coordinate system
    double SquaredLength() const noexcept;
  
    /// enforce offset of point from coordinate origin (when point is used to store a vector)
    void NormalizeLengthTo( double len=1. ) noexcept;
  
    /// return distance between current and other point
    double DistanceTo( const Point& ) const noexcept;
  
    /// checks whether points coincide within the giving tolerance
    bool CoincidesWithWithinTolerance( const Point&, double tolerance=1.0e-5 ) const noexcept;
  
    /// checks whether point lies on a straight line between the supplied to points
    bool IsBetween( const Point& pt1, const Point& pt2 ) noexcept;
  
    /// returns point coordinates into an STL vector
    std::vector<double> Coordinates() const;
    std::array<double,dim> CoordinateArray() const noexcept;
  
    /// prints point cooordinates to screen
    void  Out() const;
};


/// subtracts coordinates of point (3nd arg) from double (1st arg)
template<uint32_t dim>
Point<dim>  operator-( double, const Point<dim>& ) noexcept;

/// adds point coordinates
template<uint32_t dim>
Point<dim>  operator+( double, const Point<dim>& ) noexcept;

/// multiplies the coordinates of the 2 points
template<uint32_t dim>
Point<dim>  operator*( double, const Point<dim>& ) noexcept;

/// writes the point coordinates to an output stream
template<uint32_t dim>
std::ostream&  operator<<( std::ostream&, const Point<dim>& );

/// returns the midpoint of the 2 points
template<uint32_t dim>
Point<dim>  midPoint( const Point<dim>&, const Point<dim>& ) noexcept;

/// treating the points as vectors originating in the origin, computes their scalar product
template<uint32_t dim>
double  dotProduct( const Point<dim>&, const Point<dim>& ) noexcept;

/// treating the points as vectors originating in the origin, computes their cross product vector
template<uint32_t dim>
Point<dim>  crossProduct( const Point<dim>&, const Point<dim>& ) noexcept;

/// returns the angle in degrees between the line segments that start with the first point and terminate at the last point of thedge, ignoring edge direction
template<uint32_t dim>
double angleBetweenEdges( const std::pair<Point<dim>,Point<dim> >& edge1, const std::pair<Point<dim>,Point<dim> >& edge2 ) noexcept;

/// between two points
template<uint32_t dim>
double distance( const Point<dim>&, const Point<dim>& ) noexcept;

// specialisations

double dotProduct( const Point<1U>& p1, const Point<1U>& p2 ) noexcept;
double dotProduct( const Point<2U>& p1, const Point<2U>& p2 ) noexcept;
double dotProduct( const Point<3U>& p1, const Point<3U>& p2 ) noexcept;

Point<1U> crossProduct( const Point<1U>& p1, const Point<1U>& p2 ) noexcept;
Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 ) noexcept;
Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 ) noexcept;

template<>
class Point<1U> {
  public:
    Point( double = 0. ) noexcept; ///< explicit keyword is not required because conversion is desired
    explicit Point( const std::array<double,1>& ) noexcept;
    explicit Point( const std::vector<double>& );
    
    Point( const Point<1U>& )            noexcept = default;
    Point( Point<1U>&& )                 noexcept = default;
    Point<1U>& operator=( const Point<1U>& ) noexcept = default;
    Point<1U>& operator=( Point<1U>&& )      noexcept = default;
    ~Point()                                 noexcept = default;
    
    Point<1U>& operator=( double ) noexcept;
    
    Point<1U> operator+( const Point<1U>& ) const noexcept;
    Point<1U> operator-( const Point<1U>& ) const noexcept;
    Point<1U> operator*( const Point<1U>& ) const noexcept;
    Point<1U> operator/( const Point<1U>& ) const noexcept;
    Point<1U> operator+( double ) const noexcept;
    Point<1U> operator-( double ) const noexcept;
    Point<1U> operator*( double ) const noexcept;
    Point<1U> operator/( double ) const noexcept;
    Point<1U>& operator+=( const Point<1U>& ) noexcept;
    Point<1U>& operator-=( const Point<1U>& ) noexcept;
    Point<1U>& operator*=( const Point<1U>& ) noexcept;
    Point<1U>& operator/=( const Point<1U>& ) noexcept;
    Point<1U>& operator+=( double ) noexcept;
    Point<1U>& operator-=( double ) noexcept;
    Point<1U>& operator*=( double ) noexcept;
    Point<1U>& operator/=( double ) noexcept;
    double&  operator[]( uint32_t);
    const double&  operator[]( uint32_t) const;
    bool   operator==( const Point<1U>& ) const noexcept;
    bool   operator!=( const Point<1U>& ) const noexcept;
    bool   operator<( const Point<1U>& ) const noexcept;
    bool   operator>( const Point<1U>& ) const noexcept;
    void   Set( const std::array<double,1U>& ) noexcept;
    void   Set( const std::vector<double>& );
    double Length() const noexcept;
    double SquaredLength() const noexcept;
    void   NormalizeLengthTo( double len=1. ) noexcept;
    double DistanceTo( const Point<1U>& ) const noexcept;
    bool   CoincidesWithWithinTolerance( const Point<1U>&, double tolerance=1.0e-5 ) const noexcept;
    bool   IsBetween( const Point& pt1, const Point& pt2 ) noexcept;
    std::vector<double> Coordinates() const noexcept;
    std::array<double,1U> CoordinateArray() const  noexcept{ return std::array<double,1U>{ x_ }; }

    void   Out() const;

    friend Point<1U> operator-( double, const Point<1U>& ) noexcept;
    friend Point<1U> operator+( double, const Point<1U>& ) noexcept;
    friend Point<1U> operator*( double, const Point<1U>& ) noexcept;

  protected:
    double x_;
};

double distance( const Point<1U>&, const Point<1U>& ) noexcept;



template<>
class Point<2U> {
  public:
    explicit Point( double = 0. ) noexcept;
    Point( double, double ) noexcept;
    explicit Point( const std::array<double,2U>& ) noexcept;
    explicit Point( const std::vector<double>& );

    ///  caller can move if they wish
    explicit Point<2U>( std::vector<double>&& v )
        : x_( v[0] ), y_( v[1] ) {
        assert( v.size() == 2U );
        // v is discarded after this; no copy of the vector's buffer occurs
    }

    Point( const Point<2U>& )            noexcept = default;
    Point( Point<2U>&& )                 noexcept = default;
    Point<2U>& operator=( const Point<2U>& ) noexcept = default;
    Point<2U>& operator=( Point<2U>&& )      noexcept = default;
    ~Point()                                 noexcept = default;
    
    Point<2U>& operator=( double ) noexcept;
    
    Point<2U> operator+( const Point<2U>& ) const noexcept;
    Point<2U> operator-( const Point<2U>& ) const noexcept;
    Point<2U> operator*( const Point<2U>& ) const noexcept;
    Point<2U> operator/( const Point<2U>& ) const noexcept;
    Point<2U> operator+( double ) const noexcept;
    Point<2U> operator-( double ) const noexcept;
    Point<2U> operator*( double ) const noexcept;
    Point<2U> operator/( double ) const noexcept;
    Point<2U>& operator+=( const Point<2U>& ) noexcept;
    Point<2U>& operator-=( const Point<2U>& ) noexcept;
    Point<2U>& operator*=( const Point<2U>& ) noexcept;
    Point<2U>& operator/=( const Point<2U>& ) noexcept;
    Point<2U>& operator+=( double ) noexcept;
    Point<2U>& operator-=( double ) noexcept;
    Point<2U>& operator*=( double ) noexcept;
    Point<2U>& operator/=( double ) noexcept;
    double&  operator[]( uint32_t);
    const double&  operator[]( uint32_t) const;
    bool   operator==( const Point<2U>& ) const noexcept;
    bool   operator!=( const Point<2U>& ) const noexcept;
    bool   operator<( const Point<2U>& ) const noexcept;
    bool   operator>( const Point<2U>& ) const noexcept;
    void   Set( const std::array<double,2U>& ) noexcept;
    void   Set( const std::vector<double>& );
    void   Set( double, double ) noexcept;
    double Length() const noexcept;
    double SquaredLength() const noexcept;
    void   NormalizeLengthTo( double len=1. ) noexcept;
    double DistanceTo( const Point& ) const noexcept;
    bool   CoincidesWithWithinTolerance( const Point&, double tolerance=1.0e-5 ) const noexcept;
    bool   IsBetween( const Point& pt1, const Point& pt2 ) noexcept;
    std::vector<double> Coordinates() const noexcept;
    std::array<double,2U> CoordinateArray() const  noexcept{ return std::array<double,2U>{ x_, y_ }; }
    void   Out() const;
  
    friend Point<2U> operator-( double, const Point<2U>& ) noexcept;
    friend Point<2U> operator+( double, const Point<2U>& ) noexcept;
    friend Point<2U> operator*( double, const Point<2U>& ) noexcept;

  protected:
    double x_, y_;
};

Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 ) noexcept;

double    distance( const Point<2U>&, const Point<2U>& ) noexcept;

/// distance between point, A, and line segment BC, does not check case where projection of A is outside of segment BC
double    distanceFromLine( const Point<2U>& A, const Point<2U>& B, const Point<2U>& C ) noexcept;



template<>
class Point<3U> {
  public:
    explicit Point( double = 0. );
    Point( double, double, double ) noexcept;
    explicit Point( const std::array<double,3U>& ) noexcept;
    explicit Point( const std::vector<double>& );

    ///  caller can move if they wish
    explicit Point<3U>( std::vector<double>&& v )
        : x_( v[0] ), y_( v[1] ), z_( v[2] ) {
        assert( v.size() == 3U );
        // v is discarded after this; no copy of the vector's buffer occurs
    }

    Point( const Point<3U>& )            noexcept = default;
    Point( Point<3U>&& )                 noexcept = default;
    Point<3U>& operator=( const Point<3U>& ) noexcept = default;
    Point<3U>& operator=( Point<3U>&& )      noexcept = default;
    ~Point()                                 noexcept = default;
    
    Point<3U>& operator=( double ) noexcept;
    
    Point<3U> operator+( const Point<3U>& ) const noexcept;
    Point<3U> operator-( const Point<3U>& ) const noexcept;
    Point<3U> operator*( const Point<3U>& ) const noexcept;
    Point<3U> operator/( const Point<3U>& ) const noexcept;
    Point<3U> operator+( double ) const noexcept;
    Point<3U> operator-( double ) const noexcept;
    Point<3U> operator*( double ) const noexcept;
    Point<3U> operator/( double ) const noexcept;
    Point<3U>& operator+=( const Point<3U>& ) noexcept;
    Point<3U>& operator-=( const Point<3U>& ) noexcept;
    Point<3U>& operator*=( const Point<3U>& ) noexcept;
    Point<3U>& operator/=( const Point<3U>& ) noexcept;
    Point<3U>& operator+=( double ) noexcept;
    Point<3U>& operator-=( double ) noexcept;
    Point<3U>& operator*=( double ) noexcept;
    Point<3U>& operator/=( double ) noexcept;
    double&  operator[]( uint32_t);
    const double&  operator[]( uint32_t) const;
    bool   operator==( const Point<3U>& ) const noexcept;
    bool   operator!=( const Point<3U>& ) const noexcept;
    bool   operator<( const Point<3U>& ) const noexcept;
    bool   operator>( const Point<3U>& ) const noexcept;
    void   Set( const std::array<double,3U>& ) noexcept;
    void   Set( const std::vector<double>& );
    void   Set( double, double, double ) noexcept;
    double Length() const noexcept;
    double SquaredLength() const noexcept;
    void   NormalizeLengthTo( double len=1. ) noexcept;
    double DistanceTo( const Point& ) const noexcept;
    bool   CoincidesWithWithinTolerance( const Point&, double tolerance=1.0e-5 ) const noexcept;
    bool   IsBetween( const Point& pt1, const Point& pt2 ) noexcept;
    std::vector<double> Coordinates() const noexcept;
    std::array<double,3U> CoordinateArray() const  noexcept { return std::array<double,3U>{ x_, y_, z_ }; }
    void   Out() const;
 
    friend Point<3U> operator-( double, const Point<3U>& ) noexcept;
    friend Point<3U> operator+( double, const Point<3U>& ) noexcept;
    friend Point<3U> operator*( double, const Point<3U>& ) noexcept;

  protected:
    double x_, y_, z_;
};

/// The length of the exterior product.
template<uint32_t dim>
double exteriorProductLength( const Point<dim>& p1, const Point<dim>& p2 ) noexcept;

Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 ) noexcept;

/// distance between 2 points in 3D space
double    distance( const Point<3U>&, const Point<3U>& ) noexcept;

/// distance between point, A, and line segment BC (2D and 3D only), if projection of A on BC falls outside of this segment, ditance to nearest end-point is returned
template<uint32_t dim>
double    distanceFromLine( const Point<dim>& A, const Point<dim>& B, const Point<dim>& C ) noexcept;

std::ostream&   operator<<( std::ostream& stream, const Point<1U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<2U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<3U>& pt );


// INLINED METHODS

// 2D

inline Point<2U>::Point(double val) noexcept
 : x_(val), y_(val)
 {
 }

inline Point<2U>::Point( double px, double py ) noexcept : x_(px), y_(py)
 {
 }

inline Point<2U>::Point( const std::vector<double>& v )
 : x_(v[0]), y_(v[1])
 {
 }

inline double& Point<2U>::operator[]( uint32_t i )
 {
    return (i == 0U) ? x_ : y_;
 }

inline const double& Point<2U>::operator[]( uint32_t i ) const
 {
    return (i == 0U) ? x_ : y_;
 }

inline void Point<2U>::Set( const std::array<double,2U>& v ) noexcept
 {
    x_ = v[0];
    y_ = v[1];
 }

inline void Point<2U>::Set( const std::vector<double>& v )
 {
    assert( v.size() == 2U );
    x_ = v[0];
    y_ = v[1];
 }

inline void Point<2U>::Set( double px, double py ) noexcept
 {
    x_ = px;
    y_ = py;
 }

inline Point<2U>& Point<2U>::operator=( double val ) noexcept
 {
    x_ = val;
    y_ = val;
    return *this;
 }
    

inline Point<2U>  Point<2U>::operator+( const Point<2U>& pt ) const noexcept
 {
    return (Point<2U>(x_ + pt.x_,y_ + pt.y_));
 }

inline Point<2U>  Point<2U>::operator-( const Point<2U>& pt ) const noexcept
 {
    return (Point<2U>(x_ - pt.x_,y_ - pt.y_));
 }

inline Point<2U> Point<2U>::operator*( const Point<2U>& pt ) const noexcept
 {
    return (Point<2U>(x_ * pt.x_,y_ * pt.y_));
 }

inline Point<2U>  Point<2U>::operator/( const Point<2U>& pt ) const noexcept
 {
    return (Point<2U>(x_ / pt.x_,y_ / pt.y_));
 }
    

inline Point<2U>  Point<2U>::operator+( double val ) const noexcept
 {
    return (Point<2U>(x_ + val, y_ + val));
 }

inline Point<2U>  Point<2U>::operator-( double val ) const noexcept
 {
    return (Point<2U>(x_ - val, y_ - val));
 }

inline Point<2U>  Point<2U>::operator*( double val ) const noexcept
 {
    return (Point<2U>(x_ * val, y_ * val));
 }

inline Point<2U>  Point<2U>::operator/( double val ) const noexcept
 {
    return (Point<2U>(x_ / val, y_ / val));
 }


inline Point<2U>& Point<2U>::operator+=( const Point<2U>& pt ) noexcept
 {
    x_ += pt.x_;
    y_ += pt.y_;
    return *this;
 }

inline Point<2U>& Point<2U>::operator-=( const Point<2U>& pt ) noexcept
 {
    x_ -= pt.x_;
    y_ -= pt.y_;
    return *this;
 }

inline Point<2U>& Point<2U>::operator*=( const Point<2U>& pt ) noexcept
 {
    x_ *= pt.x_;
    y_ *= pt.y_;
    return *this;
 }

inline Point<2U>& Point<2U>::operator/=( const Point<2U>& pt ) noexcept
 {
    x_ /= pt.x_;
    y_ /= pt.y_;
    return *this;
 }
 
inline Point<2U>& Point<2U>::operator+=( double val ) noexcept
 {
    x_ += val;
    y_ += val;
    return *this;
 }

inline Point<2U>& Point<2U>::operator-=( double val ) noexcept
 {
    x_ -= val;
    y_ -= val;
    return *this;
 }

inline Point<2U>& Point<2U>::operator*=( double val ) noexcept
 {
    x_ *= val;
    y_ *= val;
    return *this;
 }

inline Point<2U>& Point<2U>::operator/=( double val ) noexcept
 {
    x_ /= val;
    y_ /= val;
    return *this;
 }

// stand-alone operator functions

inline Point<2U> operator-( double val, const Point<2U>& pt ) noexcept
 {
    return (Point<2U>( val - pt.x_, val - pt.y_ ));
 }   

inline Point<2U> operator+( double val, const Point<2U>& pt ) noexcept
 {
    return (Point<2U>( pt.x_ + val, pt.y_ + val ));
 }   

inline Point<2U> operator*( double val, const Point<2U>& pt ) noexcept
 {
    return (Point<2U>( pt.x_ * val, pt.y_ * val ));
 }   

inline bool Point<2U>::operator<( const Point<2U>& pt )  const noexcept
 {
    return (( x_ < pt.x_ ) || (( x_ == pt.x_ ) && ( y_ < pt.y_ )));
 }

inline bool Point<2U>::operator>( const Point<2U>& pt )  const noexcept
 {
    return (( x_ > pt.x_ ) || (( x_ == pt.x_ ) && ( y_ > pt.y_ )));
 }
 
// other

inline double  Point<2U>::SquaredLength() const noexcept
 {
    return x_*x_ + y_*y_;
 }

inline std::vector<double> Point<2U>::Coordinates() const noexcept
 {
    return std::vector<double>{x_,y_};
 }
 
inline double distance( const Point<2U>& a, const Point<2U>& b ) noexcept {
     return hypot( a[0]-b[0], a[1]-b[1] );
  }
  
inline double dotProduct( const Point<2U>& p1, const Point<2U>& p2 ) noexcept
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U];
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
inline Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 ) noexcept
  {
     return (Point<2U>( p1[0U] * p2[1U] - p2[0U] * p1[1U], std::numeric_limits<double>::quiet_NaN() ));
  }

template<uint32_t dim>
inline Point<dim>  midPoint( const Point<dim>& p1, const Point<dim>& p2 ) noexcept
 {
    return (Point<dim>( p1 + p2 ) / 2.);
 }


// 3D

// constructors

inline Point<3U>::Point( double val )
 : x_(val), y_(val), z_(val)
 {
 }

inline Point<3U>::Point( double px, double py, double pz ) noexcept : x_(px), y_(py), z_(pz)
 {
 }

inline Point<3U>::Point( const std::array<double,3U>& v ) noexcept
 : x_(v[0]), y_(v[1]), z_(v[2])
 {
 }

inline Point<3U>::Point( const std::vector<double>& v )
 : x_(v[0]), y_(v[1]), z_(v[2])
 {
 }

// accessors / mutators

inline double& Point<3U>::operator[]( uint32_t i )
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }

inline const double& Point<3U>::operator[]( uint32_t i ) const
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }

inline void Point<3U>::Set( const std::array<double,3U>& v ) noexcept
 {
    x_ = v[0];
    y_ = v[1];
    z_ = v[2];
 }

inline void Point<3U>::Set( const std::vector<double>& v )
 {
    assert( v.size() == 3U );
    x_ = v[0];
    y_ = v[1];
    z_ = v[2];
 }

inline void Point<3U>::Set( double px, double py, double pz ) noexcept
 {
    x_ = px;
    y_ = py;
    z_ = pz;
 }

// assignment

inline Point<3U>& Point<3U>::operator=( double val ) noexcept
 {
    x_ = val;
    y_ = val;
    z_ = val;
    return *this;
 }

// operators

inline Point<3U>  Point<3U>::operator+( const Point<3U>& pt ) const noexcept
 {
    return (Point<3U>(x_ + pt.x_, y_ + pt.y_, z_ + pt.z_));
 }

inline Point<3U>  Point<3U>::operator-( const Point<3U>& pt ) const noexcept
 {
    return (Point<3U>(x_ - pt.x_,y_ - pt.y_, z_ - pt.z_));
 }

inline Point<3U>  Point<3U>::operator*( const Point<3U>& pt ) const noexcept
 {
    return (Point<3U>(x_ * pt.x_,y_ * pt.y_, z_ * pt.z_));
 }

inline Point<3U>  Point<3U>::operator/( const Point<3U>& pt ) const noexcept
 {
    return (Point<3U>(x_ / pt.x_,y_ / pt.y_, z_ / pt.z_));
 }
    

inline Point<3U>  Point<3U>::operator+( double val ) const noexcept
 {
    return (Point<3U>(x_ + val, y_ + val, z_ + val));
 }

inline Point<3U>  Point<3U>::operator-( double val ) const noexcept
 {
    return (Point<3U>(x_ - val, y_ - val, z_ - val));
 }

inline Point<3U>  Point<3U>::operator*( double val ) const noexcept
 {
    return (Point<3U>(x_ * val, y_ * val, z_ * val));
 }

inline Point<3U>  Point<3U>::operator/( double val ) const noexcept
 {
    return (Point<3U>(x_ / val, y_ / val, z_ / val));
 }


inline Point<3U>& Point<3U>::operator+=( const Point<3U>& pt ) noexcept
 {
    x_ += pt.x_;
    y_ += pt.y_;
    z_ += pt.z_;
    return *this;
 }

inline Point<3U>& Point<3U>::operator-=( const Point<3U>& pt ) noexcept
 {
    x_ -= pt.x_;
    y_ -= pt.y_;
    z_ -= pt.z_;
    return *this;
 }

inline Point<3U>& Point<3U>::operator*=( const Point<3U>& pt ) noexcept
 {
    x_ *= pt.x_;
    y_ *= pt.y_;
    z_ *= pt.z_;
    return *this;
 }

inline Point<3U>& Point<3U>::operator/=( const Point<3U>& pt ) noexcept
 {
    x_ /= pt.x_;
    y_ /= pt.y_;
    z_ /= pt.z_;
    return *this;
 }
 

inline Point<3U>& Point<3U>::operator+=( double val ) noexcept
 {
    x_ += val;
    y_ += val;
    z_ += val;
    return *this;
 }

inline Point<3U>& Point<3U>::operator-=( double val ) noexcept
 {
    x_ -= val;
    y_ -= val;
    z_ -= val;
    return *this;
 }

inline Point<3U>& Point<3U>::operator*=( double val ) noexcept
 {
    x_ *= val;
    y_ *= val;
    z_ *= val;
    return *this;
 }

inline Point<3U>& Point<3U>::operator/=( double val ) noexcept
 {
    x_ /= val;
    y_ /= val;
    z_ /= val;
    return *this;
 }

inline Point<3U> operator-( double val, const Point<3U>& pt ) noexcept
 {
    return (Point<3U>( val - pt.x_, val - pt.y_, val - pt.z_ ));
 }   

inline Point<3U> operator+( double val, const Point<3U>& pt ) noexcept
 {
    return (Point<3U>( pt.x_ + val, pt.y_ + val, pt.z_ + val ));
 }   

inline Point<3U> operator*( double val, const Point<3U>& pt ) noexcept
 {
    return (Point<3U>( pt.x_ * val, pt.y_ * val, pt.z_ * val ));
 }   

inline double dotProduct( const Point<3U>& p1, const Point<3U>& p2 ) noexcept
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U] + p1[2U] * p2[2U];
  }

 // tested: SKM O.K.
inline Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 ) noexcept
  {
     return (Point<3U>( p1[1U] * p2[2U] - p1[2U] * p2[1U],
                                 p1[2U] * p2[0U] - p1[0U] * p2[2U],
                                 p1[0U] * p2[1U] - p1[1U] * p2[0U] ));
  }


inline std::vector<double> Point<3U>::Coordinates() const noexcept
 {
    return std::vector<double>{x_,y_,z_};
 }


inline double Point<3U>::SquaredLength() const noexcept
 {
    return x_ * x_ + y_ * y_ + z_ * z_;
 }

inline bool Point<3U>::IsBetween( const Point<3U>& pt1, const Point<3U>& pt2 ) noexcept
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (z_ > pt1.z_) &&
         (x_ < pt2.x_) && (y_ < pt2.y_) && (z_ < pt2.z_) ) return true;
    return false;
}


} // end namespace csmp

#endif

