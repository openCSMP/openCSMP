#ifndef CSMP_POINT_H
#define CSMP_POINT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <cstdlib>
#include <cassert>
#include <limits>
#include "CSMP_number_types.h"

namespace csmp {

/**
    @brief Generic class that is used in CSMP to stored Node coordinate values
    and to manipulate them efficiently.

    @author Stephan K. Matthai and Adriana Paluszny
    @date 2007
    @note refactored by SKM 20/12/2016
*/
template<size_t dim>
class Point {
  public:
    /// initialises point to default position of zero
    explicit Point( double64 = 0. );
    ~Point();
  
    /// construct point from an STL vector of coordinate values
    explicit Point( const std::vector<double64>& );
    Point( const Point& );
    Point( Point&& ) = default;
    Point& operator=( const Point& );
    Point& operator=( Point&& ) = default;
    Point& operator=( double64 );
    Point operator+( const Point& ) const;
    Point operator-( const Point& ) const;
    Point operator*( const Point& ) const;
    Point operator/( const Point& ) const;
    Point operator+( double64 ) const;
    Point operator-( double64 ) const;
    Point operator*( double64 ) const;
    Point operator/( double64 ) const;
    Point& operator+=( const Point& );
    Point& operator-=( const Point& );
    Point& operator*=( const Point& );
    Point& operator/=( const Point& );
    Point& operator+=( double64 );
    Point& operator-=( double64 );
    Point& operator*=( double64 );
    Point& operator/=( double64 );
  
    /// accessor and mutator of point (0=x coordinate, 1=y...)
    double64& operator[](size_t);
  
    /// accessor of point (0=x coordinate, 1=y...)
    const double64& operator[](size_t) const;
  
    /// compares points using epsilon from numeric_limits
    bool      operator==( const Point& ) const;
    bool      operator!=( const Point& ) const;
  
    /// comparison of position in space; ascertains spatial proximity ordering in associative containers
    bool      operator<( const Point& ) const;

    /// comparison of position in space; ascertains spatial proximity ordering in associative containers
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


// specialisations

double64 dotProduct( const Point<1U>& p1, const Point<1U>& p2 );
double64 dotProduct( const Point<2U>& p1, const Point<2U>& p2 );
double64 dotProduct( const Point<3U>& p1, const Point<3U>& p2 );

Point<1U> crossProduct( const Point<1U>& p1, const Point<1U>& p2 );
Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 );
Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 );

template<>
class Point<1U> {
  public:
    Point( double64 = 0. ); ///< explicit keyword is not required because conversion is desired
    ~Point();
    explicit Point( const std::vector<double64>& );
    Point( const Point<1U>& );
    Point( Point<1U>&& ) = default;
    Point<1U>& operator=( const Point<1U>& );
    Point<1U>& operator=( Point<1U>&& ) = default;
    Point<1U>& operator=( double64 );
    Point<1U> operator+( const Point<1U>& ) const;
    Point<1U> operator-( const Point<1U>& ) const;
    Point<1U> operator*( const Point<1U>& ) const;
    Point<1U> operator/( const Point<1U>& ) const;
    Point<1U> operator+( double64 ) const;
    Point<1U> operator-( double64 ) const;
    Point<1U> operator*( double64 ) const;
    Point<1U> operator/( double64 ) const;
    Point<1U>& operator+=( const Point<1U>& );
    Point<1U>& operator-=( const Point<1U>& );
    Point<1U>& operator*=( const Point<1U>& );
    Point<1U>& operator/=( const Point<1U>& );
    Point<1U>& operator+=( double64 );
    Point<1U>& operator-=( double64 );
    Point<1U>& operator*=( double64 );
    Point<1U>& operator/=( double64 );
    double64&  operator[](size_t);
    const double64&  operator[](size_t) const;
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
    Point( Point&& ) = default;
    Point<2U>& operator=( const Point<2U>& );
    Point<2U>& operator=( Point<2U>&& ) = default;
    Point<2U>& operator=( double64 );
    Point<2U> operator+( const Point<2U>& ) const;
    Point<2U> operator-( const Point<2U>& ) const;
    Point<2U> operator*( const Point<2U>& ) const;
    Point<2U> operator/( const Point<2U>& ) const;
    Point<2U> operator+( double64 ) const;
    Point<2U> operator-( double64 ) const;
    Point<2U> operator*( double64 ) const;
    Point<2U> operator/( double64 ) const;
    Point<2U>& operator+=( const Point<2U>& );
    Point<2U>& operator-=( const Point<2U>& );
    Point<2U>& operator*=( const Point<2U>& );
    Point<2U>& operator/=( const Point<2U>& );
    Point<2U>& operator+=( double64 );
    Point<2U>& operator-=( double64 );
    Point<2U>& operator*=( double64 );
    Point<2U>& operator/=( double64 );
    double64&  operator[](size_t);
    const double64&  operator[](size_t) const;
    bool   operator==( const Point<2U>& ) const;
    bool   operator!=( const Point<2U>& ) const;
    bool   operator<( const Point<2U>& ) const;
    bool   operator>( const Point<2U>& ) const;
    void   Set( const std::vector<double64>& );
    void   Set( double64, double64 );
    double64 Length() const;
    void     NormalizeLengthTo( double64 len=1. );
    double64 DistanceTo( const Point& ) const;
    bool     CoincidesWithWithinTolerance( const Point&, double64 tolerance=1.0e-5 ) const;
    bool     IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double64> Coordinates() const;
    void     Out() const;
  
    friend Point<2U> operator-( double64, const Point<2U>& );
    friend Point<2U> operator+( double64, const Point<2U>& );
    friend Point<2U> operator*( double64, const Point<2U>& );

  protected:
    double64 x_, y_;
};

Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 );




template<>
class Point<3U> {
  public:
    explicit Point( double64 = 0. );
    ~Point();
    Point( double64, double64, double64 );
    explicit Point( const std::vector<double64>& );
    Point( const Point& );
    Point( Point&& ) = default;
    Point<3U>& operator=( const Point<3U>& );
    Point<3U>& operator=( Point<3U>&& ) = default;
    Point<3U>& operator=( double64 );
    Point<3U> operator+( const Point<3U>& ) const;
    Point<3U> operator-( const Point<3U>& ) const;
    Point<3U> operator*( const Point<3U>& ) const;
    Point<3U> operator/( const Point<3U>& ) const;
    Point<3U> operator+( double64 ) const;
    Point<3U> operator-( double64 ) const;
    Point<3U> operator*( double64 ) const;
    Point<3U> operator/( double64 ) const;
    Point<3U>& operator+=( const Point<3U>& );
    Point<3U>& operator-=( const Point<3U>& );
    Point<3U>& operator*=( const Point<3U>& );
    Point<3U>& operator/=( const Point<3U>& );
    Point<3U>& operator+=( double64 );
    Point<3U>& operator-=( double64 );
    Point<3U>& operator*=( double64 );
    Point<3U>& operator/=( double64 );
    double64&  operator[](size_t);
    const double64&  operator[](size_t) const;
    bool     operator==( const Point<3U>& ) const;
    bool     operator!=( const Point<3U>& ) const;
    bool     operator<( const Point<3U>& ) const;
    bool     operator>( const Point<3U>& ) const;
    void     Set( const std::vector<double64>& );
    void     Set( double64, double64, double64 );
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

Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 );

std::ostream&   operator<<( std::ostream& stream, const Point<1U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<2U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<3U>& pt );

} // end namespace csmp

#endif

