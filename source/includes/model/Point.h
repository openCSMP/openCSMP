#ifndef CSMP_POINT_H
#define CSMP_POINT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <cstdlib>
#include <cassert>
#include <limits>

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
    explicit Point( double = 0. );
    ~Point();
  
    /// construct point from an STL vector of coordinate values
    explicit Point( const std::vector<double>& );
    Point( const Point& );
    Point( Point&& ) = default;
    Point& operator=( const Point& );
    Point& operator=( Point&& ) = default;
    Point& operator=( double );
    Point operator+( const Point& ) const;
    Point operator-( const Point& ) const;
    Point operator*( const Point& ) const;
    Point operator/( const Point& ) const;
    Point operator+( double ) const;
    Point operator-( double ) const;
    Point operator*( double ) const;
    Point operator/( double ) const;
    Point& operator+=( const Point& );
    Point& operator-=( const Point& );
    Point& operator*=( const Point& );
    Point& operator/=( const Point& );
    Point& operator+=( double );
    Point& operator-=( double );
    Point& operator*=( double );
    Point& operator/=( double );
  
    /// accessor and mutator of point (0=x coordinate, 1=y...)
    double& operator[]( uint32_t);
  
    /// accessor of point (0=x coordinate, 1=y...)
    const double& operator[]( uint32_t) const;
  
    /// compares points using epsilon from numeric_limits
    bool operator==( const Point& ) const;
    bool operator!=( const Point& ) const;
  
    /// comparison of position in space; ascertains spatial proximity ordering in associative containers
    bool operator<( const Point& ) const;

    /// comparison of position in space; ascertains spatial proximity ordering in associative containers
    bool operator>( const Point& ) const;
  
    /// change the coordinates of an existing point to those stored in the supplied STL vector
    void Set( const std::vector<double>& );
  
    /// returns the offset of th point from the origin of the coordinate system
    double Length() const;
  
    /// returns the square of the distance of the point from the origin of the coordinate system
    double SquaredLength() const;
  
    /// enforce offset of point from coordinate origin (when point is used to store a vector)
    void NormalizeLengthTo( double len=1. );
  
    /// return distance between current and other point
    double DistanceTo( const Point& ) const;
  
    /// checks whether points coincide within the giving tolerance
    bool CoincidesWithWithinTolerance( const Point&, double tolerance=1.0e-5 ) const;
  
    /// checks whether point lies on a straight line between the supplied to points
    bool IsBetween( const Point& pt1, const Point& pt2 );
  
    /// returns point coordinates into an STL vector
    std::vector<double> Coordinates() const;
  
    /// prints point cooordinates to screen
    void  Out() const;
};


/// subtracts coordinates of point (3nd arg) from double (1st arg)
template<uint32_t dim>
Point<dim>  operator-( double, const Point<dim>& );

/// adds point coordinates
template<uint32_t dim>
Point<dim>  operator+( double, const Point<dim>& );

/// multiplies the coordinates of the 2 points
template<uint32_t dim>
Point<dim>  operator*( double, const Point<dim>& );

/// writes the point coordinates to an output stream
template<uint32_t dim>
std::ostream&  operator<<( std::ostream&, const Point<dim>& );

/// returns the midpoint of the 2 points
template<uint32_t dim>
Point<dim>  midPoint( const Point<dim>&, const Point<dim>& );

/// treating the points as vectors originating in the origin, computes their scalar product
template<uint32_t dim>
double  dotProduct( const Point<dim>&, const Point<dim>& );

/// treating the points as vectors originating in the origin, computes their cross product vector
template<uint32_t dim>
Point<dim>  crossProduct( const Point<dim>&, const Point<dim>& );

/// returns the angle in degrees between the line segments that start with the first point and terminate at the last point of thedge, ignoring edge direction
template<uint32_t dim>
double angleBetweenEdges( const std::pair<Point<dim>,Point<dim> >& edge1, const std::pair<Point<dim>,Point<dim> >& edge2 );


// specialisations

double dotProduct( const Point<1U>& p1, const Point<1U>& p2 );
double dotProduct( const Point<2U>& p1, const Point<2U>& p2 );
double dotProduct( const Point<3U>& p1, const Point<3U>& p2 );

Point<1U> crossProduct( const Point<1U>& p1, const Point<1U>& p2 );
Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 );
Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 );

template<>
class Point<1U> {
  public:
    Point( double = 0. ); ///< explicit keyword is not required because conversion is desired
    ~Point();
    explicit Point( const std::vector<double>& );
    Point( const Point<1U>& );
    Point( Point<1U>&& ) = default;
    Point<1U>& operator=( const Point<1U>& );
    Point<1U>& operator=( Point<1U>&& ) = default;
    Point<1U>& operator=( double );
    Point<1U> operator+( const Point<1U>& ) const;
    Point<1U> operator-( const Point<1U>& ) const;
    Point<1U> operator*( const Point<1U>& ) const;
    Point<1U> operator/( const Point<1U>& ) const;
    Point<1U> operator+( double ) const;
    Point<1U> operator-( double ) const;
    Point<1U> operator*( double ) const;
    Point<1U> operator/( double ) const;
    Point<1U>& operator+=( const Point<1U>& );
    Point<1U>& operator-=( const Point<1U>& );
    Point<1U>& operator*=( const Point<1U>& );
    Point<1U>& operator/=( const Point<1U>& );
    Point<1U>& operator+=( double );
    Point<1U>& operator-=( double );
    Point<1U>& operator*=( double );
    Point<1U>& operator/=( double );
    double&  operator[]( uint32_t);
    const double&  operator[]( uint32_t) const;
    bool   operator==( const Point<1U>& ) const;
    bool   operator!=( const Point<1U>& ) const;
    bool   operator<( const Point<1U>& ) const;
    bool   operator>( const Point<1U>& ) const;
    void   Set( const std::vector<double>& );
    double Length() const;
    double SquaredLength() const;
    void   NormalizeLengthTo( double len=1. );
    double DistanceTo( const Point<1U>& ) const;
    bool   CoincidesWithWithinTolerance( const Point<1U>&, double tolerance=1.0e-5 ) const;
    bool   IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double> Coordinates() const;
    void   Out() const;

    friend Point<1U> operator-( double, const Point<1U>& );
    friend Point<1U> operator+( double, const Point<1U>& );
    friend Point<1U> operator*( double, const Point<1U>& );

  protected:
    double x_;
};



template<>
class Point<2U> {
  public:
    explicit Point( double = 0. );
    ~Point();
    Point( double, double );
    explicit Point( const std::vector<double>& );
    Point( const Point& );
    Point( Point&& ) = default;
    Point<2U>& operator=( const Point<2U>& );
    Point<2U>& operator=( Point<2U>&& ) = default;
    Point<2U>& operator=( double );
    Point<2U> operator+( const Point<2U>& ) const;
    Point<2U> operator-( const Point<2U>& ) const;
    Point<2U> operator*( const Point<2U>& ) const;
    Point<2U> operator/( const Point<2U>& ) const;
    Point<2U> operator+( double ) const;
    Point<2U> operator-( double ) const;
    Point<2U> operator*( double ) const;
    Point<2U> operator/( double ) const;
    Point<2U>& operator+=( const Point<2U>& );
    Point<2U>& operator-=( const Point<2U>& );
    Point<2U>& operator*=( const Point<2U>& );
    Point<2U>& operator/=( const Point<2U>& );
    Point<2U>& operator+=( double );
    Point<2U>& operator-=( double );
    Point<2U>& operator*=( double );
    Point<2U>& operator/=( double );
    double&  operator[]( uint32_t);
    const double&  operator[]( uint32_t) const;
    bool   operator==( const Point<2U>& ) const;
    bool   operator!=( const Point<2U>& ) const;
    bool   operator<( const Point<2U>& ) const;
    bool   operator>( const Point<2U>& ) const;
    void   Set( const std::vector<double>& );
    void   Set( double, double );
    double Length() const;
    double SquaredLength() const;
    void   NormalizeLengthTo( double len=1. );
    double DistanceTo( const Point& ) const;
    bool   CoincidesWithWithinTolerance( const Point&, double tolerance=1.0e-5 ) const;
    bool   IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double> Coordinates() const;
    void   Out() const;
  
    friend Point<2U> operator-( double, const Point<2U>& );
    friend Point<2U> operator+( double, const Point<2U>& );
    friend Point<2U> operator*( double, const Point<2U>& );

  protected:
    double x_, y_;
};

Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 );




template<>
class Point<3U> {
  public:
    explicit Point( double = 0. );
    ~Point();
    Point( double, double, double );
    explicit Point( const std::vector<double>& );
    Point( const Point& );
    Point( Point&& ) = default;
    Point<3U>& operator=( const Point<3U>& );
    Point<3U>& operator=( Point<3U>&& ) = default;
    Point<3U>& operator=( double );
    Point<3U> operator+( const Point<3U>& ) const;
    Point<3U> operator-( const Point<3U>& ) const;
    Point<3U> operator*( const Point<3U>& ) const;
    Point<3U> operator/( const Point<3U>& ) const;
    Point<3U> operator+( double ) const;
    Point<3U> operator-( double ) const;
    Point<3U> operator*( double ) const;
    Point<3U> operator/( double ) const;
    Point<3U>& operator+=( const Point<3U>& );
    Point<3U>& operator-=( const Point<3U>& );
    Point<3U>& operator*=( const Point<3U>& );
    Point<3U>& operator/=( const Point<3U>& );
    Point<3U>& operator+=( double );
    Point<3U>& operator-=( double );
    Point<3U>& operator*=( double );
    Point<3U>& operator/=( double );
    double&  operator[]( uint32_t);
    const double&  operator[]( uint32_t) const;
    bool   operator==( const Point<3U>& ) const;
    bool   operator!=( const Point<3U>& ) const;
    bool   operator<( const Point<3U>& ) const;
    bool   operator>( const Point<3U>& ) const;
    void   Set( const std::vector<double>& );
    void   Set( double, double, double );
    double Length() const;
    double SquaredLength() const;
    void   NormalizeLengthTo( double len=1. );
    double DistanceTo( const Point& ) const;
    bool   CoincidesWithWithinTolerance( const Point&, double tolerance=1.0e-5 ) const;
    bool   IsBetween( const Point& pt1, const Point& pt2 );
    std::vector<double> Coordinates() const;
    void   Out() const;
 
    friend Point<3U> operator-( double, const Point<3U>& );
    friend Point<3U> operator+( double, const Point<3U>& );
    friend Point<3U> operator*( double, const Point<3U>& );

  protected:
    double x_, y_, z_;
};

/// The length of the exterior product.
template<uint32_t dim>
double exteriorProductLength( const Point<dim>& p1, const Point<dim>& p2 );

Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 );

std::ostream&   operator<<( std::ostream& stream, const Point<1U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<2U>& pt );


std::ostream&   operator<<( std::ostream& stream, const Point<3U>& pt );

} // end namespace csmp

#endif

