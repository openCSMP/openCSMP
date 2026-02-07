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
    Point<1U>& operator=( double );
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
    Point<2U>& operator=( double );
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
    Point<3U>& operator=( double );
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
    std::vector<double> Coordinates() const;
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

} // end namespace csmp

#endif

