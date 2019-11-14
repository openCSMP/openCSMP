#include <cassert>
#include "Point.h"
#include "CSMP_definitions.h"

using namespace std;

namespace csmp {

// generic non-specialized operators for the class

template<size_t dim>
bool Point<dim>::operator<( const Point<dim>& )  const
 {
   cerr <<"\nPoint<dim>::operator<() not defined for generic case."<< endl;
   return false;
 }


template<size_t dim>
bool Point<dim>::operator>( const Point<dim>& )  const
 {
   cerr <<"\nPoint<dim>::operator>() not defined for generic case."<< endl;
   return false;
 }



template<size_t dim>
void  Point<dim>::NormalizeLengthTo( double64 len )
 {
    *this /= Length(); 
    *this *= len;
 }


ostream&  operator<<( ostream& stream, const Point<1U>& pt )
 {
     stream <<"csmp::Point<1> "<< pt[0U];
     return stream;
 }


ostream&  operator<<( ostream& stream, const Point<2U>& pt )
 {
     stream <<"csmp::Point<2> ";
     stream << pt[0U] <<" "<< pt[1U]; 
     return stream;
 }


ostream&  operator<<( ostream& stream, const Point<3U>& pt )
 {
     stream <<"csmp::Point<3> ";
     stream << pt[0U] <<" "<< pt[1U] <<" "<< pt[2U]; 
     return stream;
 }



// ---------------------------------------------------------------------------

// MEMBER DEFINITIONS

// ---------------------------------------------------------------------------
template<size_t dim>
Point<dim>::~Point()
 {
 }

/*
template<size_t dim>
double64& Point<dim>::operator[](size_t i)
 {
    return xyz_[i];
 }
 
 
template<size_t dim>
double64 Point<dim>::operator[](size_t i) const
 {
    return xyz_[i];
 }
*/



// ----------------------------------------------------------------------------
//
//  1D SPECIALIZATION
//
// ----------------------------------------------------------------------------

Point<1U>::Point( double64 val ) : x_(val)
 {
 }

Point<1U>::Point( const Point<1U>& pt ) : x_(pt.x_)
 {
 }

Point<1U>::Point( const std::vector<double64>& v )
 : x_(v[0])
 {
    assert( v.size() == 1U );
 }

Point<1U>::~Point()
 {
 }


double64& Point<1U>::operator[](size_t)
 {
    return x_;
 }

const double64& Point<1U>::operator[](size_t) const
 {
    return x_;
 }


void Point<1U>::Set( const std::vector<double64>& v )
 {
    x_ = v[0];
 }


Point<1U>& Point<1U>::operator=( const Point<1U>& pt )
 {
    if ( &pt != this ) x_ = pt.x_;
    return *this;
 }

Point<1U>& Point<1U>::operator=( double64 val )
 {
    x_ = val;
    return *this;
 }
    

Point<1U>  Point<1U>::operator+( const Point<1U>& pt ) const
 {
    return (Point<1U>(x_ + pt.x_));
 }

Point<1U> Point<1U>::operator-( const Point<1U>& pt ) const
 {
    return (Point<1U>(x_ - pt.x_));
 }

Point<1U> Point<1U>::operator*( const Point<1U>& pt ) const
 {
    return (Point<1U>(x_ * pt.x_));
 }

Point<1U>  Point<1U>::operator/( const Point<1U>& pt ) const
 {
    return (Point<1U>(x_ / pt.x_));
 }
    

Point<1U>  Point<1U>::operator+( double64 val ) const
 {
    return (Point<1U>(x_ + val));
 }

Point<1U> Point<1U>::operator-( double64 val ) const
 {
    return (Point<1U>(x_ - val));
 }

Point<1U> Point<1U>::operator*( double64 val ) const
 {
    return (Point<1U>(x_ * val));
 }

Point<1U> Point<1U>::operator/( double64 val ) const
 {
    return (Point<1U>(x_ / val));
 }



Point<1U>& Point<1U>::operator+=( const Point<1U>& pt )
 {
    x_ += pt.x_;
    return *this;
 }

Point<1U>& Point<1U>::operator-=( const Point<1U>& pt )
 {
    x_ -= pt.x_;
    return *this;
 }

Point<1U>& Point<1U>::operator*=( const Point<1U>& pt )
 {
    x_ *= pt.x_;
    return *this;
 }

Point<1U>& Point<1U>::operator/=( const Point<1U>& pt )
 {
    x_ /= pt.x_;
    return *this;
 }
 

Point<1U>& Point<1U>::operator+=( double64 val )
 {
    x_ += val;
    return *this;
 }

Point<1U>& Point<1U>::operator-=( double64 val )
 {
    x_ -= val;
    return *this;
 }

Point<1U>& Point<1U>::operator*=( double64 val )
 {
    x_ *= val;
    return *this;
 }

Point<1U>& Point<1U>::operator/=( double64 val )
 {
    x_ /= val;
    return *this;
 }
 

bool Point<1U>::operator==( const Point<1U>& pt ) const
 {
    if ( x_ != pt.x_ ) return false;
    return true;
 }

bool Point<1U>::operator!=( const Point<1U>& pt ) const
 {
    if ( x_ != pt.x_ ) return true;
    return false;
 }

bool Point<1U>::operator<( const Point<1U>& pt ) const
 {
   if ( x_ < pt.x_ ) return true;
   return false;
 }

bool Point<1U>::operator>( const Point<1U>& pt ) const
 {
   if ( x_ > pt.x_ ) return true;
   return false;
 }


double64  Point<1U>::Length() const
 {
    return std::fabs(x_);
 }

double64  Point<1U>::SquaredLength() const
 {
    return x_ * x_; 
 }


void Point<1U>::NormalizeLengthTo( double64 len ) 
 {
    x_ = len; 
 }


double64  Point<1U>::DistanceTo( const Point<1U>& pt ) const
 {
    return std::fabs(pt.x_ - x_); 
 }


bool Point<1U>::CoincidesWithWithinTolerance( const Point<1U>& pt, 
                                              double64 tolerance ) const
 {
    if ( std::fabs(pt.x_ - x_) > tolerance ) return false;
    return true;
 }

bool Point<1U>::IsBetween( const Point<1U>& pt1, const Point<1U>& pt2)
 {
    if ( (x_ > pt1.x_)  && (x_ < pt2.x_) ) return true;
    return false;
 }


std::vector<double64> Point<1U>::Coordinates() const
 {
    return std::vector<double64>(1U,x_);
 }


void Point<1U>::Out() const
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

Point<2U>::Point(double64 val)
 : x_(val), y_(val)
 {
 }

Point<2U>::Point( double64 px, double64 py ) : x_(px), y_(py)
 {
 }

Point<2U>::Point( const Point<2U>& pt ) : x_(pt.x_), y_(pt.y_)
 {
 }

Point<2U>::Point( const std::vector<double64>& v )
 : x_(v[0]), y_(v[1])
 {
 }

Point<2U>::~Point()
 {
 }


double64& Point<2U>::operator[]( size_t i )
 {
    return (i == 0U) ? x_ : y_;
 }

const double64& Point<2U>::operator[]( size_t i ) const
 {
    return (i == 0U) ? x_ : y_;
 }


void Point<2U>::Set( const std::vector<double64>& v )
 {
    assert( v.size() == 2U );
    x_ = v[0];
    y_ = v[1];
 }


void Point<2U>::Set( double64 px, double64 py )
 {
    x_ = px;
    y_ = py;
 }


Point<2U>& Point<2U>::operator=( const Point<2U>& pt )
 {
    if ( &pt != this ) {
         x_ = pt.x_;
         y_ = pt.y_;
      }
    return *this;
 }


Point<2U>& Point<2U>::operator=( double64 val )
 {
    x_ = val;
    y_ = val;
    return *this;
 }
    

Point<2U>  Point<2U>::operator+( const Point<2U>& pt ) const
 {
    return (Point<2U>(x_ + pt.x_,y_ + pt.y_));
 }

Point<2U>  Point<2U>::operator-( const Point<2U>& pt ) const
 {
    return (Point<2U>(x_ - pt.x_,y_ - pt.y_));
 }

Point<2U> Point<2U>::operator*( const Point<2U>& pt ) const
 {
    return (Point<2U>(x_ * pt.x_,y_ * pt.y_));
 }

Point<2U>  Point<2U>::operator/( const Point<2U>& pt ) const
 {
    return (Point<2U>(x_ / pt.x_,y_ / pt.y_));
 }
    

Point<2U>  Point<2U>::operator+( double64 val ) const
 {
    return (Point<2U>(x_ + val, y_ + val));
 }

Point<2U>  Point<2U>::operator-( double64 val ) const
 {
    return (Point<2U>(x_ - val, y_ - val));
 }

Point<2U>  Point<2U>::operator*( double64 val ) const
 {
    return (Point<2U>(x_ * val, y_ * val));
 }

Point<2U>  Point<2U>::operator/( double64 val ) const
 {
    return (Point<2U>(x_ / val, y_ / val));
 }



Point<2U>& Point<2U>::operator+=( const Point<2U>& pt )
 {
    x_ += pt.x_;
    y_ += pt.y_;
    return *this;
 }

Point<2U>& Point<2U>::operator-=( const Point<2U>& pt )
 {
    x_ -= pt.x_;
    y_ -= pt.y_;
    return *this;
 }

Point<2U>& Point<2U>::operator*=( const Point<2U>& pt )
 {
    x_ *= pt.x_;
    y_ *= pt.y_;
    return *this;
 }

Point<2U>& Point<2U>::operator/=( const Point<2U>& pt )
 {
    x_ /= pt.x_;
    y_ /= pt.y_;
    return *this;
 }
 

Point<2U>& Point<2U>::operator+=( double64 val )
 {
    x_ += val;
    y_ += val;
    return *this;
 }

Point<2U>& Point<2U>::operator-=( double64 val )
 {
    x_ -= val;
    y_ -= val;
    return *this;
 }

Point<2U>& Point<2U>::operator*=( double64 val )
 {
    x_ *= val;
    y_ *= val;
    return *this;
 }

Point<2U>& Point<2U>::operator/=( double64 val )
 {
    x_ /= val;
    y_ /= val;
    return *this;
 }
 

bool Point<2U>::operator==( const Point<2U>& pt ) const
 {
    return (x_ == pt.x_) && (y_ == pt.y_);
 }

bool Point<2U>::operator!=( const Point<2U>& pt ) const
 {
    return !(*this == pt);
 }

// see M. J. Lazlo, Computational Geometry & Computer Graphics in C++, 1996, p. 75
bool Point<2U>::operator<( const Point<2U>& pt )  const
 {
    return (( x_ < pt.x_ ) || (( x_ == pt.x_ ) && ( y_ < pt.y_ )));
 }

bool Point<2U>::operator>( const Point<2U>& pt )  const
 {
    return (( x_ > pt.x_ ) || (( x_ == pt.x_ ) && ( y_ > pt.y_ )));
 }


double64  Point<2U>::Length() const
 {
    return std::hypot( x_, y_ );
 }

double64  Point<2U>::SquaredLength() const
 {
    return x_*x_ + y_*y_;
 }


void Point<2U>::NormalizeLengthTo( double64 len )
 {
    *this /= Length();
    *this *= len; 
 }


double64  Point<2U>::DistanceTo( const Point<2U>& pt ) const
 {
    return Point<2U>(pt - *this).Length();
 }


bool Point<2U>::CoincidesWithWithinTolerance( const Point<2U>& pt, 
                                              double64 tolerance ) const
 {
    if ( Point<2U>(pt - *this).Length() > tolerance ) return false;
    return true;
 }

bool Point<2U>::IsBetween( const Point<2U>& pt1, const Point<2U>& pt2 )
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (x_ < pt2.x_) && (y_ < pt2.y_) ) return true;
    return false;
}


std::vector<double64> Point<2U>::Coordinates() const
 {
    return std::vector<double64>{x_,y_};
 }



void Point<2U>::Out() const
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

Point<3U>::Point( double64 val )
 : x_(val), y_(val), z_(val)
 {
 }

Point<3U>::Point( double64 px, double64 py, double64 pz ) : x_(px), y_(py), z_(pz)
 {
 }

Point<3U>::Point( const Point<3U>& pt ) : x_(pt.x_), y_(pt.y_), z_(pt.z_)
 {
 }

Point<3U>::Point( const std::vector<double64>& v )
 : x_(v[0]), y_(v[1]), z_(v[2])
 {
    assert( v.size() == 3U );
 }

Point<3U>::~Point()
 {
 }

double64& Point<3U>::operator[]( size_t i )
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }

const double64& Point<3U>::operator[]( size_t i ) const
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }


void Point<3U>::Set( const std::vector<double64>& v )
 {
    assert( v.size() == 3U );
    x_ = v[0];
    y_ = v[1];
    z_ = v[2];
 }


void Point<3U>::Set( double64 px, double64 py, double64 pz )
 {
    x_ = px;
    y_ = py;
    z_ = pz;
 }


Point<3U>& Point<3U>::operator=( const Point<3U>& pt )
 {
    if ( &pt != this ) {
         x_ = pt.x_;
         y_ = pt.y_;
         z_ = pt.z_;
      }
    return *this;
 }


Point<3U>& Point<3U>::operator=( double64 val )
 {
    x_ = val;
    y_ = val;
    z_ = val;
    return *this;
 }
    

Point<3U>  Point<3U>::operator+( const Point<3U>& pt ) const
 {
    return (Point<3U>(x_ + pt.x_, y_ + pt.y_, z_ + pt.z_));
 }

Point<3U>  Point<3U>::operator-( const Point<3U>& pt ) const
 {
    return (Point<3U>(x_ - pt.x_,y_ - pt.y_, z_ - pt.z_));
 }

Point<3U>  Point<3U>::operator*( const Point<3U>& pt ) const
 {
    return (Point<3U>(x_ * pt.x_,y_ * pt.y_, z_ * pt.z_));
 }

Point<3U>  Point<3U>::operator/( const Point<3U>& pt ) const
 {
    return (Point<3U>(x_ / pt.x_,y_ / pt.y_, z_ / pt.z_));
 }
    

Point<3U>  Point<3U>::operator+( double64 val ) const
 {
    return (Point<3U>(x_ + val, y_ + val, z_ + val));
 }

Point<3U>  Point<3U>::operator-( double64 val ) const
 {
    return (Point<3U>(x_ - val, y_ - val, z_ - val));
 }

Point<3U>  Point<3U>::operator*( double64 val ) const
 {
    return (Point<3U>(x_ * val, y_ * val, z_ * val));
 }

Point<3U>  Point<3U>::operator/( double64 val ) const
 {
    return (Point<3U>(x_ / val, y_ / val, z_ / val));
 }


Point<3U>& Point<3U>::operator+=( const Point<3U>& pt )
 {
    x_ += pt.x_;
    y_ += pt.y_;
    z_ += pt.z_;
    return *this;
 }

Point<3U>& Point<3U>::operator-=( const Point<3U>& pt )
 {
    x_ -= pt.x_;
    y_ -= pt.y_;
    z_ -= pt.z_;
    return *this;
 }

Point<3U>& Point<3U>::operator*=( const Point<3U>& pt )
 {
    x_ *= pt.x_;
    y_ *= pt.y_;
    z_ *= pt.z_;
    return *this;
 }

Point<3U>& Point<3U>::operator/=( const Point<3U>& pt )
 {
    x_ /= pt.x_;
    y_ /= pt.y_;
    z_ /= pt.z_;
    return *this;
 }
 

Point<3U>& Point<3U>::operator+=( double64 val )
 {
    x_ += val;
    y_ += val;
    z_ += val;
    return *this;
 }

Point<3U>& Point<3U>::operator-=( double64 val )
 {
    x_ -= val;
    y_ -= val;
    z_ -= val;
    return *this;
 }

Point<3U>& Point<3U>::operator*=( double64 val )
 {
    x_ *= val;
    y_ *= val;
    z_ *= val;
    return *this;
 }

Point<3U>& Point<3U>::operator/=( double64 val )
 {
    x_ /= val;
    y_ /= val;
    z_ /= val;
    return *this;
 }
 

bool Point<3U>::operator==( const Point<3U>& pt ) const
 {
    return (x_ == pt.x_) && (y_ == pt.y_) && (z_ == pt.z_);
 }

bool Point<3U>::operator!=( const Point<3U>& pt ) const
 {
    return !(*this == pt);
 }

bool Point<3U>::operator<( const Point<3U>& p )  const
 {
    return x_<p.x_ || (x_==p.x_ && y_<p.y_) || (x_==p.x_ &&  y_==p.y_ && z_<p.z_);
 }

bool Point<3U>::operator>( const Point<3U>& p )  const
 {
    return x_>p.x_ || (x_==p.x_ && y_>p.y_) || (x_==p.x_ && y_==p.y_ && z_>p.z_);
 }


double64  Point<3U>::Length() const
 {
    return std::sqrt( x_ * x_ + y_ * y_ + z_ * z_ );
 }

double64  Point<3U>::SquaredLength() const
 {
    return x_ * x_ + y_ * y_ + z_ * z_;
 }


void Point<3U>::NormalizeLengthTo( double64 len )
 {
    *this /= Length();
    *this *= len; 
 }


double64  Point<3U>::DistanceTo( const Point<3U>& pt ) const
 {
    return Point<3U>(pt - *this).Length(); 
 }


bool Point<3U>::CoincidesWithWithinTolerance( const Point<3U>& pt, 
                                                     double64 tolerance ) const
 {
    if ( Point<3U>(pt - *this).Length() > tolerance ) return false;
    return true;
 }


bool Point<3U>::IsBetween( const Point<3U>& pt1, const Point<3U>& pt2 )
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (z_ > pt1.z_) && (x_ < pt2.x_) && (y_ < pt2.y_) && (z_ < pt2.z_) ) return true;
    return false;
}


std::vector<double64> Point<3U>::Coordinates() const
 {
    return std::vector<double64>{x_,y_,z_};
 }



void Point<3U>::Out() const
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
double64  Point<dim>::DistanceTo( const Point& pt ) const
 {
    return Point<dim>( pt - *this ).Length();
 }

/*
template<size_t dim>
Point<dim> operator-( double64 val, const Point<dim>& pt )
 {
    Point<dim> temp(val);
    return (temp - pt);
 } 

template<size_t dim>
Point<dim> operator+( double64 val, const Point<dim>& pt )
 {
    return (pt + val);
 } 
  
template<size_t dim>
Point<dim> operator*( double64 val, const Point<dim>& pt )
 {
    return (pt * val);
 }   
*/

// 1D

Point<1U> operator-( double64 val, const Point<1U>& pt )
 {
    return (Point<1U>( val - pt.x_ ));
 }   

Point<1U> operator+( double64 val, const Point<1U>& pt )
 {
    return (Point<1U>( pt.x_ + val ));
 }   

Point<1U> operator*( double64 val, const Point<1U>& pt )
 {
    return (Point<1U>( pt.x_ * val ));
 }   

// 2D

Point<2U> operator-( double64 val, const Point<2U>& pt )
 {
    return (Point<2U>( val - pt.x_, val - pt.y_ ));
 }   

Point<2U> operator+( double64 val, const Point<2U>& pt )
 {
    return (Point<2U>( pt.x_ + val, pt.y_ + val ));
 }   

Point<2U> operator*( double64 val, const Point<2U>& pt )
 {
    return (Point<2U>( pt.x_ * val, pt.y_ * val ));
 }   
 
// 3D

Point<3U> operator-( double64 val, const Point<3U>& pt )
 {
    return (Point<3U>( val - pt.x_, val - pt.y_, val - pt.z_ ));
 }   

Point<3U> operator+( double64 val, const Point<3U>& pt )
 {
    return (Point<3U>( pt.x_ + val, pt.y_ + val, pt.z_ + val ));
 }   

Point<3U> operator*( double64 val, const Point<3U>& pt )
 {
    return (Point<3U>( pt.x_ * val, pt.y_ * val, pt.z_ * val ));
 }   


template<size_t dim>
Point<dim>  midPoint( const Point<dim>& p1, const Point<dim>& p2 )
 { 
    return (Point<dim>( p1 + p2 ) / 2.);
 }

template Point<1U>  midPoint( const Point<1U>&, const Point<1U>& );
template Point<2U>  midPoint( const Point<2U>&, const Point<2U>& );
template Point<3U>  midPoint( const Point<3U>&, const Point<3U>& );



// dot = scalar product
double64 dotProduct( const Point<1U>& p1, const Point<1U>& p2 )
  {
    return p1[0U] * p2[0U];
  }

double64 dotProduct( const Point<2U>& p1, const Point<2U>& p2 )
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U];
  }

double64 dotProduct( const Point<3U>& p1, const Point<3U>& p2 )
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U] + p1[2U] * p2[2U];
  }
 

Point<1U> crossProduct( const Point<1U>& p1, const Point<1U>& p2 )
  {
     return (Point<1U>( p1[0U] * p2[0U] ));
  }



template<>
double64 exteriorProductLength( const Point<1u>& p1, const Point<1u>& p2 )
{
    return 0.0;
}

template<>
double64 exteriorProductLength( const Point<2u>& p1, const Point<2u>& p2 )
{
    return p1[0] * p2[1] - p1[1] * p2[0];
}

template<>
double64 exteriorProductLength( const Point<3u>& p1, const Point<3u>& p2 )
{
    return crossProduct(p1,p2).Length();
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
Point<2U> crossProduct( const Point<2U>& p1, const Point<2U>& p2 )
  {
     return (Point<2U>( p1[0U] * p2[1U] - p2[0U] * p1[1U], std::numeric_limits<double64>::quiet_NaN() ));
  }

 // tested: SKM O.K.
Point<3U> crossProduct( const Point<3U>& p1, const Point<3U>& p2 )
  {
     return (Point<3U>( p1[1U] * p2[2U] - p1[2U] * p2[1U],
                                 p1[2U] * p2[0U] - p1[0U] * p2[2U],
                                 p1[0U] * p2[1U] - p1[1U] * p2[0U] ));
  }


// ----------------------------------------------------------------------------------------
//
// operators and other functions involving points
//
// ----------------------------------------------------------------------------------------
template<size_t dim>
Point<dim> crossProduct( const Point<dim>&, const Point<dim>& )
  {
     cerr <<"\ncrossProduct<dim>: not defined for generic case."<< endl;
     return std::numeric_limits<double64>::quiet_NaN();
  }

} // end namespace csmp
