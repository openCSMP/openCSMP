#include "Point.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

// generic non-specialized operators for the class

template<uint32_t dim>
bool Point<dim>::operator<( const Point<dim>& )  const
 {
   cerr <<"\nPoint<dim>::operator<() not defined for generic case."<< endl;
   return false;
 }


template<uint32_t dim>
bool Point<dim>::operator>( const Point<dim>& )  const
 {
   cerr <<"\nPoint<dim>::operator>() not defined for generic case."<< endl;
   return false;
 }



template<uint32_t dim>
void  Point<dim>::NormalizeLengthTo( double len )
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

/*
template<uint32_t dim>
double& Point<dim>::operator[]( uint32_t i)
 {
    return xyz_[i];
 }
 
 
template<uint32_t dim>
double Point<dim>::operator[]( uint32_t i) const
 {
    return xyz_[i];
 }
*/



// ----------------------------------------------------------------------------
//
//  1D SPECIALIZATION
//
// ----------------------------------------------------------------------------

Point<1U>::Point( double val ) : x_(val)
 {
 }

Point<1U>::Point( const array<double,1U>& v )
 : x_(v[0])
 {
 }


Point<1U>::Point( const vector<double>& v )
 : x_(v[0])
 {
    assert( v.size() == 1U );
 }


double& Point<1U>::operator[]( uint32_t)
 {
    return x_;
 }

const double& Point<1U>::operator[]( uint32_t) const
 {
    return x_;
 }


void Point<1U>::Set( const array<double,1U>& v )
 {
    x_ = v[0];
 }


void Point<1U>::Set( const vector<double>& v )
 {
    assert( v.size() >= 1 );
    x_ = v[0];
 }


Point<1U>& Point<1U>::operator=( double val )
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
    

Point<1U>  Point<1U>::operator+( double val ) const
 {
    return (Point<1U>(x_ + val));
 }

Point<1U> Point<1U>::operator-( double val ) const
 {
    return (Point<1U>(x_ - val));
 }

Point<1U> Point<1U>::operator*( double val ) const
 {
    return (Point<1U>(x_ * val));
 }

Point<1U> Point<1U>::operator/( double val ) const
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
 

Point<1U>& Point<1U>::operator+=( double val )
 {
    x_ += val;
    return *this;
 }

Point<1U>& Point<1U>::operator-=( double val )
 {
    x_ -= val;
    return *this;
 }

Point<1U>& Point<1U>::operator*=( double val )
 {
    x_ *= val;
    return *this;
 }

Point<1U>& Point<1U>::operator/=( double val )
 {
    x_ /= val;
    return *this;
 }
 

bool Point<1U>::operator==( const Point<1U>& pt ) const
 {
    if ( essentiallyEqual(x_,pt.x_) ) return true;
    return false;
 }

bool Point<1U>::operator!=( const Point<1U>& pt ) const
 {
    return !(*this == pt);
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


double  Point<1U>::Length() const
 {
    return fabs(x_);
 }

double  Point<1U>::SquaredLength() const
 {
    return x_ * x_; 
 }


void Point<1U>::NormalizeLengthTo( double len ) 
 {
    x_ = len; 
 }


double  Point<1U>::DistanceTo( const Point<1U>& pt ) const
 {
    return fabs(pt.x_ - x_);
 }


bool Point<1U>::CoincidesWithWithinTolerance( const Point<1U>& pt, 
                                              double tolerance ) const
 {
    if ( fabs(pt.x_ - x_) > tolerance ) return false;
    return true;
 }

bool Point<1U>::IsBetween( const Point<1U>& pt1, const Point<1U>& pt2)
 {
    if ( (x_ > pt1.x_)  && (x_ < pt2.x_) ) return true;
    return false;
 }


vector<double> Point<1U>::Coordinates() const
 {
    return vector<double>(1U,x_);
 }


void Point<1U>::Out() const
 {
    cout <<"\nPoint<" << 1U;
    cout <<">::Out(): coordinates: "<< x_ << endl;
    cout.flush();
 }




// ----------------------------------------------------------------------------
//
//  2D SPECIALIZATION
//
// ----------------------------------------------------------------------------

Point<2U>::Point(double val)
 : x_(val), y_(val)
 {
 }

Point<2U>::Point( double px, double py ) : x_(px), y_(py)
 {
 }

Point<2U>::Point( const vector<double>& v )
 : x_(v[0]), y_(v[1])
 {
 }

double& Point<2U>::operator[]( uint32_t i )
 {
    return (i == 0U) ? x_ : y_;
 }

const double& Point<2U>::operator[]( uint32_t i ) const
 {
    return (i == 0U) ? x_ : y_;
 }


void Point<2U>::Set( const array<double,2U>& v )
 {
    x_ = v[0];
    y_ = v[1];
 }


void Point<2U>::Set( const vector<double>& v )
 {
    assert( v.size() == 2U );
    x_ = v[0];
    y_ = v[1];
 }


void Point<2U>::Set( double px, double py )
 {
    x_ = px;
    y_ = py;
 }


Point<2U>& Point<2U>::operator=( double val )
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
    

Point<2U>  Point<2U>::operator+( double val ) const
 {
    return (Point<2U>(x_ + val, y_ + val));
 }

Point<2U>  Point<2U>::operator-( double val ) const
 {
    return (Point<2U>(x_ - val, y_ - val));
 }

Point<2U>  Point<2U>::operator*( double val ) const
 {
    return (Point<2U>(x_ * val, y_ * val));
 }

Point<2U>  Point<2U>::operator/( double val ) const
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
 

Point<2U>& Point<2U>::operator+=( double val )
 {
    x_ += val;
    y_ += val;
    return *this;
 }

Point<2U>& Point<2U>::operator-=( double val )
 {
    x_ -= val;
    y_ -= val;
    return *this;
 }

Point<2U>& Point<2U>::operator*=( double val )
 {
    x_ *= val;
    y_ *= val;
    return *this;
 }

Point<2U>& Point<2U>::operator/=( double val )
 {
    x_ /= val;
    y_ /= val;
    return *this;
 }
 

bool Point<2U>::operator==( const Point<2U>& pt ) const
 {
    return essentiallyEqual(x_,pt.x_) && essentiallyEqual(y_,pt.y_);
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


double  Point<2U>::Length() const
 {
    return hypot( x_, y_ );
 }

double  Point<2U>::SquaredLength() const
 {
    return x_*x_ + y_*y_;
 }


void Point<2U>::NormalizeLengthTo( double len )
 {
    *this /= Length();
    *this *= len; 
 }


double  Point<2U>::DistanceTo( const Point<2U>& pt ) const
 {
    return hypot(x_ - pt.x_, y_ - pt.y_);
 }


bool Point<2U>::CoincidesWithWithinTolerance( const Point<2U>& pt, 
                                              double tolerance ) const
 {
    if ( hypot(x_ - pt.x_, y_ - pt.y_) > tolerance ) return false;
    return true;
 }

bool Point<2U>::IsBetween( const Point<2U>& pt1, const Point<2U>& pt2 )
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (x_ < pt2.x_) && (y_ < pt2.y_) ) return true;
    return false;
}


vector<double> Point<2U>::Coordinates() const
 {
    return vector<double>{x_,y_};
 }



void Point<2U>::Out() const
 {
    cout <<"\nPoint<"<< 2U;
    cout <<">::Out(): coordinates: "<< x_ <<","<< y_ << endl;
    cout.flush();
 }






// ----------------------------------------------------------------------------
//
//  3D SPECIALIZATION
//
// ----------------------------------------------------------------------------

Point<3U>::Point( double val )
 : x_(val), y_(val), z_(val)
 {
 }

Point<3U>::Point( double px, double py, double pz ) : x_(px), y_(py), z_(pz)
 {
 }


Point<3U>::Point( const array<double,3U>& v )
 : x_(v[0]), y_(v[1]), z_(v[2])
 {
 }


Point<3U>::Point( const vector<double>& v )
 : x_(v[0]), y_(v[1]), z_(v[2])
 {
 }


double& Point<3U>::operator[]( uint32_t i )
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }

const double& Point<3U>::operator[]( uint32_t i ) const
 {
    return ((i==0U) ? x_ : ((i==1U) ? y_ : z_));
 }


void Point<3U>::Set( const array<double,3U>& v )
 {
    x_ = v[0];
    y_ = v[1];
    z_ = v[2];
 }


void Point<3U>::Set( const vector<double>& v )
 {
    assert( v.size() == 3U );
    x_ = v[0];
    y_ = v[1];
    z_ = v[2];
 }


void Point<3U>::Set( double px, double py, double pz )
 {
    x_ = px;
    y_ = py;
    z_ = pz;
 }



Point<3U>& Point<3U>::operator=( double val )
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
    

Point<3U>  Point<3U>::operator+( double val ) const
 {
    return (Point<3U>(x_ + val, y_ + val, z_ + val));
 }

Point<3U>  Point<3U>::operator-( double val ) const
 {
    return (Point<3U>(x_ - val, y_ - val, z_ - val));
 }

Point<3U>  Point<3U>::operator*( double val ) const
 {
    return (Point<3U>(x_ * val, y_ * val, z_ * val));
 }

Point<3U>  Point<3U>::operator/( double val ) const
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
 

Point<3U>& Point<3U>::operator+=( double val )
 {
    x_ += val;
    y_ += val;
    z_ += val;
    return *this;
 }

Point<3U>& Point<3U>::operator-=( double val )
 {
    x_ -= val;
    y_ -= val;
    z_ -= val;
    return *this;
 }

Point<3U>& Point<3U>::operator*=( double val )
 {
    x_ *= val;
    y_ *= val;
    z_ *= val;
    return *this;
 }

Point<3U>& Point<3U>::operator/=( double val )
 {
    x_ /= val;
    y_ /= val;
    z_ /= val;
    return *this;
 }
 

bool Point<3U>::operator==( const Point<3U>& pt ) const
 {
    return essentiallyEqual(x_,pt.x_) && essentiallyEqual(y_,pt.y_) && essentiallyEqual(z_,pt.z_);
//    return !( *this < pt && pt < *this );
 }

bool Point<3U>::operator!=( const Point<3U>& pt ) const
 {
    return !(*this == pt);
 }


/**
    Lexicographical compare via standard array as needed by predicate less<> to store points in std containers..
    Else one would need iterators.
    
    @code
    std::lexicographical_compare()
    @endcode

    Old code does not implement a strict weak ordering and might cause undefined behaviour if used in a map
    @code
    return x_<p.x_ || (x_==p.x_ && y_<p.y_) || (x_==p.x_ &&  y_==p.y_ && z_<p.z_);
    @endcode
*/
bool Point<3U>::operator<( const Point<3U>& p )  const
 {
    array<double,3> p0{ x_, y_, z_ }, p1{ p.x_, p.y_, p.z_ };
    return p0 < p1;
 }

bool Point<3U>::operator>( const Point<3U>& p )  const
 {
    array<double,3> p0{ x_, y_, z_ }, p1{ p.x_, p.y_, p.z_ };
    return p0 > p1;
 }


double  Point<3U>::Length() const
 {
    return hypot( x_ , y_, z_ );
 }

double  Point<3U>::SquaredLength() const
 {
    return x_ * x_ + y_ * y_ + z_ * z_;
 }


void Point<3U>::NormalizeLengthTo( double len )
 {
    *this /= Length();
    *this *= len; 
 }


double  Point<3U>::DistanceTo( const Point<3U>& pt ) const
 {
    return hypot(hypot(x_-pt.x_,y_-pt.y_),z_-pt.z_);
 }


bool Point<3U>::CoincidesWithWithinTolerance( const Point<3U>& pt, 
                                              double tolerance ) const
 {
    if ( hypot(hypot(x_-pt.x_,y_-pt.y_),z_-pt.z_) > tolerance ) return false;
    return true;
 }


bool Point<3U>::IsBetween( const Point<3U>& pt1, const Point<3U>& pt2 )
 {
    if ( (x_ > pt1.x_) && (y_ > pt1.y_) && (z_ > pt1.z_) && (x_ < pt2.x_) && (y_ < pt2.y_) && (z_ < pt2.z_) ) return true;
    return false;
}


vector<double> Point<3U>::Coordinates() const
 {
    return vector<double>{x_,y_,z_};
 }



void Point<3U>::Out() const
 {
    cout <<"\nPoint<"<< 3U;
    cout <<">::Out(): coordinates: ";
    cout << x_ <<","<< y_ <<","<< z_ << endl;
    cout.flush();
 }





// -------------------------------------------------------------------------------
//
//    GENERIC INLINE FUNCTIONS
//
// -------------------------------------------------------------------------------
/*
template<uint32_t dim>
Point<dim>::Point()
 {
    static_assert( dim <= 3, "Point<dim>::Point: default constructor: wrong template parameter value" );
 }

template<uint32_t dim>
double  Point<dim>::DistanceTo( const Point& pt ) const
 {
    return Point<dim>( pt - *this ).Length();
 }

template<uint32_t dim>
Point<dim> operator-( double val, const Point<dim>& pt )
 {
    Point<dim> temp(val);
    return (temp - pt);
 } 

template<uint32_t dim>
Point<dim> operator+( double val, const Point<dim>& pt )
 {
    return (pt + val);
 } 
  
template<uint32_t dim>
Point<dim> operator*( double val, const Point<dim>& pt )
 {
    return (pt * val);
 }   
*/

// 1D

Point<1U> operator-( double val, const Point<1U>& pt )
 {
    return (Point<1U>( val - pt.x_ ));
 }   

Point<1U> operator+( double val, const Point<1U>& pt )
 {
    return (Point<1U>( pt.x_ + val ));
 }   

Point<1U> operator*( double val, const Point<1U>& pt )
 {
    return (Point<1U>( pt.x_ * val ));
 }   

// 2D

Point<2U> operator-( double val, const Point<2U>& pt )
 {
    return (Point<2U>( val - pt.x_, val - pt.y_ ));
 }   

Point<2U> operator+( double val, const Point<2U>& pt )
 {
    return (Point<2U>( pt.x_ + val, pt.y_ + val ));
 }   

Point<2U> operator*( double val, const Point<2U>& pt )
 {
    return (Point<2U>( pt.x_ * val, pt.y_ * val ));
 }   
 
// 3D

Point<3U> operator-( double val, const Point<3U>& pt )
 {
    return (Point<3U>( val - pt.x_, val - pt.y_, val - pt.z_ ));
 }   

Point<3U> operator+( double val, const Point<3U>& pt )
 {
    return (Point<3U>( pt.x_ + val, pt.y_ + val, pt.z_ + val ));
 }   

Point<3U> operator*( double val, const Point<3U>& pt )
 {
    return (Point<3U>( pt.x_ * val, pt.y_ * val, pt.z_ * val ));
 }   


template<uint32_t dim>
Point<dim>  midPoint( const Point<dim>& p1, const Point<dim>& p2 )
 { 
    return (Point<dim>( p1 + p2 ) / 2.);
 }

template Point<1U>  midPoint( const Point<1U>&, const Point<1U>& );
template Point<2U>  midPoint( const Point<2U>&, const Point<2U>& );
template Point<3U>  midPoint( const Point<3U>&, const Point<3U>& );



// dot = scalar product
double dotProduct( const Point<1U>& p1, const Point<1U>& p2 )
  {
    return p1[0U] * p2[0U];
  }

double dotProduct( const Point<2U>& p1, const Point<2U>& p2 )
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U];
  }

double dotProduct( const Point<3U>& p1, const Point<3U>& p2 )
  {
    return p1[0U] * p2[0U] + p1[1U] * p2[1U] + p1[2U] * p2[2U];
  }
 

Point<1U> crossProduct( const Point<1U>& p1, const Point<1U>& p2 )
  {
     return (Point<1U>( p1[0U] * p2[0U] ));
  }

/// distance between the two points
double distance( const Point<1U>& a, const Point<1U>& b ) {
     return fabs( b[0] - a[0] );
  }
double distance( const Point<2U>& a, const Point<2U>& b ) {
     return hypot( a[0]-b[0], a[1]-b[1] );
  }
double distance( const Point<3U>& a, const Point<3U>& b ) {
     // hypot(hypot(x1-x2,y1-y2),z1-z2) avoids potential roundoff-related degeneracy after squaring
     return hypot(hypot(a[0]-b[0],a[1]-b[1]),a[2]-b[2]);
  }


/// distance between point, A, and line segment BC, does not check case where projection of A is outside of segment BC
template<uint32_t dim>
double distanceFromLine( const Point<dim>& A, const Point<dim>& B, const Point<dim>& C )
 {
    // find direction vector BC
    const double length_BC = distance( C, B );
    Point<dim> dirVec = (C - B) / length_BC;
    // vector connecting A with B
    Point<dim> Vec = A - B;
    // distance of projection A onto BC from B
    const double t = dotProduct( dirVec, Vec );
    
    // if the projection of A on BC hits BC beyond endpoints
    if ( t <= 0. )
      return distance( A, B );
    else if ( t > length_BC )
      return distance( A, C );

    // projection P of A onto BC
    Point<dim> P = B + (dirVec * t);
    // distance of point A from line BC
    return distance( P, A );
 }
template double distanceFromLine( const Point<3U>&, const Point<3U>&, const Point<3U>& );
template double distanceFromLine( const Point<2U>&, const Point<2U>&, const Point<2U>& );




template<>
double exteriorProductLength( const Point<1u>& p1, const Point<1u>& p2 )
{
    return 0.0;
}

template<>
double exteriorProductLength( const Point<2u>& p1, const Point<2u>& p2 )
{
    return p1[0] * p2[1] - p1[1] * p2[0];
}

template<>
double exteriorProductLength( const Point<3u>& p1, const Point<3u>& p2 )
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
     return (Point<2U>( p1[0U] * p2[1U] - p2[0U] * p1[1U], numeric_limits<double>::quiet_NaN() ));
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
template<uint32_t dim>
Point<dim> crossProduct( const Point<dim>&, const Point<dim>& )
  {
     cerr <<"\ncrossProduct<dim>: not defined for generic case."<< endl;
     return numeric_limits<double>::signaling_NaN();
  }


/**
    Calculates smallest angle in degrees between the line segments that start with the first point and terminate at the last point of thedge
    
    @return smallest angle between edges in degrees (= always less than 90o).
    
    @note, to get the acute (smallest angle) between the edges ignoring their orientation, adjust the angle as follows:
       
    @note this is equivalent to the      acute_angle  =  (angle > 90.) ? 180. -angle : angle;
*/
template<uint32_t dim>
double angleBetweenEdges( const pair<Point<dim>,Point<dim> >& edge1, const pair<Point<dim>,Point<dim> >& edge2 )
 {
    const Point<dim> a(edge1.second - edge1.first), b(edge2.second - edge2.first);
    
    // a . b
    // -----
    double ab{0.};
    for ( auto i{0U}; i<dim; ++i ) ab += a[i] * b[i];
    
    // ||a||  ||b||
    // ------------
    double a_b{numeric_limits<double>::quiet_NaN()};
    if constexpr (dim == 2) a_b = sqrt( (a[0]*a[0]+a[1]*a[1]) * (b[0]*b[0]+b[1]*b[1]) );
    else if constexpr (dim == 3 )
      a_b = sqrt( (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) * (b[0]*b[0]+b[1]*b[1]+b[2]*b[2]) );
      
    double cos_angle = ab / a_b;

    // if zero intercept
    if ( cos_angle == 0. ) return 90.;
    // if outside of range of 'acos' function
    if ( cos_angle >  1. ) return   0.;
    if ( cos_angle < -1. ) return 180.;
        
    return (180./3.14159265358979323) * acos(cos_angle);

 } // end angleBetweenEdges

template double angleBetweenEdges( const pair<Point<2>,Point<2> >&, const pair<Point<2>,Point<2> >& );
template double angleBetweenEdges( const pair<Point<3>,Point<3> >&, const pair<Point<3>,Point<3> >& );



} // end namespace csmp
