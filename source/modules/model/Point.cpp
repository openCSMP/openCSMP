#include "Point.h"
#include "CSMP_definitions.h"
#include <cassert>

using namespace std;

namespace csmp {

// generic boring non-specialized operators for the class
template<size_t dim>
Point<dim>::Point( double64 val )
 {
    for ( size_t i=0U; i<dim; i++ ) 
      xyz_[i] = val;
 }



  
template<size_t dim>
Point<dim>::Point( const vector<double64>& v )
 {
    assert( v.size()==dim );
    for ( size_t i=0U; i<dim; i++ ) 
      xyz_[i] = v[i];
 }



template<size_t dim>
void Point<dim>::Set( const vector<double64>& v )
 {
    assert( v.size()==dim );
    for ( size_t i=0U; i<dim; i++ ) 
      xyz_[i] = v[i];
 }



template<size_t dim>
vector<double64> Point<dim>::Coordinates() const
 {
    vector<double64> temp;
    temp.reserve(dim);
    for ( size_t i=0U; i<dim; i++ )  
      temp.push_back( xyz_[i] );
      
    return temp;
 }

 
template<size_t dim>
Point<dim>& Point<dim>::operator=( const Point<dim>& pt )
 {
    if ( &pt != this ) {
         for ( size_t i=0U; i<dim; i++ ) 
           xyz_[i] = pt.xyz_[i];
      }
    return *this;
 }
 
    
template<size_t dim>
Point<dim>  Point<dim>::operator+( const Point<dim>& pt ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] + pt.xyz_[i];
      
    return temp;
 }
 
    
template<size_t dim>
Point<dim>  Point<dim>::operator-( const Point<dim>& pt ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] - pt.xyz_[i];
      
    return temp;
 }
 
    
template<size_t dim>
Point<dim>  Point<dim>::operator*( const Point<dim>& pt ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] * pt.xyz_[i];
      
    return temp;
 }
 
 
template<size_t dim>
Point<dim>  Point<dim>::operator/( const Point<dim>& pt ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] / pt.xyz_[i];
      
    return temp;
 }
 
 
template<size_t dim>
Point<dim>  Point<dim>::operator+( double64 val ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] + val;
      
    return temp;
 }
 
    
template<size_t dim>
Point<dim>  Point<dim>::operator-( double64 val ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] - val;
      
    return temp;
 }
 
    
template<size_t dim>
Point<dim>  Point<dim>::operator*( double64 val ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] * val;
      
    return temp;
 }
 
 
template<size_t dim>
Point<dim>  Point<dim>::operator/( double64 val ) const
 {
    Point<dim> temp;
    for ( size_t i=0U; i<dim; i++ ) 
      temp[i] = xyz_[i] / val;
      
    return temp;
 }
 
 
 
 
 
template<size_t dim>
Point<dim>& Point<dim>::operator+=( const Point<dim>& pt )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] += pt.xyz_[i];
    return *this;
 }
 
 
template<size_t dim>
Point<dim>& Point<dim>::operator-=( const Point<dim>& pt )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] -= pt.xyz_[i];
    return *this;
 }
 
 
template<size_t dim>
Point<dim>& Point<dim>::operator*=( const Point<dim>& pt )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] *= pt.xyz_[i];
    return *this;
 }
 
 
template<size_t dim>
Point<dim>& Point<dim>::operator/=( const Point<dim>& pt )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] /= pt.xyz_[i];
    return *this;
 }
 
    
template<size_t dim>
Point<dim>& Point<dim>::operator+=( double64 val )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] += val;
    return *this;
 }


template<size_t dim>
Point<dim>& Point<dim>::operator-=( double64 val )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] -= val;
    return *this;
 }


template<size_t dim>
Point<dim>& Point<dim>::operator*=( double64 val )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] *= val;
    return *this;
 }
 
 
template<size_t dim>
Point<dim>& Point<dim>::operator/=( double64 val )
 {
    for ( size_t i=0U; i<dim; i++ ) xyz_[i] /= val;
    return *this;
 }
 
 
template<size_t dim>
bool Point<dim>::operator==( const Point<dim>& pt ) const
 {
    for ( size_t i=0U; i<dim; i++ )
      if ( xyz_[i] != pt.xyz_[i] ) return false;
    return true;
 }

 
template<size_t dim>
bool Point<dim>::operator!=( const Point<dim>& pt ) const
 {
    for ( size_t i=0U; i<dim; i++ )
      if ( xyz_[i] != pt.xyz_[i] ) return true;
    return false;
 }

 
template<size_t dim>
bool Point<dim>::operator<( const Point<dim>& )  const
 {
   cout <<"\nPoint<dim>::operator<() not defined for generic case."<< endl;
   return false;
 }


template<size_t dim>
bool Point<dim>::operator>( const Point<dim>& )  const
 {
   cout <<"\nPoint<dim>::operator>() not defined for generic case."<< endl;
   return false;
 }


// distance of point from origin
template<size_t dim>
double64  Point<dim>::Length() const
 {
    double64 sum_squares(0.);
    for ( size_t i=0U; i<dim; i++ )
      sum_squares += xyz_[i] * xyz_[i];
    
    return sqrt( sum_squares ); 
 }


template<size_t dim>
void  Point<dim>::NormalizeLengthTo( double64 len )
 {
    *this /= Length(); 
    *this *= len;
 }


template<size_t dim>
bool Point<dim>::CoincidesWithWithinTolerance( const Point<dim>& pt, 
                                               double64 tolerance ) const
 {
    for ( size_t i=0U; i<dim; i++ )
      if ( DistanceTo( pt ) > tolerance ) return false;
    return true;
 }

template<size_t dim>
bool Point<dim>::IsBetween( const Point<dim>& pt1, const Point<dim>& pt2 )
 {
    if ( (xyz_[0] > pt1.xyz_[0]) && (xyz_[1] > pt1.xyz_[1]) && (xyz_[2] > pt1.xyz_[2]) &&
         (xyz_[0] < pt2.xyz_[0]) && (xyz_[1] < pt2.xyz_[1]) && (xyz_[2] < pt2.xyz_[2]) ) return true;
    return false;
}

template<size_t dim>
void Point<dim>::Out() const
 {
    cout <<"\nPoint<"<<  dim <<">::Out(): coordinates: ";
    for ( size_t i=0U; i<dim; i++ ) cout << xyz_[i] <<" ";
    cout << endl;
    cout.flush();
 }


template<size_t dim>
ostream&  operator<<( ostream& stream, const Point<dim>& pt )
 {
     stream <<"\ncsmp::Point<"<< dim <<"> "; 
     for ( size_t i=0; i<dim; i++ ) stream << pt[i] <<" ";
     return stream;
 }

ostream&  operator<<( ostream& stream, const Point<1U>& pt )
 {
     stream <<"\ncsmp::Point<1> "<< pt[0U]; 
     return stream;
 }


ostream&  operator<<( ostream& stream, const Point<2U>& pt )
 {
     stream <<"\ncsmp::Point<2> ";
     stream << pt[0U] <<" "<< pt[1U]; 
     return stream;
 }


ostream&  operator<<( ostream& stream, const Point<3U>& pt )
 {
     stream <<"\ncsmp::Point<3> ";
     stream << pt[0U] <<" "<< pt[1U] <<" "<< pt[2U]; 
     return stream;
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


#ifndef _MSC_VER
//template class Point<4U>;
template class Point<1U>;
template class Point<2U>;
template class Point<3U>;
#endif
} // end namespace csmp
