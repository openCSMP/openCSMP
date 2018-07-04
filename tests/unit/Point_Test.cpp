/*
 *  Point_Test.cpp
 *
 *  Created by Stephan Matthai on 10/26/10.
 */

#include <limits>
#include "Point.h"
#include "Point_Test.h"

using namespace std;

namespace csmp {

Point_Test::Point_Test()
 {
 }
 
 
void Point_Test::run()
 {
    Test_1D_Point();
    Test_2D_Point();
    Test_3D_Point();
 }




void Point_Test::Test_1D_Point()
 {
   const size_t dim(1U);

   // constructing point from vector
   vector<double64>  a(1U,1.); 
   Point<dim>  p1( a );
   Point<dim>  p2( p1 );
   Point<dim>  p3; p3 = p2;
   
   // testing assignment and equal operators
   _test( p1==p3 );
   _test( !(p1!=p3) ); // operator!=( const Point& )

    // adding, multiplying and substracting points
    /* 
        Point&    operator+=( const Point& );    
        Point&    operator-=( const Point& );    
        Point&    operator*=( const Point& );    
        Point&    operator/=( const Point& );   
    */
    p2 += p1;
    p2 -= p1;
    // also testing:     double64  operator[](size_t) const;
    _equal( p1[0], p2[0], numeric_limits<double64>::epsilon() );
    _equal( p1[1], p2[1], numeric_limits<double64>::epsilon() );
    _test( p1==p2 ); // operator==( const Point& )
    p1 *= p3; 
    p1 /= p3;  
    _test( p1==p3 );

    // adding, multiplying and substracting doubles
    /* 
        Point&    operator+=( double64 );    
        Point&    operator-=( double64 );    
        Point&    operator*=( double64 );    
        Point&    operator/=( double64 );    
    */
    p1 += 10.;
    p1 -= 5.;
    p1 -= 5.;
    _test( p1==p3 );
    p2 *= 4.;
    p2 /= 2.;
    p2 /= 2.;
    _test( p1==p2 );
    
    // using temporaries
    /*
        Point     operator+( const Point& ) const;    
        Point     operator-( const Point& ) const;    
        Point     operator*( const Point& ) const;    
        Point     operator/( const Point& ) const;  
    */  
    Point<dim> res = p1 + p3;
    p2 *= 2.;
    _test( res==p2 );
    
    res = p1 - p3;
    Point<dim>  p0; 
    p0 = 0.; // assignment to double
    _test( res==p0 );
   
    p2 = p1 * p3;
    res = p2 / p3;
    _test( res==p1 );
    p2 = p2 / p3; // involving point itself
    _test( p2==p1 );

    // double64 post multiplicator
    /*
        Point     operator*( double64 ) const;    
        Point     operator/( double64 ) const;    
    */
    res = p1  * 64.5;
    p2  = res / 64.5;
    _test( p2==p1 );

    // assignments 
    p1[0] = p1[1] = p1[2];  
    _test( !(p1!=p3) );
    
    res = p1 * 5.;
    _test( p1  < res );
    _test( res > p1 );

    p1.Set( a );
    p2.Set( a );
    p3.Set( a );
    _test( p1 == p2 );
    res[0]=1.; 
    _equal( res.Length(), 1., numeric_limits<double64>::epsilon() );
    res.NormalizeLengthTo();
    _equal( res.Length(), 1., numeric_limits<double64>::epsilon() );    
    res[0]=1.; 
    p3[0]=3.;
    _equal( res.DistanceTo(p3), 2., numeric_limits<double64>::epsilon() );
    _test( !res.CoincidesWithWithinTolerance(p3) );
    vector<double64> b(p1.Coordinates());
    _test( a == b );
    p1.Out();
    cout <<"Point_Test::Test_3D_Point: ostream test: "<< p1 << endl << endl;
    
    // related non-member binary operators
    res = 1. * p1;
    _test( res == p1 );
    res = 1. + p1;
    p2  = p1 + 1.;
    _test( res == p2 );
    p1.Set( a );
    res = 3. - p1;
    _equal( res[0], 2., numeric_limits<double64>::epsilon() );
    
    // mid point
    p1[0]=0.; 
    p2[0]=1.; 
    p3[0]=2.; 
    _test( p2==midPoint(p1,p3) );
    Point<dim> pMP = midPoint(p1,p3);
    _equal( p2[0], pMP[0], numeric_limits<double64>::epsilon() );

    // dot product
    _equal( dotProduct( p2, p3 ), 2., numeric_limits<double64>::epsilon() );

    // cross product 
    p1[0]=1.; 
    p2[0]=4.;
    res = crossProduct( p1, p2 );
    _equal( res[0], 4., numeric_limits<double64>::epsilon() );

 } // end 1D test 
 


void Point_Test::Test_2D_Point()
 {
   const size_t dim(2U);

   // constructing point from vector
   vector<double64>  a(2U,1.); a[1]=2.;
   Point<dim>  p1( a );
   Point<dim>  p2( p1 );
   Point<dim>  p3; p3 = p2;
   
   // testing assignment and equal operators
   _test( p1==p3 );
   _test( !(p1!=p3) ); // operator!=( const Point& )

    // adding, multiplying and substracting points
    /* 
        Point&    operator+=( const Point& );    
        Point&    operator-=( const Point& );    
        Point&    operator*=( const Point& );    
        Point&    operator/=( const Point& );   
    */
    p2 += p1;
    p2 -= p1;
    // also testing:     double64  operator[](size_t) const;
    _equal( p1[0], p2[0], numeric_limits<double64>::epsilon() );
    _equal( p1[1], p2[1], numeric_limits<double64>::epsilon() );
    _test( p1==p2 ); // operator==( const Point& )
    p1 *= p3; 
    p1 /= p3;  
    _test( p1==p3 );

    // adding, multiplying and substracting doubles
    /* 
        Point&    operator+=( double64 );    
        Point&    operator-=( double64 );    
        Point&    operator*=( double64 );    
        Point&    operator/=( double64 );    
    */
    p1 += 10.;
    p1 -= 5.;
    p1 -= 5.;
    _test( p1==p3 );
    p2 *= 4.;
    p2 /= 2.;
    p2 /= 2.;
    _test( p1==p2 );
    
    // using temporaries
    /*
        Point     operator+( const Point& ) const;    
        Point     operator-( const Point& ) const;    
        Point     operator*( const Point& ) const;    
        Point     operator/( const Point& ) const;  
    */  
    Point<dim> res = p1 + p3;
    p2 *= 2.;
    _test( res==p2 );
    
    res = p1 - p3;
    Point<dim>  p0; 
    p0 = 0.; // assignment to double
    _test( res==p0 );
   
    p2 = p1 * p3;
    res = p2 / p3;
    _test( res==p1 );
    p2 = p2 / p3; // involving point itself
    _test( p2==p1 );

    // double64 post multiplicator
    /*
        Point     operator*( double64 ) const;    
        Point     operator/( double64 ) const;    
    */
    res = p1  * 64.5;
    p2  = res / 64.5;
    _test( p2==p1 );

    // assignments 
    p1[0] = p1[1] = p1[2];  
    _test( !(p1==p3) );
    
    res = p1 * 5.;
    _test( p1  < res );
    _test( res > p1 );

    p1.Set( a );
    p2.Set( a );
    p3.Set( a );
    _test( p1 == p2 );
    res[0]=1.; res[1]=1.;
    _equal( res.Length(), sqrt(2.), numeric_limits<double64>::epsilon() );
    res.NormalizeLengthTo();
    _equal( res.Length(), 1., numeric_limits<double64>::epsilon() );    
    res[0]=1.; res[1]=1.;
    _equal( res.DistanceTo(p3), 1., numeric_limits<double64>::epsilon() );
    _test( !res.CoincidesWithWithinTolerance(p3) );
    vector<double64> b(p1.Coordinates());
    _test( a == b );
    p1.Out();
    cout <<"Point_Test::Test_3D_Point: ostream test: "<< p1 << endl << endl;
    
    // related non-member binary operators
    res = 1. * p1;
    _test( res == p1 );
    res = 1. + p1;
    p2  = p1 + 1.;
    _test( res == p2 );
    p1.Set( a );
    res = 3. - p1;
    _equal( res[0], 2., numeric_limits<double64>::epsilon() );
    _equal( res[1], 1., numeric_limits<double64>::epsilon() );
    
    // mid point
    p1[0]=0.; p1[1]=0.;
    p2[0]=1.; p2[1]=1.;
    p3[0]=2.; p3[1]=2.;
    _test( p2==midPoint(p1,p3) );
    Point<dim> pMP = midPoint(p1,p3);
    _equal( p2[0], pMP[0], numeric_limits<double64>::epsilon() );
    _equal( p2[1], pMP[1], numeric_limits<double64>::epsilon() );

    // dot product
    _equal( dotProduct( p2, p3 ), 4., numeric_limits<double64>::epsilon() );

    // cross product 
    p1[0]=1.; p1[1]=0.;
    p2[0]=0.; p2[1]=1.;
    res = crossProduct( p1, p2 );
    _equal( res[0], 1., numeric_limits<double64>::epsilon() );
    
 } // end 2D test

 



void Point_Test::Test_3D_Point()
 {
   const size_t dim(3U);

   // constructing point from vector
   vector<double64>  a(3U,1.); a[1]=2.; a[2]=3.;
   Point<dim>  p1( a );
   Point<dim>  p2( p1 );
   Point<dim>  p3; p3 = p2;
   
   // testing assignment and equal operators
   _test( p1==p3 );
   _test( !(p1!=p3) ); // operator!=( const Point& )

    // adding, multiplying and substracting points
    /* 
        Point&    operator+=( const Point& );    
        Point&    operator-=( const Point& );    
        Point&    operator*=( const Point& );    
        Point&    operator/=( const Point& );   
    */
    p2 += p1;
    p2 -= p1;
    // also testing:     double64  operator[](size_t) const;
    _equal( p1[0], p2[0], numeric_limits<double64>::epsilon() );
    _equal( p1[1], p2[1], numeric_limits<double64>::epsilon() );
    _equal( p1[2], p2[2], numeric_limits<double64>::epsilon() );
    _test( p1==p2 ); // operator==( const Point& )
    p1 *= p3; 
    p1 /= p3;  
    _test( p1==p3 );

    // adding, multiplying and substracting doubles
    /* 
        Point&    operator+=( double64 );    
        Point&    operator-=( double64 );    
        Point&    operator*=( double64 );    
        Point&    operator/=( double64 );    
    */
    p1 += 10.;
    p1 -= 5.;
    p1 -= 5.;
    _test( p1==p3 );
    p2 *= 4.;
    p2 /= 2.;
    p2 /= 2.;
    _test( p1==p2 );
    
    // using temporaries
    /*
        Point     operator+( const Point& ) const;    
        Point     operator-( const Point& ) const;    
        Point     operator*( const Point& ) const;    
        Point     operator/( const Point& ) const;  
    */  
    Point<dim> res = p1 + p3;
    p2 *= 2.;
    _test( res==p2 );
    
    res = p1 - p3;
    Point<dim>  p0; 
    p0 = 0.; // assignment to double
    _test( res==p0 );
   
    p2 = p1 * p3;
    res = p2 / p3;
    _test( res==p1 );
    p2 = p2 / p3; // involving point itself
    _test( p2==p1 );

    // double64 post multiplicator
    /*
        Point     operator*( double64 ) const;    
        Point     operator/( double64 ) const;    
    */
    res = p1  * 64.5;
    p2  = res / 64.5;
    _test( p2==p1 );

    // assignments 
    p1[0] = p1[1] = p1[2];  
    _test( !(p1==p3) );
    
    res = p1 * 5.;
    // operators <, >
    _test( p1  < res );
    _test( res > p1 );

    p1.Set( a );
    p2.Set( a );
    p3.Set( a );
    _test( p1 == p2 );
    res[0]=1.; res[1]=1.; res[2]=0.;
    _equal( res.Length(), sqrt(2.), numeric_limits<double64>::epsilon() );
    res.NormalizeLengthTo();
    _equal( res.Length(), 1., numeric_limits<double64>::epsilon() );    
    res[0]=1.; res[1]=1.; res[2]=0.;
    _equal( res.DistanceTo(p3), sqrt(10.), numeric_limits<double64>::epsilon() );
    _test( !res.CoincidesWithWithinTolerance(p3) );
    vector<double64> b(p1.Coordinates());
    _test( a == b );
    p1.Out();
    cout <<"Point_Test::Test_3D_Point: ostream test: "<< p1 << endl << endl;
    
    // related non-member binary operators
    res = 1. * p1;
    _test( res == p1 );
    res = 1. + p1;
    p2  = p1 + 1.;
    _test( res == p2 );
    p1.Set( a );
    res = 3. - p1;
    _equal( res[0], 2., numeric_limits<double64>::epsilon() );
    _equal( res[1], 1., numeric_limits<double64>::epsilon() );
    _equal( res[2], 0., numeric_limits<double64>::epsilon() );
    
    // mid point
    p1[0]=0.; p1[1]=0.; p1[2]=0.;
    p2[0]=1.; p2[1]=1.; p2[2]=1.;
    p3[0]=2.; p3[1]=2.; p3[2]=2.;
    _test( p2==midPoint(p1,p3) );
    Point<dim> pMP = midPoint(p1,p3);
    _equal( p2[0], pMP[0], numeric_limits<double64>::epsilon() );
    _equal( p2[1], pMP[1], numeric_limits<double64>::epsilon() );
    _equal( p2[2], pMP[2], numeric_limits<double64>::epsilon() );

    // dot product
    _equal( dotProduct( p2, p3 ), 6., numeric_limits<double64>::epsilon() );

    // cross product 
    p1[0]=1.; p1[1]=0.; p1[2]=0.;
    p2[0]=0.; p2[1]=1.; p2[2]=0.;
    res = crossProduct( p1, p2 );
    _equal( res[0], 0., numeric_limits<double64>::epsilon() );
    _equal( res[1], 0., numeric_limits<double64>::epsilon() );
    _equal( res[2], 1., numeric_limits<double64>::epsilon() );

 } // end 3D test

} // csmp
