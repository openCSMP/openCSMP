/*
 *  Point_Test.cpp
 *
 *  Created by Stephan Matthai on 10/26/10.
 */

#include <limits>
#include <random>

#include "Point.h"
#include "Point_Test.h"
#include "compareFloats.h"

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
    TestPointComparitors2D();
    TestPointComparitors3D();
 }




void Point_Test::Test_1D_Point()
 {
   const uint32_t dim(1U);

   // constructing point from vector
   vector<double>  a(1U,1.); 
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
    // also testing:     double  operator[](size_t) const;
    _equal( p1[0], p2[0], numeric_limits<double>::epsilon() );
    _equal( p1[1], p2[1], numeric_limits<double>::epsilon() );
    _test( p1==p2 ); // operator==( const Point& )
    p1 *= p3; 
    p1 /= p3;  
    _test( p1==p3 );

    // adding, multiplying and substracting doubles
    /* 
        Point&    operator+=( double );    
        Point&    operator-=( double );    
        Point&    operator*=( double );    
        Point&    operator/=( double );    
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

    // double post multiplicator
    /*
        Point     operator*( double ) const;    
        Point     operator/( double ) const;    
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
    _equal( res.Length(), 1., numeric_limits<double>::epsilon() );
    res.NormalizeLengthTo();
    _equal( res.Length(), 1., numeric_limits<double>::epsilon() );    
    res[0]=1.; 
    p3[0]=3.;
    _equal( res.DistanceTo(p3), 2., numeric_limits<double>::epsilon() );
    _test( !res.CoincidesWithWithinTolerance(p3) );
    vector<double> b(p1.Coordinates());
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
    _equal( res[0], 2., numeric_limits<double>::epsilon() );
    
    // mid point
    p1[0]=0.; 
    p2[0]=1.; 
    p3[0]=2.; 
    _test( p2==midPoint(p1,p3) );
    Point<dim> pMP = midPoint(p1,p3);
    _equal( p2[0], pMP[0], numeric_limits<double>::epsilon() );

    // dot product
    _equal( dotProduct( p2, p3 ), 2., numeric_limits<double>::epsilon() );

    // cross product 
    p1[0]=1.; 
    p2[0]=4.;
    res = crossProduct( p1, p2 );
    _equal( res[0], 4., numeric_limits<double>::epsilon() );

 } // end 1D test 
 


void Point_Test::Test_2D_Point()
 {
   const uint32_t dim(2U);

   // constructing point from vector
   vector<double>  a(2U,1.); a[1]=2.;
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
    // also testing:     double  operator[](size_t) const;
    _equal( p1[0], p2[0], numeric_limits<double>::epsilon() );
    _equal( p1[1], p2[1], numeric_limits<double>::epsilon() );
    _test( p1==p2 ); // operator==( const Point& )
    p1 *= p3; 
    p1 /= p3;  
    _test( p1==p3 );

    // adding, multiplying and substracting doubles
    /* 
        Point&    operator+=( double );    
        Point&    operator-=( double );    
        Point&    operator*=( double );    
        Point&    operator/=( double );    
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

    // double post multiplicator
    /*
        Point     operator*( double ) const;    
        Point     operator/( double ) const;    
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
    _equal( res.Length(), sqrt(2.), numeric_limits<double>::epsilon() );
    res.NormalizeLengthTo();
    _equal( res.Length(), 1., numeric_limits<double>::epsilon() );    
    res[0]=1.; res[1]=1.;
    _equal( res.DistanceTo(p3), 1., numeric_limits<double>::epsilon() );
    _test( !res.CoincidesWithWithinTolerance(p3) );
    vector<double> b(p1.Coordinates());
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
    _equal( res[0], 2., numeric_limits<double>::epsilon() );
    _equal( res[1], 1., numeric_limits<double>::epsilon() );
    
    // mid point
    p1[0]=0.; p1[1]=0.;
    p2[0]=1.; p2[1]=1.;
    p3[0]=2.; p3[1]=2.;
    _test( p2==midPoint(p1,p3) );
    Point<dim> pMP = midPoint(p1,p3);
    _equal( p2[0], pMP[0], numeric_limits<double>::epsilon() );
    _equal( p2[1], pMP[1], numeric_limits<double>::epsilon() );

    // dot product
    _equal( dotProduct( p2, p3 ), 4., numeric_limits<double>::epsilon() );

    // cross product 
    p1[0]=1.; p1[1]=0.;
    p2[0]=0.; p2[1]=1.;
    res = crossProduct( p1, p2 );
    _equal( res[0], 1., numeric_limits<double>::epsilon() );
    
 } // end 2D test

 



void Point_Test::Test_3D_Point()
 {
   const uint32_t dim(3U);

   // constructing point from vector
   vector<double>  a(3U,1.); a[1]=2.; a[2]=3.;
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
    // also testing:     double  operator[](size_t) const;
    _equal( p1[0], p2[0], numeric_limits<double>::epsilon() );
    _equal( p1[1], p2[1], numeric_limits<double>::epsilon() );
    _equal( p1[2], p2[2], numeric_limits<double>::epsilon() );
    _test( p1==p2 ); // operator==( const Point& )
    p1 *= p3; 
    p1 /= p3;  
    _test( p1==p3 );

    // adding, multiplying and substracting doubles
    /* 
        Point&    operator+=( double );    
        Point&    operator-=( double );    
        Point&    operator*=( double );    
        Point&    operator/=( double );    
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

    // double post multiplicator
    /*
        Point     operator*( double ) const;    
        Point     operator/( double ) const;    
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
    _equal( res.Length(), sqrt(2.), numeric_limits<double>::epsilon() );
    res.NormalizeLengthTo();
    _equal( res.Length(), 1., numeric_limits<double>::epsilon() );    
    res[0]=1.; res[1]=1.; res[2]=0.;
    _equal( res.DistanceTo(p3), sqrt(10.), numeric_limits<double>::epsilon() );
    _test( !res.CoincidesWithWithinTolerance(p3) );
    vector<double> b(p1.Coordinates());
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
    _equal( res[0], 2., numeric_limits<double>::epsilon() );
    _equal( res[1], 1., numeric_limits<double>::epsilon() );
    _equal( res[2], 0., numeric_limits<double>::epsilon() );
    
    // mid point
    p1[0]=0.; p1[1]=0.; p1[2]=0.;
    p2[0]=1.; p2[1]=1.; p2[2]=1.;
    p3[0]=2.; p3[1]=2.; p3[2]=2.;
    _test( p2==midPoint(p1,p3) );
    Point<dim> pMP = midPoint(p1,p3);
    _equal( p2[0], pMP[0], numeric_limits<double>::epsilon() );
    _equal( p2[1], pMP[1], numeric_limits<double>::epsilon() );
    _equal( p2[2], pMP[2], numeric_limits<double>::epsilon() );

    // dot product
    _equal( dotProduct( p2, p3 ), 6., numeric_limits<double>::epsilon() );

    // cross product 
    p1[0]=1.; p1[1]=0.; p1[2]=0.;
    p2[0]=0.; p2[1]=1.; p2[2]=0.;
    res = crossProduct( p1, p2 );
    _equal( res[0], 0., numeric_limits<double>::epsilon() );
    _equal( res[1], 0., numeric_limits<double>::epsilon() );
    _equal( res[2], 1., numeric_limits<double>::epsilon() );

    // testing distanceFromLine
    {
      Point<3> A(1.,1.,0.), B(0.,0.,0.), C(2.,0.,0.);
      double dist = distanceFromLine( A, B, C );
      _test( approximatelyEqual(dist,1.) );
    }
    {
      // beyond beginning
      Point<3> A(3.,0.,0.), B(0.,0.,0.), C(2.,0.,0.);
      double dist = distanceFromLine( A, B, C );
      _test( approximatelyEqual(dist,distance(A,C)) );
    }
    {
      // beyond end-point
      Point<3> A(-1.,0.,0.), B(0.,0.,0.), C(2.,0.,0.);
      double dist = distanceFromLine( A, B, C );
      _test( approximatelyEqual(dist,distance(A,B)) );
    }
    {
      // small difference which should have no effect
      Point<3> A(1.,1.,0.), B(0.,0.,1.0e-13), C(2.,0.,-1.0e-13);
      double dist = distanceFromLine( A, B, C );
      _test( approximatelyEqual(dist,1.) );
    }
    {
      // distance vs DistanceTo
      Point<3> A(1.,1.,1.), B(3.,3.,3.0e10);
      double dist1 = A.DistanceTo( B );
      double dist2 = distance( A, B );
      _test( approximatelyEqual(dist1,dist2) );
    }




 } // end 3D test



void Point_Test::TestPointComparitors2D()
{
   const uint32_t dim{2U};
   Point<dim> p1(0.,0.), p2(1.,1.), p2b(1.,1.), p3(1.1,1.), p4(1.1,1.);
   // basic
   _test( p1 < p2 );
   _test( p2 > p1 );
   // different by one of the values
   _test( p2 < p3 );
   _test( p2 < p4 );
   _test( p3 > p2 );
   _test( p4 > p2 );
   // if points are the same comparison should fail
   _test( !(p2 < p2b) );
   _test( !(p2b > p2) );
   _test( !(p2 < p2) );
   _test( !(p2 > p2) );
   _test( !(p1 < p1) );
   _test( !(p1 > p1) );

} // end TestPointComparitors2D




void Point_Test::TestPointComparitors3D()
{
   const uint32_t dim{3U};
   Point<dim> p1(0.,0.,0.), p2(1.,1.,1.), p2b(1.,1.,1.), p3(1.,1.,1.1), p4(1.1,1.,1.), p5(1.,1.1,1.);
   // basic
   _test( p1 < p2 );
   _test( p2 > p1 );
   // different by one of the values
   _test( p2 < p3 );
   _test( p2 < p4 );
   _test( p2 < p5 );
   _test( p3 > p2 );
   _test( p4 > p2 );
   _test( p5 > p2 );
   // if points are the same comparison should fail
   _test( !(p2 < p2b) );
   _test( !(p2b > p2) );
   _test( !(p2 < p2) );
   _test( !(p2 > p2) );
   _test( !(p1 < p1) );
   _test( !(p1 > p1) );
   
   // TODO: add more challenging comparisons involving very small discrepancies
   // ---------------------------------------------------------------------------
   // generating point locations with random number generator
   random_device rd; // obtain a random number from hardware
   mt19937 gen(rd()); // seed the generator
   default_random_engine generator;
   // define the range
   uniform_real_distribution<double> double_distr( 1., 1e15 );
   lognormal_distribution<double> log_distr( 1., 10. );

    for(int n=0; n<4; ++n)
        std::cout << double_distr(gen) << ' '; // generate numbers
    cout << endl;

} // end TestPointComparitors3D



} // csmp
