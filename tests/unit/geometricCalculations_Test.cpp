//
//  geometricCalculations_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 19/04/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "geometricCalculations_Test.h"
#include "geometricCalculations_UoM.h"

using namespace std;

namespace csmp {

namespace eclipse {

void geometricCalculations_Test::run()
 {
   /// axes ( 1D )
   csmp::Point<3U> pt1(1.,1.,1.), pt2(2.,4.,6.), e1;
   getCartesianAxes( pt1, pt2, e1 );
   // testing
   // unit length
   _equal( e1.Length(), 1., numeric_limits<double64>::epsilon() );
   // does the vector point in the right direction (from origin to axes end)
   pt2 -= pt1;
   _test( dotProduct(e1,pt2) > 0. );



/*
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& pt1, csmp::Point<dim>& pt2 );

/// axes ( 2D )
template<size_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1, const csmp::Point<dim>& pt2, const csmp::Point<dim>& pt3, csmp::Point<dim>& e1, csmp::Point<dim>& e2 );
template<size_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& e2, const csmp::Point<dim>& pt1, csmp::Point<dim>& pt2 );

/// axes ( 3D )
template<size_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1, const csmp::Point<dim>& pt2, const csmp::Point<dim>& pt3, const csmp::Point<dim>& pt4, csmp::Point<dim>& e1, csmp::Point<dim>& e2, csmp::Point<dim>& e3 );
template<size_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& e2, const csmp::Point<dim>& e3, const csmp::Point<dim>& pt1, csmp::Point<dim>& pt2 );

/// bar ( length )
double unsignedLength( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2 );
double unsignedLength( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2 );
double unsignedLength( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2 );

/// point in bar
bool   isPointInsideTheBar( const csmp::Point<1U>& pt, const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2 );
bool   isPointInsideTheBar( const csmp::Point<2U>& pt, const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2 );
bool   isPointInsideTheBar( const csmp::Point<3U>& pt, const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2 );

/// triangle ( area )
double unsignedArea( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3 );
double unsignedArea( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3 );
double unsignedArea( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3 );

double signedArea( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3 );
double signedArea( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3 );
double signedArea( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3 );

/// triangle ( orientation )
bool   isCounterClockWiseOrientation( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3 );
bool   isCounterClockWiseOrientation( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3 );
bool   isCounterClockWiseOrientation( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3 );

/// point in triangle
bool   isPointInsideTheTriangle( const csmp::Point<1U>& pt, const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3 );
bool   isPointInsideTheTriangle( const csmp::Point<2U>& pt, const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3 );
bool   isPointInsideTheTriangle( const csmp::Point<2U>& pt, const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3 );

/// point with respect to line
bool   isSameSide( const csmp::Point<2U>& pt, const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3 );
bool   isSameSide( const csmp::Point<3U>& pt, const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3 );


/// quad ( area )
double unsignedArea( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 );
double unsignedArea( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 );
double unsignedArea( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

double signedArea( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 );
double signedArea( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 );
double signedArea( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

/// point in quadrilateral


/// tetra ( volume )
double unsignedVolume( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );
double signedVolume( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

/// tetra orientation
bool   isCounterClockWiseOrientation( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 );
bool   isCounterClockWiseOrientation( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 );
bool   isCounterClockWiseOrientation( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

/// dihedral angle
double dihedralDegAngle( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 );
double dihedralDegAngle( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 );
double dihedralDegAngle( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

double dihedralRadAngle( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 );
double dihedralRadAngle( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 );
double dihedralRadAngle( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

/// tetra definition
bool   isTetra( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 );
bool   isTetra( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 );
bool   isTetra( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 );

bool   isTetra( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4, std::map<size_t,size_t>& order );
bool   isTetra( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4, std::map<size_t,size_t>& order );
bool   isTetra( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4, std::map<size_t,size_t>& order );
*/
 
 } // end run

} // eclipse

} // end csmp
