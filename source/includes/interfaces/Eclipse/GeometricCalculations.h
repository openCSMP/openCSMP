#ifndef GEOMETRIC_CALCULATIONS_H
#define GEOMETRIC_CALCULATIONS_H

#include "CSMP_definitions.h"
#include "CSMP_ElementSpecifications.h"
#include "Point.h"
#include "PolygonGrid.h"

namespace csmp {

// TODO: student code: needs refactoring and elimination of duplicate functionality, used by EclipseInterface

/// axes ( 1D )
template<size_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1, const csmp::Point<dim>& pt2, csmp::Point<dim>& e1 );
template<size_t dim>
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
bool   isPointInsideTheBar( const csmp::Point<2U>& pt, const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2 );

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


/// point in tetrahedron

/**

@author R. Manasipov
@date 2015

Geometric calculations

*/

}// end namespace csmp

#endif

