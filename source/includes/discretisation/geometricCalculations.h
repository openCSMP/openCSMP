#ifndef CSMP_GEOMETRIC_CALCULATIONS_UOM_H
#define CSMP_GEOMETRIC_CALCULATIONS_UOM_H

#include "CSMP_definitions.h"

/**

@author R. Manasipov
@date 2015
@author SKM
@date modified in 2017

geometric calculations

*/

namespace csmp {

template<uint32_t> class Point;
template<uint32_t> class Element;

// CALCULATIONS ORIGINALLY DEVELOPED TO ASSESS THE DEGENERACY OF CORNER-POINT GRIDS
// TODO: useful functionality but needs documentation and refactoring.

/// axes ( 1D ) - creates vector from p1 to p2 and normalises it to a length of 1; vector is returned into e1
template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>&, const csmp::Point<dim>&, csmp::Point<dim>& e1 );
template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>&, csmp::Point<dim>& );

/// axes ( 2D )
template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>&, const csmp::Point<dim>&, const csmp::Point<dim>&, csmp::Point<dim>& e1, csmp::Point<dim>& e2 );
template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& e2, const csmp::Point<dim>&, csmp::Point<dim>& );

/// axes ( 3D )
template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>&, const csmp::Point<dim>&, const csmp::Point<dim>&, const csmp::Point<dim>&, csmp::Point<dim>& e1, csmp::Point<dim>& e2, csmp::Point<dim>& e3 );
template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& e2, const csmp::Point<dim>& e3, const csmp::Point<dim>&, csmp::Point<dim>& );

/// bar ( length )
double unsignedLength( const csmp::Point<1U>&, const csmp::Point<1U>& );
double unsignedLength( const csmp::Point<2U>&, const csmp::Point<2U>& );
double unsignedLength( const csmp::Point<3U>&, const csmp::Point<3U>& );

/// point in bar
bool   isPointInsideTheBar( const csmp::Point<1U>& pt, const csmp::Point<1U>&, const csmp::Point<1U>& );
bool   isPointInsideTheBar( const csmp::Point<2U>& pt, const csmp::Point<2U>&, const csmp::Point<2U>& );
bool   isPointInsideTheBar( const csmp::Point<3U>& pt, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// triangle ( area )
double unsignedArea( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
double unsignedArea( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
double unsignedArea( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

double signedArea( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
double signedArea( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
double signedArea( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// triangle ( orientation )
bool   isCounterClockWiseOrientation( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
bool   isCounterClockWiseOrientation( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
bool   isCounterClockWiseOrientation( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// point in triangle
bool   isPointInsideTheTriangle( const csmp::Point<1U>& pt, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
bool   isPointInsideTheTriangle( const csmp::Point<2U>& pt, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
bool   isPointInsideTheTriangle( const csmp::Point<3U>& pt, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// point with respect to line
bool   isSameSide( const csmp::Point<2U>& pt, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
bool   isSameSide( const csmp::Point<3U>& pt, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );


/// quad ( area )
double unsignedArea( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
double unsignedArea( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
double unsignedArea( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

double signedArea( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
double signedArea( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
double signedArea( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// point in quadrilateral


/// tetra ( volume )
double unsignedVolume( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );
double signedVolume( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// tetra orientation
bool   isCounterClockWiseOrientation( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
bool   isCounterClockWiseOrientation( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
bool   isCounterClockWiseOrientation( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// dihedral angle
double dihedralDegAngle( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
double dihedralDegAngle( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
double dihedralDegAngle( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

double dihedralRadAngle( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>& );
double dihedralRadAngle( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>& );
double dihedralRadAngle( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// tetra definition - checks whether provided 4 points define a space filling tetrahedron or quadrilatera (only 3D)l
bool   isTetra( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>& );

/// if the order of points is not counter-clockwise in right-hand rule coordinate system, the correct order of vertices is returned into the supplied map
bool   isTetra( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&,
                std::map<size_t,size_t>& order );

/// using equiangular skewness, evaluate whether the elements are fit for computation
bool isTetrahedron( const std::vector<Point<3U> >& vertexList );
bool isHexahedron( const std::vector<Point<3U> >& );
bool isPrism( const std::vector<Point<3U> >& );
bool isPyramid( const std::vector<Point<3U> >& );

bool isValidElement( const std::vector<Point<3U> >& vertexList );

/// using the above functionality, tests whether element is fit for computations
template<uint32_t dim>
bool isValidElement( const Element<dim>* const );

}// end namespace csmp

#endif

