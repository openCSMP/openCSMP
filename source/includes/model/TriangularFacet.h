#ifndef CSMP_TRIANGULAR_FACET_H
#define CSMP_TRIANGULAR_FACET_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

/**
@file triangularFacet.h
@addtogroup CSMPglobalFunctions
@{
*/

/// triangle area
double triangleArea( const Point<3U>&, 
                     const Point<3U>&,
                     const Point<3U>& );
/// area  in 2D
double triangleArea( const Point<2U>&,
                     const Point<2U>&,
                     const Point<2U>& );

/// returns the normal of a triangle
Point<3U>  normalOfTriangle( const Point<3U>&, 
                             const Point<3U>&, 
                             const Point<3U>& );      
                                             
/// returns zero because the normal is not captured by coordinate system
Point<2U>  normalOfTriangle( const Point<2U>&, 
                             const Point<2U>&, 
                             const Point<2U>& ); 

/**
@}
*/

} // end namespace csmp

#endif

