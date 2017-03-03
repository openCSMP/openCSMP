#ifndef CSMP_TRIANGULAR_FACET_H
#define CSMP_TRIANGULAR_FACET_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

/**
@file TriangularFacet.h
@addtogroup CSMPglobalFunctions
@{
*/

/// triangle area
double64 triangleArea( const Point<3U>&, 
                       const Point<3U>&, 
                       const Point<3U>& );

double64 triangleArea( const Point<2U>&,
                       const Point<2U>&,
                       const Point<2U>& );

double64 triangleArea( const Point<1U>&,
                       const Point<1U>&,
                       const Point<1U>& );

/// returns the normal of a triangle
Point<3U>  normalOfTriangle( const Point<3U>&, 
                             const Point<3U>&, 
                             const Point<3U>& );                      

Point<2U>  normalOfTriangle( const Point<2U>&, 
                             const Point<2U>&, 
                             const Point<2U>& ); 

Point<1U>  normalOfTriangle( const Point<1U>&,
                             const Point<1U>&,
                             const Point<1U>& );


/**
@}
*/

} // end namespace csmp

#endif

