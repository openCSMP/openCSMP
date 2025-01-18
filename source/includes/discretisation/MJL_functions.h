#ifndef MJL_FUNCTIONS_H
#define MJL_FUNCTIONS_H

#include "MJL_Triangle3D.h"
#include "MJL_Polygon.h"

namespace mjl {

/// @file CSMP_definitions.h
/**
@namespace mjl
@brief Mathematical Library(Lazlo)
*/
/**
  @addtogroup CSMPglobalFunctions
  @{
*/

/// convex hull around a point cloud
Polygon* insertionHull( const std::vector<mjl::Point>& ); // MJL, p. 114

Polygon* starShapedPolygon( const std::vector<mjl::Point>& ); // MJL, p. 110
    
Polygon* convexPolygonIntersect( Polygon& a, Polygon& b ); // MJL, p. 158    
    
Polygon* merge( Polygon&, Polygon& ); // MJL, p. 214;

void advance( const Polygon&, Polygon&, bool ); // MJL, p. 159
    
void bridge( Polygon&, Polygon&, Vertex*, Vertex*, int );  // MJL, p. 215  

void supportingLine( const mjl::Point&, Polygon*, LOCATION ); // MJL, p. 114

/// convex polygon method
bool pointInConvexPolygon( const mjl::Point&, Polygon& ); // MJL, p. 118

/// Ray shooting method for arbitrary polygons
POINT_CLASSIFICATION pointInPolygonRS( const mjl::Point, const Polygon& ); // MJL, p. 118

/// signed angle method for arbitrary polygons (NB: function has a tolerance issue
POINT_CLASSIFICATION pointInPolygonSA( const mjl::Point, const Polygon& ); // MJL, p.

bool firstContainsSecondRS( const Polygon&, const Polygon& ); // SKM

INTERSECTION  lineTriangle3DIntersect( const Edge3D& e, const Triangle3D& p, double& t );

Polygon*  project( const Triangle3D& p, int h, int v );

bool clipPolygonToEdge( const Polygon&, const Edge&, Polygon& ); // MJL, p. 127
 
/// variation of clipPolygonToEdge(), MJL, p. 127
bool clipEdgeToPolygon( const Polygon&, Edge& );



// inlined functions



/**

@warning polygon must have been built in counterclockwise way (positive) for
this to work !  

advances the current edge of polygon A clockwise and inserts this edge's 
(destination) endpoint x into intersection polygon R in counter-clockwise 
fashion, if x of A is inside R and x was not the last point inserted into 
R. 
*/
inline void advance( const Polygon& A, Polygon& R, bool inside )
 {
    // advancing A clockwise (polygon was built counter-clockwise) 
    A.Advance(CLOCKWISE);
    
    // add to new polygon in counter-clockwise fashion
    if ( inside && (R.Point() != A.Point()) )  R.Insert( A.Point() );
 }
/**
@}
*/


} // end namespace mjl

#endif

