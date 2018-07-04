#include "GeometricCalculations.h"
#include "CSMP_mathUtilities.h"
#include "ErrorHandler.h"

namespace csmp {

// AXES

/// axes ( 1D )
template<size_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1, const csmp::Point<dim>& pt2, csmp::Point<dim>& e1 )
{
    e1  = pt2;
    e1 -= pt1;
    e1 /= e1.Length();
    return;
}

template void getCartesianAxes( const csmp::Point<1U>&,const csmp::Point<1U>&,csmp::Point<1U>& );
template void getCartesianAxes( const csmp::Point<2U>&,const csmp::Point<2U>&,csmp::Point<2U>& );
template void getCartesianAxes( const csmp::Point<3U>&,const csmp::Point<3U>&,csmp::Point<3U>& );

template<size_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& pt1, csmp::Point<dim>& pt2 )
{
    pt2  = pt1;
    pt2 -= pt0;
    pt2[0] = csmp::dotProduct( pt2, e1 );
    for( size_t i = 1; i < dim; ++i )
        pt2[i] = 0.0;
}

template void getCoordinate( const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,csmp::Point<1U>& );
template void getCoordinate( const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,csmp::Point<2U>& );
template void getCoordinate( const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,csmp::Point<3U>& );


/// axes ( 2D )
template<size_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1, const csmp::Point<dim>& pt2, const csmp::Point<dim>& pt3, csmp::Point<dim>& e1, csmp::Point<dim>& e2 )
{
    /// Gram-Schmidt orthogonalization
    e1 = pt2; e1 -= pt1; /// vector e1
    e1 /= e1.Length();
    e2 = pt3; e2 -= pt1; /// vector e2
    csmp::Point<dim> e2_e1_proj = e1; e2_e1_proj *= csmp::dotProduct( e2, e1 );  /// e2 proj on e1
    e2 -= e2_e1_proj;
    e2 /= e2.Length();
}

template void getCartesianAxes( const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,csmp::Point<1U>&,csmp::Point<1U>& );
template void getCartesianAxes( const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,csmp::Point<2U>&,csmp::Point<2U>& );
template void getCartesianAxes( const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,csmp::Point<3U>&,csmp::Point<3U>& );

template<size_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& e2, const csmp::Point<dim>& pt1, csmp::Point<dim>& pt2 )
{
    pt2  = pt1;
    pt2 -= pt0;
    pt2[0] = csmp::dotProduct( pt2, e1 );
    pt2[1] = csmp::dotProduct( pt2, e2 );
    for( size_t i = 2U; i < dim; ++i )
        pt2[i] = 0.0;
}

template void getCoordinate( const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,csmp::Point<1U>& );
template void getCoordinate( const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,csmp::Point<2U>& );
template void getCoordinate( const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,csmp::Point<3U>& );

/// axes ( 3D )
template<size_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1, const csmp::Point<dim>& pt2, const csmp::Point<dim>& pt3, const csmp::Point<dim>& pt4, csmp::Point<dim>& e1, csmp::Point<dim>& e2, csmp::Point<dim>& e3 )
{
    /// Gram-Schmidt orthogonalization
    e1 = pt2; e1 -= pt1; /// vector e1
    e1 /= e1.Length();
    e2 = pt3; e2 -= pt1; /// vector e2
    csmp::Point<dim> e2_e1_proj = e1; e2_e1_proj *= csmp::dotProduct( e2, e1 );  /// e2 proj on e1
    e2 -= e2_e1_proj;
    e2 /= e2.Length();
    e3 = pt4; e3 -= pt1; /// vector e3
    csmp::Point<dim> e3_e1_proj = e1; e3_e1_proj *= csmp::dotProduct( e3, e1 );  /// e3 proj on e1
    e3 -= e3_e1_proj;
    csmp::Point<dim> e3_e2_proj = e2; e3_e2_proj *= csmp::dotProduct( e3, e2 );  /// e3 proj on e1
    e3 -= e3_e2_proj;
    e3 /= e3.Length();
}

template void getCartesianAxes( const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,csmp::Point<1U>&,csmp::Point<1U>&,csmp::Point<1U>& );
template void getCartesianAxes( const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,csmp::Point<2U>&,csmp::Point<2U>&,csmp::Point<2U>& );
template void getCartesianAxes( const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,csmp::Point<3U>&,csmp::Point<3U>&,csmp::Point<3U>& );

template<size_t dim>
void getCoordinate( const csmp::Point<dim>& pt0, const csmp::Point<dim>& e1, const csmp::Point<dim>& e2, const csmp::Point<dim>& e3, const csmp::Point<dim>& pt1, csmp::Point<dim>& pt2 )
{
    pt2  = pt1;
    pt2 -= pt0;
    pt2[0] = csmp::dotProduct( pt2, e1 );
    pt2[1] = csmp::dotProduct( pt2, e2 );
    pt2[2] = csmp::dotProduct( pt2, e3 );
}

template void getCoordinate( const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,const csmp::Point<1U>&,csmp::Point<1U>& );
template void getCoordinate( const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,const csmp::Point<2U>&,csmp::Point<2U>& );
template void getCoordinate( const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,const csmp::Point<3U>&,csmp::Point<3U>& );


// BAR

double unsignedLength( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2 )
{
    csmp::Point<1U> e = pt2; e -= pt1;
    return e.Length();
}

double unsignedLength( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2 )
{
    csmp::Point<2U> e = pt2; e -= pt1;
    return e.Length();
}

double unsignedLength( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2 )
{
    csmp::Point<3U> e = pt2; e -= pt1;
    return e.Length();
}

bool isPointInsideTheBar( const csmp::Point<1U>& pt, const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2 )
{
    if( ( pt[0] - pt1[0] )*( pt2[0] - pt[0] ) < 0.0 )
        return false;
    return true;
}

bool isPointInsideTheBar( const csmp::Point<2U>& pt, const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2 )
{
    if( unsignedArea(pt,pt1,pt2) > 0.0 )
        return false;
    csmp::Point<2U> e0 = pt2; e0 -= pt1;
    csmp::Point<2U> e1 = pt;  e1 -= pt1;
    const double len   = e0.Length();
    const double proj  = csmp::dotProduct( e0, e1 );
    if( ( proj < 0.0 ) || ( proj > len ) )
            return false;
    return true;
}

bool isPointInsideTheBar( const csmp::Point<3U>& pt, const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2 )
{
    if( unsignedArea(pt,pt1,pt2) > 0.0 )
        return false;
    csmp::Point<3U> e0 = pt2; e0 -= pt1;
    csmp::Point<3U> e1 = pt;  e1 -= pt1;
    const double len   = e0.Length();
    const double proj  = csmp::dotProduct( e0, e1 );
    if( ( proj < 0.0 ) || ( proj > len ) )
            return false;
    return true;

}

// TRIANGLES

double unsignedArea( const csmp::Point<1U>& pt1,
                     const csmp::Point<1U>& pt2,
                     const csmp::Point<1U>& pt3 )
{
    return 0.0;
}

double unsignedArea( const csmp::Point<2U>& pt1,
                     const csmp::Point<2U>& pt2,
                     const csmp::Point<2U>& pt3 )
{
    return std::abs( signedArea( pt1,pt2,pt3 ) );
}

double unsignedArea( const csmp::Point<3U>& pt1,
                     const csmp::Point<3U>& pt2,
                     const csmp::Point<3U>& pt3 )
{
    csmp::Point<3U> avec = pt2; avec -= pt1; /// vector a
    csmp::Point<3U> bvec = pt3; avec -= pt1; /// vector b
    return 0.5*csmp::crossProduct( avec, bvec ).Length();
}

/**
         | pt1[0] pt1[1] 1 |
  A = det| pt2[0] pt2[1] 1 |
         | pt3[0] pt3[1] 1 |
*/

double signedArea( const csmp::Point<1U>& pt1,
                   const csmp::Point<1U>& pt2,
                   const csmp::Point<1U>& pt3 )
{
    return 0.0;
}

double signedArea( const csmp::Point<2U>& pt1,
                   const csmp::Point<2U>& pt2,
                   const csmp::Point<2U>& pt3 )
{
    /// signed area
    double sarea = 0.0;
    sarea += (pt1[0U]-pt3[0U])*(pt2[1U]-pt3[1U]);
    sarea -= (pt2[0U]-pt3[0U])*(pt1[1U]-pt3[1U]);
    sarea *= 0.5;
    return sarea;
}

double signedArea( const csmp::Point<3U>& pt1,
                   const csmp::Point<3U>& pt2,
                   const csmp::Point<3U>& pt3 )
{
    return unsignedArea(pt1,pt2,pt3);
/*
    /// construct new coordinate system (i,j)
    csmp::Point<3U> ivec = pt2; ivec -= pt1; /// vector a, later i
    csmp::Point<3U> bvec = pt3; bvec -= pt1; /// vector b
    csmp::Point<3U> jvec = bvec;             /// vector j
    /// normalize vector a, get i direction
    const double ai_proj( ivec.Length() );
    ivec /= ai_proj;
    /// find projection of vector b on i direction
    const double bi_proj( csmp::dotProduct( bvec, ivec ) );
    /// find projection of vector b on j direction
    jvec -= ivec*bi_proj;
    jvec /= jvec.Length();
    const double bj_proj( csmp::dotProduct( bvec, jvec ) );
    /// calculate signed area
    return 0.5*ai_proj*bj_proj;
*/
}

bool isCounterClockWiseOrientation( const csmp::Point<1U>& pt1,
                                    const csmp::Point<1U>& pt2,
                                    const csmp::Point<1U>& pt3 )
{
    return true;
}

bool isCounterClockWiseOrientation( const csmp::Point<2U>& pt1,
                                    const csmp::Point<2U>& pt2,
                                    const csmp::Point<2U>& pt3 )
{
    /// signed area
    double sarea = 0.0;
    sarea += (pt1[0U]-pt3[0U])*(pt2[1U]-pt3[1U]);
    sarea -= (pt2[0U]-pt3[0U])*(pt1[1U]-pt3[1U]);
    return ( sarea > 0.0 );
}

bool isCounterClockWiseOrientation( const csmp::Point<3U>& pt1,
                                    const csmp::Point<3U>& pt2,
                                    const csmp::Point<3U>& pt3 )
{
    return true;
}

bool isPointInsideTheTriangle( const csmp::Point<1U>& pt,
                    const csmp::Point<1U>& pt1,
                    const csmp::Point<1U>& pt2,
                    const csmp::Point<1U>& pt3 )
{
    double min_pt = pt1[0];
    min_pt = std::min(min_pt,pt2[0]);
    min_pt = std::min(min_pt,pt3[0]);
    double max_pt = pt1[0];
    max_pt = std::max(max_pt,pt2[0]);
    max_pt = std::max(max_pt,pt3[0]);
    if( ( pt[0] < min_pt ) || ( pt[0] > max_pt ) )
        return false;
    return true;
}

bool isPointInsideTheTriangle( const csmp::Point<2U>& pt,
                               const csmp::Point<2U>& pt1,
                               const csmp::Point<2U>& pt2,
                               const csmp::Point<2U>& pt3 )
{
    const double area1 = signedArea( pt, pt1, pt2 );
    if( area1 == 0.0 )
        return true;
    const bool side1 = ( area1 > 0.0 );
    const double area2 = signedArea( pt, pt2, pt3 );
    if( area2 == 0.0 )
        return true;
    const bool side2 = ( area2 > 0.0 );
    if( side1 != side2 )
        return false;
    const double area3 = signedArea( pt, pt3, pt1 );
    if( area3 == 0.0 )
        return true;
    const bool side3 = ( area3 > 0.0 );
    if( side2 != side3 )
        return false;
    return true;
}

bool isPointInsideTheTriangle( const csmp::Point<3U>& pt,
                               const csmp::Point<3U>& pt1,
                               const csmp::Point<3U>& pt2,
                               const csmp::Point<3U>& pt3 )
{
    /// double area of triangle
    csmp::Point<3U> e12 = pt2; e12 -= pt1;
    csmp::Point<3U> e23 = pt3; e23 -= pt2;
    csmp::Point<3U> tri_area_vec = csmp::crossProduct( e12, e23 );

    /// process edge12:
    csmp::Point<3U> e10 = pt;  e10 -= pt1;
    csmp::Point<3U> cp120 = csmp::crossProduct( e12, e10 );

    const double area1 = csmp::dotProduct( cp120, tri_area_vec );
    if( area1 == 0.0 )
        return true;
    const bool side1 = ( area1 > 0.0 );

    /// process edge23:
    csmp::Point<3U> e20 = pt;  e20 -= pt2;
    csmp::Point<3U> cp230 = csmp::crossProduct( e23, e20 );

    const double area2 = csmp::dotProduct( cp230, tri_area_vec );
    if( area2 == 0.0 )
        return true;
    const bool side2 = ( area2 > 0.0 );
    if( side1 != side2 )
        return false;

    /// process edge31:
    csmp::Point<3U> e30 = pt;  e20 -= pt3;
    csmp::Point<3U> e31 = pt1; tri_area_vec -= pt3;
    csmp::Point<3U> cp310 = csmp::crossProduct( e31, e30 );

    const double area3 = csmp::dotProduct( cp310, tri_area_vec );
    if( area3 == 0.0 )
        return true;
    const bool side3 = ( area3 > 0.0 );
    if( side2 != side3 )
        return false;
    return true;
}

/**
  isSameSide() defines whether pt1 is on the same side as pt2 with respect to (pt3,pt4)

           pt1
          x
         /
        /
   pt3 x-------x pt4
        \
         \
          x
           pt2

 */

bool isSameSide( const csmp::Point<1U>& pt1,
                 const csmp::Point<1U>& pt2,
                 const csmp::Point<1U>& pt3,
                 const csmp::Point<1U>& pt4 )
{
    return true;
}

bool isSameSide( const csmp::Point<2U>& pt1,
                 const csmp::Point<2U>& pt2,
                 const csmp::Point<2U>& pt3,
                 const csmp::Point<2U>& pt4 )
{
    const double area1 = signedArea( pt1, pt3, pt4 );
    if( area1 == 0.0 )
        return true;
    const bool side1 = ( area1 > 0.0 );
    const double area2 = signedArea( pt2, pt3, pt4 );
    if( area2 == 0.0 )
        return true;
    const bool side2 = ( area2 > 0.0 );
    if( side1 != side2 )
        return false;
    return true;
}

bool isSameSide( const csmp::Point<3U>& pt1,
                 const csmp::Point<3U>& pt2,
                 const csmp::Point<3U>& pt3,
                 const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> e0 = pt4; e0 -= pt3;
    csmp::Point<3U> e1 = pt1; e1 -= pt3;
    csmp::Point<3U> e2 = pt2; e2 -= pt3;
    csmp::Point<3U> e01 = csmp::crossProduct( e0, e1 );
    csmp::Point<3U> e02 = csmp::crossProduct( e0, e2 );
    if( csmp::dotProduct( e01, e02 ) < 0.0)
        return false;
    return true;
}



// TETRAHEDRON AND QUADRILATERAL

double unsignedArea( const csmp::Point<1U>& pt1,
                     const csmp::Point<1U>& pt2,
                     const csmp::Point<1U>& pt3,
                     const csmp::Point<1U>& pt4 )
{
    return 0.0;
}

double unsignedArea( const csmp::Point<2U>& pt1,
                     const csmp::Point<2U>& pt2,
                     const csmp::Point<2U>& pt3,
                     const csmp::Point<2U>& pt4 )
{
    double area = 0.0;
    area += unsignedArea(pt1,pt2,pt3);
    area += unsignedArea(pt1,pt3,pt4);
    return area;
}

double unsignedArea( const csmp::Point<3U>& pt1,
                     const csmp::Point<3U>& pt2,
                     const csmp::Point<3U>& pt3,
                     const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;
    double area = 0.0;
    area += csmp::crossProduct( a, b ).Length();
    area += csmp::crossProduct( b, c ).Length();
    area *= 0.5;
    return area;
}

double signedArea( const csmp::Point<1U>& pt1,
                   const csmp::Point<1U>& pt2,
                   const csmp::Point<1U>& pt3,
                   const csmp::Point<1U>& pt4 )
{
    return 0.0;
}

double signedArea( const csmp::Point<2U>& pt1,
                   const csmp::Point<2U>& pt2,
                   const csmp::Point<2U>& pt3,
                   const csmp::Point<2U>& pt4 )
{
    double area = 0.0;
    area += signedArea(pt1,pt2,pt3);
    area += signedArea(pt1,pt3,pt4);
    return area;
}

double signedArea( const csmp::Point<3U>& pt1,
                   const csmp::Point<3U>& pt2,
                   const csmp::Point<3U>& pt3,
                   const csmp::Point<3U>& pt4 )
{
    double area = 0.0;
    area += signedArea(pt1,pt2,pt3);
    area += signedArea(pt1,pt3,pt4);
    return area;
}

double unsignedVolume( const csmp::Point<1U>& pt1,
                       const csmp::Point<1U>& pt2,
                       const csmp::Point<1U>& pt3,
                       const csmp::Point<1U>& pt4 )
{
    return 0.0;
}

double unsignedVolume( const csmp::Point<2U>& pt1,
                       const csmp::Point<2U>& pt2,
                       const csmp::Point<2U>& pt3,
                       const csmp::Point<2U>& pt4 )
{
    return 0.0;
}

double unsignedVolume( const csmp::Point<3U>& pt1,
                       const csmp::Point<3U>& pt2,
                       const csmp::Point<3U>& pt3,
                       const csmp::Point<3U>& pt4 )
{
    return std::abs( unsignedVolume( pt1,pt2, pt3, pt4 ) );
}

double signedVolume( const csmp::Point<1U>& pt1,
                     const csmp::Point<1U>& pt2,
                     const csmp::Point<1U>& pt3,
                     const csmp::Point<1U>& pt4 )
{
    return 0.0;
}

double signedVolume( const csmp::Point<2U>& pt1,
                     const csmp::Point<2U>& pt2,
                     const csmp::Point<2U>& pt3,
                     const csmp::Point<2U>& pt4 )
{
    return 0.0;
}

double signedVolume( const csmp::Point<3U>& pt1,
                     const csmp::Point<3U>& pt2,
                     const csmp::Point<3U>& pt3,
                     const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;
    csmp::Point<3U> d = csmp::crossProduct( b, c );
    return  csmp::dotProduct( a, d );
}

bool isCounterClockWiseOrientation( const csmp::Point<1U>& pt1,
                                    const csmp::Point<1U>& pt2,
                                    const csmp::Point<1U>& pt3,
                                    const csmp::Point<1U>& pt4 )
{
    return true;
}

bool isCounterClockWiseOrientation( const csmp::Point<2U>& pt1,
                                    const csmp::Point<2U>& pt2,
                                    const csmp::Point<2U>& pt3,
                                    const csmp::Point<2U>& pt4 )
{
    return isCounterClockWiseOrientation( pt2, pt3, pt4 );
}

bool isCounterClockWiseOrientation( const csmp::Point<3U>& pt1,
                                    const csmp::Point<3U>& pt2,
                                    const csmp::Point<3U>& pt3,
                                    const csmp::Point<3U>& pt4 )
{
    return ( signedVolume( pt1,pt2,pt3,pt4 ) > 0.0 );
}



/** dihedral angle between planes ( pt1,pt3,pt4 ) and ( pt2,pt3,pt4 )

             pt1
              x
             /
            /
    pt3 x-------x pt4
            \
             \
              x
               pt2

*/
double dihedralDegAngle( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 )
{
    return 0.0;
}
double dihedralRadAngle( const csmp::Point<1U>& pt1, const csmp::Point<1U>& pt2, const csmp::Point<1U>& pt3, const csmp::Point<1U>& pt4 )
{
    return 0.0;
}
double dihedralDegAngle( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 )
{
    if( isSameSide( pt1, pt2, pt3, pt4 ) )
        return 0.0;
    return 180.;
}
double dihedralRadAngle( const csmp::Point<2U>& pt1, const csmp::Point<2U>& pt2, const csmp::Point<2U>& pt3, const csmp::Point<2U>& pt4 )
{
    if( isSameSide( pt1, pt2, pt3, pt4 ) )
        return 0.0;
    return PI;
}
double dihedralDegAngle( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 )
{
    /// directional vectors
    csmp::Point<3U> a = pt4; a -= pt3;
    csmp::Point<3U> b = pt1; b -= pt3;
    csmp::Point<3U> c = pt2; c -= pt3;

    /// normals
    csmp::Point<3U> n1 = csmp::crossProduct( b, a );
    n1 /= n1.Length();
    csmp::Point<3U> n2 = csmp::crossProduct( c, a );
    n2 /= n2.Length();

    double angle = csmp::dotProduct( n1, n2 );
    if( angle <= -1.0 )
        return 180.0;
    if( angle >= 1.0 )
        return 0.0;
    return ( std::acos( angle )*180.0 )/csmp::PI;
}
double dihedralRadAngle( const csmp::Point<3U>& pt1, const csmp::Point<3U>& pt2, const csmp::Point<3U>& pt3, const csmp::Point<3U>& pt4 )
{
    /// directional vectors
    csmp::Point<3U> a = pt4; a -= pt3;
    csmp::Point<3U> b = pt1; b -= pt3;
    csmp::Point<3U> c = pt2; c -= pt3;

    /// normals
    csmp::Point<3U> n1 = csmp::crossProduct( b, a );
    n1 /= n1.Length();
    csmp::Point<3U> n2 = csmp::crossProduct( c, a );
    n2 /= n2.Length();

    double angle = csmp::dotProduct( n1, n2 );
    if( angle <= -1.0 )
        return PI;
    if( angle >= 1.0 )
        return 0.0;
    return std::acos( angle )*180.0;
}

/**
 Check whether provided 4 point create tetrahedron or quadrilateral
*/

bool isTetra( const csmp::Point<1U>& pt1,
              const csmp::Point<1U>& pt2,
              const csmp::Point<1U>& pt3,
              const csmp::Point<1U>& pt4 )
{
    return false;
}

bool isTetra( const csmp::Point<2U>& pt1,
              const csmp::Point<2U>& pt2,
              const csmp::Point<2U>& pt3,
              const csmp::Point<2U>& pt4 )
{
    return false;
}

bool isTetra( const csmp::Point<3U>& pt1,
              const csmp::Point<3U>& pt2,
              const csmp::Point<3U>& pt3,
              const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;
    csmp::Point<3U> d = csmp::crossProduct( b, c );

    double volume =  csmp::dotProduct( a, d );
    return ( std::abs( volume ) > std::numeric_limits<double>::epsilon() );
}

/**
 Check whether provided 4 point create tetrahedron or quadrilateral.
 Return order of points which will give proper orientation of element
*/

bool isTetra( const csmp::Point<1U>& pt1,
              const csmp::Point<1U>& pt2,
              const csmp::Point<1U>& pt3,
              const csmp::Point<1U>& pt4,
              std::map<size_t,size_t>& order )
{
    return false;
}

bool isTetra( const csmp::Point<2U>& pt1,
              const csmp::Point<2U>& pt2,
              const csmp::Point<2U>& pt3,
              const csmp::Point<2U>& pt4,
              std::map<size_t,size_t>& order )
{
    return false;
}

bool isTetra( const csmp::Point<3U>& pt1,
              const csmp::Point<3U>& pt2,
              const csmp::Point<3U>& pt3,
              const csmp::Point<3U>& pt4,
              std::map<size_t,size_t>& order )
{
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;
    csmp::Point<3U> d = csmp::crossProduct( b, c );

    double volume =  csmp::dotProduct( a, d );
    if( std::abs( volume ) > std::numeric_limits<double>::epsilon() )
    {
        if( volume > 0.0 )
        {
            order[ 0 ] = 0;
            order[ 1 ] = 1;
            order[ 2 ] = 2;
            order[ 3 ] = 3;
        }
        else
        {
            order[ 0 ] = 0;
            order[ 1 ] = 1;
            order[ 2 ] = 3;
            order[ 3 ] = 2;
        }
        return true;
    }
    double area;
    double area_ab = csmp::crossProduct( a, b ).Length();
    double area_bc = csmp::crossProduct( b, c ).Length();
    double area_ac = csmp::crossProduct( a, c ).Length();
    area = area_ab + area_bc - 2.0*area_ac;
    if( std::abs( area ) < std::numeric_limits<double>::epsilon() )
    {
        order[ 0 ] = 0;
        order[ 1 ] = 1;
        order[ 2 ] = 2;
        order[ 3 ] = 3;
        return false;
    }
    area = area_ac + area_bc - 2.0*area_ab;
    if( std::abs( area ) < std::numeric_limits<double>::epsilon() )
    {
        order[ 0 ] = 0;
        order[ 1 ] = 1;
        order[ 2 ] = 3;
        order[ 3 ] = 2;
        return false;
    }
    area = area_ab + area_ac - 2.0*area_bc;
    if( std::abs( area ) < std::numeric_limits<double>::epsilon() )
    {
        order[ 0 ] = 0;
        order[ 1 ] = 2;
        order[ 2 ] = 1;
        order[ 3 ] = 3;
        return false;
    }
    return false;
}

} // end namespace csmp









