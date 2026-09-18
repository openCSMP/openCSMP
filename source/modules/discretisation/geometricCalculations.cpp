// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometricCalculations.h"
#include "CSMP_mathUtilities.h"

#include "Point.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>

namespace csmp {

// -----------------------------------------------------------------------
// Normal to point cloud
// -----------------------------------------------------------------------

Point<2> normalToAveragePlaneThroughPointCloud(
    const std::vector<Point<2>>& points )
{
    assert(points.size() >= 2);

    Point<2> mean{0.0, 0.0};
    for (const auto& p : points)
        mean += p;
    mean *= 1.0 / static_cast<double>(points.size());

    double Cxx = 0.0, Cxy = 0.0, Cyy = 0.0;
    for (const auto& p : points)
    {
        const Point<2> d = p - mean;
        Cxx += d[0]*d[0];
        Cxy += d[0]*d[1];
        Cyy += d[1]*d[1];
    }
    const double invN = 1.0 / static_cast<double>(points.size());
    Cxx *= invN;  Cxy *= invN;  Cyy *= invN;

    const double delta = std::sqrt((Cxx-Cyy)*(Cxx-Cyy) + 4.0*Cxy*Cxy);
    Point<2> unrml{ 2.0*Cxy, Cyy - Cxx - delta };
    unrml.NormalizeLengthTo(1.0);
    return unrml;
}

// -----------------------------------------------------------------------
// Cartesian axes — 1D
// -----------------------------------------------------------------------

template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1,
                        const csmp::Point<dim>& pt2,
                        csmp::Point<dim>& e1 )
{
    e1  = pt2;
    e1 -= pt1;
    e1 /= e1.Length();
}

template void getCartesianAxes( const csmp::Point<1U>&, const csmp::Point<1U>&, csmp::Point<1U>& );
template void getCartesianAxes( const csmp::Point<2U>&, const csmp::Point<2U>&, csmp::Point<2U>& );
template void getCartesianAxes( const csmp::Point<3U>&, const csmp::Point<3U>&, csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Cartesian axes — 2D (Gram-Schmidt)
// -----------------------------------------------------------------------

template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1,
                        const csmp::Point<dim>& pt2,
                        const csmp::Point<dim>& pt3,
                        csmp::Point<dim>& e1,
                        csmp::Point<dim>& e2 )
{
    e1  = pt2; e1 -= pt1;
    e1 /= e1.Length();

    e2  = pt3; e2 -= pt1;
    csmp::Point<dim> proj = e1;
    proj *= csmp::dotProduct(e2, e1);
    e2 -= proj;
    e2 /= e2.Length();
}

template void getCartesianAxes( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, csmp::Point<1U>&, csmp::Point<1U>& );
template void getCartesianAxes( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, csmp::Point<2U>&, csmp::Point<2U>& );
template void getCartesianAxes( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, csmp::Point<3U>&, csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Cartesian axes — 3D (Gram-Schmidt)
// -----------------------------------------------------------------------

template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1,
                        const csmp::Point<dim>& pt2,
                        const csmp::Point<dim>& pt3,
                        const csmp::Point<dim>& pt4,
                        csmp::Point<dim>& e1,
                        csmp::Point<dim>& e2,
                        csmp::Point<dim>& e3 )
{
    e1  = pt2; e1 -= pt1;
    e1 /= e1.Length();

    e2  = pt3; e2 -= pt1;
    { csmp::Point<dim> p = e1; p *= csmp::dotProduct(e2,e1); e2 -= p; }
    e2 /= e2.Length();

    e3  = pt4; e3 -= pt1;
    { csmp::Point<dim> p = e1; p *= csmp::dotProduct(e3,e1); e3 -= p; }
    { csmp::Point<dim> p = e2; p *= csmp::dotProduct(e3,e2); e3 -= p; }
    e3 /= e3.Length();
}

template void getCartesianAxes( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, csmp::Point<1U>&, csmp::Point<1U>&, csmp::Point<1U>& );
template void getCartesianAxes( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, csmp::Point<2U>&, csmp::Point<2U>&, csmp::Point<2U>& );
template void getCartesianAxes( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, csmp::Point<3U>&, csmp::Point<3U>&, csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Coordinate projection — 1D
// -----------------------------------------------------------------------

template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0,
                     const csmp::Point<dim>& e1,
                     const csmp::Point<dim>& pt1,
                     csmp::Point<dim>& pt2 )
{
    pt2    = pt1;
    pt2   -= pt0;
    pt2[0] = csmp::dotProduct(pt2, e1);
    for (uint32_t i=1; i<dim; ++i) pt2[i] = 0.0;
}

template void getCoordinate( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, csmp::Point<1U>& );
template void getCoordinate( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, csmp::Point<2U>& );
template void getCoordinate( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Coordinate projection — 2D
// -----------------------------------------------------------------------

template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0,
                     const csmp::Point<dim>& e1,
                     const csmp::Point<dim>& e2,
                     const csmp::Point<dim>& pt1,
                     csmp::Point<dim>& pt2 )
{
    pt2    = pt1;
    pt2   -= pt0;
    pt2[0] = csmp::dotProduct(pt2, e1);
    pt2[1] = csmp::dotProduct(pt2, e2);
    for (uint32_t i=2; i<dim; ++i) pt2[i] = 0.0;
}

template void getCoordinate( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, csmp::Point<1U>& );
template void getCoordinate( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, csmp::Point<2U>& );
template void getCoordinate( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Coordinate projection — 3D
// -----------------------------------------------------------------------

template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0,
                     const csmp::Point<dim>& e1,
                     const csmp::Point<dim>& e2,
                     const csmp::Point<dim>& e3,
                     const csmp::Point<dim>& pt1,
                     csmp::Point<dim>& pt2 )
{
    pt2    = pt1;
    pt2   -= pt0;
    pt2[0] = csmp::dotProduct(pt2, e1);
    pt2[1] = csmp::dotProduct(pt2, e2);
    pt2[2] = csmp::dotProduct(pt2, e3);
}

template void getCoordinate( const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, const csmp::Point<1U>&, csmp::Point<1U>& );
template void getCoordinate( const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, const csmp::Point<2U>&, csmp::Point<2U>& );
template void getCoordinate( const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, const csmp::Point<3U>&, csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Bar
// -----------------------------------------------------------------------

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

bool isPointInsideTheBar( const csmp::Point<1U>& pt,
                           const csmp::Point<1U>& pt1,
                           const csmp::Point<1U>& pt2 )
{
    return (pt[0]-pt1[0]) * (pt2[0]-pt[0]) >= 0.0;
}

bool isPointInsideTheBar( const csmp::Point<2U>& pt,
                           const csmp::Point<2U>& pt1,
                           const csmp::Point<2U>& pt2 )
{
    if ( unsignedArea(pt, pt1, pt2) > 0.0 ) return false;
    csmp::Point<2U> e0 = pt2; e0 -= pt1;
    csmp::Point<2U> e1 = pt;  e1 -= pt1;
    const double proj = csmp::dotProduct(e0, e1);
    return (proj >= 0.0) && (proj <= e0.Length());
}

bool isPointInsideTheBar( const csmp::Point<3U>& pt,
                           const csmp::Point<3U>& pt1,
                           const csmp::Point<3U>& pt2 )
{
    if ( unsignedArea(pt, pt1, pt2) > 0.0 ) return false;
    csmp::Point<3U> e0 = pt2; e0 -= pt1;
    csmp::Point<3U> e1 = pt;  e1 -= pt1;
    const double proj = csmp::dotProduct(e0, e1);
    return (proj >= 0.0) && (proj <= e0.Length());
}

// -----------------------------------------------------------------------
// Triangle — area
// -----------------------------------------------------------------------

double unsignedArea( const csmp::Point<1U>&,
                      const csmp::Point<1U>&,
                      const csmp::Point<1U>& )
{
    return 0.0;
}

double unsignedArea( const csmp::Point<2U>& pt1,
                      const csmp::Point<2U>& pt2,
                      const csmp::Point<2U>& pt3 )
{
    return std::abs(signedArea(pt1, pt2, pt3));
}

double unsignedArea( const csmp::Point<3U>& pt1,
                      const csmp::Point<3U>& pt2,
                      const csmp::Point<3U>& pt3 )
{
    csmp::Point<3U> avec = pt2; avec -= pt1;   // a = pt2 - pt1
    csmp::Point<3U> bvec = pt3; bvec -= pt1;   // b = pt3 - pt1  (fixed)
    return 0.5 * csmp::crossProduct(avec, bvec).Length();
}

double signedArea( const csmp::Point<1U>&,
                    const csmp::Point<1U>&,
                    const csmp::Point<1U>& )
{
    return 0.0;
}

double signedArea( const csmp::Point<2U>& pt1,
                    const csmp::Point<2U>& pt2,
                    const csmp::Point<2U>& pt3 )
{
    double sarea = 0.0;
    sarea += (pt1[0]-pt3[0]) * (pt2[1]-pt3[1]);
    sarea -= (pt2[0]-pt3[0]) * (pt1[1]-pt3[1]);
    return 0.5 * sarea;
}

double signedArea( const csmp::Point<3U>& pt1,
                    const csmp::Point<3U>& pt2,
                    const csmp::Point<3U>& pt3 )
{
    return unsignedArea(pt1, pt2, pt3);
}

// -----------------------------------------------------------------------
// Triangle — orientation
// -----------------------------------------------------------------------

bool isCounterClockWiseOrientation( const csmp::Point<1U>&,
                                     const csmp::Point<1U>&,
                                     const csmp::Point<1U>& )
{
    return true;
}

bool isCounterClockWiseOrientation( const csmp::Point<2U>& pt1,
                                     const csmp::Point<2U>& pt2,
                                     const csmp::Point<2U>& pt3 )
{
    double sarea = (pt1[0]-pt3[0])*(pt2[1]-pt3[1])
                 - (pt2[0]-pt3[0])*(pt1[1]-pt3[1]);
    return sarea > 0.0;
}

bool isCounterClockWiseOrientation( const csmp::Point<3U>&,
                                     const csmp::Point<3U>&,
                                     const csmp::Point<3U>& )
{
    return true;
}

// -----------------------------------------------------------------------
// Triangle — point containment
// -----------------------------------------------------------------------

bool isPointInsideTheTriangle( const csmp::Point<1U>& pt,
                                const csmp::Point<1U>& pt1,
                                const csmp::Point<1U>& pt2,
                                const csmp::Point<1U>& pt3 )
{
    const double lo = std::min({pt1[0], pt2[0], pt3[0]});
    const double hi = std::max({pt1[0], pt2[0], pt3[0]});
    return (pt[0] >= lo) && (pt[0] <= hi);
}

bool isPointInsideTheTriangle( const csmp::Point<2U>& pt,
                                const csmp::Point<2U>& pt1,
                                const csmp::Point<2U>& pt2,
                                const csmp::Point<2U>& pt3 )
{
    const double a1 = signedArea(pt, pt1, pt2);
    if (a1 == 0.0) return true;
    const double a2 = signedArea(pt, pt2, pt3);
    if (a2 == 0.0) return true;
    if ((a1 > 0.0) != (a2 > 0.0)) return false;
    const double a3 = signedArea(pt, pt3, pt1);
    if (a3 == 0.0) return true;
    return (a2 > 0.0) == (a3 > 0.0);
}

bool isPointInsideTheTriangle( const csmp::Point<3U>& pt,
                                const csmp::Point<3U>& pt1,
                                const csmp::Point<3U>& pt2,
                                const csmp::Point<3U>& pt3 )
{
    csmp::Point<3U> e12 = pt2; e12 -= pt1;
    csmp::Point<3U> e23 = pt3; e23 -= pt2;
    const csmp::Point<3U> tri_normal = csmp::crossProduct(e12, e23);

    // Edge 12
    csmp::Point<3U> e10 = pt; e10 -= pt1;
    const double a1 = csmp::dotProduct(csmp::crossProduct(e12, e10), tri_normal);
    if (a1 == 0.0) return true;

    // Edge 23
    csmp::Point<3U> e20 = pt; e20 -= pt2;
    const double a2 = csmp::dotProduct(csmp::crossProduct(e23, e20), tri_normal);
    if (a2 == 0.0) return true;
    if ((a1 > 0.0) != (a2 > 0.0)) return false;

    // Edge 31
    csmp::Point<3U> e31 = pt1; e31 -= pt3;   // fixed: was corrupting tri_normal
    csmp::Point<3U> e30 = pt;  e30 -= pt3;   // fixed: was modifying e20
    const double a3 = csmp::dotProduct(csmp::crossProduct(e31, e30), tri_normal);
    if (a3 == 0.0) return true;
    return (a2 > 0.0) == (a3 > 0.0);
}

// -----------------------------------------------------------------------
// Same side test
// -----------------------------------------------------------------------

bool isSameSide( const csmp::Point<2U>& pt1,
                  const csmp::Point<2U>& pt2,
                  const csmp::Point<2U>& pt3,
                  const csmp::Point<2U>& pt4 )
{
    const double a1 = signedArea(pt1, pt3, pt4);
    if (a1 == 0.0) return true;
    const double a2 = signedArea(pt2, pt3, pt4);
    if (a2 == 0.0) return true;
    return (a1 > 0.0) == (a2 > 0.0);
}

bool isSameSide( const csmp::Point<3U>& pt1,
                  const csmp::Point<3U>& pt2,
                  const csmp::Point<3U>& pt3,
                  const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> e0  = pt4; e0  -= pt3;
    csmp::Point<3U> e1  = pt1; e1  -= pt3;
    csmp::Point<3U> e2  = pt2; e2  -= pt3;
    const csmp::Point<3U> e01 = csmp::crossProduct(e0, e1);
    const csmp::Point<3U> e02 = csmp::crossProduct(e0, e2);
    return csmp::dotProduct(e01, e02) >= 0.0;
}

// -----------------------------------------------------------------------
// Quadrilateral / 4-point area
// -----------------------------------------------------------------------

double unsignedArea( const csmp::Point<1U>&,
                      const csmp::Point<1U>&,
                      const csmp::Point<1U>&,
                      const csmp::Point<1U>& )
{
    return 0.0;
}

double unsignedArea( const csmp::Point<2U>& pt1,
                      const csmp::Point<2U>& pt2,
                      const csmp::Point<2U>& pt3,
                      const csmp::Point<2U>& pt4 )
{
    return unsignedArea(pt1, pt2, pt3)
         + unsignedArea(pt1, pt3, pt4);
}

double unsignedArea( const csmp::Point<3U>& pt1,
                      const csmp::Point<3U>& pt2,
                      const csmp::Point<3U>& pt3,
                      const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;
    return 0.5 * ( csmp::crossProduct(a, b).Length()
                 + csmp::crossProduct(b, c).Length() );
}

double signedArea( const csmp::Point<1U>&,
                    const csmp::Point<1U>&,
                    const csmp::Point<1U>&,
                    const csmp::Point<1U>& )
{
    return 0.0;
}

double signedArea( const csmp::Point<2U>& pt1,
                    const csmp::Point<2U>& pt2,
                    const csmp::Point<2U>& pt3,
                    const csmp::Point<2U>& pt4 )
{
    return signedArea(pt1, pt2, pt3)
         + signedArea(pt1, pt3, pt4);
}

double signedArea( const csmp::Point<3U>& pt1,
                    const csmp::Point<3U>& pt2,
                    const csmp::Point<3U>& pt3,
                    const csmp::Point<3U>& pt4 )
{
    return signedArea(pt1, pt2, pt3)
         + signedArea(pt1, pt3, pt4);
}

// -----------------------------------------------------------------------
// Tetrahedron volume
// Fixed: added factor of 1/6 to signedVolume
// -----------------------------------------------------------------------

double signedVolume( const csmp::Point<3U>& pt1,
                      const csmp::Point<3U>& pt2,
                      const csmp::Point<3U>& pt3,
                      const csmp::Point<3U>& pt4 )
{
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;
    return csmp::dotProduct(a, csmp::crossProduct(b, c)) / 6.0;
}

double unsignedVolume( const csmp::Point<3U>& pt1,
                        const csmp::Point<3U>& pt2,
                        const csmp::Point<3U>& pt3,
                        const csmp::Point<3U>& pt4 )
{
    return std::abs(signedVolume(pt1, pt2, pt3, pt4));
}

// -----------------------------------------------------------------------
// Orientation (4-point)
// -----------------------------------------------------------------------

bool isCounterClockWiseOrientation( const csmp::Point<1U>&,
                                     const csmp::Point<1U>&,
                                     const csmp::Point<1U>&,
                                     const csmp::Point<1U>& )
{
    return true;
}

bool isCounterClockWiseOrientation( const csmp::Point<2U>&,
                                     const csmp::Point<2U>& pt2,
                                     const csmp::Point<2U>& pt3,
                                     const csmp::Point<2U>& pt4 )
{
    return isCounterClockWiseOrientation(pt2, pt3, pt4);
}

bool isCounterClockWiseOrientation( const csmp::Point<3U>& pt1,
                                     const csmp::Point<3U>& pt2,
                                     const csmp::Point<3U>& pt3,
                                     const csmp::Point<3U>& pt4 )
{
    return signedVolume(pt1, pt2, pt3, pt4) > 0.0;
}

// -----------------------------------------------------------------------
// Dihedral angle
//
// Convention: angle between the two half-planes sharing edge (edge_a,
// edge_b), one containing side1 and the other containing side2.
// Fixed: dihedralRadAngle 3D was incorrectly multiplying by 180.
// -----------------------------------------------------------------------

double dihedralDegAngle( const csmp::Point<1U>&,
                          const csmp::Point<1U>&,
                          const csmp::Point<1U>&,
                          const csmp::Point<1U>& )
{
    return 0.0;
}

double dihedralRadAngle( const csmp::Point<1U>&,
                          const csmp::Point<1U>&,
                          const csmp::Point<1U>&,
                          const csmp::Point<1U>& )
{
    return 0.0;
}

double dihedralDegAngle( const csmp::Point<2U>& edge_a,
                          const csmp::Point<2U>& edge_b,
                          const csmp::Point<2U>& side1,
                          const csmp::Point<2U>& side2 )
{
    return isSameSide(side1, side2, edge_a, edge_b) ? 0.0 : 180.0;
}

double dihedralRadAngle( const csmp::Point<2U>& edge_a,
                          const csmp::Point<2U>& edge_b,
                          const csmp::Point<2U>& side1,
                          const csmp::Point<2U>& side2 )
{
    return isSameSide(side1, side2, edge_a, edge_b)
           ? 0.0
           : std::numbers::pi;
}

double dihedralDegAngle( const csmp::Point<3U>& edge_a,
                          const csmp::Point<3U>& edge_b,
                          const csmp::Point<3U>& side1,
                          const csmp::Point<3U>& side2 )
{
    csmp::Point<3U> a = edge_b; a -= edge_a;   // edge direction
    csmp::Point<3U> b = side1;  b -= edge_a;   // toward side1
    csmp::Point<3U> c = side2;  c -= edge_a;   // toward side2

    csmp::Point<3U> n1 = csmp::crossProduct(b, a);
    const double n1len = n1.Length();
    if (n1len < std::numeric_limits<double>::epsilon()) return 0.0;
    n1 /= n1len;

    csmp::Point<3U> n2 = csmp::crossProduct(c, a);
    const double n2len = n2.Length();
    if (n2len < std::numeric_limits<double>::epsilon()) return 0.0;
    n2 /= n2len;

    const double cosA = std::clamp(csmp::dotProduct(n1, n2), -1.0, 1.0);
    return std::acos(cosA) * 180.0 / csmp::PI;
}

double dihedralRadAngle( const csmp::Point<3U>& edge_a,
                          const csmp::Point<3U>& edge_b,
                          const csmp::Point<3U>& side1,
                          const csmp::Point<3U>& side2 )
{
    csmp::Point<3U> a = edge_b; a -= edge_a;
    csmp::Point<3U> b = side1;  b -= edge_a;
    csmp::Point<3U> c = side2;  c -= edge_a;

    csmp::Point<3U> n1 = csmp::crossProduct(b, a);
    const double n1len = n1.Length();
    if (n1len < std::numeric_limits<double>::epsilon()) return 0.0;
    n1 /= n1len;

    csmp::Point<3U> n2 = csmp::crossProduct(c, a);
    const double n2len = n2.Length();
    if (n2len < std::numeric_limits<double>::epsilon()) return 0.0;
    n2 /= n2len;

    const double cosA = std::clamp(csmp::dotProduct(n1, n2), -1.0, 1.0);
    return std::acos(cosA);   // fixed: was incorrectly multiplying by 180
}

// -----------------------------------------------------------------------
// Tetra identification
// -----------------------------------------------------------------------

bool isTetra( const csmp::Point<3U>& pt1,
               const csmp::Point<3U>& pt2,
               const csmp::Point<3U>& pt3,
               const csmp::Point<3U>& pt4 )
{
    // Non-zero volume => tetrahedron, zero volume => coplanar (quad)
    return std::abs(signedVolume(pt1, pt2, pt3, pt4))
           > std::numeric_limits<double>::epsilon();
}

bool isTetra( const csmp::Point<3U>& pt1,
               const csmp::Point<3U>& pt2,
               const csmp::Point<3U>& pt3,
               const csmp::Point<3U>& pt4,
               std::map<size_t,size_t>& order )
{
    const double vol = signedVolume(pt1, pt2, pt3, pt4);

    if ( std::abs(vol) > std::numeric_limits<double>::epsilon() )
    {
        if (vol > 0.0) {
            order = {{0,0},{1,1},{2,2},{3,3}};
        } else {
            // Swap vertices 2 and 3 to flip orientation
            order = {{0,0},{1,1},{2,3},{3,2}};
        }
        return true;
    }

    // Coplanar — identify which diagonal gives the quadrilateral
    csmp::Point<3U> a = pt2; a -= pt1;
    csmp::Point<3U> b = pt3; b -= pt1;
    csmp::Point<3U> c = pt4; c -= pt1;

    const double area_ab = csmp::crossProduct(a, b).Length();
    const double area_bc = csmp::crossProduct(b, c).Length();
    const double area_ac = csmp::crossProduct(a, c).Length();

    // Identify which point lies on the diagonal of the quad
    if ( std::abs(area_ab + area_bc - 2.0*area_ac)
         < std::numeric_limits<double>::epsilon() )
    {
        order = {{0,0},{1,1},{2,2},{3,3}};
        return false;
    }
    if ( std::abs(area_ac + area_bc - 2.0*area_ab)
         < std::numeric_limits<double>::epsilon() )
    {
        order = {{0,0},{1,1},{2,3},{3,2}};
        return false;
    }
    if ( std::abs(area_ab + area_ac - 2.0*area_bc)
         < std::numeric_limits<double>::epsilon() )
    {
        order = {{0,0},{1,2},{2,1},{3,3}};
        return false;
    }
    return false;
}

// -----------------------------------------------------------------------
// Element validity — equiangular skewness
// -----------------------------------------------------------------------

/**
 * @brief Compute equiangular skewness from a set of dihedral angles.
 *
 * @param angles    Vector of dihedral angles in degrees
 * @param theta_e   Equilateral reference angle in degrees
 * @return          Skewness in [0,1]
 */
static double equiangularSkewness( const std::vector<double>& angles,
                                    double theta_e )
{
    const double qmax = *std::max_element(angles.begin(), angles.end());
    const double qmin = *std::min_element(angles.begin(), angles.end());
    return std::max( (qmax - theta_e) / (180.0 - theta_e),
                     (theta_e - qmin) / theta_e );
}

// -----------------------------------------------------------------------

bool isTetrahedron( const std::vector<Point<3U>>& v )
{
    assert(v.size() == 4);

    // A regular tetrahedron has 6 edges and 6 dihedral angles.
    // Reference angle: arccos(1/3) ≈ 70.5288°
    //
    // Edge labelling (vertex indices):
    //   01, 02, 03, 12, 13, 23
    const std::vector<double> angles = {
        dihedralDegAngle(v[0], v[1], v[2], v[3]),   // edge 01
        dihedralDegAngle(v[0], v[2], v[1], v[3]),   // edge 02  (was missing)
        dihedralDegAngle(v[0], v[3], v[1], v[2]),   // edge 03
        dihedralDegAngle(v[1], v[2], v[0], v[3]),   // edge 12
        dihedralDegAngle(v[1], v[3], v[0], v[2]),   // edge 13
        dihedralDegAngle(v[2], v[3], v[0], v[1]),   // edge 23
    };

    return equiangularSkewness(angles, TETRA_EQUILATERAL_DIHEDRAL_DEG)
           <= SKEWNESS_THRESHOLD_TETRA;
}

// -----------------------------------------------------------------------

bool isPyramid( const std::vector<Point<3U>>& v )
{
    assert(v.size() == 5);

    // Vertex layout:
    //   Base quad: v[0], v[1], v[2], v[3]  (bottom face)
    //   Apex:      v[4]
    //
    // 8 edges: 4 base edges + 4 lateral edges
    // Evaluate base edges (ref 90°) and lateral edges (ref 60°) separately.

    const std::vector<double> lateral_angles = {
        dihedralDegAngle(v[0], v[1], v[3], v[4]),   // lateral edge 01
        dihedralDegAngle(v[1], v[2], v[0], v[4]),   // lateral edge 12
        dihedralDegAngle(v[2], v[3], v[1], v[4]),   // lateral edge 23
        dihedralDegAngle(v[3], v[0], v[2], v[4]),   // lateral edge 30
    };

    const std::vector<double> base_angles = {
        dihedralDegAngle(v[3], v[1], v[0], v[4]),   // base edge 31
        dihedralDegAngle(v[0], v[2], v[1], v[4]),   // base edge 02
        dihedralDegAngle(v[1], v[3], v[2], v[4]),   // base edge 13
        dihedralDegAngle(v[2], v[0], v[3], v[4]),   // base edge 20
    };

    const double q_lateral = equiangularSkewness(lateral_angles, 60.0);
    const double q_base    = equiangularSkewness(base_angles,    PYRAMID_BASE_DIHEDRAL_DEG);

    return std::max(q_lateral, q_base) <= SKEWNESS_THRESHOLD_PYRAMID;
}

// -----------------------------------------------------------------------

bool isPrism( const std::vector<Point<3U>>& v )
{
    assert(v.size() == 6);

    // Vertex layout:
    //   Bottom triangle: v[0], v[1], v[2]
    //   Top    triangle: v[3], v[4], v[5]
    //   Lateral edges:   v[0]-v[3], v[1]-v[4], v[2]-v[5]
    //
    // 9 edges total.
    // Bottom/top triangle edges: ref 60° (equilateral triangle)
    // Lateral quad edges:        ref 90° (rectangle)

    const std::vector<double> tri_angles = {
        dihedralDegAngle(v[0], v[1], v[2], v[3]),   // bottom edge 01
        dihedralDegAngle(v[1], v[2], v[0], v[4]),   // bottom edge 12
        dihedralDegAngle(v[2], v[0], v[1], v[5]),   // bottom edge 20
        dihedralDegAngle(v[3], v[4], v[5], v[0]),   // top edge 34
        dihedralDegAngle(v[4], v[5], v[3], v[1]),   // top edge 45
        dihedralDegAngle(v[5], v[3], v[4], v[2]),   // top edge 53
    };

    const std::vector<double> lateral_angles = {
        dihedralDegAngle(v[0], v[3], v[1], v[2]),   // lateral edge 03
        dihedralDegAngle(v[1], v[4], v[0], v[2]),   // lateral edge 14
        dihedralDegAngle(v[2], v[5], v[0], v[1]),   // lateral edge 25
    };

    const double q_tri     = equiangularSkewness(tri_angles,     60.0);
    const double q_lateral = equiangularSkewness(lateral_angles, PRISM_LATERAL_DIHEDRAL_DEG);

    return std::max(q_tri, q_lateral) <= SKEWNESS_THRESHOLD_PRISM;
}

// -----------------------------------------------------------------------

bool isHexahedron( const std::vector<Point<3U>>& v )
{
    assert(v.size() == 8);

    // Vertex layout (standard hex):
    //   Bottom face: v[0],v[1],v[2],v[3]
    //   Top    face: v[4],v[5],v[6],v[7]
    //   Lateral edges: v[0]-v[4], v[1]-v[5], v[2]-v[6], v[3]-v[7]
    //
    // 12 edges, all reference angle 90° (perfect cube).

    const std::vector<double> angles = {
        // Bottom face edges
        dihedralDegAngle(v[0], v[1], v[3], v[4]),
        dihedralDegAngle(v[1], v[2], v[0], v[5]),
        dihedralDegAngle(v[2], v[3], v[1], v[6]),
        dihedralDegAngle(v[3], v[0], v[2], v[7]),
        // Top face edges
        dihedralDegAngle(v[4], v[5], v[7], v[0]),
        dihedralDegAngle(v[5], v[6], v[4], v[1]),
        dihedralDegAngle(v[6], v[7], v[5], v[2]),
        dihedralDegAngle(v[7], v[4], v[6], v[3]),
        // Lateral edges
        dihedralDegAngle(v[0], v[4], v[1], v[3]),
        dihedralDegAngle(v[1], v[5], v[0], v[2]),
        dihedralDegAngle(v[2], v[6], v[1], v[3]),
        dihedralDegAngle(v[3], v[7], v[0], v[2]),
    };

    return equiangularSkewness(angles, HEX_EQUILATERAL_DIHEDRAL_DEG)
           <= SKEWNESS_THRESHOLD_HEX;
}

// -----------------------------------------------------------------------

bool isValidElement( const std::vector<Point<3U>>& vertexList )
{
    switch (vertexList.size())
    {
        case 4:  return isTetrahedron(vertexList);
        case 5:  return isPyramid(vertexList);
        case 6:  return isPrism(vertexList);
        case 8:  return isHexahedron(vertexList);
        default:
            // Unrecognised element type — cannot assess validity
            return false;
    }
}

// -----------------------------------------------------------------------

/// Complete element validity check combining multiple criteria
template<uint32_t dim>
bool isValidElement( const Element<dim>* const elmt )
{
    // 1. Positive Jacobian at all Gauss points (primary check)
    for (uint32_t i=0; i<elmt->FE()->IntegrationPoints(); ++i)
    {
        DenseMatrix<DM_MIN> B;
        if (elmt->dN_AtIntegrationPoint(B, i, dim) <= 0.0)
            return false;
    }

    // 2. Equiangular skewness (secondary check, tighter threshold)
    std::vector<Point<dim>> verts;
    verts.reserve(elmt->Nodes());
    for (auto n=elmt->NodesBegin(); n!=elmt->NodesEnd(); ++n)
        verts.push_back((*n)->Coordinate());

    return isValidElement(verts);   // uses per-type thresholds
}

template bool isValidElement( const Element<3U>* const );

} // namespace csmp
