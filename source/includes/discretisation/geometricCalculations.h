#ifndef CSMP_GEOMETRIC_CALCULATIONS_UOM_H
#define CSMP_GEOMETRIC_CALCULATIONS_UOM_H

#include "CSMP_definitions.h"
#include <numbers>
#include "Point.h"

/**
 * @file geometricCalculations.h
 * @author R. Manasipov (2015)
 * @author SKM (2017, 2024)
 *
 * @brief Geometric calculations for finite element validity checking
 *        and general mesh geometry operations.
 *
 * Element validity is assessed using equiangular skewness:
 *
 *   q = max( (theta_max - theta_e) / (180 - theta_e),
 *            (theta_e  - theta_min) / theta_e )
 *
 * where theta_e is the dihedral angle of the corresponding perfect element.
 * Elements with q > SKEWNESS_THRESHOLD are considered too distorted for
 * reliable finite element computation.
 */

namespace csmp {

template<uint32_t> class Element;

// -----------------------------------------------------------------------
// Mesh quality constants
// -----------------------------------------------------------------------

/// Equiangular skewness threshold above which an element is rejected.
/// Range [0,1]: 0 = perfect, 1 = degenerate.
static constexpr double SKEWNESS_THRESHOLD_TETRA = 0.85;
static constexpr double SKEWNESS_THRESHOLD_HEX   = 0.90;
static constexpr double SKEWNESS_THRESHOLD_PRISM  = 0.85;
static constexpr double SKEWNESS_THRESHOLD_PYRAMID= 0.80;

/// Dihedral angle of a regular tetrahedron: arccos(1/3) ≈ 70.5288°
static constexpr double TETRA_EQUILATERAL_DIHEDRAL_DEG = 70.5287793655;

/// Dihedral angle of a perfect triangular prism lateral face: 90°
static constexpr double PRISM_LATERAL_DIHEDRAL_DEG = 90.0;

/// Dihedral angle of a perfect cube: 90°
static constexpr double HEX_EQUILATERAL_DIHEDRAL_DEG = 90.0;

/// Dihedral angle of a perfect square pyramid lateral face: 90°
static constexpr double PYRAMID_BASE_DIHEDRAL_DEG = 90.0;

// -----------------------------------------------------------------------
// Normal to point cloud
// -----------------------------------------------------------------------

Point<2> normalToAveragePlaneThroughPointCloud( const std::vector<Point<2>>& );

// -----------------------------------------------------------------------
// Smallest angle between two vectors
// -----------------------------------------------------------------------

/// Returns the smallest angle (degrees) between vectors a and b
/// (both originating from the coordinate origin).
template<uint32_t dim>
inline double smallestAngleBetween( const Point<dim>& a, const Point<dim>& b )
{
    const double n1 = a.Length();
    const double n2 = b.Length();

    if ( n1 < std::numeric_limits<double>::epsilon() ||
         n2 < std::numeric_limits<double>::epsilon() )
        throw std::domain_error(
            "smallestAngleBetween: angle undefined for zero-length vector");

    double c = dotProduct(a, b) / (n1 * n2);
    c = std::clamp(c, -1.0, 1.0);
    return std::acos(c) * 180.0 / std::numbers::pi;
}

// -----------------------------------------------------------------------
// Cartesian axes (Gram-Schmidt orthogonalisation)
// -----------------------------------------------------------------------

/// 1D: unit vector from pt1 to pt2 → e1
template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1,
                        const csmp::Point<dim>& pt2,
                        csmp::Point<dim>& e1 );

/// 2D: orthonormal basis {e1, e2} from three points
template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1,
                        const csmp::Point<dim>& pt2,
                        const csmp::Point<dim>& pt3,
                        csmp::Point<dim>& e1,
                        csmp::Point<dim>& e2 );

/// 3D: orthonormal basis {e1, e2, e3} from four points
template<uint32_t dim>
void getCartesianAxes( const csmp::Point<dim>& pt1,
                        const csmp::Point<dim>& pt2,
                        const csmp::Point<dim>& pt3,
                        const csmp::Point<dim>& pt4,
                        csmp::Point<dim>& e1,
                        csmp::Point<dim>& e2,
                        csmp::Point<dim>& e3 );

// -----------------------------------------------------------------------
// Local coordinate projection
// -----------------------------------------------------------------------

/// 1D: project pt1 onto axis e1 relative to origin pt0
template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0,
                     const csmp::Point<dim>& e1,
                     const csmp::Point<dim>& pt1,
                     csmp::Point<dim>& pt2 );

/// 2D: project pt1 onto axes {e1,e2} relative to origin pt0
template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0,
                     const csmp::Point<dim>& e1,
                     const csmp::Point<dim>& e2,
                     const csmp::Point<dim>& pt1,
                     csmp::Point<dim>& pt2 );

/// 3D: project pt1 onto axes {e1,e2,e3} relative to origin pt0
template<uint32_t dim>
void getCoordinate( const csmp::Point<dim>& pt0,
                     const csmp::Point<dim>& e1,
                     const csmp::Point<dim>& e2,
                     const csmp::Point<dim>& e3,
                     const csmp::Point<dim>& pt1,
                     csmp::Point<dim>& pt2 );

// -----------------------------------------------------------------------
// Bar (line segment)
// -----------------------------------------------------------------------

double unsignedLength( const csmp::Point<1U>&, const csmp::Point<1U>& );
double unsignedLength( const csmp::Point<2U>&, const csmp::Point<2U>& );
double unsignedLength( const csmp::Point<3U>&, const csmp::Point<3U>& );

bool isPointInsideTheBar( const csmp::Point<1U>& pt,
                           const csmp::Point<1U>&,
                           const csmp::Point<1U>& );
bool isPointInsideTheBar( const csmp::Point<2U>& pt,
                           const csmp::Point<2U>&,
                           const csmp::Point<2U>& );
bool isPointInsideTheBar( const csmp::Point<3U>& pt,
                           const csmp::Point<3U>&,
                           const csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Triangle
// -----------------------------------------------------------------------

double unsignedArea( const csmp::Point<1U>&,
                      const csmp::Point<1U>&,
                      const csmp::Point<1U>& );
double unsignedArea( const csmp::Point<2U>&,
                      const csmp::Point<2U>&,
                      const csmp::Point<2U>& );
double unsignedArea( const csmp::Point<3U>&,
                      const csmp::Point<3U>&,
                      const csmp::Point<3U>& );

double signedArea( const csmp::Point<1U>&,
                    const csmp::Point<1U>&,
                    const csmp::Point<1U>& );
double signedArea( const csmp::Point<2U>&,
                    const csmp::Point<2U>&,
                    const csmp::Point<2U>& );
double signedArea( const csmp::Point<3U>&,
                    const csmp::Point<3U>&,
                    const csmp::Point<3U>& );

bool isCounterClockWiseOrientation( const csmp::Point<1U>&,
                                     const csmp::Point<1U>&,
                                     const csmp::Point<1U>& );
bool isCounterClockWiseOrientation( const csmp::Point<2U>&,
                                     const csmp::Point<2U>&,
                                     const csmp::Point<2U>& );
bool isCounterClockWiseOrientation( const csmp::Point<3U>&,
                                     const csmp::Point<3U>&,
                                     const csmp::Point<3U>& );

bool isPointInsideTheTriangle( const csmp::Point<1U>& pt,
                                const csmp::Point<1U>&,
                                const csmp::Point<1U>&,
                                const csmp::Point<1U>& );
bool isPointInsideTheTriangle( const csmp::Point<2U>& pt,
                                const csmp::Point<2U>&,
                                const csmp::Point<2U>&,
                                const csmp::Point<2U>& );
bool isPointInsideTheTriangle( const csmp::Point<3U>& pt,
                                const csmp::Point<3U>&,
                                const csmp::Point<3U>&,
                                const csmp::Point<3U>& );

bool isSameSide( const csmp::Point<2U>& pt,
                  const csmp::Point<2U>&,
                  const csmp::Point<2U>&,
                  const csmp::Point<2U>& );
bool isSameSide( const csmp::Point<3U>& pt,
                  const csmp::Point<3U>&,
                  const csmp::Point<3U>&,
                  const csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Quadrilateral / Tetrahedron (4-point)
// -----------------------------------------------------------------------

double unsignedArea( const csmp::Point<1U>&,
                      const csmp::Point<1U>&,
                      const csmp::Point<1U>&,
                      const csmp::Point<1U>& );
double unsignedArea( const csmp::Point<2U>&,
                      const csmp::Point<2U>&,
                      const csmp::Point<2U>&,
                      const csmp::Point<2U>& );
double unsignedArea( const csmp::Point<3U>&,
                      const csmp::Point<3U>&,
                      const csmp::Point<3U>&,
                      const csmp::Point<3U>& );

double signedArea( const csmp::Point<1U>&,
                    const csmp::Point<1U>&,
                    const csmp::Point<1U>&,
                    const csmp::Point<1U>& );
double signedArea( const csmp::Point<2U>&,
                    const csmp::Point<2U>&,
                    const csmp::Point<2U>&,
                    const csmp::Point<2U>& );
double signedArea( const csmp::Point<3U>&,
                    const csmp::Point<3U>&,
                    const csmp::Point<3U>&,
                    const csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Tetrahedron volume
// -----------------------------------------------------------------------

/// Unsigned volume of tetrahedron (pt1,pt2,pt3,pt4)
double unsignedVolume( const csmp::Point<3U>&,
                        const csmp::Point<3U>&,
                        const csmp::Point<3U>&,
                        const csmp::Point<3U>& );

/// Signed volume = (1/6) * det[pt2-pt1, pt3-pt1, pt4-pt1]
/// Positive when vertices are in right-hand (CCW) order.
double signedVolume( const csmp::Point<3U>&,
                      const csmp::Point<3U>&,
                      const csmp::Point<3U>&,
                      const csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Orientation
// -----------------------------------------------------------------------

bool isCounterClockWiseOrientation( const csmp::Point<1U>&,
                                     const csmp::Point<1U>&,
                                     const csmp::Point<1U>&,
                                     const csmp::Point<1U>& );
bool isCounterClockWiseOrientation( const csmp::Point<2U>&,
                                     const csmp::Point<2U>&,
                                     const csmp::Point<2U>&,
                                     const csmp::Point<2U>& );
bool isCounterClockWiseOrientation( const csmp::Point<3U>&,
                                     const csmp::Point<3U>&,
                                     const csmp::Point<3U>&,
                                     const csmp::Point<3U>& );

// -----------------------------------------------------------------------
// Dihedral angle
//
// Convention: angle between the two half-planes that share edge (edge_a,
// edge_b), one containing side1 and the other containing side2.
//
//   side1
//     x
//      \
//       \
//  edge_a x-------x edge_b
//       /
//      /
//     x
//   side2
// -----------------------------------------------------------------------

double dihedralDegAngle( const csmp::Point<1U>& edge_a,
                          const csmp::Point<1U>& edge_b,
                          const csmp::Point<1U>& side1,
                          const csmp::Point<1U>& side2 );
double dihedralDegAngle( const csmp::Point<2U>& edge_a,
                          const csmp::Point<2U>& edge_b,
                          const csmp::Point<2U>& side1,
                          const csmp::Point<2U>& side2 );
double dihedralDegAngle( const csmp::Point<3U>& edge_a,
                          const csmp::Point<3U>& edge_b,
                          const csmp::Point<3U>& side1,
                          const csmp::Point<3U>& side2 );

double dihedralRadAngle( const csmp::Point<1U>& edge_a,
                          const csmp::Point<1U>& edge_b,
                          const csmp::Point<1U>& side1,
                          const csmp::Point<1U>& side2 );
double dihedralRadAngle( const csmp::Point<2U>& edge_a,
                          const csmp::Point<2U>& edge_b,
                          const csmp::Point<2U>& side1,
                          const csmp::Point<2U>& side2 );
double dihedralRadAngle( const csmp::Point<3U>& edge_a,
                          const csmp::Point<3U>& edge_b,
                          const csmp::Point<3U>& side1,
                          const csmp::Point<3U>& side2 );

// -----------------------------------------------------------------------
// Tetra identification
// -----------------------------------------------------------------------

/// Returns true if the 4 points form a non-degenerate tetrahedron
/// (non-zero volume), false if they are coplanar (quadrilateral).
bool isTetra( const csmp::Point<3U>&,
               const csmp::Point<3U>&,
               const csmp::Point<3U>&,
               const csmp::Point<3U>& );

/// As above, and additionally returns the vertex permutation that gives
/// CCW (right-hand rule) orientation in @p order.
bool isTetra( const csmp::Point<3U>&,
               const csmp::Point<3U>&,
               const csmp::Point<3U>&,
               const csmp::Point<3U>&,
               std::map<size_t,size_t>& order );

// -----------------------------------------------------------------------
// Element validity — equiangular skewness
// -----------------------------------------------------------------------

/**
 * @brief Returns true if the tetrahedron is fit for FE computation.
 *
 * Uses equiangular skewness on all 6 dihedral angles.
 * Reference angle: arccos(1/3) ≈ 70.53° (regular tetrahedron).
 * Rejects elements with skewness > SKEWNESS_THRESHOLD.
 */
bool isTetrahedron( const std::vector<Point<3U>>& vertexList );

/**
 * @brief Returns true if the hexahedron is fit for FE computation.
 *
 * Uses equiangular skewness on all 12 dihedral angles.
 * Reference angle: 90° (perfect cube).
 */
bool isHexahedron( const std::vector<Point<3U>>& vertexList );

/**
 * @brief Returns true if the triangular prism is fit for FE computation.
 *
 * Evaluates triangular face angles (ref 60°) and lateral face angles
 * (ref 90°) separately.
 */
bool isPrism( const std::vector<Point<3U>>& vertexList );

/**
 * @brief Returns true if the square pyramid is fit for FE computation.
 *
 * Evaluates triangular face angles (ref 60°) and base face angles
 * (ref 90°) separately.
 */
bool isPyramid( const std::vector<Point<3U>>& vertexList );

/**
 * @brief Dispatches to the appropriate element validity check based on
 *        vertex count. Returns false for unrecognised element types.
 *
 * Supported: tetrahedron (4), pyramid (5), prism (6), hexahedron (8).
 */
bool isValidElement( const std::vector<Point<3U>>& vertexList );

/**
 * @brief Tests whether a finite element is fit for computation by
 *        extracting its node coordinates and calling isValidElement.
 *
 * Currently instantiated for dim=3 only.
 */
template<uint32_t dim>
bool isValidElement( const Element<dim>* const );

} // namespace csmp

#endif // CSMP_GEOMETRIC_CALCULATIONS_UOM_H

