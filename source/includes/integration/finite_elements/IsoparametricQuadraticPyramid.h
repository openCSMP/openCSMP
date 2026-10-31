// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ISOPARAMETRIC_QUADRATIC_PYRAMID_H
#define ISOPARAMETRIC_QUADRATIC_PYRAMID_H

#include "FiniteElement.h"

#include <array>
#include <cstdint>
#include <vector>

namespace csmp {

/**
 @class IsoparametricQuadraticPyramid
 @brief 13-node isoparametric quadratic pyramid finite element.

 @author S.K. Matthaei
 @author Stephen G. Roberts
 @date   1998
 @author (refactored) 2026

 @section geometry Geometry

 Node numbering:

 @code
          4 (apex)
         /|\
        / | \
       /  |  \
      3---+---2
      |   |   |
      0-------1   (base, CCW from outside)
 @endcode

 Reference coordinates:
 @code
   Corner nodes:
     Node  0: r=-1, s=-1, t=0
     Node  1: r=+1, s=-1, t=0
     Node  2: r=+1, s=+1, t=0
     Node  3: r=-1, s=+1, t=0
     Node  4: r= 0, s= 0, t=1  (apex)

   Midside nodes — base edges:
     Node  5: r= 0, s=-1, t=0   midpoint of edge 0-1
     Node  6: r=+1, s= 0, t=0   midpoint of edge 1-2
     Node  7: r= 0, s=+1, t=0   midpoint of edge 2-3
     Node  8: r=-1, s= 0, t=0   midpoint of edge 3-0

   Midside nodes — lateral edges to apex:
     Node  9: r=-0.5, s=-0.5, t=0.5   midpoint of edge 0-4
     Node 10: r=+0.5, s=-0.5, t=0.5   midpoint of edge 1-4
     Node 11: r=+0.5, s=+0.5, t=0.5   midpoint of edge 2-4
     Node 12: r=-0.5, s=+0.5, t=0.5   midpoint of edge 3-4
 @endcode

 Corner node reference coordinates are identical to those of
 IsoparametricLinearPyramid.

 @section faces Faces

 @code
   Face 0: {0,1,4,  5,10, 9}   triangular
   Face 1: {1,2,4,  6,11,10}   triangular
   Face 2: {2,3,4,  7,12,11}   triangular
   Face 3: {3,0,4,  8, 9,12}   triangular
   Face 4: {0,3,2,1, 8, 7, 6, 5}  quadrilateral (base)
 @endcode

 @section shapefunctions Shape Functions

 Bedrosian/Zganski rational formulation, corrected by Matthaei & Pichler
 (2006). Second-order complete on the pyramid domain. The rational term
 rs/(1-t) is singular at the apex (t=1); all methods guard against this
 with an explicit t==1 branch.

 @section integration Integration

 8-point conical product rule (same as IsoparametricLinearPyramid).
 Weight sum = 4/3 (volume of the reference pyramid).
*/
class IsoparametricQuadraticPyramid final : public FiniteElement
{
public:

    // -----------------------------------------------------------------------
    // Construction / destruction
    // -----------------------------------------------------------------------

    /**
     @brief Constructs the element and initialises reference coordinates,
            integration points and weights.
    */
    IsoparametricQuadraticPyramid();

    ~IsoparametricQuadraticPyramid() override;

    // -----------------------------------------------------------------------
    // Geometry
    // -----------------------------------------------------------------------

    /// @brief Computes element volume by 8-point numerical integration.
    double   Volume()                                        override final;

    /// @brief Aspect ratio: longest edge / shortest edge.
    double   AspectRatio()                                   override final;

    /// @brief Empirical inner radius estimate: V / (sum_edges / 6).
    double   InnerRadius()                                   override final;

    /**
     @brief Computes the lengths of the 8 corner edges.
     @param vec  Output: edge lengths in segment order (see NodesOfSegment).
    */
    void     EdgeLengths( std::vector<double>& vec )         override final;

    // -----------------------------------------------------------------------
    // Topology
    // -----------------------------------------------------------------------

    /**
     @brief Returns the 5 corner node local indices {0,1,2,3,4}.
     @param ids  Output vector of size 5.
    */
    void     CornerNodes( std::vector<uint32_t>& ids ) const override final;

    /// @brief Returns the number of corner nodes (always 5).
    uint32_t CornerNodes() const noexcept override final { return 5U; }

    /**
     @brief Returns the 8 midside node local indices {5,...,12}.
     @param ids  Output vector of size 8.
    */
    void     MidSideNodes( std::vector<uint32_t>& ids ) const override final;

    /// @brief Returns the number of midside nodes (always 8).
    uint32_t MidSideNodes() const noexcept override final { return 8U; }

    /**
     @brief Returns the number of nodes on a given face.
     @param face  Face index in [0,4].
     @return      8 for face 4 (base quad), 6 for faces 0–3 (triangles).
    */
    uint32_t NodesPerFace( uint32_t face ) const noexcept   override final;

    /**
     @brief Returns the element type of a given face.
     @param face  Face index in [0,4].
     @return      ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 for face 4,
                  ISOPARAMETRIC_QUADRATIC_TRIANGLE for faces 0–3.
    */
    CSMP_FEM_TYPE ElementTypeOfFace( uint32_t face ) const noexcept override final;

    /// @brief Returns ISOPARAMETRIC_QUADRATIC_BAR for all segments.
    CSMP_FEM_TYPE ElementTypeOfSegment( uint32_t ) const noexcept override final
        { return ISOPARAMETRIC_QUADRATIC_BAR; }

    /**
     @brief Returns the two endpoint nodes and the midside node of a segment.
     @param segm_id  Segment index in [0,7].
     @param snids    Output: {endpoint_a, endpoint_b, midside_node}.
    */
    void NodesOfSegment( uint32_t segm_id,
                         std::vector<uint32_t>& snids ) const override final;

    /**
     @brief Returns all nodes of a face (corner + midside), in order.
     @param face_id  Face index in [0,4].
     @return         Node index vector.
    */
    std::vector<uint32_t> NodesOfFace( uint32_t face_id )       const override final;

    /**
     @brief Returns only the corner nodes of a face.
     @param face_id  Face index in [0,4].
     @return         Corner node index vector (size 3 or 4).
    */
    std::vector<uint32_t> CornerNodesOfFace( uint32_t face_id ) const override final;

    /**
     @brief Returns the corner nodes directly connected to a given node.
     @param node_id  Local node index in [0,12].
     @return         Connected corner node index vector.
    */
    std::vector<uint32_t> NodesConnectedTo( uint32_t node_id )  const override final;

    /**
     @brief Computes the unit outward normal to a face at its centroid.
     @param face   Face index in [0,4].
     @param unrml  Output: unit normal vector of size 3.
    */
    void UnitNormalToFace( uint32_t face,
                           std::vector<double>& unrml ) const   override final;

    // -----------------------------------------------------------------------
    // Shape functions
    // -----------------------------------------------------------------------

    /**
     @brief Evaluates shape functions at a physical point xyz.
     @param N    Output: shape function values, size 13.
     @param xyz  Physical coordinates {x, y, z}.
    */
    void N( std::vector<double>& N,
            const std::vector<double>& xyz )                    override final;

    /**
     @brief Evaluates shape functions at integration point ip.
     @param ip  Integration point index in [0, gpe).
     @param N   Output: shape function values, size 13.
    */
    void N_AtIntegrationPoint( uint32_t ip,
                               std::vector<double>& N )         override final;

    /**
     @brief Evaluates shape functions at the element barycentre (0, 0, 0.25).
     @param N  Output: shape function values, size 13.
    */
    void N_AtBaryCenter( std::vector<double>& N )               override final;

    /**
     @brief Evaluates the 13 shape functions at parametric coordinates (r,s,t).

     Singular at t=1 (apex) — handled by an explicit t==1 branch.

     @param r     Parametric r in [-(1-t), 1-t].
     @param s     Parametric s in [-(1-t), 1-t].
     @param t     Parametric t in [0, 1].
     @param nrst  Output vector of size 13.
    */
    void Nrst( double r, double s, double t,
               std::vector<double>& nrst ) const                override final;

    // -----------------------------------------------------------------------
    // Shape function derivatives
    // -----------------------------------------------------------------------

    /**
     @brief Partial derivatives dN/dr at (r,s,t).
     @param r    Parametric r.
     @param s    Parametric s.
     @param t    Parametric t.
     @param dNr  Output vector of size 13.
    */
    void dNr( double r, double s, double t,
              std::vector<double>& dNr ) const                  override final;

    /**
     @brief Partial derivatives dN/ds at (r,s,t).
     @param r    Parametric r.
     @param s    Parametric s.
     @param t    Parametric t.
     @param dNs  Output vector of size 13.
    */
    void dNs( double r, double s, double t,
              std::vector<double>& dNs ) const                  override final;

    /**
     @brief Partial derivatives dN/dt at (r,s,t).

     Limiting values at t=1 are computed analytically via L'Hôpital's rule.

     @param r    Parametric r.
     @param s    Parametric s.
     @param t    Parametric t.
     @param dNt  Output vector of size 13.
    */
    void dNt( double r, double s, double t,
              std::vector<double>& dNt ) const                  override final;

    // -----------------------------------------------------------------------
    // Global shape function derivatives
    // -----------------------------------------------------------------------

    /**
     @brief Global shape function derivatives at each node's reference coords.
     @param DN14  Output: 3 x 13 matrix (row = spatial dim, col = node).
    */
    void   dN( DenseMatrix<DM_MIN>& DN14 )                      override final;

    /**
     @brief Global shape function derivatives at a physical point xyz.
     @param dn   Output: 3 x 13 B-matrix.
     @param xyz  Physical coordinates.
     @return     Jacobian determinant at xyz.
    */
    double dN_At( DenseMatrix<DM_MIN>& dn,
                  const std::vector<double>& xyz )              override final;

    /**
     @brief Global shape function derivatives at node nd.
     @param M   Output: 3 x 13 B-matrix.
     @param nd  Node index in [0, npe).
     @return    Jacobian determinant at this node.
    */
    double dN_AtNode( DenseMatrix<DM_MIN>& M,
                      uint32_t nd )                             override final;

    /**
     @brief Global shape function derivatives at integration point gauss_point.
     @param M             Output: 3 x 13 B-matrix.
     @param gauss_point   Integration point index in [0, gpe).
     @return              Jacobian determinant at this point.
    */
    double dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M,
                                  uint32_t gauss_point )        override final;

    /**
     @brief Global shape function derivatives at the element barycentre.
     @param M  Output: 3 x 13 B-matrix.
     @return   Jacobian determinant at the barycentre.
    */
    double dN_AtBarycenter( DenseMatrix<DM_MIN>& M )            override final;

    // -----------------------------------------------------------------------
    // Jacobian
    // -----------------------------------------------------------------------

    /**
     @brief Computes and stores the Jacobian matrix at integration point ip.
            Does not invert. Result is in base-class member JAC.
     @param ip  Integration point index in [0, gpe).
    */
    void JacobianAtIntegrationPoint( uint32_t ip )              override final;

    /**
     @brief Computes and stores the Jacobian matrix at arbitrary parametric
            coordinates rst. Does not invert.
     @param rst  Parametric coordinates {r, s, t}.
    */
    void JacobianAt( const std::vector<double>& rst ) noexcept override final;

    // -----------------------------------------------------------------------
    // Integration
    // -----------------------------------------------------------------------

    /**
     @brief Returns the integration weight at point i.
     @param i  Integration point index in [0, gpe).
     @return   Weight value.
    */
    double WeightAtIntegrationPoint( uint32_t i ) const noexcept override final;

    /**
     @brief Returns the physical coordinates of integration point i.
     @param i    Integration point index in [0, gpe).
     @param xyz  Output: physical coordinates {x, y, z}.
    */
    void IntegrationPoint( uint32_t i,
                           std::vector<double>& xyz ) const     override final;

    /**
     @brief Extrapolates integration point variables to element nodes.

     For each node, evaluates the shape functions at that node's reference
     coordinates and uses them as weights over the integration point values.
     Exact for fields representable by the element's shape functions.

     @param nvars  Number of variables per integration point.
     @param IVAR   Input:  size gpe * nvars, stored ip-major.
     @param NVAR   Output: size npe * nvars, stored node-major.
     @throws csmp::Exception (FATAL_ERROR) if IVAR.size() != gpe * nvars.
    */
    void ExtrapolateIntegrationPointVariableToNodes(
             uint32_t nvars,
             const std::vector<double>& IVAR,
             std::vector<double>& NVAR ) const                  override final;

    // -----------------------------------------------------------------------
    // Reference coordinates
    // -----------------------------------------------------------------------

    /**
     @brief Returns the reference (parametric) coordinates of all nodes.
     @param matCoords  Output: npe x 3 matrix. Row i = (r_i, s_i, t_i).
    */
    void ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const override final;

    // -----------------------------------------------------------------------
    // Output
    // -----------------------------------------------------------------------

    /**
     @brief Writes nodal data to a VTK unstructured grid file.

     Outputs the element as 5 VTK polylines (one per face outline) since
     VTK has no native 13-node pyramid cell type.

     @param file_name  Base filename (element index and ".vtk" are appended).
     @param var_name   Variable name written into the VTK header.
     @param DATA       1 x npe for scalars, dim x npe for vectors.
    */
    void OutputNodeDataToVTK( const char* file_name,
                              const char* var_name,
                              DenseMatrix<DM_MIN>& DATA ) const  override final;

private:

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------
    void Nrst( double r, double s, double t, double* nrst ) const noexcept;

    double assembleBMatrix( double r, double s, double t, DenseMatrix<DM_MIN>& B );

    /**
     @brief Generates the 8-point conical product integration rule.
     @param Ip  Output: 8 x 3 matrix of parametric coordinates.
     @param We  Output: vector of 8 weights. Sum = 4/3.
    */
    void GenerateIntegrationPoints( DenseMatrix<DM_MIN>& Ip,
                                    std::vector<double>& We );

    /**
     @brief Maps parametric coordinates to physical coordinates.
     @param rst  Input:  {r, s, t}.
     @param xyz  Output: {x, y, z}.
    */
    void ParametricToPhysical( std::vector<double>& rst,
                               std::vector<double>& xyz );

    /**
     @brief Maps physical coordinates to parametric coordinates via
            Newton-Raphson iteration.

     Issues a WARNING if convergence is not achieved within 20 iterations.

     @param rSt  Output: {r, s, t}.
     @param xyz  Input:  {x, y, z}.
    */
    void PhysicalToParametric( std::vector<double>& rSt,
                               const std::vector<double>& xyz );

    // -----------------------------------------------------------------------
    // Data members
    // -----------------------------------------------------------------------

    std::vector<double>  W;      ///< Integration weights [gpe]
    DenseMatrix<DM_MIN>  DN;     ///< Derivative matrix workspace [3 x 13]
    DenseMatrix<DM_MIN>  NXYZ;   ///< Reference node coordinates [13 x 3]
    DenseMatrix<DM_MIN>  IP;     ///< Integration point coordinates [8 x 3]

    /// Reference node coordinates as constexpr for zero-overhead access.
    static constexpr std::array<std::array<double,3>,13> kRefCoords = {{
        {{ -1.0, -1.0,  0.0 }},  // 0 corner
        {{  1.0, -1.0,  0.0 }},  // 1 corner
        {{  1.0,  1.0,  0.0 }},  // 2 corner
        {{ -1.0,  1.0,  0.0 }},  // 3 corner
        {{  0.0,  0.0,  1.0 }},  // 4 apex
        {{  0.0, -1.0,  0.0 }},  // 5 midside 0-1
        {{  1.0,  0.0,  0.0 }},  // 6 midside 1-2
        {{  0.0,  1.0,  0.0 }},  // 7 midside 2-3
        {{ -1.0,  0.0,  0.0 }},  // 8 midside 3-0
        {{ -0.5, -0.5,  0.5 }},  // 9 midside 0-4
        {{  0.5, -0.5,  0.5 }},  // 10 midside 1-4
        {{  0.5,  0.5,  0.5 }},  // 11 midside 2-4
        {{ -0.5,  0.5,  0.5 }}   // 12 midside 3-4
    }};

    /// Segment (edge) connectivity as constexpr.
    static constexpr std::array<std::array<uint32_t,3>,8> kSegments = {{
        {{0U,1U, 5U}}, {{1U,2U, 6U}}, {{2U,3U, 7U}}, {{3U,0U, 8U}},
        {{0U,4U, 9U}}, {{1U,4U,10U}}, {{2U,4U,11U}}, {{3U,4U,12U}}
    }};
    
    friend class IsoparametricQuadraticPyramid_Test;
};

} // namespace csmp

#endif // ISOPARAMETRIC_QUADRATIC_PYRAMID_H

