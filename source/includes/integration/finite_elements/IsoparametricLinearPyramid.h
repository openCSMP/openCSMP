#ifndef ISOPARAMETRIC_LINEAR_PYRAMID_H
#define ISOPARAMETRIC_LINEAR_PYRAMID_H

#include "FiniteElement.h"

#include <array>
#include <cstdint>
#include <vector>

namespace csmp {

/**
@brief Linear 5-node pyramid finite element.

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1998
@author (refactored) 2024

@section geometry Geometry

Node numbering:

         4 (apex)
        /|\
       / | \
      /  |  \
     3---+---2
     |   |   |
     0-------1   (base, CCW from outside)

Reference coordinates:
  Node 0: r=-1, s=-1, t=0
  Node 1: r=+1, s=-1, t=0
  Node 2: r=+1, s=+1, t=0
  Node 3: r=-1, s=+1, t=0
  Node 4: r= 0, s= 0, t=1  (apex)

@section shapefunctions Shape Functions

Rational Bedrosian/Zganski formulation — exact for the pyramid domain:

  N_i = 0.25 * ( (1-t) ± r ± s + r*s/(1-t) )   for i=0..3
  N_4 = t

These are the correct shape functions for the pyramid parametric domain
|r| <= 1-t, |s| <= 1-t, 0 <= t <= 1.
Reference: Zganski et al. 1996; Bedrosian 1992.

@section integration Integration

Supported integration point counts: 1, 5 (default), 8.

  1-point: centroid rule, exact for linear functions. Weight = 4/3.
  5-point: Stroud conical product rule, equal weights 4/15 each.
  8-point: full bi-linear rule via GenerateIntegrationPoints.

@section faces Faces

  0: {0,1,4}   triangular
  1: {1,2,4}   triangular
  2: {2,3,4}   triangular
  3: {0,4,3}   triangular
  4: {0,3,2,1} quadrilateral (base)
*/
class IsoparametricLinearPyramid : public FiniteElement
{
public:

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
    Constructs the pyramid finite element.
    @param integrationPoints  Number of Gauss points: 1, 5 (default), or 8.
    */
    explicit IsoparametricLinearPyramid( uint32_t integrationPoints = 5U );

    // -----------------------------------------------------------------------
    // Geometry
    // -----------------------------------------------------------------------

    double Volume()                                          override final;
    double AspectRatio()                                     override final;
    double InnerRadius()                                     override final;
    void   EdgeLengths( std::vector<double>& vec )           override final;

    // -----------------------------------------------------------------------
    // Topology
    // -----------------------------------------------------------------------

    void     CornerNodes( std::vector<uint32_t>& ids ) const override final;
    uint32_t CornerNodes() const noexcept override final { return 5U; }

    void CounterClockwiseNodes( std::vector<uint32_t>& ids ) const;

    void     MidSideNodes( std::vector<uint32_t>& ids ) const override final;
    uint32_t MidSideNodes() const noexcept override final { return 0U; }

    uint32_t NodesPerFace( uint32_t face ) const noexcept override final;

    CSMP_FEM_TYPE ElementTypeOfFace( uint32_t face ) const noexcept override final;
    CSMP_FEM_TYPE ElementTypeOfSegment( uint32_t ) const noexcept override final
        { return ISOPARAMETRIC_LINEAR_BAR; }

    void NodesOfSegment( uint32_t segm_id,
                         std::vector<uint32_t>& snids ) const override final;

    std::vector<uint32_t> NodesOfFace( uint32_t face_id )       const override final;
    std::vector<uint32_t> CornerNodesOfFace( uint32_t face_id ) const override final;
    std::vector<uint32_t> NodesConnectedTo( uint32_t node_id )  const override final;

    void UnitNormalToFace( uint32_t face,
                           std::vector<double>& unrml ) const override final;

    // -----------------------------------------------------------------------
    // Shape functions
    // -----------------------------------------------------------------------

    void N( std::vector<double>& N,
            const std::vector<double>& xyz )                 override final;
    void N_AtIntegrationPoint( uint32_t ip,
                               std::vector<double>& N )      override final;
    void N_AtBaryCenter( std::vector<double>& N )            override final;


    // -----------------------------------------------------------------------
    // Jacobian
    // -----------------------------------------------------------------------

    void JacobianAtIntegrationPoint( uint32_t ip )           override final;
    void JacobianAt( const std::vector<double>& rst )        override final;

    // -----------------------------------------------------------------------
    // Shape function derivatives
    // -----------------------------------------------------------------------

    double dN( DenseMatrix<DM_MIN>& dn,
               const std::vector<double>& xyz );
    void   dN( DenseMatrix<DM_MIN>& DN5 )                    override final;

    double dN_AtNode( DenseMatrix<DM_MIN>& M,
                      uint32_t node )                        override final;
    double dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M,
                                  uint32_t gauss_point )     override final;
    double dN_AtBarycenter( DenseMatrix<DM_MIN>& M )         override final;

    // --- Reference-coordinate shape functions and derivatives ---

    /**
    Rational Bedrosian/Zganski shape functions at (r,s,t).
    Singular at apex (t=1) — handled by the t!=1 guard.
    */
    void Nrst( double r, double s, double t,
               std::vector<double>& N ) const                override final;

    void dNr( double r, double s, double t,
              std::vector<double>& dNr ) const               override final;
    void dNs( double r, double s, double t,
              std::vector<double>& dNs ) const               override final;
    void dNt( double r, double s, double t,
              std::vector<double>& dNt ) const               override final;


    // -----------------------------------------------------------------------
    // Integration
    // -----------------------------------------------------------------------

    double WeightAtIntegrationPoint( uint32_t i ) const noexcept override final;

    void IntegrationPoint( uint32_t i,
                           std::vector<double>& xyz ) const  override final;

    void ExtrapolateIntegrationPointVariableToNodes(
             uint32_t nvars,
             const std::vector<double>& IVAR,
             std::vector<double>& NVAR ) const               override final;

    // -----------------------------------------------------------------------
    // Output
    // -----------------------------------------------------------------------

    void OutputNodeDataToVTK( const char* file_name,
                              const char* var_name,
                              DenseMatrix<DM_MIN>& DATA ) const override final;

    void ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const override final;

private:

    double VolumeOfTetra( uint32_t verticeIndex1, uint32_t verticeIndex2, uint32_t verticeIndex3, uint32_t verticeIndex4 );

    // -----------------------------------------------------------------------
    // Coordinate mapping
    // -----------------------------------------------------------------------

    void ParametricToPhysical( std::vector<double> &rst, std::vector<double> &xyz );
    void PhysicalToParametric( std::vector<double>& rst,const std::vector<double>& xyz );


    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /**
    Rational shape functions into a raw array — avoids vector allocation
    in tight loops.
    */
    void Nrst( double r, double s, double t, double* N ) const noexcept;

    /**
    Volume of tetrahedron defined by four node indices using the
    scalar triple product formula: V = |det[b-a, c-a, d-a]| / 6.
    */
    double VolumeOfTetra( uint32_t i, uint32_t j,
                          uint32_t k, uint32_t l ) const noexcept;

    double VolumeOfPyramid();

    /**
    Generates the 8-point integration rule for full bi-linear integration.
    */
    void GenerateIntegrationPoints( DenseMatrix<DM_MIN>& Ip,
                                    std::vector<double>& We );

    // -----------------------------------------------------------------------
    // Data members
    // -----------------------------------------------------------------------

    std::vector<double>  W;      ///< Integration weights [gpe]
    DenseMatrix<DM_MIN>  DN;     ///< Derivative matrix workspace
    DenseMatrix<DM_MIN>  NXYZ;   ///< Reference node coordinates [5 x 3]
    DenseMatrix<DM_MIN>  IP;     ///< Integration point coordinates [gpe x 3]

    /// Reference node coordinates as constexpr for zero-overhead access.
    static constexpr std::array<std::array<double,3>,5> kRefCoords = {{
        {{ -1.0, -1.0,  0.0 }},
        {{  1.0, -1.0,  0.0 }},
        {{  1.0,  1.0,  0.0 }},
        {{ -1.0,  1.0,  0.0 }},
        {{  0.0,  0.0,  1.0 }}
    }};

    /// Segment (edge) connectivity as constexpr.
    static constexpr std::array<std::array<uint32_t,2>,8> kSegments = {{
        {{0U,1U}}, {{1U,2U}}, {{2U,3U}}, {{3U,0U}},   // base edges
        {{0U,4U}}, {{1U,4U}}, {{2U,4U}}, {{3U,4U}}    // lateral edges to apex
    }};
};

} // namespace csmp

#endif // ISOPARAMETRIC_LINEAR_PYRAMID_H

