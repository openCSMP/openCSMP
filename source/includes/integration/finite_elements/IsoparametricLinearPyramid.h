#ifndef ISOPARAMETRIC_LINEAR_PYRAMID_H
#define ISOPARAMETRIC_LINEAR_PYRAMID_H

#include "FiniteElement.h"

namespace csmp {

/**
 * @class IsoparametricLinearPyramid
 * @brief Linear 5-node pyramid finite element.
 *
 * Node numbering:
 *
 *        4 (apex)
 *       /|\
 *      / | \
 *     /  |  \
 *    3---+---2
 *    |   |   |
 *    0-------1  (base, CCW from outside)
 *
 * Reference coordinates:
 *   Node 0: r=-1, s=-1, t=0
 *   Node 1: r=+1, s=-1, t=0
 *   Node 2: r=+1, s=+1, t=0
 *   Node 3: r=-1, s=+1, t=0
 *   Node 4: r= 0, s= 0, t=1  (apex)
 *
 * Shape functions (serendipity, no singularity at apex):
 *   N_i = 0.25*(1±r)*(1±s)*(1-t)  for i=0..3
 *   N_4 = t
 *
 * Faces:
 *   0: {0,1,4}   triangular
 *   1: {1,2,4}   triangular
 *   2: {2,3,4}   triangular
 *   3: {3,0,4}   triangular
 *   4: {0,3,2,1} quadrilateral (base)
 *
 * Integration points: 1, 5, or 8
 *
 * @author S.K. Matthaei (1998)
 * @author R. Manasipov (2013)
 * @author Revised 2024
 */
class IsoparametricLinearPyramid final : public FiniteElement {

  public:

    /// @param integrationPoints  Number of Gauss points: 1, 5, or 8
    explicit IsoparametricLinearPyramid( uint32_t integrationPoints = 5 );

    // ------------------------------------------------------------------
    // Geometry
    // ------------------------------------------------------------------

    double    Volume()                                        override final;
    double    AspectRatio()                                   override final;
    double    InnerRadius()                                   override final;
    void      EdgeLengths( std::vector<double>& vec )         override final;

    // ------------------------------------------------------------------
    // Topology — noexcept: return fixed data, no allocation, no throw
    // ------------------------------------------------------------------

    void      CornerNodes( std::vector<uint32_t>& ids ) const noexcept override final;
    uint32_t  CornerNodes() const noexcept override final { return 5U; }

    void      MidSideNodes( std::vector<uint32_t>& ids ) const override final;
    uint32_t  MidSideNodes() const noexcept override final { return 0U; }

    void      CounterClockwiseNodes( std::vector<uint32_t>& ids ) const noexcept;

    uint32_t  NodesPerFace( uint32_t face ) const noexcept override final;

    CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const override final;

    CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t ) const noexcept override final { return ISOPARAMETRIC_LINEAR_BAR; }

    void      NodesOfSegment( uint32_t segm_id,
                               std::vector<uint32_t>& snids ) const noexcept override final;

    std::vector<uint32_t>  NodesOfFace( uint32_t face_id )       const override final;
    std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const override final;
    std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id )  const override final;

    void      UnitNormalToFace( uint32_t face,
                                 std::vector<double>& unrml )    const override final;

    // ------------------------------------------------------------------
    // Shape functions
    // ------------------------------------------------------------------

    void      N( std::vector<double>& N,
                  const std::vector<double>& xyz )               override final;

    void      N_AtIntegrationPoint( uint32_t ip,
                                     std::vector<double>& N )    override final;

    void      N_AtBaryCenter( std::vector<double>& N )           override final;

    // ------------------------------------------------------------------
    // Jacobian
    // ------------------------------------------------------------------

    void      JacobianAtIntegrationPoint( uint32_t ip )          override final;
    void      JacobianAt( const std::vector<double>& rst )       override final;

    // ------------------------------------------------------------------
    // Shape function derivatives
    // ------------------------------------------------------------------

    double    dN( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz );
    void      dN( DenseMatrix<DM_MIN>& DN5 )                     override final;
    double    dN_AtNode( DenseMatrix<DM_MIN>& M,
                          uint32_t node )                        override final;
    double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M,
                                      uint32_t gauss_point )     override final;
    double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M )          override final;

    // ------------------------------------------------------------------
    // Reference-coordinate shape functions and derivatives
    // noexcept: pure arithmetic on inputs, no allocation, no throw
    // ------------------------------------------------------------------

    void      Nrst( double r, double s, double t,
                     std::vector<double>& N )    const noexcept override final;
    void      Nrst( double r, double s, double t, double* N ) const noexcept;
    void      dNr(  double r, double s, double t,
                     std::vector<double>& dNr )  const noexcept override final;
    void      dNs(  double r, double s, double t,
                     std::vector<double>& dNs )  const noexcept override final;
    void      dNt(  double r, double s, double t,
                     std::vector<double>& dNt )  const noexcept override final;

    // ------------------------------------------------------------------
    // Integration
    // ------------------------------------------------------------------

    /// noexcept: simple array lookup, no allocation, no throw
    double    WeightAtIntegrationPoint( uint32_t i ) const noexcept override final;

    void      IntegrationPoint( uint32_t i,
                                 std::vector<double>& xyz ) const override final;

    void      ExtrapolateIntegrationPointVariableToNodes(
                  uint32_t nvars,
                  const std::vector<double>& IVAR,
                  std::vector<double>& NVAR )            const override final;

    // ------------------------------------------------------------------
    // Output
    // ------------------------------------------------------------------

    void      OutputNodeDataToVTK( const char* file_name,
                                    const char* var_name,
                                    DenseMatrix<DM_MIN>& DATA ) const override final;

    void      ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const noexcept override final;

  private:

    std::vector<double>  W;
    DenseMatrix<DM_MIN>  DN, NXYZ, IP;

    /// Correct tetrahedron volume: V = |det[b-a, c-a, d-a]| / 6
    double VolumeOfTetra( uint32_t i, uint32_t j,
                           uint32_t k, uint32_t l ) const noexcept;

    double VolumeOfPyramid();

    void GenerateIntegrationPoints( DenseMatrix<DM_MIN>& Ip,
                                     std::vector<double>& We );

    void ParametricToPhysical( std::vector<double>& rst,
                                std::vector<double>& xyz );

    void PhysicalToParametric( std::vector<double>& rst,
                                const std::vector<double>& xyz );
};

} // namespace csmp

#endif // ISOPARAMETRIC_LINEAR_PYRAMID_H
