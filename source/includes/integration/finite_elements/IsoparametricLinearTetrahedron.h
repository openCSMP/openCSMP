// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ISOPARAMETRIC_LINEAR_TETRAHEDRON_H
#define ISOPARAMETRIC_LINEAR_TETRAHEDRON_H

#include "FiniteElement.h"
#include <vector>
#include <cmath>

namespace csmp {

/// 4-noded tetrahedron with a choice of Gaussian quadrature points.
class IsoparametricLinearTetrahedron final : public FiniteElement {
public:
    explicit IsoparametricLinearTetrahedron(uint32_t integrationPoints = 4 /* 1 or 4 */);

    // ----------------------------------------------------------------
    // Geometry
    // ----------------------------------------------------------------
    double   Volume()                                          override final;
    double   AspectRatio()                                     override final;
    double   InnerRadius()                                     override final;
    void     EdgeLengths(std::vector<double>& vec)             override final;

    void     CornerNodes(std::vector<uint32_t>& ids) const     override final;
    uint32_t CornerNodes()                           const noexcept override final { return 4U; }

    void     MidSideNodes(std::vector<uint32_t>& ids) const    override final;
    uint32_t MidSideNodes()                           const noexcept override final { return 0U; }

    void     CounterClockwiseNodes(std::vector<uint32_t>& ids) const;

    // ----------------------------------------------------------------
    // Topology
    // ----------------------------------------------------------------
    CSMP_FEM_TYPE ElementTypeOfFace(uint32_t face)          const noexcept override final;
    CSMP_FEM_TYPE ElementTypeOfSegment(uint32_t /*segment*/) const noexcept override final
    {
        return ISOPARAMETRIC_LINEAR_BAR;
    }

    void                  NodesOfSegment(uint32_t segm_id,
                                         std::vector<uint32_t>& snids) const override final;
    std::vector<uint32_t> NodesOfFace(uint32_t face_id)               const override final;
    std::vector<uint32_t> CornerNodesOfFace(uint32_t face_id)         const override final;
    std::vector<uint32_t> NodesConnectedTo(uint32_t node_id)          const override final;

    void UnitNormalToFace(uint32_t face,
                          std::vector<double>& unrml)                  const override final;

    // ----------------------------------------------------------------
    // Shape functions — local coordinates (inlined: trivially small,
    // called at every integration point in every assembly loop)
    // ----------------------------------------------------------------
    void Nrst(double r, double s, double t, std::vector<double>& N) const override final
    {
        N.resize(npe);
        N[0] = 1.0 - r - s - t;
        N[1] = r;
        N[2] = s;
        N[3] = t;
    }

    void Nrst(double r, double s, double t, double* N) const noexcept
    {
        N[0] = 1.0 - r - s - t;
        N[1] = r;
        N[2] = s;
        N[3] = t;
    }

    void dNr(double /*r*/, double /*s*/, double /*t*/,
             std::vector<double>& DNR) const override
    {
        DNR.resize(npe);
        DNR[0] = -1.0; DNR[1] = 1.0; DNR[2] = 0.0; DNR[3] = 0.0;
    }

    void dNs(double /*r*/, double /*s*/, double /*t*/,
             std::vector<double>& DNS) const override
    {
        DNS.resize(npe);
        DNS[0] = -1.0; DNS[1] = 0.0; DNS[2] = 1.0; DNS[3] = 0.0;
    }

    void dNt(double /*r*/, double /*s*/, double /*t*/,
             std::vector<double>& DNT) const override
    {
        DNT.resize(npe);
        DNT[0] = -1.0; DNT[1] = 0.0; DNT[2] = 0.0; DNT[3] = 1.0;
    }

    // ----------------------------------------------------------------
    // Shape functions — physical / integration-point interface
    // ----------------------------------------------------------------
    void   N(std::vector<double>& N, const std::vector<double>& xyz)   override final;
    void   N(std::vector<double>& N, uint32_t& iterations,
             double& distance, const std::vector<double>& xyz);

    void   N_AtIntegrationPoint(uint32_t ip, std::vector<double>& N)   override final;
    void   N_AtBaryCenter(std::vector<double>& N)                       override final;

    void   JacobianAtIntegrationPoint(uint32_t ip)                      override final;
    void   JacobianAt(const std::vector<double>& rst)                   override final;

    // ----------------------------------------------------------------
    // Shape-function derivatives — global coordinates
    // ----------------------------------------------------------------
    double dN(DenseMatrix<DM_MIN>& DN, const std::vector<double>& xyz);
    void   dN(DenseMatrix<DM_MIN>& DN4)                                  override final;
    double dN_AtNode(DenseMatrix<DM_MIN>& DN, uint32_t node)             override final;
    double dN_AtIntegrationPoint(DenseMatrix<DM_MIN>& M,
                                 uint32_t gauss_point)                   override final;
    double dN_AtBarycenter(DenseMatrix<DM_MIN>& M)                       override final;

    // ----------------------------------------------------------------
    // Integration
    // ----------------------------------------------------------------
    double WeightAtIntegrationPoint(uint32_t i) const noexcept override
    {
        return W[i];
    }

    void IntegrationPoint(uint32_t i,
                          std::vector<double>& xyz) const override; // global coords

    // ----------------------------------------------------------------
    // Post-processing
    // ----------------------------------------------------------------
    void ExtrapolateIntegrationPointVariableToNodes(
             uint32_t nvars,
             const std::vector<double>& IVAR,
             std::vector<double>& NVAR) const override;

    void OutputNodeDataToVTK(const char* file_name,
                             const char* var_name,
                             DenseMatrix<DM_MIN>& DATA) const override;

    void ReferenceCoordinates(DenseMatrix<DM_MIN>& matCoords) const override;

private:
    std::vector<double> W;
    DenseMatrix<DM_MIN> DN, NXYZ, IP;
    double   accDistance;
    uint32_t totIterations;
    uint32_t nonConvergenceOfProjections;
    uint32_t projectionCalledNTimes;

    void ParametricToPhysical(std::vector<double>& rst,
                              std::vector<double>& xyz);
    void PhysicalToParametric(std::vector<double>& rst,
                              const std::vector<double>& xyz);

    inline uint32_t n(int32_t i, int32_t a) const noexcept {
        return static_cast<uint32_t>( (i + a >= 4) ? (i + a - 4) : (i + a) );
    }
};

/**
 * @class IsoparametricLinearTetrahedron
 * @file  finite_elements/IsoparametricLinearTetrahedron.h
 * @date  1998
 * @author S.K. Matthai
 * @author Stephen G. Roberts
 */

} // namespace csmp

#endif // ISOPARAMETRIC_LINEAR_TETRAHEDRON_H

