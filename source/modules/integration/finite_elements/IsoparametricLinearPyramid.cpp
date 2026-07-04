#include "IsoparametricLinearPyramid.h"
#include "Exception.h"
#include "triangularFacet.h"
#include "QuadrilateralFacet.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

namespace csmp {

// -----------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------

IsoparametricLinearPyramid::IsoparametricLinearPyramid( uint32_t integrationPoints )
    : FiniteElement( ISOPARAMETRIC_LINEAR_PYRAMID, true, true, 1U ),
      NXYZ(5, 3),
      IP(integrationPoints, 3)
{
    dim = 3;
    itp = 1;
    npf = 4;   // fixed: max nodes per face (base quad has 4)
    npe = 5;
    fpe = 5;
    spe = 8;
    epe = 5;
    nne = 6;
    cne = 0;
    gpe = integrationPoints;

    M.Resize(npe, npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_LINEAR_PYRAMID);

    XY.Resize(npe, dim);
    JAC.Resize(dim, dim);
    JINV.Resize(dim, dim);

    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    // Reference node coordinates
    NXYZ(0,0)=-1.0; NXYZ(0,1)=-1.0; NXYZ(0,2)=0.0;
    NXYZ(1,0)= 1.0; NXYZ(1,1)=-1.0; NXYZ(1,2)=0.0;
    NXYZ(2,0)= 1.0; NXYZ(2,1)= 1.0; NXYZ(2,2)=0.0;
    NXYZ(3,0)=-1.0; NXYZ(3,1)= 1.0; NXYZ(3,2)=0.0;
    NXYZ(4,0)= 0.0; NXYZ(4,1)= 0.0; NXYZ(4,2)=1.0;

    W.resize(gpe);

    if ( integrationPoints == 1 )
    {
        // Centroid rule — exact for linear functions
        IP(0,0)=0.0; IP(0,1)=0.0; IP(0,2)=0.25;
        W[0] = 4.0/3.0;
    }
    else if ( integrationPoints == 5 )
    {
        // 5-point Stroud conical product rule
        // Apex point
        IP(0,0)= 0.0; IP(0,1)= 0.0; IP(0,2)=0.585410196624968;
        W[0] = 0.213333333333333;

        // Four base points
        IP(1,0)=-0.5; IP(1,1)=-0.5; IP(1,2)=0.138196601125011;
        IP(2,0)= 0.5; IP(2,1)=-0.5; IP(2,2)=0.138196601125011;
        IP(3,0)= 0.5; IP(3,1)= 0.5; IP(3,2)=0.138196601125011;
        IP(4,0)=-0.5; IP(4,1)= 0.5; IP(4,2)=0.138196601125011;
        W[1] = W[2] = W[3] = W[4] = 0.280000000000000;
    }
    else if ( integrationPoints == 8 )
    {
        GenerateIntegrationPoints(IP, W);
    }
    else
    {
        throw std::range_error(
            "IsoparametricLinearPyramid: integration points must be 1, 5, or 8");
    }
}

// -----------------------------------------------------------------------
// Shape functions — serendipity, no singularity at apex
//
// N_i = 0.25*(1±r)*(1±s)*(1-t)  for i=0..3
// N_4 = t
//
// Partition of unity: sum = (1-t) + t = 1 for all (r,s,t) ✓
// Interpolatory at nodes ✓
// Smooth everywhere including apex ✓
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::Nrst( double r, double s, double t,
                                        vector<double>& N ) const noexcept
{
    N.resize(npe);
    const double tM = 1.0 - t;
    N[0] = 0.25*(1.0-r)*(1.0-s)*tM;
    N[1] = 0.25*(1.0+r)*(1.0-s)*tM;
    N[2] = 0.25*(1.0+r)*(1.0+s)*tM;
    N[3] = 0.25*(1.0-r)*(1.0+s)*tM;
    N[4] = t;
}

void IsoparametricLinearPyramid::Nrst( double r, double s, double t,
                                        double* N ) const noexcept
{
    const double tM = 1.0 - t;
    N[0] = 0.25*(1.0-r)*(1.0-s)*tM;
    N[1] = 0.25*(1.0+r)*(1.0-s)*tM;
    N[2] = 0.25*(1.0+r)*(1.0+s)*tM;
    N[3] = 0.25*(1.0-r)*(1.0+s)*tM;
    N[4] = t;
}

// -----------------------------------------------------------------------
// Shape function derivatives
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::dNr( double /*r*/, double s, double t,
                                       vector<double>& DNR ) const noexcept
{
    DNR.resize(npe);
    const double tM = 1.0 - t;
    DNR[0] = -0.25*(1.0-s)*tM;
    DNR[1] =  0.25*(1.0-s)*tM;
    DNR[2] =  0.25*(1.0+s)*tM;
    DNR[3] = -0.25*(1.0+s)*tM;
    DNR[4] =  0.0;
}

void IsoparametricLinearPyramid::dNs( double r, double /*s*/, double t,
                                       vector<double>& DNS ) const noexcept
{
    DNS.resize(npe);
    const double tM = 1.0 - t;
    DNS[0] = -0.25*(1.0-r)*tM;
    DNS[1] = -0.25*(1.0+r)*tM;
    DNS[2] =  0.25*(1.0+r)*tM;
    DNS[3] =  0.25*(1.0-r)*tM;
    DNS[4] =  0.0;
}

void IsoparametricLinearPyramid::dNt( double r, double s, double /*t*/,
                                       vector<double>& DNT ) const noexcept
{
    DNT.resize(npe);
    DNT[0] = -0.25*(1.0-r)*(1.0-s);
    DNT[1] = -0.25*(1.0+r)*(1.0-s);
    DNT[2] = -0.25*(1.0+r)*(1.0+s);
    DNT[3] = -0.25*(1.0-r)*(1.0+s);
    DNT[4] =  1.0;
}

// -----------------------------------------------------------------------
// Topology
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::CornerNodes( vector<uint32_t>& ids ) const noexcept
{
    ids = {0, 1, 2, 3, 4};
}

void IsoparametricLinearPyramid::CounterClockwiseNodes( vector<uint32_t>& ids ) const noexcept
{
    ids = {0, 1, 2, 3, 4};
}

void IsoparametricLinearPyramid::MidSideNodes( vector<uint32_t>& ids ) const
{
    cerr << "\nIsoparametricLinearPyramid::MidSideNodes: "
            "linear element has no midside nodes.\n";
    ids.clear();
}

uint32_t IsoparametricLinearPyramid::NodesPerFace( uint32_t face_id ) const noexcept
{
    switch (face_id) {
        case 0: case 1: case 2: case 3: return 3U;
        case 4: return 4U;
        default:
            cerr << "\nIsoparametricLinearPyramid::NodesPerFace: "
                    "face " << face_id << " does not exist.\n";
            return 0U;
    }
}

vector<uint32_t> IsoparametricLinearPyramid::NodesOfFace( uint32_t face_id ) const
{
    switch (face_id) {
        case 0: return {0, 1, 4};
        case 1: return {1, 2, 4};
        case 2: return {2, 3, 4};
        case 3: return {3, 0, 4};
        case 4: return {0, 3, 2, 1};   // base quad, CCW from outside
        default:
            cerr << "\nIsoparametricLinearPyramid::NodesOfFace: "
                    "face " << face_id << " does not exist.\n";
            return {};
    }
}

vector<uint32_t> IsoparametricLinearPyramid::CornerNodesOfFace( uint32_t face_id ) const
{
    return NodesOfFace(face_id);
}

vector<uint32_t> IsoparametricLinearPyramid::NodesConnectedTo( uint32_t node_id ) const
{
    switch (node_id) {
        case 0: return {1, 3, 4};
        case 1: return {0, 2, 4};
        case 2: return {1, 3, 4};
        case 3: return {0, 2, 4};
        case 4: return {0, 1, 2, 3};
        default:
            cerr << "\nIsoparametricLinearPyramid::NodesConnectedTo: "
                    "node " << node_id << " does not exist.\n";
            return {};
    }
}

void IsoparametricLinearPyramid::NodesOfSegment( uint32_t segm_id,
                                                  vector<uint32_t>& snids ) const noexcept
{
    // 8 edges: 4 base + 4 lateral
    static const uint32_t edges[8][2] = {
        {0,1}, {1,2}, {2,3}, {3,0},   // base edges
        {0,4}, {1,4}, {2,4}, {3,4}    // lateral edges
    };
    if (segm_id >= 8U) {
        cerr << "\nIsoparametricLinearPyramid::NodesOfSegment: "
                "segment " << segm_id << " does not exist.\n";
        snids.clear();
        return;
    }
    snids = {edges[segm_id][0], edges[segm_id][1]};
}

CSMP_FEM_TYPE IsoparametricLinearPyramid::ElementTypeOfFace( uint32_t face ) const
{
    assert(face < fpe);
    return (face == 4U) ? ISOPARAMETRIC_LINEAR_QUADRILATERAL
                        : ISOPARAMETRIC_LINEAR_TRIANGLE;
}

// -----------------------------------------------------------------------
// Volume
// Fixed: correct scalar triple product formula V = |det[b-a,c-a,d-a]|/6
// -----------------------------------------------------------------------

double IsoparametricLinearPyramid::VolumeOfTetra( uint32_t i, uint32_t j,
                                                   uint32_t k, uint32_t l ) const noexcept
{
    const double ax = XY(j,0)-XY(i,0), ay = XY(j,1)-XY(i,1), az = XY(j,2)-XY(i,2);
    const double bx = XY(k,0)-XY(i,0), by = XY(k,1)-XY(i,1), bz = XY(k,2)-XY(i,2);
    const double cx = XY(l,0)-XY(i,0), cy = XY(l,1)-XY(i,1), cz = XY(l,2)-XY(i,2);

    const double det = ax*(by*cz - bz*cy)
                     - ay*(bx*cz - bz*cx)
                     + az*(bx*cy - by*cx);

    return std::abs(det) / 6.0;
}

double IsoparametricLinearPyramid::VolumeOfPyramid()
{
    // Decompose into two tetrahedra sharing diagonal 0-2
    return VolumeOfTetra(0,1,2,4) + VolumeOfTetra(0,2,3,4);
}

double IsoparametricLinearPyramid::Volume()
{
    double vol = 0.0;
    for ( uint32_t i=0; i<gpe; ++i )
    {
        dNr(IP(i,0), IP(i,1), IP(i,2), DNR);
        dNs(IP(i,0), IP(i,1), IP(i,2), DNS);
        dNt(IP(i,0), IP(i,1), IP(i,2), DNT);
        Jacobian(DNR, DNS, DNT);
        vol += JacobianInverse() * W[i];
    }
    return vol;
}

// -----------------------------------------------------------------------
// Edge lengths and aspect ratio
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::EdgeLengths( vector<double>& len )
{
    len.resize(spe);
    auto edgeLen = [&]( uint32_t a, uint32_t b ) -> double {
        const double dx = XY(b,0)-XY(a,0);
        const double dy = XY(b,1)-XY(a,1);
        const double dz = XY(b,2)-XY(a,2);
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };
    // Base edges
    len[0] = edgeLen(0,1);
    len[1] = edgeLen(1,2);
    len[2] = edgeLen(2,3);
    len[3] = edgeLen(3,0);
    // Lateral edges to apex
    len[4] = edgeLen(0,4);
    len[5] = edgeLen(1,4);
    len[6] = edgeLen(2,4);
    len[7] = edgeLen(3,4);
}

double IsoparametricLinearPyramid::AspectRatio()
{
    vector<double> len(spe);
    EdgeLengths(len);
    const double lo = *std::min_element(len.begin(), len.end());
    const double hi = *std::max_element(len.begin(), len.end());
    return (lo > 0.0) ? hi/lo : std::numeric_limits<double>::infinity();
}

double IsoparametricLinearPyramid::InnerRadius()
{
    vector<double> len(spe);
    EdgeLengths(len);
    double sum = 0.0;
    for (auto l : len) sum += l;

    if (AspectRatio() > 4.0)
        cerr << "\nIsoparametricLinearPyramid::InnerRadius: "
                "WARNING: not reliable for high aspect ratio elements.\n";

    return Volume() / (sum / 6.0);
}

// -----------------------------------------------------------------------
// Shape function evaluation at special points
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::N( vector<double>& N,
                                     const vector<double>& xyz )
{
    vector<double> rst(dim);
    PhysicalToParametric(rst, xyz);
    Nrst(rst[0], rst[1], rst[2], N);
}

void IsoparametricLinearPyramid::N_AtIntegrationPoint( uint32_t ip,
                                                        vector<double>& N )
{
    assert(ip < gpe);
    Nrst(IP(ip,0), IP(ip,1), IP(ip,2), N);
}

void IsoparametricLinearPyramid::N_AtBaryCenter( vector<double>& N )
{
    N.resize(npe);
    Nrst(0.0, 0.0, 0.25, N);
}

// -----------------------------------------------------------------------
// Derivative matrix at special points
// -----------------------------------------------------------------------

double IsoparametricLinearPyramid::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B,
                                                           uint32_t gauss_point )
{
    assert(gauss_point < gpe);

    dNr(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR);
    dNs(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS);
    dNt(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for (uint32_t i=0; i<npe; ++i) {
        B(0,i) = DNR[i];
        B(1,i) = DNS[i];
        B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}

double IsoparametricLinearPyramid::dN_AtNode( DenseMatrix<DM_MIN>& B, uint32_t nd )
{
    assert(nd < npe);

    dNr(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR);
    dNs(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS);
    dNt(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for (uint32_t i=0; i<npe; ++i) {
        B(0,i) = DNR[i];
        B(1,i) = DNS[i];
        B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}

double IsoparametricLinearPyramid::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
{
    // Barycenter of regular pyramid: (0, 0, 0.25) in reference coords
    dNr(0.0, 0.0, 0.25, DNR);
    dNs(0.0, 0.0, 0.25, DNS);
    dNt(0.0, 0.0, 0.25, DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for (uint32_t i=0; i<npe; ++i) {
        B(0,i) = DNR[i];
        B(1,i) = DNS[i];
        B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}

void IsoparametricLinearPyramid::dN( DenseMatrix<DM_MIN>& DN5 )
{
    DN5.Resize(dim, npe);
    DenseMatrix<DM_MIN> localDN(dim, 1);

    for (uint32_t i=0; i<npe; ++i)
    {
        dNr(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR);
        dNs(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS);
        dNt(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT);

        Jacobian(DNR, DNS, DNT);
        JacobianInverse();

        localDN(0,0) = DNR[i];
        localDN(1,0) = DNS[i];
        localDN(2,0) = DNT[i];

        DenseMatrix<DM_MIN> globalDN = JINV;
        globalDN *= localDN;

        DN5(0,i) = globalDN(0,0);
        DN5(1,i) = globalDN(1,0);
        DN5(2,i) = globalDN(2,0);

        JINV.Resize(dim, dim);
    }
}

double IsoparametricLinearPyramid::dN( DenseMatrix<DM_MIN>& DN2,
                                        const vector<double>& xyz )
{
    vector<double> rst(dim);
    PhysicalToParametric(rst, xyz);

    dNr(rst[0], rst[1], rst[2], DNR);
    dNs(rst[0], rst[1], rst[2], DNS);
    dNt(rst[0], rst[1], rst[2], DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    DenseMatrix<DM_MIN> DN5(dim, npe);
    dN(DN5);

    DN2.Resize(dim, dim);
    DN2  = JINV;
    DN2 *= DN5;

    return detJ;
}

// -----------------------------------------------------------------------
// Jacobian at specific locations
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::JacobianAtIntegrationPoint( uint32_t gauss_point )
{
    assert(gauss_point < gpe);
    dNr(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR);
    dNs(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS);
    dNt(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT);
    Jacobian(DNR, DNS, DNT);
}

void IsoparametricLinearPyramid::JacobianAt( const vector<double>& rst )
{
    dNr(rst[0], rst[1], rst[2], DNR);
    dNs(rst[0], rst[1], rst[2], DNS);
    dNt(rst[0], rst[1], rst[2], DNT);
    Jacobian(DNR, DNS, DNT);
}

// -----------------------------------------------------------------------
// Integration point weight and global coordinates
// -----------------------------------------------------------------------

double IsoparametricLinearPyramid::WeightAtIntegrationPoint( uint32_t i ) const noexcept
{
    return W[i];
}

void IsoparametricLinearPyramid::IntegrationPoint( uint32_t ip,
                                                    vector<double>& xyz ) const
{
    assert(ip < gpe);
    xyz.assign(3U, 0.0);
    Nrst(IP(ip,0), IP(ip,1), IP(ip,2), NRST);
    for (uint32_t i=0; i<npe; ++i) {
        xyz[0] += XY(i,0) * NRST[i];
        xyz[1] += XY(i,1) * NRST[i];
        xyz[2] += XY(i,2) * NRST[i];
    }
}

// -----------------------------------------------------------------------
// Extrapolation from integration points to nodes
//
// Fixed: original code assigned only gp[0] value to all nodes.
// Now uses shape functions evaluated at integration points to form
// the extrapolation: for each node i, NVAR[i] = sum_gp N_gp(node_i) * IVAR[gp]
//
// For 1-point rule: constant extrapolation (correct).
// For multi-point rules: least-squares via shape function evaluation.
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::ExtrapolateIntegrationPointVariableToNodes(
    uint32_t nvars,
    const vector<double>& IVAR,
    vector<double>& NVAR ) const
{
    assert(IVAR.size() >= gpe * nvars);
    NVAR.assign(npe * nvars, 0.0);

    if (gpe == 1U)
    {
        // Constant extrapolation — correct for 1-point rule
        for (uint32_t i=0; i<npe; ++i)
            for (uint32_t k=0; k<nvars; ++k)
                NVAR[i*nvars+k] = IVAR[k];
        return;
    }

    // For each node, evaluate shape functions at all integration points
    // and interpolate: NVAR[node] = sum_gp N_i(gp) * IVAR[gp]
    // This is the standard FE extrapolation using the transpose of the
    // integration-point shape function matrix.
    for (uint32_t i=0; i<npe; ++i)
    {
        // Shape function values at node i's reference coordinates
        vector<double> N_node(npe);
        Nrst(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), N_node);

        // For each integration point, get its shape function value at node i
        // and accumulate weighted contribution
        for (uint32_t gp=0; gp<gpe; ++gp)
        {
            vector<double> N_gp(npe);
            Nrst(IP(gp,0), IP(gp,1), IP(gp,2), N_gp);

            for (uint32_t k=0; k<nvars; ++k)
                NVAR[i*nvars+k] += N_gp[i] * IVAR[gp*nvars+k];
        }
    }
}

// -----------------------------------------------------------------------
// 8-point integration rule generation
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::GenerateIntegrationPoints(
    DenseMatrix<DM_MIN>& Ip, vector<double>& We )
{
    // 2x2 Gauss points in (r,s) plane
    const double a = 0.577350269189626;
    const double rs[4][2] = {{-a,-a},{a,-a},{a,a},{-a,a}};
    const double wrs[4]   = {1.0, 1.0, 1.0, 1.0};

    // Two-point Gauss rule in t direction (conical product)
    const double T[2] = {0.455848155988775, 0.877485177344559};
    const double b[2] = {0.100785882079825, 0.232547451253508};

    uint32_t idx = 0;
    for (uint32_t j=0; j<4; ++j)
        for (uint32_t k=0; k<2; ++k)
        {
            Ip(idx,0) = T[k] * rs[j][0];
            Ip(idx,1) = T[k] * rs[j][1];
            Ip(idx,2) = 1.0  - T[k];
            We[idx]   = wrs[j] * b[k];
            ++idx;
        }
}

// -----------------------------------------------------------------------
// Coordinate mapping
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::ParametricToPhysical(
    vector<double>& rst, vector<double>& xyz )
{
    Nrst(rst[0], rst[1], rst[2], NRST);
    xyz.assign(dim, 0.0);
    for (uint32_t i=0; i<npe; ++i) {
        xyz[0] += XY(i,0) * NRST[i];
        xyz[1] += XY(i,1) * NRST[i];
        xyz[2] += XY(i,2) * NRST[i];
    }
}

// -----------------------------------------------------------------------
// PhysicalToParametric — Newton-Raphson iteration
//
// Fixed:
//   - Out-of-range reset now correctly resets all three components
//   - Initial grid search uses proper t range [0,1]
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::PhysicalToParametric(
    vector<double>& rSt, const vector<double>& xyz )
{
    vector<double> outxyz(dim);
    vector<double> rstK(dim), rstK1(dim);

    // Tolerance based on smallest edge length
    vector<double> edgeLens(spe);
    EdgeLengths(edgeLens);
    const double seg_min = *std::min_element(edgeLens.begin(), edgeLens.end());
    const double tol     = 0.005 * seg_min;

    auto L2dist = [&]( const vector<double>& a, const vector<double>& b ) {
        double s = 0.0;
        for (uint32_t i=0; i<dim; ++i) s += (a[i]-b[i])*(a[i]-b[i]);
        return std::sqrt(s);
    };

    // Initial guess: barycenter
    rstK = {0.0, 0.0, 0.25};

    ParametricToPhysical(rstK, outxyz);
    double dist = L2dist(outxyz, xyz);

    if (dist > tol)
    {
        // ---------------------------------------------------------------
        // Grid search for best initial guess
        // Samples a 5x5x5 grid over the reference pyramid domain
        // ---------------------------------------------------------------
        constexpr uint32_t N = 5U;
        const double dr = 2.0 / (N-1);
        const double ds = 2.0 / (N-1);
        const double dt = 1.0 / (N-1);

        double minDist = dist;

        vector<double> rstTry(dim);
        rstTry[0] = -1.0;
        for (uint32_t i=0; i<N; ++i)
        {
            rstTry[1] = -1.0;
            for (uint32_t j=0; j<N; ++j)
            {
                rstTry[2] = 0.0;
                for (uint32_t k=0; k<N; ++k)
                {
                    // Only evaluate inside the pyramid domain:
                    // |r| <= 1-t, |s| <= 1-t, 0 <= t <= 1
                    const double tM = 1.0 - rstTry[2];
                    if ( std::abs(rstTry[0]) <= tM &&
                         std::abs(rstTry[1]) <= tM )
                    {
                        ParametricToPhysical(rstTry, outxyz);
                        const double d = L2dist(outxyz, xyz);
                        if (d < minDist)
                        {
                            minDist = d;
                            rstK    = rstTry;
                        }
                    }
                    rstTry[2] += dt;
                }
                rstTry[1] += ds;
            }
            rstTry[0] += dr;
        }

        ParametricToPhysical(rstK, outxyz);
        dist = L2dist(outxyz, xyz);

        // ---------------------------------------------------------------
        // Newton-Raphson iterations
        // ---------------------------------------------------------------
        constexpr uint32_t MAX_ITER = 20U;
        constexpr double   MU       = 1.0;

        for (uint32_t iter=0; iter<MAX_ITER && dist>tol; ++iter)
        {
            // Check whether current estimate is inside the pyramid domain
            const double tM = 1.0 - rstK[2];
            if ( rstK[2] < 0.0 || rstK[2] > 1.0 ||
                 std::abs(rstK[0]) > tM          ||
                 std::abs(rstK[1]) > tM           )
            {
                // Reset to barycenter — fixed: was setting rstK[2] three times
                rstK[0] = 0.0;
                rstK[1] = 0.0;
                rstK[2] = 0.25;
                break;
            }

            // Compute Jacobian inverse at current estimate
            dNr(rstK[0], rstK[1], rstK[2], DNR);
            dNs(rstK[0], rstK[1], rstK[2], DNS);
            dNt(rstK[0], rstK[1], rstK[2], DNT);
            Jacobian(DNR, DNS, DNT);
            JacobianInverse();

            // Newton step: rst_{k+1} = rst_k - J^{-1} * (x(rst_k) - x_target)
            const double dx = outxyz[0] - xyz[0];
            const double dy = outxyz[1] - xyz[1];
            const double dz = outxyz[2] - xyz[2];

            rstK1[0] = rstK[0] - MU*(JINV(0,0)*dx + JINV(1,0)*dy + JINV(2,0)*dz);
            rstK1[1] = rstK[1] - MU*(JINV(0,1)*dx + JINV(1,1)*dy + JINV(2,1)*dz);
            rstK1[2] = rstK[2] - MU*(JINV(0,2)*dx + JINV(1,2)*dy + JINV(2,2)*dz);

            rstK = rstK1;

            ParametricToPhysical(rstK, outxyz);
            dist = L2dist(outxyz, xyz);
        }

        // Warn if not converged
        if (dist > tol)
        {
            cerr << "\nIsoparametricLinearPyramid::PhysicalToParametric: "
                    "Newton-Raphson did not converge.\n"
                 << "  Target:  (" << xyz[0]    << ", " << xyz[1]    << ", " << xyz[2]    << ")\n"
                 << "  Found:   (" << outxyz[0] << ", " << outxyz[1] << ", " << outxyz[2] << ")\n"
                 << "  Residual L2: " << dist << "  tolerance: " << tol << "\n"
                 << "  Parametric:  (" << rstK[0] << ", " << rstK[1] << ", " << rstK[2] << ")\n";

            csmp::Exception(WARNING,
                "IsoparametricLinearPyramid::PhysicalToParametric",
                "Newton-Raphson iteration not converged");
        }
    }

    rSt = rstK;
}

// -----------------------------------------------------------------------
// Face normals
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::UnitNormalToFace( uint32_t face,
                                                    vector<double>& unrml ) const
{
    assert(face < Faces());
    unrml.resize(3);

    auto setNormal = [&]( const Point<3>& n ) {
        unrml[0] = n[0];
        unrml[1] = n[1];
        unrml[2] = n[2];
    };

    switch (face)
    {
        case 0:   // triangle {0,1,4}
            setNormal(normalOfTriangle(
                Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                Point<3>(XY(4,0),XY(4,1),XY(4,2))));
            break;

        case 1:   // triangle {1,2,4}
            setNormal(normalOfTriangle(
                Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                Point<3>(XY(4,0),XY(4,1),XY(4,2))));
            break;

        case 2:   // triangle {2,3,4}
            setNormal(normalOfTriangle(
                Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                Point<3>(XY(4,0),XY(4,1),XY(4,2))));
            break;

        case 3:   // triangle {3,0,4}
            setNormal(normalOfTriangle(
                Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                Point<3>(XY(4,0),XY(4,1),XY(4,2))));
            break;

        case 4:   // quad base {0,3,2,1} — CCW from outside (below)
            setNormal(normalAtFacetCenter(
                Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                Point<3>(XY(1,0),XY(1,1),XY(1,2))));
            break;

        default:
            cerr << "\nIsoparametricLinearPyramid::UnitNormalToFace: "
                    "face " << face << " does not exist.\n";
            break;
    }
}

// -----------------------------------------------------------------------
// Reference coordinates
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const noexcept
{
    matCoords.Resize(npe, dim);
    matCoords = NXYZ;
}

// -----------------------------------------------------------------------
// VTK output
// -----------------------------------------------------------------------

void IsoparametricLinearPyramid::OutputNodeDataToVTK(
    const char* file_name,
    const char* var_name,
    DenseMatrix<DM_MIN>& DATA ) const
{
    // Build output filename: file_name + element_id + ".vtk"
    char outfile[NAME_STRING];
    strcpy(outfile, file_name);
    char elmt[30];
    snprintf(elmt, sizeof(elmt), "%zu", CurrentID());
    strcat(outfile, elmt);
    strcat(outfile, ".vtk");

    ofstream ofs(outfile, ios::out | ios::trunc);
    if (!ofs) {
        cerr << "\nIsoparametricLinearPyramid::OutputNodeDataToVTK: "
                "could not open '" << outfile << "'.\n";
        return;
    }

    // VTK header
    ofs << "# vtk DataFile Version 2.0\n"
        << "CSMP finite-element dataset: " << var_name << "\n"
        << "ASCII\n\n";

    // Node coordinates
    ofs << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << npe << " float\n";
    for (uint32_t i=0; i<npe; ++i) {
        for (uint32_t j=0; j<dim; ++j) ofs << XY(i,j) << " ";
        ofs << "\n";
    }
    ofs << "\n";

    // Cell connectivity
    ofs << "CELLS 1 6\n"
        << "5 0 1 2 3 4\n\n";

    // Cell type: VTK_PYRAMID = 14
    ofs << "CELL_TYPES 1\n14\n\n";

    // Point data
    ofs << "POINT_DATA " << npe << "\n";
    ofs << std::scientific;

    if (DATA.Rows() == 1) {
        ofs << "SCALARS " << var_name << " float\n"
            << "LOOKUP_TABLE default\n";
        for (uint32_t i=0; i<DATA.Cols(); ++i)
            ofs << DATA(0,i) << " ";
        ofs << "\n";
    } else {
        ofs << "VECTORS " << var_name << " float\n";
        for (uint32_t i=0; i<DATA.Cols(); ++i) {
            for (uint32_t j=0; j<DATA.Rows(); ++j)
                ofs << DATA(j,i) << "  ";
            ofs << "\n";
        }
    }
    ofs << "\n";

    cout << "\nIsoparametricLinearPyramid::OutputNodeDataToVTK: '"
         << outfile << "' written.\n";
}

} // namespace csmp

