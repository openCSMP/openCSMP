#include "IsoparametricLinearTetrahedron.h"
#include "Exception.h"
#include "TriangularFacet.h"

#include <cassert>
#include <cstring>   // strcpy, strcat
#include <cstdio>    // snprintf
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace csmp {

// ============================================================
// Constructor
// ============================================================
IsoparametricLinearTetrahedron::IsoparametricLinearTetrahedron(
        uint32_t integrationPoints)
    : FiniteElement(ISOPARAMETRIC_LINEAR_TETRAHEDRON, true, true, 1U)
    , NXYZ(4, 3)
    , IP(integrationPoints, 3)
    , projectionCalledNTimes(0)
    , accDistance(0.0)
    , totIterations(0)
    , nonConvergenceOfProjections(0)
{
    dim = 3;
    itp = 1;
    npf = 3;
    npe = 4;
    fpe = 4;
    spe = 6;
    epe = 4;
    nne = 8;
    cne = 0;
    gpe = integrationPoints;

    M.Resize(npe, npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_LINEAR_TETRAHEDRON);

    XY.Resize(npe, dim);
    JAC.Resize(dim, dim);
    JINV.Resize(dim, dim);

    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    // Local node coordinates: r==ksi, s==nu, t==mu
    NXYZ(0,0) = 0.0; NXYZ(0,1) = 0.0; NXYZ(0,2) = 0.0;
    NXYZ(1,0) = 1.0; NXYZ(1,1) = 0.0; NXYZ(1,2) = 0.0;
    NXYZ(2,0) = 0.0; NXYZ(2,1) = 1.0; NXYZ(2,2) = 0.0;
    NXYZ(3,0) = 0.0; NXYZ(3,1) = 0.0; NXYZ(3,2) = 1.0;

    W.resize(gpe);

    if (integrationPoints == 1U) {
        IP(0,0) = 1./4.; IP(0,1) = 1./4.; IP(0,2) = 1./4.;
        W[0] = 1./6.;
    }
    else if (integrationPoints == 4U) {
        IP(0,0) = 0.585410196624969; IP(0,1) = 0.138196601125011; IP(0,2) = 0.138196601125011;
        W[0] = 1./24.;
        IP(1,0) = 0.138196601125011; IP(1,1) = 0.585410196624969; IP(1,2) = 0.138196601125011;
        W[1] = 1./24.;
        IP(2,0) = 0.138196601125011; IP(2,1) = 0.138196601125011; IP(2,2) = 0.585410196624969;
        W[2] = 1./24.;
        IP(3,0) = 0.138196601125011; IP(3,1) = 0.138196601125011; IP(3,2) = 0.138196601125011;
        W[3] = 1./24.;
    }
    else {
        throw std::range_error(
            "IsoparametricLinearTetrahedron: number of integration points must be 1 or 4");
    }
}

// ============================================================
// Geometry
// ============================================================
double IsoparametricLinearTetrahedron::Volume()
{
    double area{0.};
    for (uint32_t i{0U}; i < gpe; ++i) {
        dNr(IP(i,0), IP(i,1), IP(i,2), DNR);
        dNs(IP(i,0), IP(i,1), IP(i,2), DNS);
        dNt(IP(i,0), IP(i,1), IP(i,2), DNT);
        Jacobian(DNR, DNS, DNT);
        area += JacobianInverse() * W[i];
    }
    return area;
}

double IsoparametricLinearTetrahedron::AspectRatio()
{
    NRST.resize(spe);
    EdgeLengths(NRST);

    double seg_max{NRST[0]}, seg_min{NRST[0]};
    for (uint32_t i{1U}; i < spe; ++i) {
        if (NRST[i] > seg_max) seg_max = NRST[i];
        if (NRST[i] < seg_min) seg_min = NRST[i];
    }
    return seg_max / seg_min;
}

double IsoparametricLinearTetrahedron::InnerRadius()
{
    NRST.resize(spe);
    double sum{0.};
    EdgeLengths(NRST);
    for (uint32_t i{0U}; i < spe; ++i) sum += NRST[i];

    if (AspectRatio() > 4.)
        cerr << "\nIsoparametricLinearTetrahedron::InnerRadius: "
                "WARNING: function not applicable for this high aspect ratio.\n";

    return Volume() / (sum / 2.);
}

void IsoparametricLinearTetrahedron::EdgeLengths(std::vector<double>& len)
{
    len.resize(spe);
    auto seg = [&](uint32_t a, uint32_t b) {
        double dx = XY(a,0)-XY(b,0), dy = XY(a,1)-XY(b,1), dz = XY(a,2)-XY(b,2);
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };
    len[0] = seg(1,0);
    len[1] = seg(2,1);
    len[2] = seg(0,2);
    len[3] = seg(3,0);
    len[4] = seg(3,1);
    len[5] = seg(3,2);
}

// ============================================================
// Topology
// ============================================================
void IsoparametricLinearTetrahedron::CornerNodes(std::vector<uint32_t>& ids) const
{
    ids = {0U, 1U, 2U, 3U};
}

void IsoparametricLinearTetrahedron::CounterClockwiseNodes(std::vector<uint32_t>& ids) const
{
    ids = {0U, 1U, 2U, 3U};
}

void IsoparametricLinearTetrahedron::MidSideNodes(std::vector<uint32_t>&) const
{
    throw csmp::Exception(WARNING,
        "IsoparametricLinearTetrahedron::MidSideNodes:",
        "MidSideNodes not present",
        "Probably unintended use of function.");
}

CSMP_FEM_TYPE
IsoparametricLinearTetrahedron::ElementTypeOfFace(uint32_t) const noexcept
{
    return ISOPARAMETRIC_LINEAR_TRIANGLE;
}

void IsoparametricLinearTetrahedron::NodesOfSegment(
        uint32_t segm_id, std::vector<uint32_t>& snids) const
{
    // clang-format off
    static constexpr uint32_t table[6][2] = {
        {0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}
    };
    // clang-format on
    if (segm_id >= spe) {
        cerr << "\nIsoparametricLinearTetrahedron::NodesOfSegment: "
                "erratic segment id: " << segm_id << '\n';
        return;
    }
    snids = {table[segm_id][0], table[segm_id][1]};
}

vector<uint32_t>
IsoparametricLinearTetrahedron::NodesOfFace(uint32_t face_id) const
{
    switch (face_id) {
        case 0: return {1U,2U,3U};
        case 1: return {0U,3U,2U};
        case 2: return {0U,1U,3U};
        case 3: return {0U,2U,1U};
        default:
            cerr << "\nIsoparametricLinearTetrahedron::NodesOfFace: "
                    "face " << face_id << " does not exist.\n";
            return {};
    }
}

vector<uint32_t>
IsoparametricLinearTetrahedron::CornerNodesOfFace(uint32_t face_id) const
{
    // For a linear tetrahedron all face nodes are corner nodes.
    return NodesOfFace(face_id);
}

vector<uint32_t>
IsoparametricLinearTetrahedron::NodesConnectedTo(uint32_t node_id) const
{
    switch (node_id) {
        case 0: return {1U,2U,3U};
        case 1: return {0U,2U,3U};
        case 2: return {0U,1U,3U};
        case 3: return {0U,1U,2U};
        default:
            cerr << "\nIsoparametricLinearTetrahedron::NodesConnectedTo: "
                    "node " << node_id << " does not exist.\n";
            return {};
    }
}

void IsoparametricLinearTetrahedron::UnitNormalToFace(
        uint32_t face, std::vector<double>& unrml) const
{
    assert(face < Faces());
    unrml.resize(3);

    // Helper: copy Point<3> into unrml
    auto store = [&](const Point<3>& p) {
        unrml[0] = p[0]; unrml[1] = p[1]; unrml[2] = p[2];
    };

    switch (face) {
        case 0:
            store(normalOfTriangle(Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                   Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                   Point<3>(XY(3,0),XY(3,1),XY(3,2))));
            return;
        case 1:
            store(normalOfTriangle(Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                   Point<3>(XY(3,0),XY(3,1),XY(3,2)),
                                   Point<3>(XY(2,0),XY(2,1),XY(2,2))));
            return;
        case 2:
            store(normalOfTriangle(Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                   Point<3>(XY(1,0),XY(1,1),XY(1,2)),
                                   Point<3>(XY(3,0),XY(3,1),XY(3,2))));
            return;
        case 3:
            store(normalOfTriangle(Point<3>(XY(0,0),XY(0,1),XY(0,2)),
                                   Point<3>(XY(2,0),XY(2,1),XY(2,2)),
                                   Point<3>(XY(1,0),XY(1,1),XY(1,2))));
            return;
        default:
            assert(false && "UnitNormalToFace: invalid face index");
    }
}

// ============================================================
// Shape functions — physical-coordinate interface
// ============================================================
void IsoparametricLinearTetrahedron::N(
        vector<double>& N, const vector<double>& xyz)
{
    DNR.resize(dim);
    PhysicalToParametric(DNR, xyz);
    Nrst(DNR[0], DNR[1], DNR[2], N);
}

void IsoparametricLinearTetrahedron::N(
        vector<double>& N, uint32_t& /*iterations*/,
        double& /*distance*/, const vector<double>& xyz)
{
    // Delegate to the standard physical-coordinate overload.
    this->N(N, xyz);
}

void IsoparametricLinearTetrahedron::N_AtIntegrationPoint(
        uint32_t ip, std::vector<double>& N)
{
    assert(ip < gpe);
    Nrst(IP(ip,0), IP(ip,1), IP(ip,2), N);
}

void IsoparametricLinearTetrahedron::N_AtBaryCenter(std::vector<double>& N)
{
    constexpr double kQuarter{0.25};
    Nrst(kQuarter, kQuarter, kQuarter, N);
}

// ============================================================
// Jacobian
// ============================================================
void IsoparametricLinearTetrahedron::JacobianAtIntegrationPoint(uint32_t gauss_point)
{
    assert(gauss_point < gpe);
    dNr(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR);
    dNs(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS);
    dNt(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT);
    Jacobian(DNR, DNS, DNT);
}

void IsoparametricLinearTetrahedron::JacobianAt(const std::vector<double>& rst)
{
    dNr(rst[0], rst[1], rst[2], DNR);
    dNs(rst[0], rst[1], rst[2], DNS);
    dNt(rst[0], rst[1], rst[2], DNT);
    Jacobian(DNR, DNS, DNT);
}

// ============================================================
// Shape-function derivatives — global coordinates
// ============================================================
void IsoparametricLinearTetrahedron::dN(DenseMatrix<DM_MIN>& DN4)
{
    DN4.Resize(dim, npe);
    M.Resize(dim, 1);

    for (uint32_t i{0U}; i < npe; ++i) {
        dNr(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR);
        dNs(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS);
        dNt(NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT);

        Jacobian(DNR, DNS, DNT);
        JacobianInverse();

        M(0,0) = DNR[i]; M(1,0) = DNS[i]; M(2,0) = DNT[i];
        JINV *= M;
        DN4(0,i) = JINV(0,0);
        DN4(1,i) = JINV(1,0);
        DN4(2,i) = JINV(2,0);
        JINV.Resize(dim, dim);
    }
}

double IsoparametricLinearTetrahedron::dN(
        DenseMatrix<DM_MIN>& DN2, const vector<double>& xyz)
{
    vector<double> rst(dim);
    PhysicalToParametric(rst, xyz);

    dNr(rst[0], rst[1], rst[2], DNR);
    dNs(rst[0], rst[1], rst[2], DNS);
    dNt(rst[0], rst[1], rst[2], DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    DN2.Resize(dim, dim);
    DN2 = JINV;

    DenseMatrix<DM_MIN> DN_TEMP(dim, npe);
    for (uint32_t i{0U}; i < npe; ++i) {
        DN_TEMP(0,i) = DNR[i];
        DN_TEMP(1,i) = DNS[i];
        DN_TEMP(2,i) = DNT[i];
    }
    DN2 *= DN_TEMP;

    return detJ;
}

double IsoparametricLinearTetrahedron::dN_AtNode(
        DenseMatrix<DM_MIN>& B, uint32_t nd)
{
    assert(nd < npe);
    dNr(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR);
    dNs(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS);
    dNt(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for (uint32_t i{0U}; i < npe; ++i) {
        B(0,i) = DNR[i]; B(1,i) = DNS[i]; B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}

double IsoparametricLinearTetrahedron::dN_AtIntegrationPoint(
        DenseMatrix<DM_MIN>& B, uint32_t gauss_point)
{
    assert(gauss_point < gpe);
    dNr(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR);
    dNs(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS);
    dNt(IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for (uint32_t i{0U}; i < npe; ++i) {
        B(0,i) = DNR[i]; B(1,i) = DNS[i]; B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}

double IsoparametricLinearTetrahedron::dN_AtBarycenter(DenseMatrix<DM_MIN>& B)
{
    constexpr double kQ{0.25};
    dNr(kQ, kQ, kQ, DNR);
    dNs(kQ, kQ, kQ, DNS);
    dNt(kQ, kQ, kQ, DNT);

    Jacobian(DNR, DNS, DNT);
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for (uint32_t i{0U}; i < npe; ++i) {
        B(0,i) = DNR[i]; B(1,i) = DNS[i]; B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}

// ============================================================
// Integration point location (global coordinates)
// ============================================================
void IsoparametricLinearTetrahedron::IntegrationPoint(
        uint32_t ip, vector<double>& xyz) const
{
    assert(ip < gpe);
    xyz.assign(3U, 0.);
    Nrst(IP(ip,0), IP(ip,1), IP(ip,2), NRST);
    for (uint32_t i{0U}; i < npe; ++i) {
        xyz[0] += XY(i,0) * NRST[i];
        xyz[1] += XY(i,1) * NRST[i];
        xyz[2] += XY(i,2) * NRST[i];
    }
}

// ============================================================
// Coordinate transformations
// ============================================================
void IsoparametricLinearTetrahedron::ParametricToPhysical(
        vector<double>& rst, vector<double>& xyz)
{
    Nrst(rst[0], rst[1], rst[2], DNR);
    xyz.assign(dim, 0.);
    for (uint32_t i{0U}; i < npe; ++i) {
        xyz[0] += XY(i,0) * DNR[i];
        xyz[1] += XY(i,1) * DNR[i];
        xyz[2] += XY(i,2) * DNR[i];
    }
}

void IsoparametricLinearTetrahedron::PhysicalToParametric(
        vector<double>& rSt, const vector<double>& xyz)
{
    const double x12(XY(0,0)-XY(1,0)), x13(XY(0,0)-XY(2,0)), x14(XY(0,0)-XY(3,0));
    const double x23(XY(1,0)-XY(2,0)), x24(XY(1,0)-XY(3,0)), x34(XY(2,0)-XY(3,0));
    const double x21(-x12), x31(-x13), x32(-x23), x43(-x34);

    const double y12(XY(0,1)-XY(1,1)), y13(XY(0,1)-XY(2,1)), y14(XY(0,1)-XY(3,1));
    const double y23(XY(1,1)-XY(2,1)), y24(XY(1,1)-XY(3,1)), y34(XY(2,1)-XY(3,1));
    const double y21(-y12), y31(-y13), y43(-y34);

    const double z12(XY(0,2)-XY(1,2)), z13(XY(0,2)-XY(2,2)), z14(XY(0,2)-XY(3,2));
    const double z23(XY(1,2)-XY(2,2)), z24(XY(1,2)-XY(3,2)), z34(XY(2,2)-XY(3,2));
    const double z21(-z12), z31(-z13), z43(-z34);

    const double a2(y31*z43 - y34*z13);
    const double a3(y24*z14 - y14*z24);
    const double a4(y13*z21 - y12*z31);

    const double b2(x43*z31 - x13*z34);
    const double b3(x14*z24 - x24*z14);
    const double b4(x21*z13 - x31*z12);

    const double c2(x31*y43 - x34*y13);
    const double c3(x24*y14 - x14*y24);
    const double c4(x13*y21 - x12*y31);

    const double V00(x21*(y23*z34 - y34*z23)
                   + x32*(y34*z12 - y12*z34)
                   + x43*(y12*z23 - y23*z12));

    const double V02(XY(0,0)*(XY(3,1)*XY(2,2)-XY(2,1)*XY(3,2))
                   + XY(2,0)*(XY(0,1)*XY(3,2)-XY(3,1)*XY(0,2))
                   + XY(3,0)*(XY(2,1)*XY(0,2)-XY(0,1)*XY(2,2)));
    const double V03(XY(0,0)*(XY(1,1)*XY(3,2)-XY(3,1)*XY(1,2))
                   + XY(1,0)*(XY(3,1)*XY(0,2)-XY(0,1)*XY(3,2))
                   + XY(3,0)*(XY(0,1)*XY(1,2)-XY(1,1)*XY(0,2)));
    const double V04(XY(0,0)*(XY(2,1)*XY(1,2)-XY(1,1)*XY(2,2))
                   + XY(1,0)*(XY(0,1)*XY(2,2)-XY(2,1)*XY(0,2))
                   + XY(2,0)*(XY(1,1)*XY(0,2)-XY(0,1)*XY(1,2)));

    rSt[0] = (V02 + a2*xyz[0] + b2*xyz[1] + c2*xyz[2]) / V00;
    rSt[1] = (V03 + a3*xyz[0] + b3*xyz[1] + c3*xyz[2]) / V00;
    rSt[2] = (V04 + a4*xyz[0] + b4*xyz[1] + c4*xyz[2]) / V00;
}

// ============================================================
// Extrapolation
// ============================================================
// WARNING: the static locals below are NOT thread-safe.
// Consider replacing with member variables or a std::once_flag
// if this class is used in a multi-threaded context.
void IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes(
        uint32_t nvars,
        const vector<double>& IVAR,
        vector<double>& NVAR) const
{
    assert(IVAR.size() >= (gpe * nvars));
    NVAR.resize(npe * nvars);

    if (gpe == 1U) {
        for (uint32_t i{0U}; i < npe; ++i)
            for (uint32_t k{0U}; k < nvars; ++k)
                NVAR[i*nvars + k] = IVAR[k];
        return;
    }

    static double a[4], b[4], c[4], d[4], intpol[4], volume6;
    static bool   first_call{true};

    if (IVAR.size() != gpe * nvars)
        throw csmp::Exception(FATAL_ERROR,
            "IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
            "Input vector must have 'nvars' x 4 entries");

    if (NVAR.size() != npe * nvars)
        throw csmp::Exception(FATAL_ERROR,
            "IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
            "Output vector must have 'nvars' x nodes entries");

    if (first_call) {
        if (gpe != 4U)
            throw csmp::Exception(FATAL_ERROR,
                "IsoparametricLinearTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
                "This method expects four integration points");

        for ( int32_t i = 0, j = 1; i < 4; ++i) {
            a[i]  = -IP(n(i,1),0)*(IP(n(i,3),1)*IP(n(i,2),2)-IP(n(i,3),2)*IP(n(i,2),1));
            a[i] -=  IP(n(i,2),0)*(IP(n(i,1),1)*IP(n(i,3),2)-IP(n(i,1),2)*IP(n(i,3),1));
            a[i] -=  IP(n(i,3),0)*(IP(n(i,2),1)*IP(n(i,1),2)-IP(n(i,2),2)*IP(n(i,1),1));

            b[i]  = IP(n(i,3),1)*IP(n(i,2),2) - IP(n(i,3),2)*IP(n(i,2),1);
            b[i] += IP(n(i,1),1)*IP(n(i,3),2) - IP(n(i,1),2)*IP(n(i,3),1);
            b[i] += IP(n(i,2),1)*IP(n(i,1),2) - IP(n(i,2),2)*IP(n(i,1),1);

            c[i]  = IP(n(i,3),2)*IP(n(i,2),0) - IP(n(i,3),0)*IP(n(i,2),2);
            c[i] += IP(n(i,1),2)*IP(n(i,3),0) - IP(n(i,1),0)*IP(n(i,3),2);
            c[i] += IP(n(i,2),2)*IP(n(i,1),0) - IP(n(i,2),0)*IP(n(i,1),2);

            d[i]  = IP(n(i,3),0)*IP(n(i,2),1) - IP(n(i,3),1)*IP(n(i,2),0);
            d[i] += IP(n(i,1),0)*IP(n(i,3),1) - IP(n(i,1),1)*IP(n(i,3),0);
            d[i] += IP(n(i,2),0)*IP(n(i,1),1) - IP(n(i,2),1)*IP(n(i,1),0);

            a[i] *= static_cast<double>(j);
            b[i] *= static_cast<double>(j);
            c[i] *= static_cast<double>(j);
            d[i] *= static_cast<double>(j);
            j    *= -1;
        }

        volume6    = a[0] + a[1] + a[2] + a[3];
        first_call = false;
    }

    // 2. Extrapolate from integration points to nodes.
    vector<double> sum(nvars);
    for (uint32_t i{0U}; i < npe; ++i) {
        intpol[0] = (a[0] + b[0]*NXYZ(i,0) + c[0]*NXYZ(i,1) + d[0]*NXYZ(i,2)) / volume6;
        intpol[1] = (a[1] + b[1]*NXYZ(i,0) + c[1]*NXYZ(i,1) + d[1]*NXYZ(i,2)) / volume6;
        intpol[2] = (a[2] + b[2]*NXYZ(i,0) + c[2]*NXYZ(i,1) + d[2]*NXYZ(i,2)) / volume6;
        intpol[3] = (a[3] + b[3]*NXYZ(i,0) + c[3]*NXYZ(i,1) + d[3]*NXYZ(i,2)) / volume6;

        fill(sum.begin(), sum.end(), 0.0);
        for (uint32_t j{0U}; j < gpe; ++j)
            for (uint32_t k{0U}; k < nvars; ++k)
                sum[k] += intpol[j] * IVAR[j*nvars + k];

        for (uint32_t k{0U}; k < nvars; ++k)
            NVAR[i*nvars + k] = sum[k];
    }

} // end ExtrapolateIntegrationPointVariableToNodes

// ============================================================
// Reference coordinates
// ============================================================
void IsoparametricLinearTetrahedron::ReferenceCoordinates(
        DenseMatrix<DM_MIN>& matCoords) const
{
    matCoords.Resize(npe, dim);
    matCoords = NXYZ;
}

// ============================================================
// VTK output
// ============================================================
void IsoparametricLinearTetrahedron::OutputNodeDataToVTK(
        const char* file_name,
        const char* var_name,
        DenseMatrix<DM_MIN>& DATA) const
{
    char outfile[NAME_STRING], elmt[30];
    strcpy(outfile, file_name);
    snprintf(elmt, sizeof(elmt), "%zu", CurrentID());
    strcat(outfile, elmt);
    strcat(outfile, ".vtk");

    ofstream ofs;
    ofs.open(outfile, ios::out | ios::trunc);
    if (!ofs) {
        cerr << "\nIsoparametricLinearTetrahedron::OutputNodeDataToVTK: "
                "output file could not be opened.\n";
        return;
    }

    // Header
    ofs << "# vtk DataFile Version 2.0\n";
    ofs << "Finite-element dataset (CSMP): variable: " << var_name << '\n';
    ofs << "ASCII\n\n";

    // Node coordinates
    ofs << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << npe << " double\n";
    for (uint32_t i{0U}; i < npe; ++i) {
        for (uint32_t j{0U}; j < dim; ++j) ofs << XY(i,j) << ' ';
        ofs << '\n';
    }
    ofs << '\n';

    // Cells
    ofs << "CELLS 1 5\n";
    ofs << "4 0 1 2 3\n\n";

    // Cell types — 10 == VTK_TETRA
    ofs << "CELL_TYPES 1\n";
    ofs << "10\n\n";

    // Point data
    ofs << "POINT_DATA " << npe << '\n';
    ofs.setf(ios::scientific);

    if (DATA.Rows() == 1) {
        ofs << "SCALARS " << var_name << " double\n";
        ofs << "LOOKUP_TABLE default\n";
        for (uint32_t i{0U}; i < DATA.Cols(); ++i) ofs << DATA(0,i) << ' ';
        ofs << '\n';
    }
    else {
        ofs << "VECTORS " << var_name << " double\n";
        for (uint32_t i{0U}; i < DATA.Cols(); ++i) {
            for (uint32_t j{0U}; j < DATA.Rows(); ++j) ofs << DATA(j,i) << "  ";
            ofs << '\n';
        }
    }

    ofs << '\n';
    ofs.close();

    cout << "\nIsoparametricLinearTetrahedron::OutputNodeDataToVTK: file '"
         << outfile << "' written successfully.\n";

} // end OutputNodeDataToVTK

} // namespace csmp

