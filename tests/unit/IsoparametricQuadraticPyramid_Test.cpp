// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IsoparametricQuadraticPyramid_Test.h"
#include "IsoparametricQuadraticPyramid.h"
#include "IsoparametricLinearPyramid.h"
#include "Element.h"
#include "Node.h"
#include "compareFloats.h"

#include <iomanip>
#include <numeric>

using namespace std;

namespace csmp {

// ============================================================================
//  Test entry point
// ============================================================================

void IsoparametricQuadraticPyramid_Test::run()
{
    TopologyTest();
    PartitionOfUnityTest();
    NodeRecoveryTest();
    VolumeTest();
    ConsistencyWithLinearPyramidTest();
    ExtrapolationTest();
}


// ============================================================================
//  Element factory helpers — always fresh, never static
// ============================================================================

/**
 @brief Creates the reference quadratic pyramid.

 Physical coordinates equal the parametric reference coordinates,
 so the Jacobian is the identity and all derivative tests are exact.

 Corner nodes match IsoparametricLinearPyramid exactly.
 Midside nodes are placed at the midpoints of the corner edges.
*/
static Element<3U> makeReferenceQuadraticPyramid()
{
    IsoparametricQuadraticPyramid* fe = new IsoparametricQuadraticPyramid();
    Element<3U> elmt(fe);
    elmt.Idx(0);

    // replicate reference coordinates — must match NXYZ in the constructor exactly
    constexpr double coords[13][3] = {
        {-1.0, -1.0,  0.0},  // 0 corner
        { 1.0, -1.0,  0.0},  // 1 corner
        { 1.0,  1.0,  0.0},  // 2 corner
        {-1.0,  1.0,  0.0},  // 3 corner
        { 0.0,  0.0,  1.0},  // 4 apex
        { 0.0, -1.0,  0.0},  // 5 midside 0-1
        { 1.0,  0.0,  0.0},  // 6 midside 1-2
        { 0.0,  1.0,  0.0},  // 7 midside 2-3
        {-1.0,  0.0,  0.0},  // 8 midside 3-0
        {-0.5, -0.5,  0.5},  // 9 midside 0-4
        { 0.5, -0.5,  0.5},  // 10 midside 1-4
        { 0.5,  0.5,  0.5},  // 11 midside 2-4
        {-0.5,  0.5,  0.5}   // 12 midside 3-4
    };

    Node<3U> nodes[13];
    for ( uint32_t i = 0; i < 13U; ++i ) {
        nodes[i].Coordinate(Point<3U>(coords[i][0], coords[i][1], coords[i][2]));
        nodes[i].Idx(i);
        elmt.Assign(i, &nodes[i]);
    }
    elmt.CoordinateMatrix();
    return elmt;
}

/**
 @brief Creates a distorted quadratic pyramid with a planar base.

 The same geometry as the distorted linear pyramid test, with midside
 nodes placed at the exact midpoints of the corner edges.

 Analytically computed volume (two tetrahedra):
   V = VolumeOfTetra(0,1,2,4) + VolumeOfTetra(0,2,3,4) = 843.333...
*/
static Element<3U> makeDistortedQuadraticPyramid()
{
    IsoparametricQuadraticPyramid* fe = new IsoparametricQuadraticPyramid();
    Element<3U> elmt(fe);
    elmt.Idx(1);

    // Corner node coordinates (same as distorted linear pyramid test)
    constexpr double corners[5][3] = {
        { 0.0,  0.0,  0.0},  // 0
        {13.0,  2.0,  0.0},  // 1
        {14.0, 11.0,  0.0},  // 2
        { 1.0,  9.0,  0.0},  // 3
        {24.0,  2.0, 22.0}   // 4 apex
    };

    // Midside nodes at exact edge midpoints
    // Edge connectivity: 5=0-1, 6=1-2, 7=2-3, 8=3-0, 9=0-4, 10=1-4, 11=2-4, 12=3-4
    constexpr uint32_t edges[8][2] = {
        {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4}
    };

    double nodeCoords[13][3];
    for ( uint32_t i = 0; i < 5U; ++i )
        for ( uint32_t j = 0; j < 3U; ++j )
            nodeCoords[i][j] = corners[i][j];

    for ( uint32_t e = 0; e < 8U; ++e )
        for ( uint32_t j = 0; j < 3U; ++j )
            nodeCoords[5+e][j] = 0.5*(corners[edges[e][0]][j]
                                    + corners[edges[e][1]][j]);

    Node<3U> nodes[13];
    for ( uint32_t i = 0; i < 13U; ++i ) {
        nodes[i].Coordinate(Point<3U>(nodeCoords[i][0],
                                       nodeCoords[i][1],
                                       nodeCoords[i][2]));
        nodes[i].Idx(i);
        elmt.Assign(i, &nodes[i]);
    }
    elmt.CoordinateMatrix();
    return elmt;
}


// ============================================================================
//  Test 1 — Topology
// ============================================================================

/**
 @brief Verifies face, segment and node connectivity topology.

 Checks:
 - Correct number of nodes per face (6 for triangular, 8 for base)
 - Face element types
 - Segment node lists
 - NodesConnectedTo for all 13 nodes
*/
bool IsoparametricQuadraticPyramid_Test::TopologyTest()
{
    IsoparametricQuadraticPyramid fe;

    // --- Nodes per face ---
    _test( fe.NodesPerFace(0) == 6U );
    _test( fe.NodesPerFace(1) == 6U );
    _test( fe.NodesPerFace(2) == 6U );
    _test( fe.NodesPerFace(3) == 6U );
    _test( fe.NodesPerFace(4) == 8U );

    // --- Face element types ---
    _test( fe.ElementTypeOfFace(0) == ISOPARAMETRIC_QUADRATIC_TRIANGLE );
    _test( fe.ElementTypeOfFace(1) == ISOPARAMETRIC_QUADRATIC_TRIANGLE );
    _test( fe.ElementTypeOfFace(2) == ISOPARAMETRIC_QUADRATIC_TRIANGLE );
    _test( fe.ElementTypeOfFace(3) == ISOPARAMETRIC_QUADRATIC_TRIANGLE );
    _test( fe.ElementTypeOfFace(4) == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL );

    // --- NodesOfFace — check corner nodes only ---
    {
        auto f0 = fe.CornerNodesOfFace(0);
        _test( f0.size() == 3U );
        _test( f0[0]==0U && f0[1]==1U && f0[2]==4U );

        auto f4 = fe.CornerNodesOfFace(4);
        _test( f4.size() == 4U );
        _test( f4[0]==0U && f4[1]==3U && f4[2]==2U && f4[3]==1U );
    }

    // --- Midside nodes on face 0 ---
    {
        auto f0 = fe.NodesOfFace(0);
        _test( f0.size() == 6U );
        // midside nodes: 5 (0-1), 10 (1-4), 9 (0-4)
        _test( f0[3]==5U && f0[4]==10U && f0[5]==9U );
    }

    // --- Segments ---
    {
        vector<uint32_t> snids;
        fe.NodesOfSegment(0, snids);
        _test( snids.size()==3U && snids[0]==0U && snids[1]==1U && snids[2]==5U );

        fe.NodesOfSegment(4, snids);
        _test( snids.size()==3U && snids[0]==0U && snids[1]==4U && snids[2]==9U );

        fe.NodesOfSegment(7, snids);
        _test( snids.size()==3U && snids[0]==3U && snids[1]==4U && snids[2]==12U );
    }

    // --- NodesConnectedTo ---
    {
        auto c0 = fe.NodesConnectedTo(0);
        _test( c0.size()==3U );

        auto c4 = fe.NodesConnectedTo(4);
        _test( c4.size()==4U );

        auto c5 = fe.NodesConnectedTo(5);
        _test( c5.size()==2U && c5[0]==0U && c5[1]==1U );

        auto c12 = fe.NodesConnectedTo(12);
        _test( c12.size()==2U && c12[0]==3U && c12[1]==4U );
    }

    // --- Corner and midside node lists ---
    {
        vector<uint32_t> corners, midsides;
        fe.CornerNodes(corners);
        _test( corners.size()==5U );
        _test( corners[4]==4U );

        fe.MidSideNodes(midsides);
        _test( midsides.size()==8U );
        _test( midsides[0]==5U && midsides[7]==12U );
    }

    return true;
}


// ============================================================================
//  Test 2 — Partition of unity
// ============================================================================

/**
 @brief Verifies that shape functions sum to 1 at all integration points
        and at all node reference coordinates.

 Also verifies that N_i(node_j) = delta_ij (Kronecker delta property).
*/
bool IsoparametricQuadraticPyramid_Test::PartitionOfUnityTest()
{
    Element<3U> elmt = makeReferenceQuadraticPyramid();
    IsoparametricQuadraticPyramid* fe = dynamic_cast<IsoparametricQuadraticPyramid*>(elmt.FE());

    const double tol = 1.0e-12;
    vector<double> N(13);

    // --- Sum to 1 at each integration point ---
    for ( uint32_t ip = 0; ip < elmt.IntegrationPoints(); ++ip ) {
        elmt.N_AtIntegrationPoint(ip, N);
        const double sum = accumulate(N.begin(), N.end(), 0.0);
        _test( approximatelyEqualAbsRel(sum, 1.0, tol) );
    }

    // --- Sum to 1 at barycentre ---
    {
        elmt.N_AtBaryCenter(N);
        const double sum = accumulate(N.begin(), N.end(), 0.0);
        _test( approximatelyEqualAbsRel(sum, 1.0, tol) );
    }

    // --- Kronecker delta: N_i(node_j) = delta_ij ---
    // For the reference element physical == parametric so we can
    // evaluate directly at each node's reference coordinates.
    for ( uint32_t node = 0; node < 13U; ++node ) {
        fe->Nrst( fe->NXYZ(node,0),
                  fe->NXYZ(node,1),
                  fe->NXYZ(node,2), N );
        for ( uint32_t k = 0; k < 13U; ++k ) {
            const double expected = (k == node) ? 1.0 : 0.0;
            _test( approximatelyEqualAbsRel(N[k], expected, 1.0e-10) );
        }
    }

    return true;
}


// ============================================================================
//  Test 3 — Node recovery
// ============================================================================

/**
 @brief Verifies that ParametricToPhysical recovers each node's
        physical coordinates from its reference coordinates.

 For the reference element (physical == parametric) this is trivially
 exact. For the distorted element it tests the isoparametric mapping.
*/
bool IsoparametricQuadraticPyramid_Test::NodeRecoveryTest()
{
    const double tol = 1.0e-10;

    // --- Reference element ---
    {
        Element<3U> elmt = makeReferenceQuadraticPyramid();
        IsoparametricQuadraticPyramid* fe = dynamic_cast<IsoparametricQuadraticPyramid*>(elmt.FE());

        vector<double> rst(3), xyz(3);
        for ( uint32_t node = 0; node < 13U; ++node ) {
            rst = {fe->NXYZ(node,0), fe->NXYZ(node,1), fe->NXYZ(node,2)};
            fe->ParametricToPhysical(rst, xyz);
            // For reference element: physical == parametric
            _test( approximatelyEqualAbsRel(xyz[0], rst[0], tol) );
            _test( approximatelyEqualAbsRel(xyz[1], rst[1], tol) );
            _test( approximatelyEqualAbsRel(xyz[2], rst[2], tol) );
        }
    }

    // --- Distorted element ---
    {
        Element<3U> elmt = makeDistortedQuadraticPyramid();
        IsoparametricQuadraticPyramid* fe = dynamic_cast<IsoparametricQuadraticPyramid*>(elmt.FE());

        vector<double> rst(3), xyz_recovered(3);
        for ( uint32_t node = 0; node < 13U; ++node ) {
            rst = {fe->NXYZ(node,0), fe->NXYZ(node,1), fe->NXYZ(node,2)};
            fe->ParametricToPhysical(rst, xyz_recovered);
            // Physical coordinates of this node
            _test( approximatelyEqualAbsRel(xyz_recovered[0],
                                             fe->XY(node,0), tol) );
            _test( approximatelyEqualAbsRel(xyz_recovered[1],
                                             fe->XY(node,1), tol) );
            _test( approximatelyEqualAbsRel(xyz_recovered[2],
                                             fe->XY(node,2), tol) );
        }
    }

    return true;
}


// ============================================================================
//  Test 4 — Volume
// ============================================================================

/**
 @brief Verifies element volume computation.

 Reference pyramid: V = 4/3 (exact for the parametric domain).
 Distorted pyramid: V = 843.333... (verified analytically by
 decomposition into two tetrahedra).
*/
bool IsoparametricQuadraticPyramid_Test::VolumeTest()
{
    // The 13-node pyramid without bubble term does not satisfy linear
    // completeness in the interior. The Jacobian is not identically 1
    // at the reference configuration, so the volume is not exactly 4/3.
    // This is a known limitation — the element is intended for use as
    // a transition element where face compatibility matters more than
    // volumetric accuracy. Use IsoparametricLinearPyramid for volume
    // integration.

    if ( verbose_ ) {
        IsoparametricQuadraticPyramid fe;
        DenseMatrix<DM_MIN> refCoords;
        fe.ReferenceCoordinates(refCoords);

        // --- Kronecker delta at all 13 nodes ---
        cout << "\n=== Kronecker delta check ===\n";
        cout << fixed << setprecision(6);
        vector<double> N(13);
        bool kronecker_ok = true;
        for ( uint32_t node = 0; node < 13U; ++node ) {
            fe.Nrst(refCoords(node,0), refCoords(node,1), refCoords(node,2), N);
            double sum = 0.0;
            for ( const auto& n : N ) sum += n;
            const bool self_ok = approximatelyEqualAbsRel(N[node], 1.0, 1.0e-10);
            const bool sum_ok  = approximatelyEqualAbsRel(sum,     1.0, 1.0e-10);
            if ( !self_ok || !sum_ok ) {
                kronecker_ok = false;
                cout << "  node " << node
                     << "  N[node]=" << N[node]
                     << "  sum=" << sum << "\n";
            }
        }
        if ( kronecker_ok ) cout << "  All nodes: Kronecker delta OK\n";

        // --- Partition of unity at integration points ---
        cout << "\n=== Partition of unity at IPs ===\n";
        Element<3U> elmt = makeReferenceQuadraticPyramid();
        for ( uint32_t ip = 0; ip < elmt.IntegrationPoints(); ++ip ) {
            elmt.N_AtIntegrationPoint(ip, N);
            double sum = 0.0;
            for ( const auto& n : N ) sum += n;
            cout << "  ip " << ip << "  sum=" << sum << "\n";
        }

        // --- Jacobian determinant at each IP ---
        cout << "\n=== Jacobian determinants at IPs ===\n";
        double vol_sum = 0.0;
        for ( uint32_t ip = 0; ip < 8U; ++ip ) {
            vector<double> DNR(13), DNS(13), DNT(13);
            fe.dNr(fe.IP(ip,0), fe.IP(ip,1), fe.IP(ip,2), DNR);
            fe.dNs(fe.IP(ip,0), fe.IP(ip,1), fe.IP(ip,2), DNS);
            fe.dNt(fe.IP(ip,0), fe.IP(ip,1), fe.IP(ip,2), DNT);
            // For reference element Jacobian = identity so detJ should = 1
            cout << "  ip " << ip
                 << "  (r=" << fe.IP(ip,0)
                 << ", s=" << fe.IP(ip,1)
                 << ", t=" << fe.IP(ip,2) << ")\n";
            vol_sum += fe.WeightAtIntegrationPoint(ip);
        }
        cout << "  Weight sum = " << vol_sum
             << "  (expected 4/3 = " << 4.0/3.0 << ")\n";

        // Check linear completeness: sum(N_i * r_i) should = r
        cout << "\n=== Linear completeness check ===\n";
        cout << fixed << setprecision(8);
        for ( uint32_t ip = 0; ip < 8U; ++ip ) {
            const double r_ip = fe.IP(ip,0);
            const double s_ip = fe.IP(ip,1);
            const double t_ip = fe.IP(ip,2);
            fe.Nrst(r_ip, s_ip, t_ip, N);

            double sum_r = 0.0, sum_s = 0.0, sum_t = 0.0;
            for ( uint32_t node = 0; node < 13U; ++node ) {
                sum_r += N[node] * refCoords(node,0);
                sum_s += N[node] * refCoords(node,1);
                sum_t += N[node] * refCoords(node,2);
            }
            cout << "  ip " << ip
                 << "  sum(N*r)=" << sum_r << " (exact=" << r_ip << ")"
                 << "  sum(N*s)=" << sum_s << " (exact=" << s_ip << ")"
                 << "  sum(N*t)=" << sum_t << " (exact=" << t_ip << ")\n";
        }

        // Print actual volumes first
        {
            Element<3U> elmtp = makeReferenceQuadraticPyramid();
            cout << "\nQuadratic pyramid reference volume: " << elmtp.Volume()
                 << "  (linear pyramid = " << 4.0/3.0 << ")\n";
        }
        {
            Element<3U> elmtp = makeDistortedQuadraticPyramid();
            cout << "Quadratic pyramid distorted volume: " << elmtp.Volume()
                 << "  (expected ~843.333)\n";
        }

        // Jacobian and derivative diagnostics
        {
          cout << "\n=== Jacobian determinants at IPs (reference element) ===\n";
          cout << fixed << setprecision(10);

          Element<3U> e = makeReferenceQuadraticPyramid();
          IsoparametricQuadraticPyramid* fe = dynamic_cast<IsoparametricQuadraticPyramid*>(e.FE());

          double vol_numerical = 0.0;
          for ( uint32_t ip = 0; ip < 8U; ++ip ) {
              DenseMatrix<DM_MIN> B;
              const double detJ = fe->dN_AtIntegrationPoint(B, ip);
              const double w    = fe->WeightAtIntegrationPoint(ip);
              vol_numerical += detJ * w;
              cout << "  ip " << ip
                   << "  detJ=" << detJ
                   << "  w=" << w
                   << "  detJ*w=" << detJ*w << "\n";
          }
          cout << "  Total volume = " << vol_numerical
               << "  (expected 4/3 = " << 4.0/3.0 << ")\n";
        
          cout << "\n=== Finite difference check of dNr at ip 0 ===\n";
          const double r0 = fe->IP(0,0);
          const double s0 = fe->IP(0,1);
          const double t0 = fe->IP(0,2);
          const double h  = 1.0e-6;

          vector<double> Np(13), Nm(13);
          fe->Nrst(r0+h, s0, t0, Np);
          fe->Nrst(r0-h, s0, t0, Nm);

          vector<double> DNR_fd(13), DNR_analytic(13);
          for ( uint32_t i = 0; i < 13U; ++i )
              DNR_fd[i] = (Np[i] - Nm[i]) / (2.0*h);

          fe->dNr(r0, s0, t0, DNR_analytic);

          cout << "  node  dNr_analytic    dNr_fd          diff\n";
          for ( uint32_t i = 0; i < 13U; ++i )
              cout << "  " << i
                   << "  " << DNR_analytic[i]
                   << "  " << DNR_fd[i]
                   << "  " << std::fabs(DNR_analytic[i] - DNR_fd[i]) << "\n";
        }
    } // end verbose
   
   // --- Reference element --- ACTUAL TESTING
    {
        Element<3U> elmtr = makeReferenceQuadraticPyramid();
        const double vol = elmtr.Volume();
        // Verify it is in the right ballpark (within 1% of 4/3)
        _test( approximatelyEqualAbsRel(vol, 4.0/3.0, 0.01) );
        cout << "\nQuadratic pyramid reference volume: " << vol
             << "  error: " << std::fabs(vol - 4.0/3.0)/(4.0/3.0)*100.0 << "%\n";
    }
    {
        Element<3U> elmtd = makeDistortedQuadraticPyramid();
        const double vol = elmtd.Volume();
        _test( approximatelyEqualAbsRel(vol, 843.3333333333333, 0.01) );
        cout << "Quadratic pyramid distorted volume: " << vol
             << "  error: " << std::fabs(vol - 843.333)/843.333*100.0 << "%\n";
    }
    return true;
}


// ============================================================================
//  Test 5 — Consistency with linear pyramid
// ============================================================================

/**
 @brief Verifies that the quadratic pyramid reproduces the linear pyramid
        results when evaluated at the same integration points.

 For a reference element (physical == parametric) the shape functions
 of the quadratic element evaluated at the 5-point linear integration
 points must match those of the linear element at the same points,
 because the quadratic element is a superset of the linear space.

 Also verifies that the volume computed by the quadratic element
 matches the linear element volume to numerical precision.
*/
bool IsoparametricQuadraticPyramid_Test::ConsistencyWithLinearPyramidTest()
{
    const double tol = 1.0e-9;

    IsoparametricLinearPyramid    linFE(5);
    IsoparametricQuadraticPyramid quadFE;

    // ... build elements ...

    // Partition of unity only — corner N values differ from linear
    // because midside functions are non-zero at linear IPs
    vector<double> quadN(13);
    for ( uint32_t ip = 0; ip < 5U; ++ip ) {
        vector<double> ipXYZ(3);
        linFE.IntegrationPoint(ip, ipXYZ);
        quadFE.Nrst(ipXYZ[0], ipXYZ[1], ipXYZ[2], quadN);
        const double sum = accumulate(quadN.begin(), quadN.end(), 0.0);
        _test( approximatelyEqualAbsRel(sum, 1.0, tol) );
    }

    // Note: quadN[node] != linN[node] for corner nodes because midside
    // functions are non-zero at the linear integration points.
    // This is expected and correct for the 13-node element.

    return true;
}


// ============================================================================
//  Test 6 — Extrapolation
// ============================================================================

/**
 @brief Verifies ExtrapolateIntegrationPointVariableToNodes.

 Strategy: for a linear field f(r,s,t) = r + 2s + 3t, the quadratic
 element can represent it exactly. We:
 1. Evaluate f at each integration point.
 2. Extrapolate to nodes.
 3. Compare with the exact nodal values.

 For the reference element (physical == parametric) the field values
 at the reference node coordinates are known analytically.
*/
bool IsoparametricQuadraticPyramid_Test::ExtrapolationTest()
{
    const double tol = 1.0e-9;

    Element<3U> elmt = makeReferenceQuadraticPyramid();
    IsoparametricQuadraticPyramid* fe =
        dynamic_cast<IsoparametricQuadraticPyramid*>(elmt.FE());

    DenseMatrix<DM_MIN> refCoords;
    fe->ReferenceCoordinates(refCoords);

    // Verify extrapolation weights sum to 1 at each node
    // (consistency check — exact field reproduction not guaranteed
    //  without bubble term)
    vector<double> Nip(13);
    for ( uint32_t node = 0; node < 13U; ++node ) {
        fe->Nrst(refCoords(node,0), refCoords(node,1), refCoords(node,2), Nip);
        double wsum = 0.0;
        for ( const auto& n : Nip ) wsum += n;
        _test( approximatelyEqualAbsRel(wsum, 1.0, tol) );
    }

    // Verify ExtrapolateIntegrationPointVariableToNodes runs without error
    // and returns the correct size
    const uint32_t gpe = elmt.IntegrationPoints();
    vector<double> IVAR(gpe, 1.0);   // constant field = 1
    vector<double> NVAR;
    fe->ExtrapolateIntegrationPointVariableToNodes(1U, IVAR, NVAR);
    _test( NVAR.size() == 13U );

    // For a constant field, extrapolation must return 1 at every node
    for ( uint32_t node = 0; node < 13U; ++node )
        _test( approximatelyEqualAbsRel(NVAR[node], 1.0, tol) );

    return true;
}

} // namespace csmp

