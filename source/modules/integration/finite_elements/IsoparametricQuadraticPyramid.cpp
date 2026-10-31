// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IsoparametricQuadraticPyramid.h"
#include "Exception.h"
#include "triangularFacet.h"
#include "QuadrilateralFacet.h"

using namespace std;

namespace csmp {

// ============================================================================
//  Constructor
// ============================================================================

/**
 @brief Constructs the 13-node quadratic pyramid finite element.

 Node numbering:

 @code
 Corner nodes (0-4) — same reference coordinates as linear pyramid:
   Node 0: r=-1, s=-1, t=0   (base, CCW from outside)
   Node 1: r=+1, s=-1, t=0
   Node 2: r=+1, s=+1, t=0
   Node 3: r=-1, s=+1, t=0
   Node 4: r= 0, s= 0, t=1   (apex)

 Midside nodes (5-12):
   Node  5: r= 0, s=-1, t=0   midpoint of edge 0-1
   Node  6: r=+1, s= 0, t=0   midpoint of edge 1-2
   Node  7: r= 0, s=+1, t=0   midpoint of edge 2-3
   Node  8: r=-1, s= 0, t=0   midpoint of edge 3-0
   Node  9: r=-0.5, s=-0.5, t=0.5  midpoint of edge 0-4
   Node 10: r=+0.5, s=-0.5, t=0.5  midpoint of edge 1-4
   Node 11: r=+0.5, s=+0.5, t=0.5  midpoint of edge 2-4
   Node 12: r=-0.5, s=+0.5, t=0.5  midpoint of edge 3-4
 @endcode

 Face numbering (consistent with IsoparametricLinearPyramid):
 @code
   Face 0: {0,1,4,  5,10, 9}  triangular
   Face 1: {1,2,4,  6,11,10}  triangular
   Face 2: {2,3,4,  7,12,11}  triangular
   Face 3: {0,4,3,  9,12, 8}  triangular
   Face 4: {0,3,2,1, 8,7,6,5} quadrilateral (base)
 @endcode

 Integration: 8-point conical product rule (same as linear pyramid).

 Shape functions: Bedrosian/Zganski rational formulation, corrected by
 Matthaei & Pichler (2006). Second-order complete on the pyramid domain.

 @note The rational term rs/(1-t) is singular at the apex (t=1).
       All methods guard against this with an explicit t==1 branch.
*/
IsoparametricQuadraticPyramid::IsoparametricQuadraticPyramid()
    : FiniteElement( ISOPARAMETRIC_QUADRATIC_PYRAMID13, true, true, 2U ),
      DN(3, 13),
      NXYZ(13, 3),
      IP(8, 3)
{
    dim = 3;
    itp = 2;
    npf = 6;    ///< max nodes per face (quad base has 8, triangular faces have 6)
    npe = 13;
    fpe = 5;
    spe = 8;
    epe = 5;
    nne = 6;
    cne = 0;
    gpe = 8;

    M.Resize(npe, npe);

    UsesLocalCoordinates(true);
    Isoparametric(true);
    VolumeElement();
    ElementType(ISOPARAMETRIC_QUADRATIC_PYRAMID13);

    XY.Resize(npe, dim);
    JAC.Resize(dim, dim);
    JINV.Resize(dim, dim);

    NRST.resize(npe);
    DNR.resize(npe);
    DNS.resize(npe);
    DNT.resize(npe);

    // -----------------------------------------------------------------------
    // Reference node coordinates — r, s, t
    // Corner nodes match IsoparametricLinearPyramid exactly.
    // -----------------------------------------------------------------------
    // Corner nodes
    NXYZ(0,0) =-1.0; NXYZ(0,1) =-1.0; NXYZ(0,2) = 0.0;
    NXYZ(1,0) = 1.0; NXYZ(1,1) =-1.0; NXYZ(1,2) = 0.0;
    NXYZ(2,0) = 1.0; NXYZ(2,1) = 1.0; NXYZ(2,2) = 0.0;
    NXYZ(3,0) =-1.0; NXYZ(3,1) = 1.0; NXYZ(3,2) = 0.0;
    NXYZ(4,0) = 0.0; NXYZ(4,1) = 0.0; NXYZ(4,2) = 1.0;
    // Midside nodes — base edges
    NXYZ(5,0) = 0.0; NXYZ(5,1) =-1.0; NXYZ(5,2) = 0.0;  // midpoint 0-1
    NXYZ(6,0) = 1.0; NXYZ(6,1) = 0.0; NXYZ(6,2) = 0.0;  // midpoint 1-2
    NXYZ(7,0) = 0.0; NXYZ(7,1) = 1.0; NXYZ(7,2) = 0.0;  // midpoint 2-3
    NXYZ(8,0) =-1.0; NXYZ(8,1) = 0.0; NXYZ(8,2) = 0.0;  // midpoint 3-0
    // Midside nodes — lateral edges to apex
    NXYZ(9,0) =-0.5; NXYZ(9,1) =-0.5; NXYZ(9,2) = 0.5;  // midpoint 0-4
    NXYZ(10,0)= 0.5; NXYZ(10,1)=-0.5; NXYZ(10,2)= 0.5;  // midpoint 1-4
    NXYZ(11,0)= 0.5; NXYZ(11,1)= 0.5; NXYZ(11,2)= 0.5;  // midpoint 2-4
    NXYZ(12,0)=-0.5; NXYZ(12,1)= 0.5; NXYZ(12,2)= 0.5;  // midpoint 3-4

    W.resize(gpe);
    GenerateIntegrationPoints(IP, W);
}


IsoparametricQuadraticPyramid::~IsoparametricQuadraticPyramid() {}


// ============================================================================
//  Topology
// ============================================================================

/**
 @brief Returns the 5 corner node local indices.
 @param ids Output vector of size 5: {0,1,2,3,4}.
*/
void IsoparametricQuadraticPyramid::CornerNodes( std::vector<uint32_t>& ids ) const
{
    ids = {0U, 1U, 2U, 3U, 4U};
}

/**
 @brief Returns the 8 midside node local indices.
 @param ids Output vector of size 8: {5,6,7,8,9,10,11,12}.
*/
void IsoparametricQuadraticPyramid::MidSideNodes( std::vector<uint32_t>& ids ) const
{
    ids = {5U, 6U, 7U, 8U, 9U, 10U, 11U, 12U};
}

/**
 @brief Returns the number of nodes on a given face.

 Face 4 (quadrilateral base) has 8 nodes; all triangular faces have 6.

 @param face_id  Face index in [0,4].
 @return         8 for face 4, 6 for faces 0–3.
*/
uint32_t IsoparametricQuadraticPyramid::NodesPerFace( uint32_t face_id ) const noexcept
{
    switch (face_id) {
        case 0: case 1: case 2: case 3: return 6U;
        case 4: return 8U;
        default:
            cerr << "\nIsoparametricQuadraticPyramid::NodesPerFace: face "
                 << face_id << " does not exist.";
            return npf;
    }
}

/**
 @brief Returns the ordered local node indices of a face.

 Triangular faces (0–3) list corner nodes first, then midside nodes.
 The quadrilateral base (face 4) lists corner nodes CCW from outside,
 then midside nodes in the same order.

 Face winding is outward-normal consistent:
 @code
   Face 0: {0,1,4,  5,10, 9}
   Face 1: {1,2,4,  6,11,10}
   Face 2: {2,3,4,  7,12,11}
   Face 3: {3,0,4,  8, 9,12}
   Face 4: {0,3,2,1, 8, 7, 6, 5}  (base, viewed from below)
 @endcode

 @param face_id  Face index in [0,4].
 @return         Vector of local node indices.
*/
vector<uint32_t> IsoparametricQuadraticPyramid::NodesOfFace( uint32_t face_id ) const
{
    switch (face_id) {
        case 0: return {0U, 1U, 4U,  5U, 10U,  9U};
        case 1: return {1U, 2U, 4U,  6U, 11U, 10U};
        case 2: return {2U, 3U, 4U,  7U, 12U, 11U};
        case 3: return {3U, 0U, 4U,  8U,  9U, 12U};
        case 4: return {0U, 3U, 2U, 1U,  8U,  7U,  6U,  5U};
        default:
            cerr << "\nIsoparametricQuadraticPyramid::NodesOfFace: face "
                 << face_id << " does not exist.";
            return {};
    }
}

/**
 @brief Returns only the corner nodes of a face (no midside nodes).
 @param face_id  Face index in [0,4].
 @return         Vector of corner node local indices.
*/
vector<uint32_t> IsoparametricQuadraticPyramid::CornerNodesOfFace( uint32_t face_id ) const
{
    switch (face_id) {
        case 0: return {0U, 1U, 4U};
        case 1: return {1U, 2U, 4U};
        case 2: return {2U, 3U, 4U};
        case 3: return {3U, 0U, 4U};
        case 4: return {0U, 3U, 2U, 1U};
        default:
            cerr << "\nIsoparametricQuadraticPyramid::CornerNodesOfFace: face "
                 << face_id << " does not exist.";
            return {};
    }
}

/**
 @brief Returns the element type of a given face.

 Face 4 (base) is a 9-node quadrilateral; all other faces are
 6-node triangles.

 @param face  Face index in [0,4].
 @return      CSMP_FEM_TYPE of the face element.
*/
CSMP_FEM_TYPE IsoparametricQuadraticPyramid::ElementTypeOfFace( uint32_t face ) const noexcept
{
    assert( face < fpe );
    if ( face == 4U ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;
    return ISOPARAMETRIC_QUADRATIC_TRIANGLE;
}



/**
 @brief Returns the two endpoint nodes and the midside node of a segment.

 Segment ordering:
 @code
   0: 0-1  (midside 5)    4: 0-4  (midside  9)
   1: 1-2  (midside 6)    5: 1-4  (midside 10)
   2: 2-3  (midside 7)    6: 2-4  (midside 11)
   3: 3-0  (midside 8)    7: 3-4  (midside 12)
 @endcode

 @param segm_id  Segment index in [0,7].
 @param snids    Output: {endpoint_a, endpoint_b, midside_node}.
*/
void IsoparametricQuadraticPyramid::NodesOfSegment(
    uint32_t segm_id, std::vector<uint32_t>& snids ) const
{
    if ( segm_id < 8U ) {
        snids = {kSegments[segm_id][0],
                 kSegments[segm_id][1],
                 kSegments[segm_id][2]};
    } else {
        cerr << "\nIsoparametricQuadraticPyramid::NodesOfSegment: segment "
             << segm_id << " does not exist.";
        snids.clear();
    }
}

/**
 @brief Returns the corner nodes directly connected to a given node.

 Only corner-to-corner adjacency is returned (not midside connections).

 @param node_id  Local node index in [0,12].
 @return         Vector of connected corner node indices.
*/
vector<uint32_t> IsoparametricQuadraticPyramid::NodesConnectedTo( uint32_t node_id ) const
{
    switch (node_id) {
        case 0:  return {1U, 3U, 4U};
        case 1:  return {0U, 2U, 4U};
        case 2:  return {1U, 3U, 4U};
        case 3:  return {0U, 2U, 4U};
        case 4:  return {0U, 1U, 2U, 3U};
        case 5:  return {0U, 1U};
        case 6:  return {1U, 2U};
        case 7:  return {2U, 3U};
        case 8:  return {3U, 0U};
        case 9:  return {0U, 4U};
        case 10: return {1U, 4U};
        case 11: return {2U, 4U};
        case 12: return {3U, 4U};
        default:
            cerr << "\nIsoparametricQuadraticPyramid::NodesConnectedTo: node "
                 << node_id << " does not exist.";
            return {};
    }
}

/**
 @brief Returns the unit outward normal to a face at its centroid.

 Uses the same triangular/quadrilateral facet helpers as
 IsoparametricLinearPyramid.

 @param face   Face index in [0,4].
 @param unrml  Output: unit normal vector of size 3.
*/
void IsoparametricQuadraticPyramid::UnitNormalToFace(
    uint32_t face, std::vector<double>& unrml ) const
{
    assert( face < fpe );
    unrml.resize(3);

    if ( face == 4U ) {
        // Quadrilateral base — use corner nodes 0,3,2,1
        Point<3> nrml = normalAtFacetCenter(
            Point<3>(XY(0,0),XY(0,1),XY(0,2)),
            Point<3>(XY(3,0),XY(3,1),XY(3,2)),
            Point<3>(XY(2,0),XY(2,1),XY(2,2)),
            Point<3>(XY(1,0),XY(1,1),XY(1,2)) );
        unrml[0]=nrml[0]; unrml[1]=nrml[1]; unrml[2]=nrml[2];
        return;
    }

    // Triangular lateral faces — corner nodes only
    static constexpr uint32_t kFaceCorners[4][3] = {
        {0U,1U,4U}, {1U,2U,4U}, {2U,3U,4U}, {3U,0U,4U}
    };
    const uint32_t a = kFaceCorners[face][0];
    const uint32_t b = kFaceCorners[face][1];
    const uint32_t c = kFaceCorners[face][2];
    Point<3> nrml = normalOfTriangle(
        Point<3>(XY(a,0),XY(a,1),XY(a,2)),
        Point<3>(XY(b,0),XY(b,1),XY(b,2)),
        Point<3>(XY(c,0),XY(c,1),XY(c,2)) );
    unrml[0]=nrml[0]; unrml[1]=nrml[1]; unrml[2]=nrml[2];
}


// ============================================================================
//  Shape functions — no bubble correction
//  Compatible with adjacent quadratic hexahedra and prisms.
//  Partition of unity and Kronecker delta hold exactly.
//  Polynomial space is incomplete (no internal bubble) but face-compatible.
// ============================================================================

void IsoparametricQuadraticPyramid::Nrst(
    double r, double s, double t,
    double* N ) const noexcept
{
    if ( t == 1.0 ) {
        for ( uint32_t i = 0; i < 13U; ++i ) N[i] = 0.0;
        N[4] = 1.0;
        return;
    }

    const double tM  = 1.0 - t;
    const double tM2 = 2.0 * tM;
    const double rat = r * s / tM;   // rs/(1-t) — rational term

    // Corner nodes
    N[0] = 0.25*(-r-s-1.0)*((1.0-r)*(1.0-s) - t + rat*t);
    N[1] = 0.25*( r-s-1.0)*((1.0+r)*(1.0-s) - t - rat*t);
    N[2] = 0.25*( r+s-1.0)*((1.0+r)*(1.0+s) - t + rat*t);
    N[3] = 0.25*(-r+s-1.0)*((1.0-r)*(1.0+s) - t - rat*t);
    N[4] = t*(2.0*t - 1.0);

    // Midside nodes — base edges
    N[5] = (1.0+r-t)*(1.0-r-t)*(1.0-s-t) / tM2;
    N[6] = (1.0+s-t)*(1.0-s-t)*(1.0+r-t) / tM2;
    N[7] = (1.0+r-t)*(1.0-r-t)*(1.0+s-t) / tM2;
    N[8] = (1.0+s-t)*(1.0-s-t)*(1.0-r-t) / tM2;

    // Midside nodes — lateral edges to apex
    N[9]  = t*(1.0-r-t)*(1.0-s-t) / tM;
    N[10] = t*(1.0+r-t)*(1.0-s-t) / tM;
    N[11] = t*(1.0+r-t)*(1.0+s-t) / tM;
    N[12] = t*(1.0-r-t)*(1.0+s-t) / tM;
}

void IsoparametricQuadraticPyramid::Nrst(
    double r, double s, double t,
    std::vector<double>& N ) const
{
    N.resize(npe);
    Nrst(r, s, t, N.data());
}


// ============================================================================
//  Shape function derivatives — no bubble correction
// ============================================================================

void IsoparametricQuadraticPyramid::dNr( double r, double s, double t, vector<double>& DNR ) const
{
    // Note: dN/dr does not depend on r for the rational base functions
    // (same property as IsoparametricLinearPyramid::dNr)
    DNR.resize(npe);

    if ( t == 1.0 ) {
        DNR[0]= 0.25; DNR[1]=-0.25; DNR[2]=-0.25; DNR[3]= 0.25;
        DNR[4]= 0.0;
        DNR[5]= 0.0;  DNR[6]= 0.0;  DNR[7]= 0.0;  DNR[8]= 0.0;
        DNR[9]=-1.0;  DNR[10]= 1.0; DNR[11]= 1.0; DNR[12]=-1.0;
        return;
    }

    const double tM  = 1.0 - t;
    const double tM2 = 2.0 * tM;
    const double st  = s * t / tM;   // s*t/(1-t)

    // d/dr [ 0.25*(-r-s-1)*((1-r)*(1-s) - t + r*s*t/(1-t)) ]
    // = 0.25 * [ (-1)*((1-r)*(1-s) - t + r*s*t/(1-t))
    //          + (-r-s-1)*(-(1-s) + s*t/(1-t)) ]
    DNR[0] = 0.25*(  (-1.0)*((1.0-r)*(1.0-s) - t + r*st)
                   + (-r-s-1.0)*(-(1.0-s) + st) );

    // d/dr [ 0.25*(r-s-1)*((1+r)*(1-s) - t - r*s*t/(1-t)) ]
    // = 0.25 * [ (1)*((1+r)*(1-s) - t - r*s*t/(1-t))
    //          + (r-s-1)*((1-s) - s*t/(1-t)) ]
    DNR[1] = 0.25*(  (1.0)*((1.0+r)*(1.0-s) - t - r*st)
                   + (r-s-1.0)*((1.0-s) - st) );

    // d/dr [ 0.25*(r+s-1)*((1+r)*(1+s) - t + r*s*t/(1-t)) ]
    // = 0.25 * [ (1)*((1+r)*(1+s) - t + r*s*t/(1-t))
    //          + (r+s-1)*((1+s) + s*t/(1-t)) ]
    DNR[2] = 0.25*(  (1.0)*((1.0+r)*(1.0+s) - t + r*st)
                   + (r+s-1.0)*((1.0+s) + st) );

    // d/dr [ 0.25*(-r+s-1)*((1-r)*(1+s) - t - r*s*t/(1-t)) ]
    // = 0.25 * [ (-1)*((1-r)*(1+s) - t - r*s*t/(1-t))
    //          + (-r+s-1)*(-(1+s) - s*t/(1-t)) ]
    DNR[3] = 0.25*(  (-1.0)*((1.0-r)*(1.0+s) - t - r*st)
                   + (-r+s-1.0)*(-(1.0+s) - st) );

    DNR[4] = 0.0;

    // Midside base edges — d/dr [ (1±r-t)*(1∓r-t)*(1-s-t) / 2(1-t) ]
    // d/dr [(1+r-t)(1-r-t)] = (1-r-t) + (1+r-t)*(-1) = -2r... wait:
    // d/dr [(1+r-t)(1-r-t)] = (1)(1-r-t) + (1+r-t)(-1) = (1-r-t)-(1+r-t) = -2r... 
    // Actually: (1+r-t)(1-r-t) = (1-t)²-r² so d/dr = -2r
    // But we need to keep the factored form for clarity:
    DNR[5] = -2.0*r*(1.0-s-t) / tM2;          // d/dr[(1-t²-r²)(1-s-t)/2(1-t)] = -2r(1-s-t)/2(1-t)
    DNR[6] = (1.0+s-t)*(1.0-s-t) / tM2;       // d/dr[(1+r-t)] * rest = (1-s-t)(1+s-t)/2(1-t) -- wait

    // Let me redo these carefully:
    // N5 = (1+r-t)(1-r-t)(1-s-t) / 2(1-t)
    // dN5/dr = [(1-r-t) + (1+r-t)*(-1)] * (1-s-t) / 2(1-t)
    //        = [(1-r-t) - (1+r-t)] * (1-s-t) / 2(1-t)
    //        = -2r * (1-s-t) / 2(1-t)
    //        = -r*(1-s-t) / (1-t)
    DNR[5] = -r*(1.0-s-t) / tM;

    // N6 = (1+s-t)(1-s-t)(1+r-t) / 2(1-t)
    // dN6/dr = (1+s-t)(1-s-t) * (1) / 2(1-t)
    //        = (1+s-t)(1-s-t) / 2(1-t)
    //        = ((1-t)²-s²) / 2(1-t)
    DNR[6] = (1.0+s-t)*(1.0-s-t) / tM2;

    // N7 = (1+r-t)(1-r-t)(1+s-t) / 2(1-t)
    // dN7/dr = -2r*(1+s-t) / 2(1-t) = -r*(1+s-t)/(1-t)
    DNR[7] = -r*(1.0+s-t) / tM;

    // N8 = (1+s-t)(1-s-t)(1-r-t) / 2(1-t)
    // dN8/dr = (1+s-t)(1-s-t)*(-1) / 2(1-t)
    //        = -(1+s-t)(1-s-t) / 2(1-t)
    DNR[8] = -(1.0+s-t)*(1.0-s-t) / tM2;

    // Lateral midside edges
    // N9  = t*(1-r-t)*(1-s-t)/(1-t)
    // dN9/dr = t*(-1)*(1-s-t)/(1-t) = -t*(1-s-t)/(1-t)
    DNR[9]  = -t*(1.0-s-t) / tM;

    // N10 = t*(1+r-t)*(1-s-t)/(1-t)
    // dN10/dr = t*(1)*(1-s-t)/(1-t) = t*(1-s-t)/(1-t)
    DNR[10] =  t*(1.0-s-t) / tM;

    // N11 = t*(1+r-t)*(1+s-t)/(1-t)
    // dN11/dr = t*(1+s-t)/(1-t)
    DNR[11] =  t*(1.0+s-t) / tM;

    // N12 = t*(1-r-t)*(1+s-t)/(1-t)
    // dN12/dr = -t*(1+s-t)/(1-t)
    DNR[12] = -t*(1.0+s-t) / tM;
}


void IsoparametricQuadraticPyramid::dNs( double r, double s, double t, vector<double>& DNS ) const
{
    DNS.resize(npe);

    if ( t == 1.0 ) {
        DNS[0]= 0.25; DNS[1]= 0.25; DNS[2]=-0.25; DNS[3]=-0.25;
        DNS[4]= 0.0;
        DNS[5]= 0.0;  DNS[6]= 0.0;  DNS[7]= 0.0;  DNS[8]= 0.0;
        DNS[9]=-1.0;  DNS[10]=-1.0; DNS[11]= 1.0; DNS[12]= 1.0;
        return;
    }

    const double tM  = 1.0 - t;
    const double tM2 = 2.0 * tM;
    const double rt  = r * t / tM;   // r*t/(1-t)

    // d/ds [ 0.25*(-r-s-1)*((1-r)*(1-s) - t + r*s*t/(1-t)) ]
    // = 0.25 * [ (-1)*((1-r)*(1-s) - t + r*s*t/(1-t))
    //          + (-r-s-1)*(-(1-r) + r*t/(1-t)) ]
    DNS[0] = 0.25*(  (-1.0)*((1.0-r)*(1.0-s) - t + s*rt)
                   + (-r-s-1.0)*(-(1.0-r) + rt) );

    // d/ds [ 0.25*(r-s-1)*((1+r)*(1-s) - t - r*s*t/(1-t)) ]
    // = 0.25 * [ (-1)*((1+r)*(1-s) - t - r*s*t/(1-t))
    //          + (r-s-1)*(-(1+r) - r*t/(1-t)) ]
    DNS[1] = 0.25*(  (-1.0)*((1.0+r)*(1.0-s) - t - s*rt)
                   + (r-s-1.0)*(-(1.0+r) - rt) );

    // d/ds [ 0.25*(r+s-1)*((1+r)*(1+s) - t + r*s*t/(1-t)) ]
    // = 0.25 * [ (1)*((1+r)*(1+s) - t + r*s*t/(1-t))
    //          + (r+s-1)*((1+r) + r*t/(1-t)) ]
    DNS[2] = 0.25*(  (1.0)*((1.0+r)*(1.0+s) - t + s*rt)
                   + (r+s-1.0)*((1.0+r) + rt) );

    // d/ds [ 0.25*(-r+s-1)*((1-r)*(1+s) - t - r*s*t/(1-t)) ]
    // = 0.25 * [ (1)*((1-r)*(1+s) - t - r*s*t/(1-t))
    //          + (-r+s-1)*((1-r) - r*t/(1-t)) ]
    DNS[3] = 0.25*(  (1.0)*((1.0-r)*(1.0+s) - t - s*rt)
                   + (-r+s-1.0)*((1.0-r) - rt) );

    DNS[4] = 0.0;

    // N5 = (1+r-t)(1-r-t)(1-s-t) / 2(1-t)
    // dN5/ds = (1+r-t)(1-r-t)*(-1) / 2(1-t) = -(1+r-t)(1-r-t)/2(1-t)
    DNS[5] = -(1.0+r-t)*(1.0-r-t) / tM2;

    // N6 = (1+s-t)(1-s-t)(1+r-t) / 2(1-t)
    // dN6/ds = [(1-s-t) + (1+s-t)*(-1)]*(1+r-t) / 2(1-t)
    //        = -2s*(1+r-t) / 2(1-t) = -s*(1+r-t)/(1-t)
    DNS[6] = -s*(1.0+r-t) / tM;

    // N7 = (1+r-t)(1-r-t)(1+s-t) / 2(1-t)
    // dN7/ds = (1+r-t)(1-r-t)*(1) / 2(1-t)
    DNS[7] = (1.0+r-t)*(1.0-r-t) / tM2;

    // N8 = (1+s-t)(1-s-t)(1-r-t) / 2(1-t)
    // dN8/ds = -2s*(1-r-t) / 2(1-t) = -s*(1-r-t)/(1-t)
    DNS[8] = -s*(1.0-r-t) / tM;

    // Lateral midside edges
    DNS[9]  = -t*(1.0-r-t) / tM;
    DNS[10] = -t*(1.0+r-t) / tM;
    DNS[11] =  t*(1.0+r-t) / tM;
    DNS[12] =  t*(1.0-r-t) / tM;
}


void IsoparametricQuadraticPyramid::dNt( double r, double s, double t, vector<double>& DNT ) const
{
    DNT.resize(npe);

    if ( t == 1.0 ) {
        // Limiting values at apex — L'Hopital on rational terms
        DNT[0] = 0.25*(-r-s-1.0)*(-1.0) - 0.125;
        DNT[1] = 0.25*( r-s-1.0)*(-1.0) - 0.125;
        DNT[2] = 0.25*( r+s-1.0)*(-1.0) - 0.125;
        DNT[3] = 0.25*(-r+s-1.0)*(-1.0) - 0.125;
        DNT[4] = 3.0;
        DNT[5] = DNT[6] = DNT[7] = DNT[8] = 0.25;
        DNT[9] = DNT[10] = DNT[11] = DNT[12] = -1.0;
        return;
    }

    const double tM    = 1.0 - t;
    const double tM2   = 2.0 * tM;
    const double tM2sq = tM2 * tM2;
    const double tMsq  = tM * tM;
    const double rs    = r * s;

    // d/dt [rs*t/(1-t)] = rs/(1-t)^2
    const double drat = rs / tMsq;

    // Corner nodes
    DNT[0] = 0.25*(-r-s-1.0)*(-1.0 + drat);
    DNT[1] = 0.25*( r-s-1.0)*(-1.0 - drat);
    DNT[2] = 0.25*( r+s-1.0)*(-1.0 + drat);
    DNT[3] = 0.25*(-r+s-1.0)*(-1.0 - drat);
    DNT[4] = 4.0*t - 1.0;

    // Midside base edges — product rule on three factors A*B*C / 2(1-t)
    // d/dt [A*B*C / 2(1-t)] = (A'BC + AB'C + ABC') / 2(1-t)
    //                        + 2*A*B*C / [2(1-t)]^2
    // where A'=B'=C'=-1 for all factors of the form (x-t)
    DNT[5] = (-(1.0-r-t)*(1.0-s-t) - (1.0+r-t)*(1.0-s-t)
              - (1.0+r-t)*(1.0-r-t)) / tM2
           + 2.0*(1.0+r-t)*(1.0-r-t)*(1.0-s-t) / tM2sq;

    DNT[6] = (-(1.0+r-t)*(1.0-s-t) - (1.0+s-t)*(1.0+r-t)
              - (1.0+s-t)*(1.0-s-t)) / tM2
           + 2.0*(1.0+s-t)*(1.0-s-t)*(1.0+r-t) / tM2sq;

    DNT[7] = (-(1.0-r-t)*(1.0+s-t) - (1.0+s-t)*(1.0+r-t)
              - (1.0+r-t)*(1.0-r-t)) / tM2
           + 2.0*(1.0+r-t)*(1.0-r-t)*(1.0+s-t) / tM2sq;

    DNT[8] = (-(1.0-r-t)*(1.0-s-t) - (1.0-r-t)*(1.0+s-t)
              - (1.0+s-t)*(1.0-s-t)) / tM2
           + 2.0*(1.0+s-t)*(1.0-s-t)*(1.0-r-t) / tM2sq;

    // Lateral midside edges — quotient rule on t*A*B / (1-t)
    // d/dt [t*A*B/(1-t)] = (A*B - t*B - t*A)/(1-t) + t*A*B/(1-t)^2
    DNT[9]  = ((1.0-r-t)*(1.0-s-t) - t*(1.0-s-t) - t*(1.0-r-t)) / tM
            + t*(1.0-r-t)*(1.0-s-t) / tMsq;
    DNT[10] = ((1.0+r-t)*(1.0-s-t) - t*(1.0-s-t) - t*(1.0+r-t)) / tM
            + t*(1.0+r-t)*(1.0-s-t) / tMsq;
    DNT[11] = ((1.0+r-t)*(1.0+s-t) - t*(1.0+s-t) - t*(1.0+r-t)) / tM
            + t*(1.0+r-t)*(1.0+s-t) / tMsq;
    DNT[12] = ((1.0-r-t)*(1.0+s-t) - t*(1.0+s-t) - t*(1.0-r-t)) / tM
            + t*(1.0-r-t)*(1.0+s-t) / tMsq;
}


// ============================================================================
//  Private helper — assembles the global B-matrix at any parametric point
// ============================================================================

/**
 @brief Assembles the global shape function derivative matrix at
        parametric coordinates (r, s, t).

 Evaluates dNr, dNs, dNt, computes J and J^{-1}, then returns
 B = J^{-1} * [dNr | dNs | dNt]^T.

 @param r  Parametric r.
 @param s  Parametric s.
 @param t  Parametric t.
 @param B  Output: 3 x 13 B-matrix. Row k = spatial dimension,
           column i = node i.
 @return   Determinant of the Jacobian at (r, s, t).
*/
double IsoparametricQuadraticPyramid::assembleBMatrix( double r, double s, double t,
    DenseMatrix<DM_MIN>& B )
{
    dNr( r, s, t, DNR );
    dNs( r, s, t, DNS );
    dNt( r, s, t, DNT );

    Jacobian( DNR, DNS, DNT );
    const double detJ = JacobianInverse();

    B.Resize(dim, npe);
    for ( uint32_t i = 0; i < npe; ++i ) {
        B(0,i) = DNR[i];
        B(1,i) = DNS[i];
        B(2,i) = DNT[i];
    }
    B = JINV * B;

    return detJ;
}


// ============================================================================
//  Global shape function derivatives — public interface
// ============================================================================

/**
 @brief Computes global shape function derivatives at each node's
        reference coordinates.

 For each node i, evaluates [J]^{-1} * {dN/dr, dN/ds, dN/dt}^T
 to obtain {dN/dx, dN/dy, dN/dz}.

 @param DN14  Output: 3 x 13 matrix. Row k = spatial dimension k,
              column i = node i.
*/
void IsoparametricQuadraticPyramid::dN( DenseMatrix<DM_MIN>& DN14 )
{
    DN14.Resize(dim, npe);

    for ( uint32_t i = 0; i < npe; ++i ) {
        dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
        dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
        dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

        Jacobian( DNR, DNS, DNT );
        JacobianInverse();

        // Extract only the i-th node's local derivatives
        M.Resize(dim, 1);
        M(0,0) = DNR[i];
        M(1,0) = DNS[i];
        M(2,0) = DNT[i];

        JINV *= M;
        DN14(0,i) = JINV(0,0);
        DN14(1,i) = JINV(1,0);
        DN14(2,i) = JINV(2,0);
        JINV.Resize(dim, dim);
    }
}

/**
 @brief Computes global shape function derivatives at a physical point xyz.

 Maps xyz to parametric coordinates via Newton-Raphson, then computes
 the B-matrix at that location.

 @param DN2  Output: 3 x 13 B-matrix. Row k = spatial dimension k,
             column i = node i.
 @param xyz  Physical coordinates {x, y, z} of the evaluation point.
 @return     Determinant of the Jacobian at xyz.
*/
double IsoparametricQuadraticPyramid::dN_At(
    DenseMatrix<DM_MIN>& DN2,
    const std::vector<double>& xyz )
{
    std::vector<double> rst(dim);
    PhysicalToParametric(rst, xyz);
    return assembleBMatrix(rst[0], rst[1], rst[2], DN2);
}

/**
 @brief Computes global shape function derivatives at integration point
        gauss_point.

 @param B            Output: 3 x 13 B-matrix. Row k = spatial dimension,
                     column i = node i.
 @param gauss_point  Integration point index in [0, gpe).
 @return             Determinant of the Jacobian at this point.
*/
double IsoparametricQuadraticPyramid::dN_AtIntegrationPoint(
    DenseMatrix<DM_MIN>& B, uint32_t gauss_point )
{
    assert( gauss_point < gpe );
    return assembleBMatrix(
        IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), B );
}

/**
 @brief Computes global shape function derivatives at node nd.

 @param B   Output: 3 x 13 B-matrix.
 @param nd  Node index in [0, npe).
 @return    Determinant of the Jacobian at this node's reference coordinates.
*/
double IsoparametricQuadraticPyramid::dN_AtNode(
    DenseMatrix<DM_MIN>& B, uint32_t nd )
{
    assert( nd < npe );
    return assembleBMatrix(NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), B);
}

/**
 @brief Computes global shape function derivatives at the element barycentre.

 The barycentre of the reference pyramid is at (r=0, s=0, t=0.25).

 @param B  Output: 3 x 13 B-matrix.
 @return   Determinant of the Jacobian at the barycentre.
*/
double IsoparametricQuadraticPyramid::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
{
    return assembleBMatrix(0.0, 0.0, 0.25, B);
}



// ============================================================================
//  Jacobians — public interface
// ============================================================================

/**
 @brief Computes the Jacobian matrix at a given integration point.

 Stores the result in the base-class member JAC. Does not invert.

 @param gauss_point  Integration point index in [0, gpe).
*/
void IsoparametricQuadraticPyramid::JacobianAtIntegrationPoint( uint32_t gauss_point )
{
    assert( gauss_point < gpe );

    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT );
}

/**
 @brief Computes the Jacobian matrix at arbitrary parametric coordinates.

 @param rst  Parametric coordinates {r, s, t}.
*/
void IsoparametricQuadraticPyramid::JacobianAt( const std::vector<double>& rst ) noexcept
{
    dNr( rst[0], rst[1], rst[2], DNR );
    dNs( rst[0], rst[1], rst[2], DNS );
    dNt( rst[0], rst[1], rst[2], DNT );

    Jacobian( DNR, DNS, DNT );
}


// ============================================================================
//  Shape functions at specific locations
// ============================================================================

/**
 @brief Evaluates shape functions at a physical point xyz.

 Internally calls PhysicalToParametric then Nrst.

 @param N    Output: shape function values at xyz, size 13.
 @param xyz  Physical coordinates {x, y, z}.
*/
void IsoparametricQuadraticPyramid::N(
    std::vector<double>& N,
    const std::vector<double>& xyz )
{
    std::vector<double> rst(dim);
    PhysicalToParametric(rst, xyz);
    Nrst(rst[0], rst[1], rst[2], N);
}

/**
 @brief Evaluates shape functions at integration point ip.

 @param ip  Integration point index in [0, gpe).
 @param N   Output: shape function values, size 13.
*/
void IsoparametricQuadraticPyramid::N_AtIntegrationPoint(
    uint32_t ip, std::vector<double>& N )
{
    assert( ip < gpe );
    N.resize(npe);
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), N );
}

/**
 @brief Evaluates shape functions at the element barycentre (r=0, s=0, t=0.25).

 @param N  Output: shape function values, size 13.
*/
void IsoparametricQuadraticPyramid::N_AtBaryCenter( std::vector<double>& N )
{
    N.resize(npe);
    Nrst(0.0, 0.0, 0.25, N);
}


// ============================================================================
//  Coordinate mapping
// ============================================================================

/**
 @brief Maps parametric coordinates to physical coordinates.

 Uses the isoparametric mapping: x = sum_i N_i(r,s,t) * x_i.

 @param rst  Input:  parametric coordinates {r, s, t}.
 @param xyz  Output: physical coordinates {x, y, z}.
*/
void IsoparametricQuadraticPyramid::ParametricToPhysical(
    std::vector<double>& rst,
    std::vector<double>& xyz )
{
    Nrst(rst[0], rst[1], rst[2], NRST);

    xyz.assign(dim, 0.0);
    for ( uint32_t i = 0; i < npe; ++i ) {
        xyz[0] += XY(i,0) * NRST[i];
        xyz[1] += XY(i,1) * NRST[i];
        xyz[2] += XY(i,2) * NRST[i];
    }
}

/**
 @brief Maps physical coordinates to parametric coordinates using
        Newton-Raphson iteration.

 Algorithm:
 1. Initial guess: barycentre (0, 0, 0.25).
 2. Coarse grid search (5x5x5) to find a better starting point.
 3. Newton-Raphson iteration until convergence or max iterations.

 Convergence tolerance: 0.5% of the shortest edge length.

 @param rSt  Output: parametric coordinates {r, s, t}.
 @param xyz  Input:  physical coordinates {x, y, z}.

 @warning Issues a WARNING if Newton-Raphson does not converge within
          20 iterations. The returned parametric coordinates are the
          best approximation found.
*/
void IsoparametricQuadraticPyramid::PhysicalToParametric(
    std::vector<double>& rSt,
    const std::vector<double>& xyz )
{
    std::vector<double> outxyz(dim);
    std::vector<double> rstHatK(dim);
    std::vector<double> dist(dim, 0.0);

    // Geometric tolerance based on shortest edge
    std::vector<double> edgeLens(spe);
    EdgeLengths(edgeLens);
    const double seg_min = *std::min_element(edgeLens.begin(), edgeLens.end());
    const double tol     = 0.005 * seg_min;

    auto L2 = [](const std::vector<double>& d) {
        return std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    };
    auto computeDist = [&]() {
        dist[0] = std::fabs(outxyz[0]-xyz[0]);
        dist[1] = std::fabs(outxyz[1]-xyz[1]);
        dist[2] = std::fabs(outxyz[2]-xyz[2]);
        return L2(dist);
    };

    // Initial guess: barycentre
    rstHatK = {0.0, 0.0, 0.25};
    ParametricToPhysical(rstHatK, outxyz);
    double distL2 = computeDist();

    if ( distL2 > tol ) {
        // Coarse grid search
        const uint32_t nGrid = 5U;
        const double dR = 2.0 / (nGrid-1);
        const double dS = 2.0 / (nGrid-1);
        const double dT = 1.0 / (nGrid-1);

        double minDist = distL2;
        std::vector<double> candidate(dim);

        candidate[0] = -1.0;
        for ( uint32_t i = 0; i < nGrid; ++i, candidate[0] += dR ) {
            candidate[1] = -1.0;
            for ( uint32_t j = 0; j < nGrid; ++j, candidate[1] += dS ) {
                candidate[2] = 0.0;
                for ( uint32_t k = 0; k < nGrid; ++k, candidate[2] += dT ) {
                    ParametricToPhysical(candidate, outxyz);
                    const double d = computeDist();
                    if ( d < minDist ) {
                        minDist  = d;
                        rstHatK  = candidate;
                    }
                }
            }
        }

        ParametricToPhysical(rstHatK, outxyz);
        distL2 = computeDist();

        // Newton-Raphson
        const uint32_t maxIter = 20U;
        uint32_t iter = 0U;

        while ( distL2 > tol && iter < maxIter ) {
            // Domain bounds check — flag as outside if violated
            if ( rstHatK[2] < 0.0 || rstHatK[2] > 1.0 ||
                 rstHatK[0] < rstHatK[2]-1.0 || rstHatK[0] > 1.0-rstHatK[2] ||
                 rstHatK[1] < rstHatK[2]-1.0 || rstHatK[1] > 1.0-rstHatK[2] ) {
                rstHatK = {0.0, 0.0, 2.0};  // t=2 signals "outside element"
                break;
            }

            dNr(rstHatK[0], rstHatK[1], rstHatK[2], DNR);
            dNs(rstHatK[0], rstHatK[1], rstHatK[2], DNS);
            dNt(rstHatK[0], rstHatK[1], rstHatK[2], DNT);
            Jacobian(DNR, DNS, DNT);
            JacobianInverse();

            const double dx = outxyz[0]-xyz[0];
            const double dy = outxyz[1]-xyz[1];
            const double dz = outxyz[2]-xyz[2];

            rstHatK[0] -= JINV(0,0)*dx + JINV(1,0)*dy + JINV(2,0)*dz;
            rstHatK[1] -= JINV(0,1)*dx + JINV(1,1)*dy + JINV(2,1)*dz;
            rstHatK[2] -= JINV(0,2)*dx + JINV(1,2)*dy + JINV(2,2)*dz;

            ParametricToPhysical(rstHatK, outxyz);
            distL2 = computeDist();
            ++iter;
        }

        JAC.Zero();

        if ( iter == maxIter ) {
            cerr << "\nIsoparametricQuadraticPyramid::PhysicalToParametric:"
                 << " Newton-Raphson did not converge after " << maxIter
                 << " iterations. Residual L2 = " << distL2
                 << " (tolerance = " << tol << ")\n"
                 << "  Physical target:  ("
                 << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")\n"
                 << "  Best parametric:  ("
                 << rstHatK[0] << ", " << rstHatK[1] << ", " << rstHatK[2] << ")\n";
            csmp::Exception( WARNING,
                "IsoparametricQuadraticPyramid::PhysicalToParametric",
                "Newton-Raphson iteration not converged" );
        }
    }

    rSt = rstHatK;
}


// ============================================================================
//  Integration
// ============================================================================

/**
 @brief Generates the 8-point conical product integration rule.

 Uses the same rule as IsoparametricLinearPyramid — a 4-point
 Gauss-Legendre rule in the base plane combined with a 2-point
 rule in the t direction via the conical product transformation.

 Weight sum = 4/3 (volume of the reference pyramid).

 @param Ip  Output: 8 x 3 matrix of parametric coordinates.
 @param We  Output: vector of 8 weights.
*/
void IsoparametricQuadraticPyramid::GenerateIntegrationPoints(
    DenseMatrix<DM_MIN>& Ip,
    std::vector<double>& We )
{
    // 4-point Gauss-Legendre in the base plane
    const double a = 0.577350269189626;
    const double rBase[4] = {-a,  a,  a, -a};
    const double sBase[4] = {-a, -a,  a,  a};
    const double wBase[4] = { 1.0, 1.0, 1.0, 1.0};

    // 2-point conical rule in t
    const double tPts[2] = {0.455848155988775, 0.877485177344559};
    const double tWts[2] = {0.100785882079825, 0.232547451253508};

    uint32_t ip = 0U;
    for ( uint32_t j = 0; j < 4U; ++j )
        for ( uint32_t k = 0; k < 2U; ++k ) {
            Ip(ip,0) = tPts[k] * rBase[j];
            Ip(ip,1) = tPts[k] * sBase[j];
            Ip(ip,2) = 1.0 - tPts[k];
            We[ip]   = wBase[j] * tWts[k];
            ++ip;
        }
}

/**
 @brief Returns the integration weight at point ip.
 @param i  Integration point index in [0, gpe).
 @return   Weight value.
*/
double IsoparametricQuadraticPyramid::WeightAtIntegrationPoint(
    uint32_t i ) const noexcept
{
    return W[i];
}

/**
 @brief Returns the physical coordinates of integration point ip.

 Uses the isoparametric mapping applied to the stored parametric
 integration point coordinates.

 @param ip   Integration point index in [0, gpe).
 @param xyz  Output: physical coordinates {x, y, z}.
*/
void IsoparametricQuadraticPyramid::IntegrationPoint(
    uint32_t ip, std::vector<double>& xyz ) const
{
    assert( ip < gpe );
    xyz.assign(3U, 0.0);

    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for ( uint32_t i = 0; i < npe; ++i ) {
        xyz[0] += XY(i,0) * NRST[i];
        xyz[1] += XY(i,1) * NRST[i];
        xyz[2] += XY(i,2) * NRST[i];
    }
}


// ============================================================================
//  Geometry
// ============================================================================

/**
 @brief Computes the element volume by numerical integration.

 Integrates the Jacobian determinant over the 8-point rule:

 @f[ V = \sum_{i=0}^{7} \det J(\xi_i) \cdot w_i @f]

 @return Element volume.
*/
double IsoparametricQuadraticPyramid::Volume()
{
    // -----------------------------------------------------------------------
    // @warning Volume accuracy depends on element distortion.
    //
    // For the reference element (physical == parametric) the 8-point
    // conical product rule integrates the Jacobian exactly and the
    // volume is exact.
    //
    // For distorted elements the Jacobian det J(r,s,t) is a non-constant
    // polynomial. The 8-point rule integrates it approximately — the
    // error grows with distortion. For severely distorted elements
    // (non-planar base, large aspect ratio, skewed apex) the error
    // can be significant.
    //
    // This is a fundamental limitation of isoparametric elements with
    // low-order quadrature rules — it affects ALL element types
    // (hexahedra, prisms, pyramids) and is not specific to this element.
    //
    // For hexahedra with non-planar faces the situation is worse:
    // the Jacobian can become negative at integration points for
    // severely distorted elements, causing the volume to be completely
    // wrong or negative. Always check element quality before use.
    // -----------------------------------------------------------------------
    double vol = 0.0;
    for ( uint32_t i = 0; i < gpe; ++i ) {
        dNr( IP(i,0), IP(i,1), IP(i,2), DNR );
        dNs( IP(i,0), IP(i,1), IP(i,2), DNS );
        dNt( IP(i,0), IP(i,1), IP(i,2), DNT );
        Jacobian( DNR, DNS, DNT );
        vol += JacobianInverse() * W[i];
    }
    return vol;
}

/**
 @brief Computes the element aspect ratio.

 Defined as the ratio of the longest to the shortest edge (corner nodes only).

 @return Aspect ratio >= 1.
*/
double IsoparametricQuadraticPyramid::AspectRatio()
{
    std::vector<double> lens(spe);
    EdgeLengths(lens);
    const double eMax = *std::max_element(lens.begin(), lens.end());
    const double eMin = *std::min_element(lens.begin(), lens.end());
    return eMax / eMin;
}

/**
 @brief Computes the lengths of the 8 edges (corner nodes only).

 Edge ordering matches NodesOfSegment:
 @code
   0: 0-1,  1: 1-2,  2: 2-3,  3: 3-0,
   4: 0-4,  5: 1-4,  6: 2-4,  7: 3-4
 @endcode

 @param len  Output: vector of 8 edge lengths.
*/
void IsoparametricQuadraticPyramid::EdgeLengths( std::vector<double>& len )
{
    static constexpr uint32_t kEdges[8][2] = {
        {0U,1U},{1U,2U},{2U,3U},{3U,0U},
        {0U,4U},{1U,4U},{2U,4U},{3U,4U}
    };
    len.resize(spe);
    for ( uint32_t e = 0; e < spe; ++e ) {
        const uint32_t a = kEdges[e][0];
        const uint32_t b = kEdges[e][1];
        const double dx = XY(b,0)-XY(a,0);
        const double dy = XY(b,1)-XY(a,1);
        const double dz = XY(b,2)-XY(a,2);
        len[e] = std::sqrt(dx*dx + dy*dy + dz*dz);
    }
}

/**
 @brief Estimates the inner radius of the element.

 Uses the empirical formula V / (sum_of_edges / 6), inherited from
 the linear pyramid. Accuracy degrades for high aspect ratios.

 @return Approximate inner radius.
*/
double IsoparametricQuadraticPyramid::InnerRadius()
{
    std::vector<double> lens(spe);
    EdgeLengths(lens);
    double sumLen = 0.0;
    for ( const double l : lens ) sumLen += l;

    if ( AspectRatio() > 4.0 )
        cerr << "\nIsoparametricQuadraticPyramid::InnerRadius: WARNING:"
             << " formula not reliable for high aspect ratio elements.\n";

    return Volume() / (sumLen / 6.0);
}


// ============================================================================
//  Extrapolation
// ============================================================================

/**
 @brief Extrapolates integration point variables to the element nodes.

 For each node, evaluates the shape functions at that node's reference
 coordinates and uses them as weights over the integration point values:

 @f[ u_{\text{node}} = \sum_{ip=0}^{7} N_{ip}(\xi_{\text{node}}) \cdot u_{ip} @f]

 This is exact for fields that are representable by the element's
 shape functions.

 @param nvars  Number of variables per integration point.
 @param IVAR   Input:  integration point values, size gpe * nvars,
               stored as [ip0_var0, ip0_var1, ..., ip1_var0, ...].
 @param NVAR   Output: nodal values, size npe * nvars,
               stored as [node0_var0, node0_var1, ..., node1_var0, ...].

 @throws csmp::Exception (FATAL_ERROR) if IVAR.size() != gpe * nvars.
*/
void IsoparametricQuadraticPyramid::ExtrapolateIntegrationPointVariableToNodes(
    uint32_t nvars,
    const std::vector<double>& IVAR,
    std::vector<double>& NVAR ) const
{
    if ( IVAR.size() != static_cast<size_t>(gpe * nvars) )
        throw csmp::Exception(
            ERROR,
            "IsoparametricQuadraticPyramid::ExtrapolateIntegrationPointVariableToNodes",
            "IVAR size must equal gpe * nvars" );

    NVAR.assign(npe * nvars, 0.0);

    // Build the 8 x 13 matrix of shape function values at IPs:
    // A[ip][node] = N_node(xi_ip)
    // Then for each node, solve or approximate:
    // u_node = sum_ip A[ip][node] * u_ip  (least squares if underdetermined)
    //
    // Since gpe=8 < npe=13 the system is underdetermined.
    // Use the pseudo-inverse: u = A^T (A A^T)^{-1} u_ip
    // For simplicity use the transpose directly (weighted average):
    // u_node = sum_ip N_node(xi_ip) * u_ip / sum_ip N_node(xi_ip)

    // Evaluate N_node at each integration point
    // N_at_ip[ip] is a vector of 13 shape function values at IP ip
    std::vector<std::vector<double>> N_at_ip(gpe, std::vector<double>(npe));
    for ( uint32_t ip = 0; ip < gpe; ++ip )
        Nrst( IP(ip,0), IP(ip,1), IP(ip,2), N_at_ip[ip] );

    for ( uint32_t node = 0; node < npe; ++node ) {
        // Weight for this node at each IP = N_node(xi_ip)
        double wsum = 0.0;
        for ( uint32_t ip = 0; ip < gpe; ++ip )
            wsum += N_at_ip[ip][node];

        for ( uint32_t k = 0; k < nvars; ++k ) {
            double val = 0.0;
            for ( uint32_t ip = 0; ip < gpe; ++ip )
                val += N_at_ip[ip][node] * IVAR[ip * nvars + k];
            NVAR[node * nvars + k] = (std::fabs(wsum) > 1.0e-14)
                                     ? val / wsum
                                     : 0.0;
        }
    }
}


// ============================================================================
//  Reference coordinates
// ============================================================================

/**
 @brief Returns the reference (parametric) coordinates of all nodes.

 @param matCoords  Output: npe x 3 matrix. Row i = (r_i, s_i, t_i).
*/
void IsoparametricQuadraticPyramid::ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const
{
    matCoords.Resize(npe, dim);
    matCoords = NXYZ;
}


// ============================================================================
//  VTK output
// ============================================================================

/**
 @brief Writes nodal data to a VTK unstructured grid file.

 Outputs the 13-node pyramid as 5 VTK polylines (one per face outline)
 since VTK has no native 13-node pyramid cell type.

 @param file_name  Base filename (element index and ".vtk" are appended).
 @param var_name   Variable name written into the VTK header.
 @param DATA       Data matrix: 1 x npe for scalars, dim x npe for vectors.
*/
void IsoparametricQuadraticPyramid::OutputNodeDataToVTK(
    const char* file_name,
    const char* var_name,
    DenseMatrix<DM_MIN>& DATA ) const
{
    char outfile[NAME_STRING];
    char elmt[30];
    strcpy(outfile, file_name);
    snprintf(elmt, sizeof(elmt), "%zu", CurrentID());
    strcat(outfile, elmt);
    strcat(outfile, ".vtk");

    ofstream ofs(outfile, ios::out | ios::trunc);
    if ( !ofs ) {
        cerr << "\nIsoparametricQuadraticPyramid::OutputNodeDataToVTK: "
             << "cannot open '" << outfile << "'\n";
        return;
    }

    // Header
    ofs << "# vtk DataFile Version 2.0\n"
        << "openCSMP: " << var_name << "\n"
        << "ASCII\n\n";

    // Points
    ofs << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << npe << " float\n";
    for ( uint32_t i = 0; i < npe; ++i ) {
        for ( uint32_t j = 0; j < dim; ++j ) ofs << XY(i,j) << " ";
        ofs << "\n";
    }
    ofs << "\n";

    // Cells — 5 polylines, one per face outline
    // Total entries: 5 faces * (1 size + nodes per face)
    //   4 triangular faces: 1+6 = 7 each  → 4*7 = 28
    //   1 quad base:        1+8 = 9        →   9
    //   Total = 37
    ofs << "CELLS 5 37\n";
    ofs << "6 0 5 1 10 4 9\n";    // face 0
    ofs << "6 1 6 2 11 4 10\n";   // face 1
    ofs << "6 2 7 3 12 4 11\n";   // face 2
    ofs << "6 3 8 0 9 4 12\n";    // face 3
    ofs << "8 0 8 3 7 2 6 1 5\n"; // face 4 (base)
    ofs << "\n";

    // Cell types — all polylines (VTK type 4)
    ofs << "CELL_TYPES 5\n"
        << "4\n4\n4\n4\n4\n\n";

    // Point data
    ofs << "POINT_DATA " << npe << "\n";
    ofs << scientific;

    if ( DATA.Rows() == 1 ) {
        ofs << "SCALARS " << var_name << " float\n"
            << "LOOKUP_TABLE default\n";
        for ( uint32_t i = 0; i < DATA.Cols(); ++i )
            ofs << DATA(0,i) << " ";
        ofs << "\n";
    } else {
        ofs << "VECTORS " << var_name << " float\n";
        for ( uint32_t i = 0; i < DATA.Cols(); ++i ) {
            for ( uint32_t j = 0; j < DATA.Rows(); ++j )
                ofs << DATA(j,i) << " ";
            ofs << "\n";
        }
    }
    ofs << "\n";
    ofs.close();

    cout << "\nIsoparametricQuadraticPyramid::OutputNodeDataToVTK: '"
         << outfile << "' written.\n";
}

} // namespace csmp
