// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IsoparametricQuadraticPrism_Test.h"
#include "Point.h"
#include "Node.h"
#include "Element.h"
#include "IsoparametricQuadraticPrism.h"
#include "DenseMatrix.h"
#include "LocalVariables.h"
#include "compareFloats.h"
#include <array>
#include <cmath>
#include <deque>
#include <set>
#include <vector>

using namespace std;

namespace csmp {

namespace {

constexpr uint32_t dim{3};
constexpr uint32_t n_nodes{18};

/// coordinates of the 18 nodes of a quadratic prism
using PrismCoordinates = std::array<std::array<double,3>,n_nodes>;

/**
   the prism in local (r,s,t) coordinates: triangle (0,0)-(1,0)-(0,1) in the r-s plane, -1 <= t <= 1; its volume is 1

   nodes 0-5: corners (bottom triangle, top triangle); 6-8: edges of the bottom triangle; 9-11: middle of the vertical edges;
   12-14: edges of the top triangle; 15-17: centres of the quadrilateral faces
*/
const PrismCoordinates reference_prism{{ {0.,0.,-1.},  {1.,0.,-1.},  {0.,1.,-1.},  {0.,0.,1.},  {1.,0.,1.},  {0.,1.,1.},
                                         {0.5,0.,-1.}, {0.5,0.5,-1.}, {0.,0.5,-1.},
                                         {0.,0.,0.},   {1.,0.,0.},   {0.,1.,0.},
                                         {0.5,0.,1.},  {0.5,0.5,1.},  {0.,0.5,1.},
                                         {0.5,0.,0.},  {0.5,0.5,0.},  {0.,0.5,0.} }};

/// reference prism stretched by the factors (2,3,0.5) and moved by (10,-5,7); volume = 2*3*0.5 = 3
PrismCoordinates stretchedAndMovedPrism()
 {
    PrismCoordinates xyz( reference_prism );
    for ( auto& node : xyz ) {
         node[0] = 2.0*node[0] + 10.;
         node[1] = 3.0*node[1] -  5.;
         node[2] = 0.5*node[2] +  7.;
      }
    return xyz;
 }

/// reference prism sheared parallel to its triangular faces (x += 0.3 z, y -= 0.2 z); the volume remains 1
PrismCoordinates shearedPrism()
 {
    PrismCoordinates xyz( reference_prism );
    for ( auto& node : xyz ) {
         node[0] += 0.3*node[2];
         node[1] -= 0.2*node[2];
      }
    return xyz;
 }

/// a prism finite element together with the nodes and the Element that make it usable
struct PrismFixture {
    LocalVariables               lvs;
    IsoparametricQuadraticPrism  prism;
    std::deque<Node<dim> >       nodes;    // deque: the addresses of the nodes must not change
    Element<dim>                 element;

    explicit PrismFixture( const PrismCoordinates& xyz )
     : element( &prism )
     {
        element.Idx(0);
        for ( uint32_t n=0; n<n_nodes; ++n ) {
             nodes.emplace_back( n, Point<dim>( xyz[n][0], xyz[n][1], xyz[n][2] ), lvs );
             element.Assign( n, &nodes.back() );
          }
        element.CoordinateMatrix(); // hands the node coordinates to the finite element
     }
 };

/// some points inside of the reference prism, in local coordinates
const std::vector<std::array<double,3> > interior_points{ {1./3.,1./3.,0.}, {0.1,0.2,-0.7}, {0.6,0.3,0.9}, {0.05,0.9,0.25}, {0.25,0.25,-1.} };

} // anonymous namespace


void IsoparametricQuadraticPrism_Test::run()
 {
    Topology_Test();
    InterpolationFunctions_Test();
    NumericalIntegration_Test();
    Volume_Test();
    GlobalDerivatives_Test();

    cout <<"\n"<<"IsoparametricQuadraticPrism_Test::run: complete."<< endl;
 }


/// numbers of nodes, faces and segments; corner nodes of the faces
void IsoparametricQuadraticPrism_Test::Topology_Test()
 {
    IsoparametricQuadraticPrism prism;

    _test( prism.ElementType() == ISOPARAMETRIC_QUADRATIC_PRISM18 );
    _test( prism.Nodes() == n_nodes );
    _test( prism.Faces() == 5U );
    _test( prism.Segments() == 9U );
    _test( prism.IntegrationPoints() == 6U );

    // 2 triangular and 3 quadrilateral faces, identified by their corner nodes
    uint32_t triangles{0U}, quadrilaterals{0U};
    std::multiset<std::pair<uint32_t,uint32_t> > edges_of_faces;
    for ( uint32_t face=0; face<prism.Faces(); ++face ) {
         const std::vector<uint32_t> corners = prism.CornerNodesOfFace( face );
         if      ( corners.size() == 3U ) ++triangles;
         else if ( corners.size() == 4U ) ++quadrilaterals;
         // the corner nodes of a face are corner nodes of the element, and distinct
         const std::set<uint32_t> distinct( corners.begin(), corners.end() );
         _test( distinct.size() == corners.size() );
         for ( const uint32_t node : corners ) _test( node < 6U );
         // consecutive corner nodes of a face form its edges
         for ( size_t i=0; i<corners.size(); ++i ) {
              const uint32_t a = corners[i], b = corners[(i+1)%corners.size()];
              edges_of_faces.insert( std::make_pair( std::min(a,b), std::max(a,b) ) );
           }
         // all nodes of the face are valid and distinct, and include its corner nodes
         const std::vector<uint32_t> face_nodes = prism.NodesOfFace( face );
         const std::set<uint32_t> distinct_face_nodes( face_nodes.begin(), face_nodes.end() );
         _test( distinct_face_nodes.size() == face_nodes.size() );
         for ( const uint32_t node : face_nodes ) _test( node < n_nodes );
         for ( const uint32_t node : corners )    _test( distinct_face_nodes.count( node ) == 1U );
      }
    _test( triangles == 2U );
    _test( quadrilaterals == 3U );

    // closed surface: there are 9 edges and each of them is shared by exactly 2 faces
    const std::set<std::pair<uint32_t,uint32_t> > distinct_edges( edges_of_faces.begin(), edges_of_faces.end() );
    _test( distinct_edges.size() == 9U );
    bool every_edge_shared_by_two_faces{true};
    for ( const auto& edge : distinct_edges )
      if ( edges_of_faces.count( edge ) != 2U ) every_edge_shared_by_two_faces = false;
    _test( every_edge_shared_by_two_faces );
 }


/// interpolation functions and their derivatives in local coordinates
void IsoparametricQuadraticPrism_Test::InterpolationFunctions_Test()
 {
    IsoparametricQuadraticPrism prism;
    const double tol{1.0e-13};

    // local coordinates of the nodes as stored by the element
    DenseMatrix<DM_MIN> rst;
    prism.ReferenceCoordinates( rst );
    _test( rst.Rows() == n_nodes && rst.Cols() == 3U );
    for ( uint32_t node=0; node<n_nodes; ++node )
      for ( uint32_t d=0; d<dim; ++d )
        _equal( rst(node,d), reference_prism[node][d], tol );

    // N_i equals 1 at node i and 0 at all other nodes
    std::vector<double> N;
    for ( uint32_t node=0; node<n_nodes; ++node ) {
         prism.Nrst( reference_prism[node][0], reference_prism[node][1], reference_prism[node][2], N );
         _test( N.size() == n_nodes );
         for ( uint32_t i=0; i<n_nodes; ++i )
           _equal( N[i], ( i == node ? 1. : 0. ), tol );
      }

    // a complete quadratic polynomial must be interpolated exactly
    auto quadratic = []( double r, double s, double t ) {
        return 1.5 - 2.*r + 0.5*s + 3.*t + 0.7*r*r - 1.1*s*s + 0.4*t*t + 2.*r*s - 0.6*r*t + 0.9*s*t; };
    std::array<double,n_nodes> q;
    for ( uint32_t n=0; n<n_nodes; ++n ) q[n] = quadratic( reference_prism[n][0], reference_prism[n][1], reference_prism[n][2] );

    std::vector<double> dNr, dNs, dNt;
    for ( const auto& point : interior_points ) {
         // partition of unity
         prism.Nrst( point[0], point[1], point[2], N );
         double sum{0.}, interpolated{0.};
         for ( uint32_t i=0; i<n_nodes; ++i ) {
              sum          += N[i];
              interpolated += N[i] * q[i];
           }
         _equal( sum, 1., tol );
         _equal( interpolated, quadratic( point[0], point[1], point[2] ), tol );

         // the pointer version returns the same values
         double N_array[n_nodes];
         prism.Nrst( point[0], point[1], point[2], N_array );
         for ( uint32_t i=0; i<n_nodes; ++i ) _equal( N_array[i], N[i], tol );

         // the derivatives of a partition of unity add up to zero
         prism.dNr( point[0], point[1], point[2], dNr );
         prism.dNs( point[0], point[1], point[2], dNs );
         prism.dNt( point[0], point[1], point[2], dNt );
         double sum_r{0.}, sum_s{0.}, sum_t{0.};
         for ( uint32_t i=0; i<n_nodes; ++i ) { sum_r += dNr[i]; sum_s += dNs[i]; sum_t += dNt[i]; }
         _equal( sum_r, 0., tol );
         _equal( sum_s, 0., tol );
         _equal( sum_t, 0., tol );

         // the derivatives agree with central finite differences of the interpolation functions
         const double h{1.0e-6};
         std::vector<double> Np, Nm;
         prism.Nrst( point[0]+h, point[1], point[2], Np ); prism.Nrst( point[0]-h, point[1], point[2], Nm );
         for ( uint32_t i=0; i<n_nodes; ++i ) _equal( dNr[i], (Np[i]-Nm[i])/(2.*h), 1.0e-7 );
         prism.Nrst( point[0], point[1]+h, point[2], Np ); prism.Nrst( point[0], point[1]-h, point[2], Nm );
         for ( uint32_t i=0; i<n_nodes; ++i ) _equal( dNs[i], (Np[i]-Nm[i])/(2.*h), 1.0e-7 );
         prism.Nrst( point[0], point[1], point[2]+h, Np ); prism.Nrst( point[0], point[1], point[2]-h, Nm );
         for ( uint32_t i=0; i<n_nodes; ++i ) _equal( dNt[i], (Np[i]-Nm[i])/(2.*h), 1.0e-7 );
      }
 }


/// integration points and weights of the 6-point rule
void IsoparametricQuadraticPrism_Test::NumericalIntegration_Test()
 {
    const double tol{1.0e-13};
    PrismFixture fixture( reference_prism );
    IsoparametricQuadraticPrism& prism = fixture.prism;

    double sum_of_weights{0.};
    std::vector<double> rst, N;
    for ( uint32_t ip=0; ip<prism.IntegrationPoints(); ++ip ) {
         _test( prism.WeightAtIntegrationPoint( ip ) > 0. );
         sum_of_weights += prism.WeightAtIntegrationPoint( ip );

         // the integration points are located inside of the reference prism
         prism.IntegrationPoint( ip, rst );
         _test( rst.size() == 3U );
         _test( rst[0] > 0. && rst[1] > 0. && rst[0] + rst[1] < 1. );
         _test( rst[2] > -1. && rst[2] < 1. );

         // partition of unity at the integration points
         prism.N_AtIntegrationPoint( ip, N );
         double sum{0.};
         for ( const double value : N ) sum += value;
         _equal( sum, 1., tol );

         // ... where the interpolation functions are those in local coordinates
         std::vector<double> Nrst;
         prism.Nrst( rst[0], rst[1], rst[2], Nrst );
         for ( uint32_t i=0; i<n_nodes; ++i ) _equal( N[i], Nrst[i], tol );
      }
    // the weights add up to the volume of the reference prism (triangle area 1/2 x height 2)
    _equal( sum_of_weights, 1., 1.0e-14 );
 }


/// volumes of the reference prism and of affinely transformed prisms, where the integration is exact
void IsoparametricQuadraticPrism_Test::Volume_Test()
 {
    const double tol{1.0e-13};
    { PrismFixture fixture( reference_prism );
      _equal( fixture.prism.Volume(), 1., tol );
      _equal( fixture.element.Volume(), 1., tol ); }

    { PrismFixture fixture( stretchedAndMovedPrism() );
      _equal( fixture.prism.Volume(), 3., tol );
      _equal( fixture.element.Volume(), 3., tol ); }

    { PrismFixture fixture( shearedPrism() );
      _equal( fixture.prism.Volume(), 1., tol );
      _equal( fixture.element.Volume(), 1., tol ); }
 }


/// derivatives of the interpolation functions in global coordinates
void IsoparametricQuadraticPrism_Test::GlobalDerivatives_Test()
 {
    const double tol{1.0e-11};
    // a linear field f = a + b x + c y + d z must be differentiated exactly
    const double a{2.5}, b{-1.5}, c{0.75}, d{4.0};

    for ( const PrismCoordinates& xyz : { reference_prism, stretchedAndMovedPrism(), shearedPrism() } ) {
         PrismFixture fixture( xyz );
         IsoparametricQuadraticPrism& prism = fixture.prism;

         std::array<double,n_nodes> f;
         for ( uint32_t n=0; n<n_nodes; ++n ) f[n] = a + b*xyz[n][0] + c*xyz[n][1] + d*xyz[n][2];

         DenseMatrix<DM_MIN> DN;
         double volume{0.};
         for ( uint32_t ip=0; ip<prism.IntegrationPoints(); ++ip ) {
              const double detJ = prism.dN_AtIntegrationPoint( DN, ip );
              _test( detJ > 0. );
              volume += detJ * prism.WeightAtIntegrationPoint( ip );
              _test( DN.Rows() == 3U && DN.Cols() == n_nodes );

              double gradient[3] = { 0., 0., 0. }, sum[3] = { 0., 0., 0. };
              for ( uint32_t n=0; n<n_nodes; ++n )
                for ( uint32_t k=0; k<dim; ++k ) {
                     gradient[k] += DN(k,n) * f[n];
                     sum[k]      += DN(k,n);
                  }
              _equal( gradient[0], b, tol );
              _equal( gradient[1], c, tol );
              _equal( gradient[2], d, tol );
              _equal( sum[0], 0., tol );
              _equal( sum[1], 0., tol );
              _equal( sum[2], 0., tol );
           }
         _equal( volume, prism.Volume(), tol );

         // the same at the barycentre
         const double detJ = prism.dN_AtBarycenter( DN );
         _test( detJ > 0. );
         double gradient[3] = { 0., 0., 0. };
         for ( uint32_t n=0; n<n_nodes; ++n )
           for ( uint32_t k=0; k<dim; ++k ) gradient[k] += DN(k,n) * f[n];
         _equal( gradient[0], b, tol );
         _equal( gradient[1], c, tol );
         _equal( gradient[2], d, tol );

         // interpolation at the barycentre returns the value of the linear field there
         std::vector<double> N;
         prism.N_AtBaryCenter( N );
         double value{0.}, centre[3] = { 0., 0., 0. };
         for ( uint32_t n=0; n<n_nodes; ++n ) value += N[n] * f[n];
         for ( uint32_t n=0; n<6U; ++n )
           for ( uint32_t k=0; k<dim; ++k ) centre[k] += xyz[n][k] / 6.;
         _equal( value, a + b*centre[0] + c*centre[1] + d*centre[2], tol );
      }
 }


} // csmp
