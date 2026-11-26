// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  MeshManagementUtilities_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//

#include "MeshManagementUtilities_Test.h"
#include "meshManagementUtilities.h"
#include "CSMP_definitions.h"
#include "ErrorHandler.h"
#include "vsetMakers.h"
#include "VTU_Interface.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
//#include "NodeManifold.h"
#include "compareFloats.h"

#include "VTK_Interface.h"
#include "Node.h"
#include "Point.h"
#include "LocalVariables.h"
#include "ModelTopology.h"
#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearTetrahedron.h"
#include <array>
#include <cmath>
#include <deque>
#include <map>
#include <set>
#include <vector>

// #define CSMP_MESH_MANAGER_TEST_DEBUG

/*
   Open questions (functions that are deliberately not, or only partly, asserted below; SKM to decide):

   1. connectivityCheck: reports every hexahedron or pyramid with <= 3 neighbors and every prism with <= 2 neighbors
      as an error. The 8 corner cells of any box-shaped hexahedral mesh have exactly 3 neighbors, so that the regular
      3x3x3 mesh used below yields 8 'errors'. Should the limits be '< 3' and '< 2'?
   2. hasNonManifoldVertices: returns true for a healthy tetrahedron (positive determinant) and false for a flat or an
      inverted one, i.e. the opposite of what its name and comment say.
   3. findInterconnectedNodes: always returns nodes.size(), whatever the traversal visited.
   4. parentElementsSharingMultipleEdgeNodes: overwrites node and element indices, needs pointer-sorted input, records
      only the first volumetric parent of a node (which need not contain the second edge node) and always prints the
      'Segment ID could not be found' message.
   5. findContiguousMeshPatches returns 1 for a contiguous mesh, but (number of patches + 1) otherwise;
      findPointersToStandAloneMeshPatches returns 0 for a contiguous mesh, (number of patches + 1) otherwise, and does
      not record the first patch in its map.
*/

using namespace std;

namespace csmp {

namespace {

/// single element with its own nodes and finite element, built in place (no Model needed)
template<class FE>
struct ElementFixture {
    LocalVariables        lvs;
    FE                    fe;
    std::deque<Node<3> >  nodes;    // deque: the addresses of the nodes must not change
    Element<3>            element;

    template<class... FE_Args>
    explicit ElementFixture( const std::vector<std::array<double,3> >& xyz, FE_Args... fe_args )
     : fe( fe_args... ), element( &fe )
     {
        element.Idx(0);
        uint32_t n{0U};
        for ( const auto& p : xyz ) {
             nodes.emplace_back( n, Point<3>( p[0], p[1], p[2] ), lvs );
             element.Assign( n, &nodes.back() );
             ++n;
          }
        element.CoordinateMatrix(); // hands the node coordinates to the finite element
     }
};

typedef ElementFixture<IsoparametricLinearLineElement>  LineFixture;
typedef ElementFixture<IsoparametricLinearTriangle>     TriangleFixture;
typedef ElementFixture<IsoparametricLinearTetrahedron>  TetFixture;

/// finds the cell with the given barycentre by linear search; nullptr if there is none
template<uint32_t dim>
Element<dim>* cellWithBaryCenter( const Region<dim>& region, const Point<dim>& bctr )
 {
    for ( auto c : region.CellVector() )
      if ( bctr.DistanceTo( c->BaryCenter() ) < 1.0e-9 ) return c;
    return nullptr;
 }

/// finds the node with the given coordinate by linear search; nullptr if there is none
template<uint32_t dim>
Node<dim>* nodeAt( const Region<dim>& region, const Point<dim>& pt )
 {
    for ( auto nit=region.NodesBegin(); nit!=region.NodesEnd(); ++nit )
      if ( pt.DistanceTo( (*nit)->Coordinate() ) < 1.0e-9 ) return (*nit);
    return nullptr;
 }

const double tol{ 1.0e-12 };

} // anonymous namespace


void MeshManagementUtilities_Test::run()
 {
    Test_transformMeshIntoRefinedLinearAndQuadraticMeshes();
    Test_StandAloneElements();
    Test_HexahedraMesh();
    Test_TetrahedraMesh();
    Test_DisconnectedMesh2D();
    Test_MeshWithLineElements2D();
    Test_FracturedBox();
 }


/// functions that take individual cells; all inputs are constructed here so that the results are known exactly
void MeshManagementUtilities_Test::Test_StandAloneElements()
 {
    // --- angle between surface cells = angle between their unit normals [degrees]
    {
      TriangleFixture  xy(      { {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.} }, 3U, 3U ); // normal +z
      TriangleFixture  xz(      { {0.,0.,0.}, {0.,0.,1.}, {1.,0.,0.} }, 3U, 3U ); // normal +y
      TriangleFixture  flipped( { {0.,0.,0.}, {0.,1.,0.}, {1.,0.,0.} }, 3U, 3U ); // normal -z
      TriangleFixture  tilted(  { {0.,0.,0.}, {1.,0.,0.}, {0.,1.,1.} }, 3U, 3U ); // rotated by 45 degrees about x

      _equal( angleBetweenSurfaceCells<Element>( &xy.element, &xy.element ),        0., 1.0e-9 );
      _equal( angleBetweenSurfaceCells<Element>( &xy.element, &xz.element ),       90., 1.0e-9 );
      _equal( angleBetweenSurfaceCells<Element>( &xy.element, &tilted.element ),   45., 1.0e-9 );
      _equal( angleBetweenSurfaceCells<Element>( &xy.element, &flipped.element ), 180., 1.0e-9 );
      // symmetry
      _equal( angleBetweenSurfaceCells<Element>( &tilted.element, &xy.element ),   45., 1.0e-9 );
    }

    // --- angle between line cells = angle between the vectors from node 0 to node 1 [degrees]
    {
      LineFixture  along_x(      { {0.,0.,0.}, {1.,0.,0.} }, 3U );
      LineFixture  continuation( { {1.,0.,0.}, {2.,0.,0.} }, 3U );
      LineFixture  along_y(      { {1.,0.,0.}, {1.,1.,0.} }, 3U );
      LineFixture  reversed(     { {1.,0.,0.}, {0.,0.,0.} }, 3U );
      LineFixture  diagonal(     { {0.,0.,0.}, {1.,1.,1.} }, 3U );

      _equal( (angleBetweenLineCells<3,Element>( &along_x.element, &continuation.element )),   0., 1.0e-9 );
      _equal( (angleBetweenLineCells<3,Element>( &along_x.element, &along_y.element )),       90., 1.0e-9 );
      _equal( (angleBetweenLineCells<3,Element>( &along_x.element, &reversed.element )),     180., 1.0e-9 );
      // angle between cube edge and space diagonal: acos(1/sqrt(3))
      _equal( (angleBetweenLineCells<3,Element>( &along_x.element, &diagonal.element )), acos( 1./sqrt(3.) ) * 180. / acos(-1.), 1.0e-9 );
    }

    // --- single tetrahedra
    {
      TetFixture  good(     { {0.,0.,0.},    {1.,0.,0.},    {0.,1.,0.},    {0.,0.,1.}    }, 4U );
      TetFixture  twin(     { {0.,0.,0.},    {1.,0.,0.},    {0.,1.,0.},    {0.,0.,1.}    }, 4U ); // same location as 'good'
      TetFixture  shifted(  { {0.1,0.1,0.1}, {1.1,0.1,0.1}, {0.1,1.1,0.1}, {0.1,0.1,1.1} }, 4U ); // overlaps 'good'
      TetFixture  far_away( { {5.,5.,5.},    {6.,5.,5.},    {5.,6.,5.},    {5.,5.,6.}    }, 4U );
      TetFixture  sliver(   { {0.,0.,0.},    {1.,0.,0.},    {0.,1.,0.},    {0.3,0.3,1.0e-3} }, 4U ); // nearly flat
      TetFixture  repeated( { {0.,0.,0.},    {1.,0.,0.},    {0.,1.,0.},    {0.,1.,0.}    }, 4U ); // nodes 2 and 3 collocated
      TriangleFixture  triangle( { {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.} }, 3U, 3U );

      // collocatedNodes: number of nodes of the cell that share their coordinate with another one of its nodes
      _test( collocatedNodes<3>( &good.element ) == 0 );
      _test( collocatedNodes<3>( &repeated.element ) == 1 );

      // normals to the faces
      _test( areUnitNormalsToFacesAreOutwardPointing<3>( &good.element ) );
      _test( areUnitNormalsToFacesAreOutwardPointing<3>( &far_away.element ) );

      // interPenetrating: true if one element contains the barycentre of the other one
      _test( interPenetrating<3>( &good.element, &shifted.element ) );
      _test( interPenetrating<3>( &shifted.element, &good.element ) );
      _test( interPenetrating<3>( &good.element, &good.element ) );
      _test( !interPenetrating<3>( &good.element, &far_away.element ) );
      _test( !interPenetrating<3>( &good.element, &triangle.element ) ); // only defined for volume elements

      // testCellSkewing: returns the cells whose stiffness matrix kernel is (nearly) singular
      {
        const vector<Element<3>*> cells{ &good.element, &sliver.element, &far_away.element };
        const auto skewed = testCellSkewing<3,Element>( cells, 50., false );
        _test( skewed.size() == 1 );
        _test( skewed.size() == 1 && skewed[0] == &sliver.element );
        const vector<Element<3>*> good_cells{ &good.element, &far_away.element };
        _test( (testCellSkewing<3,Element>( good_cells, 50., false )).empty() );
      }

      // findCollocatedCells: cells with identical barycentres; all but the first one of each group are returned
      {
        vector<Element<3>*> collocated;
        const vector<Element<3>*> distinct{ &good.element, &shifted.element, &far_away.element };
        _test( !(findCollocatedCells<3,Element>( distinct.cbegin(), distinct.cend(), collocated )) );
        _test( collocated.empty() );

        const vector<Element<3>*> with_twin{ &good.element, &shifted.element, &twin.element, &far_away.element };
        _test( (findCollocatedCells<3,Element>( with_twin.cbegin(), with_twin.cend(), collocated )) );
        _test( collocated.size() == 1 );
        _test( collocated.size() == 1 && ( collocated[0] == &good.element || collocated[0] == &twin.element ) );

        // the same pointer listed twice is a duplicate entry, but not a collocated cell
        const vector<Element<3>*> repeated_pointer{ &good.element, &far_away.element, &good.element };
        _test( !(findCollocatedCells<3,Element>( repeated_pointer.cbegin(), repeated_pointer.cend(), collocated )) );
        _test( (duplicatesCheck<3,Element>( repeated_pointer.cbegin(), repeated_pointer.cend() )) == 1 );
        _test( (duplicatesCheck<3,Element>( distinct.cbegin(), distinct.cend() )) == 0 );
      }

      // eraseElementPointerFromVector: removes all entries equal to the pointer from a sorted vector
      {
        vector<Element<3>*> cells{ &good.element, &shifted.element, &far_away.element, &sliver.element, &good.element };
        sort( cells.begin(), cells.end() );
        eraseElementPointerFromVector<3>( cells, &good.element );
        _test( cells.size() == 3 );
        _test( find( cells.begin(), cells.end(), &good.element ) == cells.end() );
        // pointer that is not in the vector: nothing happens
        eraseElementPointerFromVector<3>( cells, &twin.element );
        _test( cells.size() == 3 );
        // remaining pointers are intact
        _test( find( cells.begin(), cells.end(), &shifted.element )  != cells.end() );
        _test( find( cells.begin(), cells.end(), &far_away.element ) != cells.end() );
        _test( find( cells.begin(), cells.end(), &sliver.element )   != cells.end() );
      }
    }
 } // end Test_StandAloneElements



/// regular mesh of 3x3x3 unit hexahedra on the node grid {0,1,2,3}^3
void MeshManagementUtilities_Test::Test_HexahedraMesh()
 {
    VSet<3> vset;
    create_Hexahedra_VSet( vset, false );
    Model<3>            model( vset, "MeshManager_Test-variables.txt" );
    MeshManager<3>&     mesh( model.Mesh() );
    Region<3>&          region( model.Region("Model") );
    const Region<3>&    mdomain( region );

    _test( mesh.Nodes() == 64 );
    _test( mesh.Elements() == 27 );
    _test( mdomain.Cells() == 27 );
    _test( mdomain.Nodes() == 64 );

    // --- cell census
    size_t volume_cells{0}, surface_cells{0}, line_cells{0};
    _test( currentCellTypes( mesh, ELEMENT, volume_cells, surface_cells, line_cells ) == 27 + mesh.Faces() );
    _test( volume_cells == 27 );
    _test( surface_cells == 0 );
    _test( line_cells == 0 );
    _test( isoparametricElementMesh( model ) );
    _test( countCornerNodes<3>( mesh.ElementsBegin(), mesh.ElementsEnd() ) == 64 );

    // every cell contains at least one of the 8 interior nodes
    set<size_t> boundary_only_elmts;
    _test( detectElementsWithAllNodesOnBoundary( mesh, boundary_only_elmts ) == 0 );
    _test( boundary_only_elmts.empty() );

    // --- traversal via the neighbor pointers
    Element<3>* const first_elmt{ &(*mesh.ElementsBegin()) };
    {
      set<Element<3>*> filled;
      floodFill( first_elmt, filled );
      _test( filled.size() == 27 );

      set<Element<3>*> patch;
      _test( findContiguousMeshPatch( first_elmt, patch ) == 27 );
      _test( patch.size() == 27 );
      _test( patch == filled );

      map<string,vector<Element<3>*> > patches;
      findContiguousMeshPatches<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd(), patches ); // return value: see open question 5
      _test( patches.size() == 1 );
      _test( patches.size() == 1 && patches.begin()->second.size() == 27 );
    }

    // --- nodes
    {
      set<Node<3>*> orphan_nodes;
      _test( detectOrphanNodes( mesh, orphan_nodes ) == 0 );
      _test( orphan_nodes.empty() );

      set<Node<3>*> cluster;
      _test( findInterconnectedNodeCluster( &(*mesh.NodesBegin()), cluster ) == 64 );
      _test( cluster.size() == 64 );

      // node neighbors along cell edges: 8 corners (3), 24 edge nodes (4), 24 face nodes (5), 8 interior nodes (6)
      vector<set<size_t> > node_neighbors;
      nodeNeighbors( mdomain, node_neighbors );
      _test( node_neighbors.size() == 64 );
      map<size_t,size_t> histogram;
      for ( const auto& nbors : node_neighbors ) histogram[ nbors.size() ]++;
      _test( histogram.size() == 4 );
      _test( histogram[3] == 8 );
      _test( histogram[4] == 24 );
      _test( histogram[5] == 24 );
      _test( histogram[6] == 8 );

      // all node neighbors are one cell edge apart
      _equal( minimumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() ), 1., tol );
      _equal( maximumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() ), 1., tol );
      _equal( averageNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() ), 1., tol );

      const auto corners = boundingBox<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
      _test( corners.first.DistanceTo(  Point<3>{0.,0.,0.} ) < tol );
      _test( corners.second.DistanceTo( Point<3>{3.,3.,3.} ) < tol );
      const auto corners1 = boundingBox1<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
      for ( uint32_t i{0U}; i<3U; ++i ) {
           _equal( corners1.first[i],  0., tol );
           _equal( corners1.second[i], 3., tol );
        }
      // empty range
      const auto no_corners = boundingBox<3>( mdomain.NodesEnd(), mdomain.NodesEnd() );
      _test( no_corners.first.DistanceTo( no_corners.second ) == 0. );

      // findNode
      Node<3>* const target{ nodeAt( mdomain, Point<3>{1.,2.,3.} ) };
      _test( target != nullptr );
      if ( target ) {
           _test( findNode( model, Point<3>{1.,2.,3.}, 1.0e-9 ) == static_cast<long>(target->Idx()) );
           _test( findNode( model, 1., 2., 3., 1.0e-9 ) == target->Idx() );
           // within the tolerance
           _test( findNode( model, Point<3>{1.,2.+1.0e-4,3.}, 1.0e-3 ) == static_cast<long>(target->Idx()) );
        }
      // there is no such node: the Point version returns -1, the coordinate version throws
      _test( findNode( model, Point<3>{1.5,2.,3.}, 1.0e-9 ) == -1 );
      bool exception_thrown{false};
      try { findNode( model, 1.5, 2., 3., 1.0e-9 ); }
      catch ( const csmp::Exception& ) { exception_thrown = true; }
      _test( exception_thrown );

      // distances from the nodes to the barycentres of their parent cells; the last entry is the sum of the distances
      vector<vector<double> > distances_and_weights;
      distancesAndWeights<3>( mdomain.NodesBegin(), mdomain.NodesEnd(), distances_and_weights );
      _test( distances_and_weights.size() == 64 );
      const double half_diagonal{ 0.5 * sqrt(3.) };
      size_t node_index{0};
      for ( auto nit=mdomain.NodesBegin(); nit!=mdomain.NodesEnd(); ++nit, ++node_index ) {
           const auto&  dw{ distances_and_weights[node_index] };
           const size_t parents{ (*nit)->Parents() };
           _test( dw.size() == parents + 1 );
           if ( dw.size() != parents + 1 ) continue;
           for ( size_t i{0}; i<parents; ++i ) _equal( dw[i], half_diagonal, tol );
           _equal( dw[parents], static_cast<double>(parents) * half_diagonal, 1.0e-10 );
        }
      Node<3>* const interior_node{ nodeAt( mdomain, Point<3>{1.,1.,1.} ) };
      Node<3>* const corner_node{   nodeAt( mdomain, Point<3>{3.,3.,3.} ) };
      _test( interior_node && interior_node->Parents() == 8 );
      _test( corner_node && corner_node->Parents() == 1 );
    }

    // --- individual cells
    Element<3>* const corner_elmt{ cellWithBaryCenter( mdomain, Point<3>{0.5,0.5,0.5} ) };
    Element<3>* const next_elmt{   cellWithBaryCenter( mdomain, Point<3>{1.5,0.5,0.5} ) };
    Element<3>* const third_elmt{  cellWithBaryCenter( mdomain, Point<3>{2.5,0.5,0.5} ) };
    Element<3>* const centre_elmt{ cellWithBaryCenter( mdomain, Point<3>{1.5,1.5,1.5} ) };
    _test( corner_elmt && next_elmt && third_elmt && centre_elmt );
    if ( !( corner_elmt && next_elmt && third_elmt && centre_elmt ) ) return;

    // only the central cell has no face on the perimeter of the region
    _test( distance( mdomain.CellsBegin(), mdomain.PerimeterCellsBegin() ) == 1 );
    _test( (*mdomain.CellsBegin()) == centre_elmt );

    {
      // shared face found via neighbor pointers and via nodes
      const auto via_neighbors = findAdjacentFacesFromNeighbors( corner_elmt, next_elmt );
      const auto via_nodes     = findAdjacentFacesFromNodes( corner_elmt, next_elmt );
      _test( via_neighbors == via_nodes );
      _test( via_neighbors.first < corner_elmt->Neighbors() && via_neighbors.second < next_elmt->Neighbors() );
      if ( via_neighbors.first < corner_elmt->Neighbors() && via_neighbors.second < next_elmt->Neighbors() ) {
           _test( corner_elmt->Neighbor( via_neighbors.first ) == next_elmt );
           _test( next_elmt->Neighbor( via_neighbors.second ) == corner_elmt );
           // the shared face lies in the plane x=1
           _equal( corner_elmt->FaceBaryCenter( via_neighbors.first )[0], 1., tol );
           _equal( next_elmt->FaceBaryCenter( via_neighbors.second )[0],  1., tol );
        }
      // swapped arguments: swapped result
      const auto swapped = findAdjacentFacesFromNodes( next_elmt, corner_elmt );
      _test( swapped.first == via_nodes.second && swapped.second == via_nodes.first );
      // cells that do not share a face
      const uint32_t unspecified{ numeric_limits<uint32_t>::max() };
      const auto none_via_neighbors = findAdjacentFacesFromNeighbors( corner_elmt, third_elmt );
      const auto none_via_nodes     = findAdjacentFacesFromNodes( corner_elmt, third_elmt );
      _test( none_via_neighbors.first == unspecified && none_via_neighbors.second == unspecified );
      _test( none_via_nodes.first == unspecified && none_via_nodes.second == unspecified );
    }

    _test( !interPenetrating<3>( corner_elmt, next_elmt ) );
    _test( !interPenetrating<3>( corner_elmt, centre_elmt ) );

    for ( auto c : mdomain.CellVector() ) {
         _test( collocatedNodes<3>( c ) == 0 );
         _test( areUnitNormalsToFacesAreOutwardPointing<3>( c ) );
      }
    _test( (testCellSkewing<3,Element>( mdomain.CellVector(), 50., false )).empty() );

    // --- point location
    _test( pointInVolumeElement( region, Point<3>{1.4,1.6,1.5} ) == centre_elmt );
    _test( pointInVolumeElement( region, Point<3>{0.2,0.3,0.4} ) == corner_elmt );
    _test( pointInVolumeElement( region, Point<3>{5.,5.,5.} ) == nullptr );
    _test( isContainedIn<3>( mdomain, array<double,3>{1.4,1.6,1.5} ) == centre_elmt );
    _test( isContainedIn<3>( mdomain, array<double,3>{2.7,0.1,0.9} ) == third_elmt );
    _test( isContainedIn<3>( mdomain, array<double,3>{5.,5.,5.} ) == nullptr );

    // --- checks of the mesh as a whole
    _test( (integrityCheck<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() )) );
    _test( !(integrityCheck<3,Element>( mesh.ElementsEnd(), mesh.ElementsEnd() )) ); // empty range
    _test( (setNeighborsWithInvalidFE_PointersTo_nullptr<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() )) == 0 );
    _test( (duplicatesCheck<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd() )) == 0 );
    _test( (detectDuplicateCells<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), false )) == 0 );
    {
      vector<Element<3>*> collocated;
      _test( !(findCollocatedCells<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), collocated )) );
      _test( collocated.empty() );
    }
    {
      // two cells listed twice, and a nullptr entry which must be skipped
      vector<Element<3>*> cells( mdomain.CellVector() );
      cells.push_back( corner_elmt );
      cells.push_back( nullptr );
      cells.push_back( centre_elmt );
      _test( (detectDuplicateCells<3,Element>( cells.cbegin(), cells.cend(), false )) == 2 );
    }
    {
      // backup of the neighbor pointers: entry i of a cell is its neighbor i, nullptr on the model boundary
      vector<vector<Element<3>*> > nbor_pointers;
      backupNeighborConnectivity<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), nbor_pointers );
      _test( nbor_pointers.size() == 27 );
      size_t valid_neighbors{0}, mismatches{0}, cell_index{0};
      for ( auto cit=mdomain.CellsBegin(); cit!=mdomain.CellsEnd() && cell_index<nbor_pointers.size(); ++cit, ++cell_index ) {
           const auto& nbors{ nbor_pointers[cell_index] };
           if ( nbors.size() != (*cit)->Neighbors() ) { ++mismatches; continue; }
           for ( uint32_t i{0U}; i<(*cit)->Neighbors(); ++i ) {
                if ( nbors[i] != (*cit)->Neighbor(i) ) ++mismatches;
                if ( nbors[i] != nullptr ) ++valid_neighbors;
             }
        }
      _test( mismatches == 0 );
      // 27 cells with 6 faces each, of which 6x9 lie on the surface of the cube
      _test( valid_neighbors == 27*6 - 54 );
    }
 } // end Test_HexahedraMesh



/// cube [-1,1]^3 subdivided into 6 tetrahedra which all share the space diagonal from node 3 to node 5
void MeshManagementUtilities_Test::Test_TetrahedraMesh()
 {
    VSet<3> vset;
    create_Tetra_VSet( vset );
    Model<3>          model( vset, "MeshManager_Test-variables.txt" );
    MeshManager<3>&   mesh( model.Mesh() );
    const Region<3>&  mdomain( model.Region("Model") );

    _test( mesh.Nodes() == 8 );
    _test( mesh.Elements() == 6 );

    size_t volume_cells{0}, surface_cells{0}, line_cells{0};
    currentCellTypes( mesh, ELEMENT, volume_cells, surface_cells, line_cells );
    _test( volume_cells == 6 && surface_cells == 0 && line_cells == 0 );
    _test( isoparametricElementMesh( model ) );
    _test( countCornerNodes<3>( mesh.ElementsBegin(), mesh.ElementsEnd() ) == 8 );

    // there is no interior node: all cells have all of their nodes on the boundary
    set<size_t> boundary_only_elmts;
    _test( detectElementsWithAllNodesOnBoundary( mesh, boundary_only_elmts ) == 6 );
    _test( boundary_only_elmts.size() == 6 );

    set<Element<3>*> filled;
    floodFill( &(*mesh.ElementsBegin()), filled );
    _test( filled.size() == 6 );

    set<Node<3>*> cluster;
    _test( findInterconnectedNodeCluster( &(*mesh.NodesBegin()), cluster ) == 8 );

    // shortest edges: cube edges (2), longest edge: the space diagonal (2 sqrt(3))
    _equal( minimumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() ), 2., tol );
    _equal( maximumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() ), 2.*sqrt(3.), tol );
    const double average_spacing{ averageNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() ) };
    _test( average_spacing > 2. && average_spacing < 2.*sqrt(3.) );

    const auto corners = boundingBox<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
    _test( corners.first.DistanceTo(  Point<3>{-1.,-1.,-1.} ) < tol );
    _test( corners.second.DistanceTo( Point<3>{ 1., 1., 1.} ) < tol );

    for ( auto c : mdomain.CellVector() ) {
         _test( collocatedNodes<3>( c ) == 0 );
         _test( areUnitNormalsToFacesAreOutwardPointing<3>( c ) );
      }
    for ( size_t i{0}; i<mdomain.Cells(); ++i )
      for ( size_t j{i+1}; j<mdomain.Cells(); ++j )
        _test( !interPenetrating<3>( mdomain.E(i), mdomain.E(j) ) );

    _test( (testCellSkewing<3,Element>( mdomain.CellVector(), 50., false )).empty() );
    _test( (integrityCheck<3,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() )) );
    _test( (duplicatesCheck<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd() )) == 0 );
    _test( (detectDuplicateCells<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), false )) == 0 );
    _test( connectivityCheck<3>( mdomain.CellsBegin(), mdomain.CellsEnd() ) == 0 );

    // every point of the cube lies in one of the tetrahedra
    _test( isContainedIn<3>( mdomain, array<double,3>{ 0.3,-0.2, 0.5} ) != nullptr );
    _test( isContainedIn<3>( mdomain, array<double,3>{-0.7, 0.6,-0.1} ) != nullptr );
    _test( isContainedIn<3>( mdomain, array<double,3>{ 1.5, 0., 0.} ) == nullptr );
 } // end Test_TetrahedraMesh



/// two quadrilaterals (nodes with x in [1,4]) and four triangles (x in [5,7]) which do not share any node
void MeshManagementUtilities_Test::Test_DisconnectedMesh2D()
 {
    VSet<2> vset;
    create_Disconnected2D_VSet( vset );
    Model<2>          model( vset, "CSMP-variables.txt" );
    MeshManager<2>&   mesh( model.Mesh() );
    const Region<2>&  mdomain( model.Region("Model") );

    _test( mesh.Nodes() == 12 );
    _test( mesh.Elements() == 6 );

    size_t volume_cells{9}, surface_cells{9}, line_cells{9};
    _test( currentCellTypes( mesh, ELEMENT, volume_cells, surface_cells, line_cells ) == 6 + mesh.Faces() );
    _test( volume_cells == 0 && surface_cells == 6 && line_cells == 0 );
    _test( countCornerNodes<2>( mesh.ElementsBegin(), mesh.ElementsEnd() ) == 12 );

    Element<2>* const quad{     cellWithBaryCenter( mdomain, Point<2>{2.5,4.} ) };
    Element<2>* const triangle{ cellWithBaryCenter( mdomain, Point<2>{17./3.,11./3.} ) };
    _test( quad && triangle );
    if ( !( quad && triangle ) ) return;

    // the flood fill stays within the patch that it was started in
    set<Element<2>*> quads, triangles;
    floodFill( quad, quads );
    floodFill( triangle, triangles );
    _test( quads.size() == 2 );
    _test( triangles.size() == 4 );
    for ( auto c : quads )     _test( c->Nodes() == 4 );
    for ( auto c : triangles ) _test( c->Nodes() == 3 );
    // ... whichever cell of the patch is used as the seed
    for ( auto seed : triangles ) {
         set<Element<2>*> patch;
         _test( findContiguousMeshPatch( seed, patch ) == 4 );
         _test( patch == triangles );
      }

    map<string,vector<Element<2>*> > patches;
    findContiguousMeshPatches<2,Element>( mesh.ElementsBegin(), mesh.ElementsEnd(), patches ); // return value: see open question 5
    _test( patches.size() == 2 );
    size_t cells_in_patches{0};
    for ( const auto& p : patches ) {
         _test( p.second.size() == 2 || p.second.size() == 4 );
         cells_in_patches += p.second.size();
      }
    _test( cells_in_patches == 6 );

    // node clusters: 6 nodes per patch
    Node<2>* const left_node{  nodeAt( mdomain, Point<2>{1.,5.} ) };
    Node<2>* const right_node{ nodeAt( mdomain, Point<2>{7.,1.} ) };
    _test( left_node && right_node );
    if ( left_node && right_node ) {
         set<Node<2>*> left_cluster, right_cluster;
         _test( findInterconnectedNodeCluster( left_node, left_cluster ) == 6 );
         _test( findInterconnectedNodeCluster( right_node, right_cluster ) == 6 );
         _test( left_cluster.count( right_node ) == 0 );
         _test( right_cluster.count( left_node ) == 0 );
      }
    set<Node<2>*> orphan_nodes;
    _test( detectOrphanNodes( mesh, orphan_nodes ) == 0 );

    // findNode in 2D
    Node<2>* const target{ nodeAt( mdomain, Point<2>{4.,3.} ) };
    _test( target != nullptr );
    if ( target ) {
         _test( findNode( model, Point<2>{4.,3.}, 1.0e-9 ) == static_cast<long>(target->Idx()) );
         _test( findNode( model, 4., 3., 1.0e-9 ) == target->Idx() );
      }
    _test( findNode( model, Point<2>{4.5,3.}, 1.0e-9 ) == -1 ); // in the gap between the patches

    const auto corners = boundingBox<2>( mdomain.NodesBegin(), mdomain.NodesEnd() );
    _test( corners.first.DistanceTo(  Point<2>{1.,1.} ) < tol );
    _test( corners.second.DistanceTo( Point<2>{7.,5.} ) < tol );

    // point location
    _test( isContainedIn<2>( mdomain, array<double,2>{2.,4.} ) == quad );
    _test( isContainedIn<2>( mdomain, array<double,2>{5.6,3.7} ) == triangle );
    _test( isContainedIn<2>( mdomain, array<double,2>{4.5,4.} ) == nullptr ); // in the gap between the patches

    for ( auto c : mdomain.CellVector() ) _test( areUnitNormalsToFacesAreOutwardPointing<2>( c ) );
    _test( !interPenetrating<2>( quad, triangle ) );
    _test( (testCellSkewing<2,Element>( mdomain.CellVector(), 50., false )).empty() );
    _test( (integrityCheck<2,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() )) );
    _test( (detectDuplicateCells<2,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), false )) == 0 );
 } // end Test_DisconnectedMesh2D



/// 25 surface elements (triangles and quadrilaterals) with line elements along some of their edges
void MeshManagementUtilities_Test::Test_MeshWithLineElements2D()
 {
    VSet<2>       vset;
    ModelTopology topo = create_MeshPatchWithLineElements_VSet( vset );
    const bool    treat_domains_as_regions{true}; // as in ANSYS_Model2D_Test
    Model<2>          model( topo, vset, "CSMP-variables.txt", treat_domains_as_regions );
    MeshManager<2>&   mesh( model.Mesh() );
    const Region<2>&  mdomain( model.Region("Model") );

    _test( mesh.Nodes() == 22 );

    // elements: the line elements that were not turned into Face objects remain in the element container
    size_t volume_cells{9}, surface_cells{9}, line_cells{9};
    const size_t n_cells{ currentCellTypes( mesh, ELEMENT, volume_cells, surface_cells, line_cells ) };
    _test( n_cells == mesh.Elements() + mesh.Faces() + mesh.Interfaces() );
    _test( n_cells == 45 ); // as many cells as the VSet has
    _test( volume_cells == 0 );
    _test( surface_cells == 25 );
    _test( surface_cells + line_cells == mesh.Elements() );
    // faces of a 2D model are lines
    const size_t line_elmts{ line_cells };
    currentCellTypes( mesh, FACE, volume_cells, surface_cells, line_cells );
    _test( volume_cells == 0 && surface_cells == 0 );
    _test( line_cells == mesh.Faces() );
    _test( line_elmts + line_cells == 20 ); // number of line elements in the VSet

    _test( isoparametricElementMesh( model ) );

    // angles between line elements of the model [degrees]
    vector<Element<2>*> lines;
    for ( auto c : mdomain.CellVector() ) if ( c->IsLine() ) lines.push_back( c );
    _test( lines.size() == line_elmts );
    for ( auto line1 : lines )
      for ( auto line2 : lines ) {
           const double angle{ angleBetweenLineCells<2,Element>( line1, line2 ) };
           _test( angle >= 0. && angle <= 180. );
           if ( line1 == line2 ) _equal( angle, 0., 1.0e-6 );
           _equal( (angleBetweenLineCells<2,Element>( line2, line1 )), angle, 1.0e-9 );
           // chained line elements: the angle is that between the two direction vectors, which can be computed by hand
           const Point<2> v1( line1->N(1)->Coordinate() - line1->N(0)->Coordinate() );
           const Point<2> v2( line2->N(1)->Coordinate() - line2->N(0)->Coordinate() );
           const double cosine{ ( v1[0]*v2[0] + v1[1]*v2[1] ) / ( sqrt( v1[0]*v1[0] + v1[1]*v1[1] ) * sqrt( v2[0]*v2[0] + v2[1]*v2[1] ) ) };
           _equal( angle, acos( max( -1., min( 1., cosine ) ) ) * 180. / acos(-1.), 1.0e-4 ); // acos is ill-conditioned near 0 and 180
        }

    // the surface elements form a single patch
    Element<2>* surface_elmt{ nullptr };
    for ( auto c : mdomain.CellVector() ) if ( c->IsSurface() ) { surface_elmt = c; break; }
    _test( surface_elmt != nullptr );
    if ( surface_elmt ) {
         set<Element<2>*> filled;
         floodFill( surface_elmt, filled );
         size_t surfaces_in_patch{0};
         for ( auto c : filled ) if ( c->IsSurface() ) ++surfaces_in_patch;
         _test( surfaces_in_patch == 25 );
      }

    set<Node<2>*> orphan_nodes;
    _test( detectOrphanNodes( mesh, orphan_nodes ) == 0 );
    _test( (integrityCheck<2,Element>( mesh.ElementsBegin(), mesh.ElementsEnd() )) );
    _test( (detectDuplicateCells<2,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), false )) == 0 );

    // print functions: must run through; printLineElementRegion returns the number of line elements in the region
    _test( printLineElementRegion( model, "FRAC1", false ) == 3 );
    printNodeCoordinates<2>( mdomain.NodesBegin(), mdomain.NodesBegin()+3 );
    printNodeAttibutes<2>( mdomain.NodesBegin(), mdomain.NodesBegin()+3 );
    printNodes<2,Element>( *mdomain.E(0) );
 } // end Test_MeshWithLineElements2D



/// unrefined 'FracBox': tetrahedra, the triangles of a planar fracture and line elements in a single region
void MeshManagementUtilities_Test::Test_FracturedBox()
 {
    VSet<3>       vset;
    ModelTopology topo = create_FracBox( vset );
    // the model name is the prefix of the optional '-regions.txt' file, which would remove cells from the model;
    // with a name of its own, this test always sees the complete mesh
    topo.ModelName("FracBox_complete");
    Model<3>          model( topo, vset, "CSMP-1phase-variables.txt", true );
    MeshManager<3>&   mesh( model.Mesh() );
    const Region<3>&  mdomain( model.Region("Model") );

    size_t volume_cells{0}, surface_cells{0}, line_cells{0};
    const size_t n_cells{ currentCellTypes( mesh, ELEMENT, volume_cells, surface_cells, line_cells ) };
    _test( n_cells == 1872 ); // cells in the VSet
    _test( n_cells == mesh.Elements() + mesh.Faces() );
    _test( volume_cells == 1363 );
    _test( surface_cells == 38 );
    _test( line_cells == 89 );
    _test( volume_cells + surface_cells + line_cells == mesh.Elements() );
    // the triangles on the 6 box boundaries have become Face objects
    size_t surface_faces{0}, line_faces{0}, no_volumes{0};
    currentCellTypes( mesh, FACE, no_volumes, surface_faces, line_faces );
    _test( surface_faces == 382 && line_faces == 0 && no_volumes == 0 );
    _test( surface_faces == mesh.Faces() );

    vector<Element<3>*> volumes, surfaces;
    for ( auto c : mdomain.CellVector() ) {
         if ( c->IsVolume() ) volumes.push_back( c );
         else if ( c->IsSurface() ) surfaces.push_back( c );
      }
    _test( volumes.size() == volume_cells );
    _test( surfaces.size() == surface_cells );
    if ( volumes.empty() || surfaces.empty() ) return;

    // cells of different dimension are not neighbors: the flood fill finds the volumes, or the fracture
    {
      set<Element<3>*> filled;
      floodFill( volumes.front(), filled );
      _test( filled.size() == volume_cells );
      filled.clear();
      floodFill( surfaces.front(), filled );
      _test( filled.size() == surface_cells );
      for ( auto c : filled ) _test( c->IsSurface() );
    }

    // the fracture is planar and its triangles are consistently oriented
    size_t pairs{0};
    for ( auto c : surfaces )
      for ( uint32_t i{0U}; i<c->Neighbors(); ++i )
        if ( c->Neighbor(i) && c->Neighbor(i)->IsSurface() ) {
             _equal( angleBetweenSurfaceCells<Element>( c, c->Neighbor(i) ), 0., 1.0e-4 );
             ++pairs;
          }
    _test( pairs > 0 );

    _test( connectivityCheck<3>( mdomain.CellsBegin(), mdomain.CellsEnd() ) == 0 );
    _test( (detectDuplicateCells<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), false )) == 0 );
    for ( auto c : volumes ) _test( areUnitNormalsToFacesAreOutwardPointing<3>( c ) );
 } // end Test_FracturedBox
 
 

void MeshManagementUtilities_Test::Test_transformMeshIntoRefinedLinearAndQuadraticMeshes()
  {
      VSet<3>       mesh;
      ModelTopology topo = create_FracBox( mesh );
      
      VSet<3>       refined_lin_mesh, quadratic_mesh;
      ModelTopology refined_linear_topo, quadratic_topo;
      
      // removing any variables with a placement other than ELEMENT or FACE since they cannot handled by the tested method yet
      set<string> props_to_delete;
      for ( auto pit=mesh.PropertyValuesBegin(); pit!=mesh.PropertyValuesEnd(); ++pit )
        if ( (*pit).second.Placement() != ELEMENT && (*pit).second.Placement() != FACE )
          props_to_delete.insert( (*pit).first );
      for ( const auto& pit : props_to_delete )
        mesh.RemoveData( pit.c_str() );
  
      transformMeshIntoRefinedLinearAndQuadraticMeshes( mesh, topo, refined_lin_mesh, refined_linear_topo,
                                                        quadratic_mesh, quadratic_topo );

      // testing the transformed meshes
      {
        const bool treat_domains_as_regions_and_use_regions_file_if_any{true}; // TODO: change wording
        Model<3> lin_model( refined_linear_topo, refined_lin_mesh, "CSMP-1phase-variables.txt", treat_domains_as_regions_and_use_regions_file_if_any );
        const Region<3>& mdomain = lin_model.Region("Model");

        _test( distance(lin_model.BoundariesBegin(), lin_model.BoundariesEnd()) == 6 );

        // Testing using MeshManagementUtilities computation
        size_t volume_cells{0}, surface_cells{0}, line_cells{0};
        size_t n_total_cells = currentCellTypes( lin_model.Mesh(), ELEMENT, volume_cells, surface_cells, line_cells );
        size_t n_model_elmts = lin_model.Mesh().Elements();
        _test( n_total_cells == n_model_elmts + lin_model.Mesh().Faces() );
        _test( volume_cells == 10904 ); // after conversion into Face objects
        _test( volume_cells < refined_lin_mesh.Elements() ); // after conversion into Face objects
        _test( line_cells == 0 ); // all removed since not mentioned in regions file
        
        _test( isoparametricElementMesh(lin_model) == true );

        set<size_t> detached_elmts;
        _test( detectElementsWithAllNodesOnBoundary( lin_model.Mesh(), detached_elmts ) == 0 );

        const auto duplicates = detectDuplicateCells<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), verbose_ );
        _test( duplicates == 0 );
    
        // Compute parent element barycentre-to-node distances for range of nodes;
        // returns them into vector [e1,e2...e_n,e_sum] with a length of parent elements+1
        vector<vector<double> > distances_and_weight;
        distancesAndWeights<3>( mdomain.NodesBegin(), mdomain.NodesEnd(), distances_and_weight );
        _test( distances_and_weight[25][2] > 0. );

        // trying corners:
        //for supplied edge nodes, find their volumetric parent elements; if find segment ids is on, their local numbers are assigned to Idx of the parent elements
        for ( auto& it : mdomain.CellVector() )
           if ( isCorner( it->AtBoundary(1) /*face=1*/ ) ) {
                 const bool find_segm_ids{ true };
                 vector<Node<3>*>                   edge_nodes;
                 map<Element<3>*,vector<Node<3>*> > segm_parents;
                 size_t n_ele = parentElementsSharingMultipleEdgeNodes( edge_nodes, segm_parents, find_segm_ids );
                 _test( n_ele >= 2 );
                 _test( edge_nodes.size() == 2 );
                 _test( segm_parents.size() == 2 );
             }

        // finds node by point coordinate; returns -1 if not found; @attention tolerance needs to account for single-precision of CAD tools
        auto test_node = mdomain.N(25);
        long idx = findNode( lin_model, test_node->Coordinate(), 1.0e-7, verbose_ );
        _test( idx >= 0 ); // node was found
        _test( test_node->Idx() == static_cast<size_t>(idx) );
        size_t id = findNode( lin_model, test_node->x(), test_node->y(), test_node->z(), 1.0e-7 );
        _test( test_node->Idx() == id );

         set<Node<3>*> nodes;
        _test( detectOrphanNodes( lin_model.Mesh(), nodes ) == 0 );
        _test( findInterconnectedNodeCluster( test_node, nodes ) == lin_model.Mesh().Nodes() ); // all nodes should be connected

        // Relying on the node-to-node connectivity, uses a recursive depth first traversal to discover
        // the interconnected nodes in the supplied range; set of node pairs defining the cel segments must be input
        set<pair<const Node<3>*,const Node<3>*>> validEdges;
        size_t n_inds = findInterconnectedNodes<3>( lin_model.Mesh().Nodes(), validEdges );
        _test( n_inds == lin_model.Mesh().Nodes() );

        double min_dn = minimumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        double max_dn = maximumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        double avg_dn = averageNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        _test( min_dn > 0.4 );
        _test( min_dn < max_dn );
        _test( max_dn < 50. );
        _test( avg_dn > min_dn );

        // must be a volumetric element 
        Point<3> bctr = (mdomain.E(1000)->IsEquidimensional()) ? mdomain.E(1000)->BaryCenter() : mdomain.E(2000)->BaryCenter();
        _test( isContainedIn<3>( mdomain, array<double,3>{ bctr[0], bctr[1], bctr[2] } ) == mdomain.E(1000) );

        size_t duplicated_elmts = duplicatesCheck<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd() );
        _test( duplicated_elmts == 0 );
                          
        size_t connectivity_errors = connectivityCheck<3>( mdomain.CellsBegin(), mdomain.CellsEnd() );
        _test( connectivity_errors == 0 );

        pair<Point<3>,Point<3>> cnrs = boundingBox<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        _test( cnrs.first < cnrs.second );
      }
  }
  
} // end csmp
