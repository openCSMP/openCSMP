// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only
//
//  VData_Test.cpp
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 26/1/2022.
//

#include "VData_Test.h"
#include "VData.h"
#include "VSet.h"
#include "VSetConverter.h"
#include "Box.h"
#include "ScalarVariable.h"
#include "ModelTopology.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "ANSYS_Interface.h"
#include "VTK_Interface.h"
#include "vsetMakers.h"
#include "Exception.h"
#include <random>

using namespace std;

namespace csmp {

  
/** TESTING DATA MEMBERS

    std::vector<double>             px, py, pz; (1)
    std::deque<std::vector<uint32_t> >  plist;      (2)
    std::deque<std::vector<int32_t> >   pfverts;    (3)
    std::map<size_t,int32_t>            bflags;     (4)
    std::vector<int32_t>                pelmt;      (5)

    (1) ... size of node# goes from 0..(n-1): node coordinates
    (2) ... plist: deque size of element#. each entry contains
            a vector with size of node count for that element.
            the vector contains the node ids wrt px,py,pz index
            that make up the given element in the csmp specific order.
    (3) ... deque size of element#, with a vector that each holds all
            neighbour elmt ids for given elememt wrt plist deque index.
    (5) ... type of elmts in int, indicating number of nodes, vertices


  */
void VData_Test::run()
{
  // TESTING BASICS
  // ==============
  
  // .)CONSTRUCTORS
  VData defaultVData;
  std::deque<uint32_t> npes;
  std::deque<uint32_t> epes;
  const size_t nodeCount1( 12 );
  // nodes per element
  npes.push_back( 8 );
  npes.push_back( 8 );
  npes.push_back( 5 );
  // number of neighbours
  epes.push_back( 6 );
  epes.push_back( 6 );
  epes.push_back( 5 );
  VData dequeVData( npes, epes, nodeCount1 );
  VData countVData( 8, 6, nodeCount1, 2 );
  VData copyDequeVData( dequeVData );
  VData copyCountVData = countVData;
  _test( countVData.HybridElementTypeMesh() == false );
  _test( dequeVData.HybridElementTypeMesh() == true );
  _test( countVData.Elements() == 2 );
  _test( countVData.Vertices() == nodeCount1 );
  _test( countVData.PfvertsSize( 1 ) == 6 );
  _test( countVData.PlistSize( 1 ) == 8 );
  _test( dequeVData.Elements() == 3 );
  _test( dequeVData.Vertices() == nodeCount1 );
  _test( dequeVData.PfvertsSize( 1 ) == 6 );
  _test( dequeVData.PfvertsSize( 2 ) == 5 );
  _test( dequeVData.PlistSize( 1 ) == 8 );
  _test( dequeVData.PlistSize( 2 ) == 5 );


  // .)COMPARISON OP
  _test( copyDequeVData == dequeVData );
  _test( copyCountVData == countVData );
  _test( !( copyDequeVData == countVData ) );
  _test( !( copyCountVData == dequeVData ) );

  // .)ELEMENT TYPE OPS
  std::vector<int8_t> elmtType;
  elmtType.push_back( 8 );
  elmtType.push_back( 8 );
  elmtType.push_back( 8 );
  elmtType.push_back( 6 );
  elmtType.push_back( 6 );
  VData dummyVData( 2, 2, 8, 5 );
  dummyVData.ElementTypes( elmtType );
  _test( dummyVData.ElementType( 0 ) == 8 );
  _test( dummyVData.ElementType( 2 ) == 8 );
  _test( dummyVData.ElementType( 4 ) == 6 );
  auto elmtStart = elmtType.begin();
  auto elmtEnd = elmtType.end();
  dummyVData.AddElementTypes( elmtStart, elmtEnd );
  _test( dummyVData.ElementType( 0 ) == 8 );
  _test( dummyVData.ElementType( 2 ) == 8 );
  _test( dummyVData.ElementType( 4 ) == 6 );
  dummyVData.ElementType( 0, 2 );
  dummyVData.ElementType( 2, 2 );
  dummyVData.ElementType( 4, 2 );
  _test( dummyVData.ElementType( 0 ) == 2 );
  _test( dummyVData.ElementType( 2 ) == 2 );
  _test( dummyVData.ElementType( 4 ) == 2 );
  dummyVData.SingleElementType( 5 );
  _test( dummyVData.ElementType( 0 ) == 5 );

  // .)RESET/RESIZE
  countVData.Resize( 4, 4, 8, 6, 3 );
  _test( countVData.Elements() == 3 );
  _test( countVData.Vertices() == 6 );
  std::deque<int8_t> etypes;
  etypes.push_back( 3 );
  etypes.push_back( 3 );
  etypes.push_back( 3 );
  std::deque<uint32_t> npes2;
  npes2.push_back( 5 );
  npes2.push_back( 5 );
  npes2.push_back( 5 );
  std::deque<uint32_t> epes2;
  epes2.push_back( 4 );
  epes2.push_back( 4 );
  epes2.push_back( 4 );
  countVData.Resize( etypes, npes2, epes2, 6, 0, 0 );
  _test( countVData.Elements() == 3 );
  _test( countVData.ElementType( 1 ) == 3 );
  _test( countVData.Vertices() == 6 );
  _test( countVData.PfvertsSize( 1 ) == 4 );
  _test( countVData.PlistSize( 1 ) == 5 );
  dequeVData.ResizeNodes( 10 );
  _test( dequeVData.Vertices() == 10 );
  std::deque<uint32_t> plistSize;
  plistSize.push_back( 9 );
  plistSize.push_back( 9 );
  plistSize.push_back( 9 );
  dequeVData.ResizePlist( plistSize );
  _test( dequeVData.PlistSize( 1 ) == 9 );
  std::deque<uint32_t> pfvertsDeque;
  pfvertsDeque.push_back( 10 );
  pfvertsDeque.push_back( 10 );
  pfvertsDeque.push_back( 10 );
  dequeVData.ResizePfverts( pfvertsDeque );

  // .)PLIST/NUMBERING
  VData vdata1( 3, 2, 3, 1 );
  vdata1.Plist( 0, 0, 1 );
  vdata1.Plist( 0, 1, 2 );
  vdata1.Plist( 0, 2, 3 );
  _test( vdata1.Plist( 0, 0 ) == 1 );
  _test( vdata1.Plist( 0, 1 ) == 2 );
  _test( vdata1.Plist( 0, 2 ) == 3 );
  vdata1.EstablishZeroBasedNumbering();
  _test( vdata1.Plist( 0, 0 ) == 0 );
  _test( vdata1.Plist( 0, 1 ) == 1 );
  _test( vdata1.Plist( 0, 2 ) == 2 );

  // .)PFVERTS
  vdata1.Pfvert( 0, 0, 2 );
  vdata1.Pfvert( 0, 1, 1 );
  _test( vdata1.Pfvert( 0, 0 ) == 2 );
  _test( vdata1.Pfvert( 0, 1 ) == 1 );

  // .)COORDINATES
  vdata1.Px( 0, 10. );
  vdata1.Py( 0, 10. );
  vdata1.Pz( 0, 10. );
  vdata1.Px( 1, 20. );
  vdata1.Py( 1, 20. );
  vdata1.Pz( 1, 20. );
  vdata1.Px( 2, 30. );
  vdata1.Py( 2, 30. );
  vdata1.Pz( 2, 30. );
  _test( vdata1.Px( 0 ) == 10. );
  _test( vdata1.Py( 1 ) == 20. );
  _test( vdata1.Pz( 2 ) == 30. );
  double min, max;
  vdata1.CoordinateRange( 'x', min, max );
  _equal( min, 10., 1.E-10 );
  _equal( max, 30., 1.E-10 );
  vdata1.CoordinateRange( 'y', min, max );
  _equal( min, 10., 1.E-10 );
  _equal( max, 30., 1.E-10 );
  vdata1.CoordinateRange( 'z', min, max );
  _equal( min, 10., 1.E-10 );
  _equal( max, 30., 1.E-10 );
  vdata1.ScaleCoordinateToRange( 'x', 15., 25. );
  vdata1.ScaleCoordinateToRange( 'y', 15., 25. );
  vdata1.ScaleCoordinateToRange( 'z', 15., 25. );
  vdata1.CoordinateRange( 'x', min, max );
  _equal( min, 15., 1.E-5 );
  _equal( max, 25., 1.E-5 );
  vdata1.CoordinateRange( 'y', min, max );
  _equal( min, 15., 1.E-5 );
  _equal( max, 25., 1.E-5 );
  vdata1.CoordinateRange( 'z', min, max );
  _equal( min, 15., 1.E-5 );
  _equal( max, 25., 1.E-5 );

  // .)CHECKS
  VData vdata2( 3, 2, 3, 1 );
  vdata2.Plist( 0, 0, 1 );
  vdata2.Plist( 0, 1, 2 );
  vdata2.Plist( 0, 2, 3 );
  _test( vdata2.CheckFix() );
  vdata2.Plist( 0, 0 , 9999 );
  vdata2.BFlag( 1, 100 );
  _test( !( vdata2.CheckFix() ) );

  // .)ORPHANS
  VData orphanVData( 8, 6, 13, 2 );
  orphanVData.Px( 0, 0. ); orphanVData.Py( 0, 0. ); orphanVData.Pz( 0, 0. );
  orphanVData.Px( 1, 1. ); orphanVData.Py( 1, 1. ); orphanVData.Pz( 1, 1. );
  orphanVData.Px( 2, 2. ); orphanVData.Py( 2, 2. ); orphanVData.Pz( 2, 2. );
  orphanVData.Px( 3, 3. ); orphanVData.Py( 3, 3. ); orphanVData.Pz( 3, 3. );
  orphanVData.Px( 4, 4. ); orphanVData.Py( 4, 4. ); orphanVData.Pz( 4, 4. );
  orphanVData.Px( 5, 5. ); orphanVData.Py( 5, 5. ); orphanVData.Pz( 5, 5. );
  orphanVData.Px( 6, 6. ); orphanVData.Py( 6, 6. ); orphanVData.Pz( 6, 6. );
  orphanVData.Px( 7, 7. ); orphanVData.Py( 7, 7. ); orphanVData.Pz( 7, 7. );
  orphanVData.Px( 8, 8. ); orphanVData.Py( 8, 8. ); orphanVData.Pz( 8, 8. );
  orphanVData.Px( 9, 9. ); orphanVData.Py( 9, 9. ); orphanVData.Pz( 9, 9. );
  orphanVData.Px( 10, 10. ); orphanVData.Py( 10, 10. ); orphanVData.Pz( 10, 10. );
  orphanVData.Px( 11, 11. ); orphanVData.Py( 11, 11. ); orphanVData.Pz( 11, 11. );
  orphanVData.Px( 12, 12. ); orphanVData.Py( 12, 12. ); orphanVData.Pz( 12, 12. );
  // ..) first hexahedron
  orphanVData.Plist( 0, 0, 0 );
  orphanVData.Plist( 0, 1, 1 );
  orphanVData.Plist( 0, 2, 2 );
  orphanVData.Plist( 0, 3, 3 );
  orphanVData.Plist( 0, 4, 4 );
  orphanVData.Plist( 0, 5, 6 );
  orphanVData.Plist( 0, 6, 7 );
  orphanVData.Plist( 0, 7, 5 );
  // ..) second hexahedron
  orphanVData.Plist( 1, 0, 1 );
  orphanVData.Plist( 1, 1, 8 );
  orphanVData.Plist( 1, 2, 9 );
  orphanVData.Plist( 1, 3, 2 );
  orphanVData.Plist( 1, 4, 6 );
  orphanVData.Plist( 1, 5, 10 );
  orphanVData.Plist( 1, 6, 11);
  orphanVData.Plist( 1, 7, 7 );
  _test( orphanVData.Vertices() == 13 );
  orphanVData.DetectAndEliminateOrphanNodes( true );
  _test( orphanVData.Vertices() == 12 );

  // .)REDUCETO
  VData reducedVData = orphanVData;
  std::map<size_t,size_t> reducedMap;
  std::map<size_t,size_t> reducedNodeMap;
  reducedMap.insert( std::make_pair( 1,0 ) );
  reducedVData.ReduceTo( reducedMap, reducedNodeMap );
  _test( reducedVData.Elements() == 1 );
  _test( reducedVData.Vertices() == 8 );
  for( unsigned int i = 0; i < reducedVData.Vertices(); ++i )
  {
    _test( reducedVData.Px( i ) != 0. ) ;  _test( reducedVData.Py( i ) != 0. ) ;  _test( reducedVData.Pz( i ) != 0. ) ;
    _test( reducedVData.Px( i ) != 3. ) ;  _test( reducedVData.Py( i ) != 3. ) ;  _test( reducedVData.Pz( i ) != 3. ) ;
    _test( reducedVData.Px( i ) != 4. ) ;  _test( reducedVData.Py( i ) != 4. ) ;  _test( reducedVData.Pz( i ) != 4. ) ;
    _test( reducedVData.Px( i ) != 5. ) ;  _test( reducedVData.Py( i ) != 5. ) ;  _test( reducedVData.Pz( i ) != 5. ) ;
  }

  // .)ITERATORS
  _test( reducedVData.PlistBegin() == reducedVData.plist.begin() );
  _test( reducedVData.PlistEnd() == reducedVData.plist.end() );
  _test( reducedVData.PfvertsBegin() == reducedVData.pfverts.begin() );
  _test( reducedVData.PfvertsEnd() == reducedVData.pfverts.end() );
  _test( reducedVData.PelmtBegin() == reducedVData.pelmt.begin() );
  _test( reducedVData.PelmtEnd() == reducedVData.pelmt.end() );
  _test( reducedVData.BFlagsBegin() == reducedVData.bflags.begin() );
  _test( reducedVData.BFlagsEnd() == reducedVData.bflags.end() );

  // .)ERASE
  VData orphanVDataCopy = orphanVData;
  orphanVDataCopy.Erase();
  _test( orphanVDataCopy.Elements() == 0 );
  _test( orphanVDataCopy.Vertices() == 0 );

  // .)IN/OUT
  // ..)binary
  VData orphanVDataCopy2 = orphanVData;
  std::fstream fp( "VData", std::ios::out | std::ios::binary );
  orphanVDataCopy2.OutBinary( fp );
  fp.close();
  fp.open( "VData", std::ios::in | std::ios::binary);
  dummyVData.InBinary( fp );
  fp.close();
  _test( dummyVData == orphanVDataCopy2 );
  // ..)ascii
  orphanVDataCopy2.OutASCII( "ASCII" );
  std::ifstream inVDataText;
  inVDataText.open( "ASCII-vdata.txt" );
  _test( inVDataText.is_open() );
  
  // TESTING ADVANCED FUNCTIONALITY
  // ==============================
  TestBinaryIO();
  
  // using handcrafted model 'MeshPatchWithLineElements' as input
  Test_EstablishElementConnectivity2D();
  
  // using handcrafted model 'MeshPatchWithLineElements' as input
  Test_CreateConsistentLineElementOrientations2D();

  // closed loops, chains along boundaries and around holes, crossing fractures, corner triangles with line elements
  Test_LineElementOrientationCornerCases();

  // 2D
  // ======================================================
  Test_RecreateConnectivityOfTriangularMesh();
  Test_RecreateConnectivityOfQuadrilateralMesh();
  
  // 3D
  // ======================================================
  // testing whether the consistent numbering of faces/neighbors can be reproduced for all-tetra model 'Tetra' (6 tets from hexahedron)
  Test_RecreateConnectivityOfTetrahedralMesh();

  // testing whether the consistent numbering of faces/neighbors can be reproduced for all-hexahedral model Rubik cube
  Test_RecreateConnectivityOfHexahedralMesh();

  Test_RecreateConnectivityOfPrismMesh();

  // SKM tests of VData mesh-fix functions
  // TODO: does not work yet:  TestReplacementOfCornerTetrahedra();

  // extra tests (stand-alone functions etc.)
  Test_InitialiseNodeTopologyIdentifiers();
  
  // global functions that modify VData
  // ==================================
  Test_refineSimplexMesh();

  Test_linearToQuadraticMeshConversion();

} // run
  


/* test read-write complete VSet, including boundaries and interfaces
*/
void VData_Test::TestBinaryIO()
 {
    VSet<2U> vset1, vset2;
    ModelTopology topo = create_BoundarySplitBoundaryPatch( vset1 );
    fstream fp( "VData_complete_model", ios::out | ios::binary );
    vset1.OutBinary( fp );
    fp.close();
    fp.open( "VData_complete_model", std::ios::in | std::ios::binary);
    vset2.InBinary( fp );
    
    _test( vset1 == vset2 );

} // end TestBinaryIO

 
 






// checks whether mesh is still intact after corner elements were split
// test model FracBox
bool VData_Test::TestReplacementOfCornerTetrahedra()
  {
     VSet<3U>      vset;
     ModelTopology topology = create_FracBox( vset );
     
     const size_t n_cnr_spanning_tets{4}; // 385, 687, 760, 821
     _test( vset.RemeshCornerSpanningTetrahedra() == n_cnr_spanning_tets );
     // adjusting pmtrl and property arrays to the new mesh size
     vset.UpdatePropertyStorage();
     
     // checking whether equidimensional mesh connectivity is still intact
     // to get isoparametric elements
     VSetConverter().ConvertElementTypesToOnesUsingLocalCoordinateSystem( vset );
     Model<3U> model( vset, "CSMP-1phase-variables.txt" );
     
     // making a Region of volumetric elements only
     Region<3U>&          model_domain = model.Region("Model");
     vector<Element<3U>*> volume_elmts, surface_elmts;
     volume_elmts.reserve( model_domain.Cells() );
     for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it ) {
          if ( (*it)->IsVolume() )
            volume_elmts.push_back( (*it) );
          else if ( (*it)->IsSurface() )
            surface_elmts.push_back( (*it) );
       }
        
     model.FormRegionFrom( "MATRIX", volume_elmts.begin(), volume_elmts.end(), false );
     
     VTK_Interface<3U>  vtk_output;
     vtk_output.OutputDataToVTK( model, "MATRIX", "FracBox_corners_fixed_elmt", "element number", 1, true );
     vtk_output.OutputDataToVTK( model, "MATRIX", "FracBox_corners_fixed_node", "node number", 1, true );

     
     // checking whether surface elements no longer match

     return true;
     
  } // end
  



/// for vsetMaker vset 'MeshPatchWithLineElements' tests whether the handcoded neighbor connectivity is reproduced
void VData_Test::Test_EstablishElementConnectivity2D()
 {
    VSet<2U>  vset1;
    // 2 fractures that are crossing so that there is a line-element manifold that needs to be disambiguated
    ModelTopology topo = create_MeshPatchWithLineElements_VSet( vset1 );
 
    // 1. making sure that EstablishElementConnectivity2D performs as expected
    VSet<2U>  vset2( vset1 );
    vset2.EstablishElementConnectivity2D();

    // TESTING
    _test( vset2 == vset1 );

 } // end Test_EstablishElementConnectivity2D






/* checks whether a neighbor-based traversal of line elements is possible after calling this method
   inside of EstablishElementConnectivity2D()
*/
void VData_Test::Test_CreateConsistentLineElementOrientations2D()
 {
    // test with predefined verified but simple model
    {
      VSet<2U>  vset1;
      // 2 fractures that are crossing so that there is a line-element manifold that needs to be disambiguated
      ModelTopology topo = create_MeshPatchWithLineElements_VSet( vset1 );
   
      // 1. making sure that EstablishElementConnectivity2D performs as expected
      VSet<2U>  vset2( vset1 );
      vset2.EstablishElementConnectivity2D();
      _test( vset2 == vset1 );
      
      // 2. checking whether the line elements can be traversed from neighbor to neighbor
      //    domain( "FRAC1", {25,29,31}
      // checking that the node0 of the next element is equal to node1 of the current one
      _test( vset2.Plist( 25, 1) ==  vset2.Plist( 29, 0) );
      _test( vset2.Plist( 29, 1) ==  vset2.Plist( 31, 0) );
      //    domain( "FRAC2", {26,27,28}
      _test( vset2.Plist( 26, 1) ==  vset2.Plist( 27, 0) );
      _test( vset2.Plist( 27, 1) ==  vset2.Plist( 28, 0) );

      // testing that the endpoints of the interfaces are correct
      // FRAC1
      _test( vset2.Plist( 25, 0) ==  5 );
      _test( vset2.Plist( 31, 1) == 15 );
      // FRAC2
      _test( vset2.Plist( 26, 0) ==  8 );
      _test( vset2.Plist( 28, 1) == 21 );
    }
    
    // testing with the faulty ANSYS model 'box2d_fault'
    {
      VSet<2U>         mesh;
      constexpr bool   isoparametric_elements( true );

      ModelTopology    mesh_topology( isoparametric_elements );
      ANSYS_Interface  mesh_interface( isoparametric_elements );

      // 0. reading mesh from ANSYS CSMP-input .asc and .dat files,
      //    eliminating unwanted line/surface element regions
      constexpr bool reassign_boundary_flags{true}, binary_input_file{true};
      mesh_interface.Read_ANSYS_Mesh( string{"fractured_slate"}, mesh, mesh_topology,
                                      binary_input_file, reassign_boundary_flags );
   
      // 1. recreating 'pfverts' information because ANSYS ICEM CFD does not get the neighbor connectivity right
      mesh.RemovePfverts();
      mesh.EstablishElementConnectivity2D(); // includes line elements
      
      // 2. Build model so that a specific region can be targeted
      Model<2U> model( mesh_topology, mesh, "example1.txt", true /* use regions file */ );
      
      // 3. test line element traversal of line-element region "Fault"
      // ---------------------------------------------------------------
      // The fault is expected to be either:
      //   (a) an open chain  — exactly 2 perimeter cells (one at each end),
      //                        each with exactly one null neighbour
      //   (b) a closed loop  — 0 perimeter cells, every element has 2 valid neighbours
      //
      // Convention (from CreateConsistentLineElementOrientations2D):
      //   Neighbor(0) = neighbour on the nd1 side (right-hand / forward neighbour)
      //   Neighbor(1) = neighbour on the nd0 side (left-hand  / backward neighbour)
      // ---------------------------------------------------------------

      auto&        fault             = model.Region("WELL_FRACTURE");
      const size_t n_fault_elements  = fault.Cells();
      const auto   region_id         = fault.DomainIndex();

      _test( n_fault_elements > 0U );

      const bool is_open_chain = ( fault.PerimeterCellsBegin() != fault.CellsEnd() );

      if ( is_open_chain )
      {
          // ---------------------------------------------------------------
          // Open-chain checks
          // ---------------------------------------------------------------

          // There must be exactly 2 perimeter cells (the two endpoints)
          {
              size_t n_perimeter{0U};
              for ( auto it = fault.PerimeterCellsBegin(); it != fault.CellsEnd(); ++it )
                  ++n_perimeter;
              _test( n_perimeter == 2U );
          }

          // Each perimeter cell must have exactly one null neighbour and one valid neighbour
          for ( auto it = fault.PerimeterCellsBegin(); it != fault.CellsEnd(); ++it )
          {
              auto* cell = *it;
              const bool nbor0_missing = ( cell->Neighbor(0) == nullptr );
              const bool nbor1_missing = ( cell->Neighbor(1) == nullptr );
              _test(  nbor0_missing != nbor1_missing );   // exactly one missing
          }

          // Traverse the chain from the first perimeter cell following Neighbor(0)
          // strictly — no fallback to Neighbor(1).
          //
          // The orientation convention is:
          //   Neighbor(0) = forward neighbour (nd1 side)
          //   Neighbor(1) = backward neighbour (nd0 side)
          //
          // Starting from a perimeter cell whose missing neighbour is on side 1
          // (i.e. Neighbor(1)==nullptr, Neighbor(0) is valid), we must be able to
          // walk the entire chain by always following Neighbor(0), and at every
          // step the cell we arrive at must know the cell we came from as its
          // Neighbor(1).  Any deviation means an element is incorrectly oriented.
          {
              // Find the endpoint whose missing neighbour is on side 1,
              // i.e. the element from which forward traversal via Neighbor(0) begins.
              auto* start = *fault.PerimeterCellsBegin();
              if ( start->Neighbor(1) != nullptr )
              {
                  // The first perimeter cell has its missing neighbour on side 0;
                  // use the second perimeter cell instead.
                  auto it = fault.PerimeterCellsBegin();
                  ++it;
                  start = *it;
              }

              // After the above, start must have Neighbor(1)==nullptr and Neighbor(0)!=nullptr
              _test( start->Neighbor(1) == nullptr );
              _test( start->Neighbor(0) != nullptr );

              auto* current  = start;
              auto* previous = static_cast<decltype(current)>( nullptr );

              size_t counter{0U};
              bool   back_links_consistent{ true };
              bool   all_cells_in_region{ true };

              while ( current != nullptr )
              {
                  ++counter;

                  if ( current->Region_ID() != region_id )
                      all_cells_in_region = false;

                  // Strict back-link check: current must know previous as Neighbor(1).
                  // No fallback — a failure here means an element is wrongly oriented.
                  if ( previous != nullptr )
                  {
                      if ( current->Neighbor(1) != previous )
                          back_links_consistent = false;
                  }

                  previous = current;
                  current  = current->Neighbor(0);  // strictly forward only
              }

              _test( all_cells_in_region );
              _test( back_links_consistent );
              _test( counter == n_fault_elements );
          }
      }
      else {
          // ---------------------------------------------------------------
          // Closed-loop checks
          // ---------------------------------------------------------------

          // Every cell in a closed loop must have exactly 2 valid neighbours
          {
              bool all_have_two_neighbours{ true };
              for ( auto it = fault.CellsBegin(); it != fault.CellsEnd(); ++it )
              {
                  auto* cell = *it;
                  if ( cell->Neighbor(0) == nullptr || cell->Neighbor(1) == nullptr )
                      all_have_two_neighbours = false;
              }
              _test( all_have_two_neighbours );
          }

          // Traverse the loop from an arbitrary seed, following Neighbor(0) strictly.
          //
          // Orientation convention:
          //   Neighbor(0) = forward neighbour (nd1 side)
          //   Neighbor(1) = backward neighbour (nd0 side)
          //
          // At every step the cell we arrive at must know the cell we came from
          // as its Neighbor(1).  Any deviation means an element is incorrectly
          // oriented.  After exactly fault.Cells() forward steps we must have
          // returned to the seed — no more, no fewer.
          {
              auto* seed     = *fault.CellsBegin();
              auto* current  = seed->Neighbor(0);  // first forward step from seed
              auto* previous = seed;

              size_t counter{1U};                   // seed counts as step 1
              bool   back_links_consistent{ true };

              // Check seed membership
              _test( seed->Region_ID() == region_id );

              // Guard: first forward step must not immediately be null
              _test( current != nullptr );

              while ( current != seed )
              {
                  ++counter;

                  // The current cell must belong to the fault region
                  _test( current->Region_ID() == region_id );

                  // Strict back-link check: current must know previous as Neighbor(1)
                  if ( current->Neighbor(1) != previous )
                      back_links_consistent = false;

                  // Guard against corrupt connectivity causing an infinite loop
                  _test( counter <= n_fault_elements );

                  previous = current;
                  current  = current->Neighbor(0);  // strictly forward only

                  // Guard: Neighbor(0) must never be null in a well-formed closed loop
                  _test( current != nullptr );
              }

              // After returning to seed: verify the closing back-link.
              // 'previous' is now the last element before the seed; the seed's
              // Neighbor(1) must be that last element.
              if ( seed->Neighbor(1) != previous )
                  back_links_consistent = false;

              _test( back_links_consistent );
              _test( counter == n_fault_elements );
          }
      }
   }
   
 } // end Test_CreateConsistentLineElementOrientations2D



// ---------------------------------------------------------------------------------------------
// Helpers for Test_LineElementOrientationCornerCases()
// ---------------------------------------------------------------------------------------------
namespace {

/// node of the test grid, addressed by its column and row
using GridNode = std::pair<int,int>;

/// index of grid node (i,j) in a grid with nx cells in x-direction
size_t gridNodeIndex( int nx, int i, int j )
 { return static_cast<size_t>( j*(nx+1) + i ); }

/**
    Creates a (nx x ny) grid of unit cells, each split into two counter-clockwise numbered linear triangles,
    and adds one linear line element for every segment of the supplied polylines. 'pfverts' is left empty.

    @param hole_cells cells that are left out, creating holes whose nodes receive the flag 'hole_flag'
    @param polylines sequences of grid nodes; successive nodes must be connected by an edge of a triangle
    @param scramble_seed if > 0, the order of the line elements and their orientations are scrambled (reproducibly)
    @param corner_triangles if false, the cell diagonals emanate from the model corners (nx, ny even) so that no triangle
           has two faces on the model boundary; if true, the diagonals cut the corners off, so that there is such a
           triangle in every corner of the model
*/
void makeTriangleGridWithLineElements( VSet<2U>& vset, int nx, int ny,
                                       const std::set<GridNode>& hole_cells, BOX_BOUNDARY hole_flag,
                                       const std::vector<std::vector<GridNode> >& polylines,
                                       unsigned scramble_seed, bool corner_triangles = false )
 {
    // nodes and their boundary flags
    const size_t nodes = static_cast<size_t>( (nx+1)*(ny+1) );
    std::deque<double>  px(nodes), py(nodes), pz(nodes,0.);
    std::vector<int8_t> bflags(nodes,NOT);
    for ( int j=0; j<=ny; ++j )
      for ( int i=0; i<=nx; ++i ) {
           const size_t n = gridNodeIndex(nx,i,j);
           px[n] = i;
           py[n] = j;
           BOX_BOUNDARY flag(NOT);
           if ( i == 0  ) flag = LEFT;
           if ( i == nx ) flag = RIGHT;
           if ( j == 0  ) flag = BOTTOM;
           if ( j == ny ) flag = TOP;
           if ( i == 0  && j == 0  ) flag = CNR1;
           if ( i == nx && j == 0  ) flag = CNR2;
           if ( i == nx && j == ny ) flag = CNR3;
           if ( i == 0  && j == ny ) flag = CNR4;
           bflags[n] = flag;
        }

    // surface elements
    std::deque<std::vector<size_t> > cells;
    for ( int j=0; j<ny; ++j )
      for ( int i=0; i<nx; ++i ) {
           const size_t a = gridNodeIndex(nx,i,j),     b = gridNodeIndex(nx,i+1,j),
                        c = gridNodeIndex(nx,i+1,j+1), d = gridNodeIndex(nx,i,j+1);
           if ( hole_cells.count( GridNode(i,j) ) > 0 ) {
                for ( const size_t n : { a, b, c, d } )
                  if ( bflags[n] == NOT ) bflags[n] = hole_flag;
                continue;
             }
           const bool diagonal_ac = ( ((i+j)%2 == 0) != corner_triangles );
           if ( diagonal_ac ) { cells.push_back( { a, b, c } ); cells.push_back( { a, c, d } ); }
           else               { cells.push_back( { a, b, d } ); cells.push_back( { b, c, d } ); }
        }
    const size_t n_triangles = cells.size();

    // line elements
    std::vector<std::pair<size_t,size_t> > segments;
    for ( const auto& polyline : polylines )
      for ( size_t k=0; k+1<polyline.size(); ++k )
        segments.push_back( std::make_pair( gridNodeIndex(nx,polyline[k].first,polyline[k].second),
                                            gridNodeIndex(nx,polyline[k+1].first,polyline[k+1].second) ) );
    if ( scramble_seed > 0 ) {
         std::mt19937 rng( scramble_seed );
         for ( size_t k=segments.size(); k>1; --k )  // Fisher-Yates shuffle, independent of the standard library
           std::swap( segments[k-1], segments[ rng() % k ] );
         for ( auto& segment : segments )
           if ( rng() & 1U ) std::swap( segment.first, segment.second );
      }
    for ( const auto& segment : segments )
      cells.push_back( { segment.first, segment.second } );

    // assembling the VSet
    std::deque<int8_t>   etypes( cells.size(), ISOPARAMETRIC_LINEAR_TRIANGLE );
    std::deque<uint32_t> npes( cells.size(), 3U ), epes( cells.size(), 3U );
    for ( size_t e=n_triangles; e<cells.size(); ++e ) {
         etypes[e] = ISOPARAMETRIC_LINEAR_BAR;
         npes[e]   = 2U;
         epes[e]   = 2U;
      }
    vset.Erase();
    vset.Resize( etypes, npes, epes, nodes, 0, 0 );
    vset.AddElementTypes( etypes.begin(), etypes.end() );
    vset.AddXYZ( px, py, pz );
    for ( size_t n=0; n<nodes; ++n ) vset.BFlag( n, bflags[n] );
    vset.AddPlist( cells.begin(), cells.end() );
    vset.RemovePfverts();
 }


/// what can be wrong with the line elements of a mesh of linear triangles and linear line elements
struct LineElementDefects {
    size_t non_reciprocal{0U};       ///< neighbor entries that are not returned by the neighbor, or neighbors that do not share the node
    size_t not_head_to_tail{0U};     ///< neighbor entries across which node 1 of one element is not node 0 of the next
    size_t reversed_on_boundary{0U}; ///< line elements numbered against the surface-element face (without neighbor) they coincide with
    size_t off_mesh{0U};             ///< line elements that do not coincide with an edge of a surface element
 };

LineElementDefects findLineElementDefects( const VSet<2U>& vset )
 {
    LineElementDefects defects;
    std::set<std::pair<size_t,size_t> > boundary_faces, edges;
    const size_t cells = vset.Cells();
    for ( size_t e=0; e<cells; ++e )
      if ( vset.PlistSize(e) == 3U )
        for ( uint32_t face=0; face<3U; ++face ) {
             // the face of a triangle is located opposite to the node with the same number
             const size_t n0 = vset.Plist( e, (face+1U)%3U ), n1 = vset.Plist( e, (face+2U)%3U );
             edges.insert( std::make_pair( std::min(n0,n1), std::max(n0,n1) ) );
             if ( vset.Pfvert( e, face ) < 0 ) boundary_faces.insert( std::make_pair( n0, n1 ) );
          }
    for ( size_t e=0; e<cells; ++e ) {
         if ( vset.PlistSize(e) != 2U ) continue;
         const size_t nd0 = vset.Plist(e,0), nd1 = vset.Plist(e,1);
         for ( uint32_t side=0; side<2U; ++side ) {
              const int64_t nbor = vset.Pfvert( e, side );
              if ( nbor < 0 ) continue;
              const size_t f = static_cast<size_t>(nbor);
              // neighbor 0 is located at node 1 and vice versa
              const size_t shared = ( side == 0U ) ? nd1 : nd0;
              if ( vset.PlistSize(f) != 2U || ( vset.Plist(f,0) != shared && vset.Plist(f,1) != shared ) ) {
                   defects.non_reciprocal++;
                   continue;
                }
              const uint32_t side_of_nbor = ( vset.Plist(f,1) == shared ) ? 0U : 1U;
              if ( vset.Pfvert( f, side_of_nbor ) != static_cast<int64_t>(e) ) defects.non_reciprocal++;
              if ( side == 0U && vset.Plist(f,0) != shared ) defects.not_head_to_tail++;
              if ( side == 1U && vset.Plist(f,1) != shared ) defects.not_head_to_tail++;
           }
         if ( boundary_faces.count( std::make_pair( nd1, nd0 ) ) > 0 ) defects.reversed_on_boundary++;
         if ( edges.count( std::make_pair( std::min(nd0,nd1), std::max(nd0,nd1) ) ) == 0 ) defects.off_mesh++;
      }
    return defects;
 }

/// sums up the defects found for 'repetitions' differently scrambled versions of the same line elements
LineElementDefects defectsOfScrambledLineElements( int nx, int ny,
                                                   const std::set<GridNode>& hole_cells, BOX_BOUNDARY hole_flag,
                                                   const std::vector<std::vector<GridNode> >& polylines,
                                                   unsigned repetitions = 20U )
 {
    LineElementDefects total;
    for ( unsigned seed=1U; seed<=repetitions; ++seed ) {
         VSet<2U> vset;
         makeTriangleGridWithLineElements( vset, nx, ny, hole_cells, hole_flag, polylines, seed );
         vset.EstablishElementConnectivity2D(); // calls CreateConsistentLineElementOrientations2D()
         const LineElementDefects defects = findLineElementDefects( vset );
         total.non_reciprocal       += defects.non_reciprocal;
         total.not_head_to_tail     += defects.not_head_to_tail;
         total.reversed_on_boundary += defects.reversed_on_boundary;
         total.off_mesh             += defects.off_mesh;
      }
    return total;
 }

} // anonymous namespace



/**
    Corner cases of CreateConsistentLineElementOrientations2D(), which is called at the end of EstablishElementConnectivity2D(),
    of the disambiguation of line-element manifolds by AngleBetweenLineElements2D(), and of SwitchCornerTriangles2D().

    All meshes are grids of triangles with line elements along selected triangle edges. The line elements are
    supplied in scrambled order and orientation (20 reproducible variants per case); afterwards
     - the neighbor entries of the line elements must be reciprocal,
     - successive line elements must be numbered head-to-tail,
     - line elements on the exterior boundary or the boundary of a hole must be numbered like the face
       of the adjacent surface element.
*/
void VData_Test::Test_LineElementOrientationCornerCases()
 {
    typedef std::vector<GridNode>  Polyline;
    const std::set<GridNode> no_holes;
    const std::set<GridNode> hole_in_centre{ {2,2}, {3,2}, {2,3}, {3,3} };
    const Polyline around_hole{ {2,2}, {3,2}, {4,2}, {4,3}, {4,4}, {3,4}, {2,4}, {2,3}, {2,2} };
    Polyline perimeter;
    for ( int i=0; i<=6; ++i ) perimeter.push_back( {i,0} );
    for ( int j=1; j<=6; ++j ) perimeter.push_back( {6,j} );
    for ( int i=5; i>=0; --i ) perimeter.push_back( {i,6} );
    for ( int j=5; j>=0; --j ) perimeter.push_back( {0,j} );

    // 1. closed loop of line elements inside of the model (has no end to start a traversal from)
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, no_holes, NOT, { around_hole } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
    }
    // 2. fracture crossing the model from the left to the right boundary
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, no_holes, NOT,
                                { Polyline{ {0,3}, {1,3}, {2,3}, {3,3}, {4,3}, {5,3}, {6,3} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
    }
    // 3. chain that follows the bottom boundary and then turns into the model, where it ends
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, no_holes, NOT,
                                { Polyline{ {0,0}, {1,0}, {2,0}, {3,0}, {3,1}, {3,2}, {3,3} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
      _test( defects.reversed_on_boundary == 0U );
    }
    // 4. chain that starts and ends with elements on the boundary, connected through the inside of the model;
    //    counter-clockwise sense on both boundaries: bottom left-to-right, top right-to-left
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, no_holes, NOT,
                                { Polyline{ {0,0}, {1,0}, {2,0}, {2,1}, {2,2}, {2,3}, {2,4}, {2,5}, {2,6}, {1,6}, {0,6} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
      _test( defects.reversed_on_boundary == 0U );
    }
    // 5. line elements all around the model plus a fracture that joins them on the left and on the right
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, no_holes, NOT,
                                { perimeter, Polyline{ {0,3}, {1,3}, {2,3}, {3,3}, {4,3}, {5,3}, {6,3} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
      _test( defects.reversed_on_boundary == 0U );
    }
    // 5b. line elements all around the model plus a fracture from corner to corner: in the corners, the fracture makes a smaller
    //     angle (45 degrees) with the boundary line elements than these make with each other (90 degrees); nevertheless, the
    //     boundary line elements must be connected with one another, so that they form one closed, counter-clockwise loop
    {
      Polyline perimeter4;
      for ( int i=0; i<=4; ++i ) perimeter4.push_back( {i,0} );
      for ( int j=1; j<=4; ++j ) perimeter4.push_back( {4,j} );
      for ( int i=3; i>=0; --i ) perimeter4.push_back( {i,4} );
      for ( int j=3; j>=0; --j ) perimeter4.push_back( {0,j} );
      const auto defects = defectsOfScrambledLineElements( 4, 4, no_holes, NOT,
                                { perimeter4, Polyline{ {0,0}, {1,1}, {2,2}, {3,3}, {4,4} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
      _test( defects.reversed_on_boundary == 0U );
    }
    // 6. loop of line elements around a hole; hole nodes flagged INTERNAL or IRREGULAR
    for ( const BOX_BOUNDARY hole_flag : { INTERNAL, IRREGULAR } )
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, hole_in_centre, hole_flag, { around_hole } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
      _test( defects.reversed_on_boundary == 0U );
    }
    // 7. chain along the left boundary that ends with an element connecting the boundary with a hole: both nodes
    //    of this element are flagged as boundary nodes, but it is not located on a boundary of the surface mesh
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, std::set<GridNode>{ {1,2} }, IRREGULAR,
                                { Polyline{ {0,4}, {0,3}, {0,2}, {1,2} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.not_head_to_tail == 0U );
      _test( defects.reversed_on_boundary == 0U );
    }
    // 8. chain that follows the bottom and the top boundary in the same direction: there is only one correct sense along
    //    a boundary, so that the boundary elements are aligned at the expense of the head-to-tail order (a warning is issued)
    {
      const auto defects = defectsOfScrambledLineElements( 6, 6, no_holes, NOT,
                                { Polyline{ {0,0}, {1,0}, {2,0}, {2,1}, {2,2}, {2,3}, {2,4}, {2,5}, {2,6}, {3,6}, {4,6} } } );
      _test( defects.non_reciprocal == 0U );
      _test( defects.reversed_on_boundary == 0U );
      _test( defects.not_head_to_tail > 0U );
    }

    // 9. angle between line elements
    {
      VSet<2U> vset;                                      // elements: n,  n+1,          n+2,          n+3
      makeTriangleGridWithLineElements( vset, 2, 2, no_holes, NOT, { Polyline{ {0,1}, {1,1}, {2,1} }, Polyline{ {1,1}, {1,2} }, Polyline{ {2,1}, {1,1} } }, 0U );
      const size_t first_line = vset.Cells() - 4U;
      VData& vdata( vset );
      _equal( vdata.AngleBetweenLineElements2D( first_line,    first_line+1U ),   0., 1.0e-10 ); // straight continuation
      _equal( vdata.AngleBetweenLineElements2D( first_line,    first_line+3U ), 180., 1.0e-10 ); // same segment direction reversed
      _equal( vdata.AngleBetweenLineElements2D( first_line,    first_line+2U ),  90., 1.0e-10 );
    }
    // 10. two straight fractures crossing at 45 degrees: at the crossing the straight continuations must be connected
    for ( unsigned seed=1U; seed<=20U; ++seed )
    {
      VSet<2U> vset;
      makeTriangleGridWithLineElements( vset, 6, 6, no_holes, NOT,
                                        { Polyline{ {1,3}, {2,3}, {3,3}, {4,3}, {5,3} }, Polyline{ {2,2}, {3,3}, {4,4} } }, seed );
      vset.EstablishElementConnectivity2D();
      const size_t crossing = gridNodeIndex( 6, 3, 3 );
      bool straight_continuations{true};
      size_t elements_at_crossing{0U};
      for ( size_t e=0; e<vset.Cells(); ++e ) {
           if ( vset.PlistSize(e) != 2U ) continue;
           for ( uint32_t nd=0; nd<2U; ++nd ) {
                if ( vset.Plist(e,nd) != crossing ) continue;
                elements_at_crossing++;
                const size_t  far_node = vset.Plist( e, 1U-nd );
                const int64_t nbor     = vset.Pfvert( e, 1U-nd ); // neighbor 0 is located at node 1 and vice versa
                if ( nbor < 0 ) { straight_continuations = false; continue; }
                const size_t f = static_cast<size_t>(nbor);
                const size_t far_node_of_nbor = ( vset.Plist(f,0) == crossing ) ? vset.Plist(f,1) : vset.Plist(f,0);
                // the crossing is the midpoint between the far nodes of two elements that continue each other
                // (the grid coordinates are integers, so that the comparison is exact)
                if ( vset.Px(far_node) + vset.Px(far_node_of_nbor) != 2.*vset.Px(crossing) ||
                     vset.Py(far_node) + vset.Py(far_node_of_nbor) != 2.*vset.Py(crossing) )
                  straight_continuations = false;
             }
        }
      _test( elements_at_crossing == 4U );
      _test( straight_continuations );
    }

    // 11. SwitchCornerTriangles2D(): triangles with two faces on the model boundary are reconnected with their neighbor ...
    {
      VSet<2U> vset;
      makeTriangleGridWithLineElements( vset, 4, 4, no_holes, NOT, { Polyline{ {1,2}, {2,2}, {3,2} } }, 0U, true /* corner triangles */ );
      bool exception_thrown{false};
      try { vset.EstablishElementConnectivity2D(); }
      catch ( const csmp::Exception& ) { exception_thrown = true; }
      _test( !exception_thrown );
      size_t corner_triangles{0U};
      for ( size_t e=0; e<vset.Cells(); ++e )
        if ( vset.PlistSize(e) == 3U ) {
             uint32_t boundary_faces{0U};
             for ( uint32_t face=0; face<3U; ++face )
               if ( vset.Pfvert( e, face ) < 0 ) boundary_faces++;
             if ( boundary_faces > 1U ) corner_triangles++;
          }
      _test( corner_triangles == 0U );
      _test( findLineElementDefects( vset ).off_mesh == 0U );
    }
    // ... unless a line element is co-located with the face that the corner triangle shares with its neighbor
    {
      VSet<2U> vset;
      makeTriangleGridWithLineElements( vset, 4, 4, no_holes, NOT, { Polyline{ {1,0}, {0,1} } }, 0U, true /* corner triangles */ );
      bool exception_thrown{false};
      try { vset.EstablishElementConnectivity2D(); }
      catch ( const csmp::Exception& ) { exception_thrown = true; }
      _test( exception_thrown );
    }

 } // end Test_LineElementOrientationCornerCases




void VData_Test::Test_RecreateConnectivityOfTriangularMesh()
 {
    VSet<2U> vset, vset_backup;
    {
      // single square test to check boundary conditions if there are no sides
      create_2_Triangle_VSet( vset );
      vset_backup = vset;
      vset.EstablishElementConnectivity2D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfTriangularMesh: 2-triangle model: testing 'bflags'..."<< endl;
      _test( vset == vset_backup );
    }

    // triangular mesh, OK 10/7/24
    create_TrianglePatch_VSet( vset );
    vset_backup = vset;
    vset.EstablishElementConnectivity2D();
    if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfTriangularMesh: 2-triangle model: testing connectivity..."<< endl;
    _test( vset == vset_backup );
 }



/**
     Mesh patch for inter-element connectivity
     Single element test for neighbor connectivity.
*/
void VData_Test::Test_RecreateConnectivityOfQuadrilateralMesh()
 {
    // 2D test cases
    {
      VSet<2U> vset;
      // quadrilateral mesh
      const int    size_sides{1}; // one square VSet
      const double dimension{5.};
      const bool   skewed=false;
      create_Square_VSet( vset, size_sides, dimension, skewed );
      VSet<2U> vset_backup = vset;
      vset.EstablishElementConnectivity2D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: 1-square model: testing model 'bflags'..."<< endl;
      _test( vset == vset_backup );
      // 4 quads
      vset = create_Quadrilateral_VSet();
      vset_backup = vset;
      vset.EstablishElementConnectivity2D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: 4 quad model: testing connectivity..."<< endl;
      _test( vset == vset_backup );
    }
 }



    /// reestablishes neighbor connectivity for 3D model Rubik cube
void VData_Test::Test_RecreateConnectivityOfTetrahedralMesh()
 {
    VSet<3U> vset;
    create_Tetra_VSet( vset );
    VSet<3U> vset_backup = vset;
    vset.EstablishElementConnectivity3D();
    if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfTetrahedralMesh: testing model 'Tetra' (6 from hexahedron)..."<< endl;
    _test( vset == vset_backup );
 }



    /// reestablishes neighbor connectivity for 3D model Rubik cube
void VData_Test::Test_RecreateConnectivityOfHexahedralMesh()
 {
    // single hexahedron
    {
      VSet<3U> vset;
      create_1Hexahedron_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: testing model '1Hexahedron_VSet'..."<< endl;
      _test( vset == vset_backup );
    }
    // all hexahedral model
    {
      VSet<3U> vset;
      create_Hexahedra_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: testing model 'Hexahedra_VSet'..."<< endl;
      _test( vset == vset_backup );
    }
    // pyramids forming a cube
    {
      VSet<3U> vset;
      create_Pyramid_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: testing model 'Pyramid_VSet'..."<< endl;
      _test( vset == vset_backup );
    }
    // mixed pyramids and hexahedra
    {
      VSet<3U> vset;
      create_Pyramid_Hexa_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: testing model 'Pyramid_Hexa_VSet'..."<< endl;
      _test( vset == vset_backup );
    }
    // mixed prisms and hexahedra
    {
      VSet<3U> vset;
      create_Prism_Hexa_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfHexahedralMesh: testing model 'create_Prism_Hexa_VSet'..."<< endl;
      _test( vset == vset_backup );
    }
 }



    /// reestablishes neighbor connectivity for 3D model Rubik cube
void VData_Test::Test_RecreateConnectivityOfPrismMesh()
 {
    // single prism
    {
      VSet<3U> vset;
      create_1Prism_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfPrismMesh: testing model '1Prism_VSet'..."<< endl;
      _test( vset == vset_backup );
    }
    // prism stack
    {
      VSet<3U> vset;
      create_Prism_VSet( vset );
      VSet<3U> vset_backup = vset;
      vset.EstablishElementConnectivity3D();
      if ( verbose_ ) cout <<"\n"<<"VData_Test::Test_RecreateConnectivityOfPrismMesh: testing model 'Prism_VSet' (54 prisms)..."<< endl;
      _test( vset == vset_backup );
    }
 }





void VData_Test::Test_InitialiseNodeTopologyIdentifiers()
 {
    // 1. Testing a VSet<2> that contains Faces and Interfaces already
    // ---------------------------------------------------------------
    {
      VSet<2U> vset1, vset2;
      ModelTopology topo = create_BoundarySplitBoundaryPatch( vset1 );
      vset2 = vset1;
      
      // recreating the node geometry identifiers, to see whether same results are obtained
      vset2.InitialiseNodeTopologyIdentifiers();
      
      _test( vset2 == vset1 );
    }

   // 2. Testing a VSet with elements only, but with line element regions
   // -------------------------------------------------------------------
   {
      VSet<2U> vset1, vset2;
      create_MeshPatchWithLineElements_VSet( vset1 );
      vset2 = vset1;
      vset2.InitialiseNodeTopologyIdentifiers();
      _test( vset2 == vset1 );
   }

   // 2. Testing the VSet FRAC_BOX with intersecting fractures & regions only
   // -----------------------------------------------------------------------
   {
      VSet<3U> vset1, vset2;
      ModelTopology topo = create_FracBox( vset1 );
      // using model topology to assign bounndary flags
      topo.AssignBoxShapedModelFlags( vset1 ); // includes FlagNodesUsingBoundaryDomains(vset1);
      
      vset2 = vset1;
      vset2.InitialiseNodeTopologyIdentifiers();
      _test( vset2 == vset1 );
   }

 } // end Test_InitialiseNodeTopologyIdentifiers






/*
    Refines Simplex (line,triangle,tetrahera) mesh by introducing midside side nodes
    Checks: that new number of elements is correct
*/
void VData_Test::Test_refineSimplexMesh()
 {
    VSet<3> vset;
    ModelTopology topo = create_FracBox( vset );
    
    // counting line, triangle and tetra element to check whether new numbers are correct
    size_t n_lines{0}, n_tria{0}, n_tets{0};
    for ( const auto& pit : vset.plist )
      switch( pit.size() ) {
        case 2: n_lines++; break;
        case 3: n_tria++;  break;
        case 4: n_tets++;
      }
    size_t n_nodes_before{ vset.Vertices() };
    
    // caculating new numbers
    const size_t n_elmts_new = n_lines * 2 + n_tria * 4 + n_tets * 8;
    
    // refining mesh
    vector<size_t> old_to_new_idx = refineSimplexMesh( vset );
    
    size_t n_nodes_after{ vset.Vertices() };
    size_t n_elmts_after{ vset.Elements() };
    
    // test
    _test ( n_nodes_after > n_nodes_before );
    _test ( n_elmts_after == n_elmts_new );
    _test ( old_to_new_idx.size() == n_elmts_new );
    
    // consistency check
    vset.EstablishElementConnectivity3D();
 }





// helper that finds midsize nodes by averaging corner node coordinates
static uint32_t midSideNode( const VSet<3>& vs, size_t elmt, uint32_t n1, uint32_t n2 ) {
   Point<3> nd1pt( vs.Px(vs.Plist(elmt,n1)), vs.Py(vs.Plist(elmt,n1)), vs.Pz(vs.Plist(elmt,n1)) ),
            nd2pt( vs.Px(vs.Plist(elmt,n2)), vs.Py(vs.Plist(elmt,n2)), vs.Pz(vs.Plist(elmt,n2)) );
   Point<3> midpt = (nd1pt + nd2pt) / 2.;
   // loop over the nodes of the element until the midpoint is found
   uint32_t n_midside_node{0};
   for ( uint32_t i{0}; i<vs.PlistSize(elmt); ++ i ) {
         Point<3> ndpt( vs.Px(vs.Plist(elmt,i)), vs.Py(vs.Plist(elmt,i)), vs.Pz(vs.Plist(elmt,i)) );
         if ( ndpt == midpt ) {
              n_midside_node = i;
              break;
           }
     }
   return n_midside_node;
}


void VData_Test::Test_linearToQuadraticMeshConversion()
 {
    // testing a simplex mesh
    // ----------------------
    {
      VSet<3> vset;
      create_FracBox( vset );
      convertLinearToQuadraticSimplexElementMesh( vset );
      
      _test( vset.Vertices() == 2233 ); // as we know from refinement test
      
      // element types
      for ( auto cit=vset.PelmtBegin(); cit!=vset.PelmtEnd(); ++cit ) {
            // _test( fe.IsSimplex(*cit) == true );
            const auto et = static_cast<CSMP_FEM_TYPE>(*cit);
            _test( isLineElement(et) || isTriangular(et) || isTetrahedral(et) );
        }
    }
    // testing polyhedral prism mesh
    {
      const bool          isoparametric(true);
      ANSYS_Interface     mesh_interface(isoparametric);
      VSet<3>             vset;
      ModelTopology       mesh_topology(isoparametric);
      const bool binary_file( true );
      mesh_interface.Read_ANSYS_Mesh( "prism_test", vset, mesh_topology, binary_file, true );
      // TODO: replace ANSYS mesh; it has no hexa and is too big for a unit test
      
      convertLinearToQuadraticPolyElementTypeMesh( vset );
 
      // checking that the nodes are assigned to the elements in the expected way
      size_t elmt{0};
      auto it=vset.PlistBegin(), itend=vset.PlistEnd();
      for ( auto et=vset.PelmtBegin(); it!=itend && et!=vset.PelmtEnd(); ++it, ++et ) {
            const auto etype = static_cast<CSMP_FEM_TYPE>(*et);
           switch( etype ) {
              case ISOPARAMETRIC_QUADRATIC_BAR:
                   _test( midSideNode( vset, elmt, 0, 1 ) == 2 ); // testing midSideNode()
                break;
              case ISOPARAMETRIC_QUADRATIC_TRIANGLE:
                   _test( midSideNode( vset, elmt, 0, 1 ) == 3 );
                   _test( midSideNode( vset, elmt, 1, 2 ) == 4 );
                   _test( midSideNode( vset, elmt, 2, 0 ) == 5 );
                break;
             case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON:
                   _test( midSideNode( vset, elmt, 0, 1 ) == 4 );
                   _test( midSideNode( vset, elmt, 1, 2 ) == 5 );
                   _test( midSideNode( vset, elmt, 2, 0 ) == 6 );
                   _test( midSideNode( vset, elmt, 3, 0 ) == 7 );
                   _test( midSideNode( vset, elmt, 3, 1 ) == 8 );
                   _test( midSideNode( vset, elmt, 3, 2 ) == 9 );
                break;
              case ISOPARAMETRIC_QUADRATIC_PRISM15:
                   _test( midSideNode( vset, elmt, 0, 1 ) == 6 );
                   _test( midSideNode( vset, elmt, 1, 2 ) == 7 );
                   _test( midSideNode( vset, elmt, 2, 0 ) == 8 );
                   _test( midSideNode( vset, elmt, 3, 0 ) == 9 );
                   _test( midSideNode( vset, elmt, 4, 1 ) == 10 );
                   _test( midSideNode( vset, elmt, 5, 2 ) == 11 );
                   _test( midSideNode( vset, elmt, 3, 4 ) == 12 );
                   _test( midSideNode( vset, elmt, 4, 5 ) == 13 );
                   _test( midSideNode( vset, elmt, 5, 3 ) == 14 );
                break;
              case ISOPARAMETRIC_QUADRATIC_PYRAMID13:
                   _test( midSideNode( vset, elmt, 0, 1 ) == 5 );
                   _test( midSideNode( vset, elmt, 1, 2 ) == 6 );
                   _test( midSideNode( vset, elmt, 2, 3 ) == 7 );
                   _test( midSideNode( vset, elmt, 3, 0 ) == 8 );
                   _test( midSideNode( vset, elmt, 4, 0 ) == 9 );
                   _test( midSideNode( vset, elmt, 4, 1 ) == 10 );
                   _test( midSideNode( vset, elmt, 4, 2 ) == 11 );
                   _test( midSideNode( vset, elmt, 4, 3 ) == 12 );
                break;
              case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20:
                   _test( midSideNode( vset, elmt, 0, 1 ) == 8 );
                   _test( midSideNode( vset, elmt, 1, 2 ) == 9 );
                   _test( midSideNode( vset, elmt, 2, 3 ) == 10 );
                   _test( midSideNode( vset, elmt, 3, 0 ) == 11 );
                   _test( midSideNode( vset, elmt, 4, 0 ) == 12 );
                   _test( midSideNode( vset, elmt, 5, 1 ) == 13 );
                   _test( midSideNode( vset, elmt, 6, 2 ) == 14 );
                   _test( midSideNode( vset, elmt, 7, 3 ) == 15 );
                   _test( midSideNode( vset, elmt, 4, 5 ) == 16 );
                   _test( midSideNode( vset, elmt, 5, 6 ) == 17 );
                   _test( midSideNode( vset, elmt, 6, 7 ) == 18 );
                   _test( midSideNode( vset, elmt, 7, 4 ) == 19 );
                break;
              default:
                throw csmp::Exception( ERROR, "VData_Test::Test_linearToQuadraticMeshConversion",
                                       parseFiniteElementType(etype), "Element type not handled." );
           }
          elmt++;
        }
    }

 } // end Test_linearToQuadraticMeshConversion



} // end csmp
