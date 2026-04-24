//
//  VData_Test.cpp
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 26/1/2022.
//  Copyright © 2022 Stephan Matthai. All rights reserved.
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
#include "Element.h"
#include "ANSYS_Interface.h"
#include "VTK_Interface.h"
#include "vsetMakers.h"
#include "Exception.h"

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
    
 } // end Test_CreateConsistentLineElementOrientations2D




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
      // TODO: Get this test to pass
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
