

#ifndef VDATA_TEST_H
#define VDATA_TEST_H

#include "Test.h"
#include "VData.h"

namespace csmp{

/// PL Nov 2010
class VData_Test : public Test
{
public:
  virtual void run();
};




/** TESTING DATA MEMBERS

    std::vector<double64>             px, py, pz; (1)
    std::deque<std::vector<size_t> >  plist;      (2)
    std::deque<std::vector<int32> >   pfverts;    (3)
    std::map<size_t,int32>            bflags;     (4)
    std::vector<int32>                pelmt;      (5)

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
  // .)CONSTRUCTORS
  VData defaultVData;
  std::deque<size_t> npes;
  std::deque<size_t> epes;
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
  std::deque<size_t> npes2;
  npes2.push_back( 5 );
  npes2.push_back( 5 );
  npes2.push_back( 5 );
  std::deque<size_t> epes2;
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
  std::deque<size_t> plistSize;
  plistSize.push_back( 9 );
  plistSize.push_back( 9 );
  plistSize.push_back( 9 );
  dequeVData.ResizePlist( plistSize );
  _test( dequeVData.PlistSize( 1 ) == 9 );
  std::deque<size_t> pfvertsDeque;
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
  vdata2.AddBFlag( 1, 100 );
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
  for( size_t i = 0; i < reducedVData.Vertices(); ++i )
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
  /* TODO: DOES NOT WORK(ASSERTION), PROBABLY DUE TO INCOMPATIBILITY OF OutASCII and InText formatting
  VData orphanVDataCopy3 = orphanVData;
  orphanVDataCopy3.Erase();
  orphanVDataCopy3.InText( inVDataText );
  inVDataText.close();
  _test( orphanVDataCopy3 == orphanVDataCopy2 );
  */

} // run

} // csmp

#endif // VDATA_TEST_H
