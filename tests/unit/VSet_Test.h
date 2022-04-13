#ifndef VSET_TEST_H
#define VSET_TEST_H

#include "Test.h"
#include "VSet.h"
#include "vsetMakers.h"
#include "FEM_Data.h"

namespace csmp {

/// PL Nov 2010
class VSet_Test : public Test
{
public:
  virtual void run();
};

void VSet_Test::run()
{
  // .)CONSTRUCTORS
  VSet<3U> vset1( 8, 6, LINEAR_TETRAHEDRON, 12, 2 );
  _test( vset1.Elements() == 2 );
  _test( vset1.Vertices() == 12 );
  _test( vset1.HybridElementTypeMesh() == false );
  _test( vset1.PfvertsSize( 0 ) == 6 );
  _test( vset1.PfvertsSize( 1 ) == 6 );
  _test( vset1.PlistSize( 0 ) == 8 );
  _test( vset1.PlistSize( 1 ) == 8 );
  _test( vset1.ElementType( 0 ) == LINEAR_TETRAHEDRON );
  _test( vset1.MeshDimension() == 3 );

  // vset2: poly-element mesh
  std::deque<uint32_t> npe;
  npe.push_back( 8 ); npe.push_back( 8 ); npe.push_back( 5 );
  std::deque<uint32_t> epe;
  epe.push_back( 6 ); epe.push_back( 6 ); epe.push_back( 5 );
  vector<int8_t> fem_types = { ISOPARAMETRIC_LINEAR_HEXAHEDRON, ISOPARAMETRIC_LINEAR_HEXAHEDRON, ISOPARAMETRIC_LINEAR_PYRAMID };
  
  VSet<3U> vset2( fem_types, npe, epe, 13 );
  _test( vset2.Elements() == 3 );
  _test( vset2.Vertices() == 13 );
  _test( vset2.HybridElementTypeMesh() == true );
  _test( vset2.PfvertsSize( 0 ) == 6 );
  _test( vset2.PfvertsSize( 1 ) == 6 );
  _test( vset2.PfvertsSize( 2 ) == 5 );
  _test( vset2.PlistSize( 0 ) == 8 );
  _test( vset2.PlistSize( 1 ) == 8 );
  _test( vset2.PlistSize( 2 ) == 5 );
  
  VSet<3U> vset3( vset1 );
  _test( vset3.Elements() == 2 );
  _test( vset3.Vertices() == 12 );
  _test( vset3.HybridElementTypeMesh() == false );
  _test( vset3.PfvertsSize( 0 ) == 6 );
  _test( vset3.PfvertsSize( 1 ) == 6 );
  _test( vset3.PlistSize( 0 ) == 8 );
  _test( vset3.PlistSize( 1 ) == 8 );
  _test( vset3.ElementType( 0 ) == LINEAR_TETRAHEDRON );
  
  VSet<3U> vset4 = vset1;
  _test( vset4.Elements() == 2 );
  _test( vset4.Vertices() == 12 );
  _test( vset4.PfvertsSize( 0 ) == 6 );
  _test( vset4.PfvertsSize( 1 ) == 6 );
  _test( vset4.PlistSize( 0 ) == 8 );
  _test( vset4.PlistSize( 1 ) == 8 );
  _test( vset1.HybridElementTypeMesh() == false );

  // .)RESIZE
  vset4.Resize(  8, 6, LINEAR_TETRAHEDRON, 12, 2 );
  _test( vset4.Elements() == 2 );
  _test( vset4.Vertices() == 12 );
  _test( vset4.HybridElementTypeMesh() == false );
  _test( vset4.PfvertsSize( 0 ) == 6 );
  _test( vset4.PfvertsSize( 1 ) == 6 );
  _test( vset4.PlistSize( 0 ) == 8 );
  _test( vset4.PlistSize( 1 ) == 8 );
  _test( vset4.ElementType( 0 ) == LINEAR_TETRAHEDRON );
  std::deque<int8_t> etypes;
  etypes.push_back( ISOPARAMETRIC_LINEAR_HEXAHEDRON );
  etypes.push_back( ISOPARAMETRIC_LINEAR_HEXAHEDRON );
  etypes.push_back( ISOPARAMETRIC_LINEAR_PYRAMID );
  vset4.Resize( etypes, npe, epe, 13, 0, 0 );
  _test( vset4.Elements() == 3 );
  _test( vset4.Vertices() == 13 );
  _test( vset4.HybridElementTypeMesh() == true );
  _test( vset4.PfvertsSize( 0 ) == 6 );
  _test( vset4.PfvertsSize( 1 ) == 6 );
  _test( vset4.PfvertsSize( 2 ) == 5 );
  _test( vset4.PlistSize( 0 ) == 8 );
  _test( vset4.PlistSize( 1 ) == 8 );
  _test( vset4.PlistSize( 2 ) == 5 );
  _test( vset4.ElementType( 0 ) == ISOPARAMETRIC_LINEAR_HEXAHEDRON );
  _test( vset4.ElementType( 1 ) == ISOPARAMETRIC_LINEAR_HEXAHEDRON );
  _test( vset4.ElementType( 2 ) == ISOPARAMETRIC_LINEAR_PYRAMID );

  // .)DATA OPS
  // TODO: use PropertyData interface rather than the deprecated FEM_Data interface
  std::deque<double> px, py, pz;
  for( auto i = 0; i < 13; ++i )
  {
    px.push_back( (double)i*10 );
    py.push_back( (double)i*10 );
    pz.push_back( (double)i*10 );
  }
  vset4.AddXYZ( px, py, pz );
  for( auto i = 0; i < 13; ++i )
  {
    _test( vset4.Px( i ) == (double)i*10 );
    _test( vset4.Py( i ) == (double)i*10 );
    _test( vset4.Pz( i ) == (double)i*10 );
  }
  std::map<size_t,std::vector<int64_t> > plist;
  std::vector<int64_t> plist1;
  plist1.push_back( 0 ); plist1.push_back( 1 ); plist1.push_back( 4 ); plist1.push_back( 3 );
  plist1.push_back( 6 ); plist1.push_back( 7 ); plist1.push_back( 10 ); plist1.push_back( 9 );
  std::vector<int64_t> plist2;
  plist2.push_back( 1 ); plist2.push_back( 2 ); plist2.push_back( 5 ); plist2.push_back( 4 );
  plist2.push_back( 7 ); plist2.push_back( 8 ); plist2.push_back( 11 ); plist2.push_back( 10 );
  std::vector<int64_t> plist3;
  plist3.push_back( 2 ); plist3.push_back( 5 ); plist3.push_back( 11 ); plist3.push_back( 8 ); plist3.push_back( 12 );
  plist.insert( make_pair( 0, plist1 ) );
  plist.insert( make_pair( 1, plist2 ) );
  plist.insert( make_pair( 2, plist3 ) );
  vset4.AddPlist( plist.begin(), plist.end() );
  _test( vset4.PlistSize( 0 ) == 8 ); _test( vset4.PlistSize( 1 ) == 8 ); _test( vset4.PlistSize( 2 ) == 5 );
  _test( vset4.Plist( 0, 0 ) == 0 ); _test( vset4.Plist( 0, 1 ) == 1 ); _test( vset4.Plist( 0, 2 ) == 4 ); _test( vset4.Plist( 0, 3 ) == 3 );
  _test( vset4.Plist( 0, 4 ) == 6 ); _test( vset4.Plist( 0, 5 ) == 7 ); _test( vset4.Plist( 0, 6 ) == 10 ); _test( vset4.Plist( 0, 7 ) == 9 );
  _test( vset4.Plist( 1, 0 ) == 1 ); _test( vset4.Plist( 1, 1 ) == 2 ); _test( vset4.Plist( 1, 2 ) == 5 ); _test( vset4.Plist( 1, 3 ) == 4 );
  _test( vset4.Plist( 1, 4 ) == 7 ); _test( vset4.Plist( 1, 5 ) == 8 ); _test( vset4.Plist( 1, 6 ) == 11 ); _test( vset4.Plist( 1, 7 ) == 10 );
  _test( vset4.Plist( 2, 0 ) == 2 ); _test( vset4.Plist( 2, 1 ) == 5 ); _test( vset4.Plist( 2, 2 ) == 11 ); _test( vset4.Plist( 2, 3 ) == 8 ); _test( vset4.Plist( 2, 4 ) == 12 );
  std::deque<std::vector<int64_t> > plistDeque;
  plistDeque.push_back( plist1 );
  plistDeque.push_back( plist2 );
  plistDeque.push_back( plist3 );
  vset4.AddPlist( plistDeque.begin(), plistDeque.end() );
  _test( vset4.PlistSize( 0 ) == 8 ); _test( vset4.PlistSize( 1 ) == 8 ); _test( vset4.PlistSize( 2 ) == 5 );
  _test( vset4.Plist( 0, 0 ) == 0 ); _test( vset4.Plist( 0, 1 ) == 1 ); _test( vset4.Plist( 0, 2 ) == 4 ); _test( vset4.Plist( 0, 3 ) == 3 );
  _test( vset4.Plist( 0, 4 ) == 6 ); _test( vset4.Plist( 0, 5 ) == 7 ); _test( vset4.Plist( 0, 6 ) == 10 ); _test( vset4.Plist( 0, 7 ) == 9 );
  _test( vset4.Plist( 1, 0 ) == 1 ); _test( vset4.Plist( 1, 1 ) == 2 ); _test( vset4.Plist( 1, 2 ) == 5 ); _test( vset4.Plist( 1, 3 ) == 4 );
  _test( vset4.Plist( 1, 4 ) == 7 ); _test( vset4.Plist( 1, 5 ) == 8 ); _test( vset4.Plist( 1, 6 ) == 11 ); _test( vset4.Plist( 1, 7 ) == 10 );
  _test( vset4.Plist( 2, 0 ) == 2 ); _test( vset4.Plist( 2, 1 ) == 5 ); _test( vset4.Plist( 2, 2 ) == 11 ); _test( vset4.Plist( 2, 3 ) == 8 ); _test( vset4.Plist( 2, 4 ) == 12 );
  std::vector<int64_t> pfverts1;
  pfverts1.push_back( 99 ); pfverts1.push_back( 1 ); pfverts1.push_back( 99 ); pfverts1.push_back( 99 ); pfverts1.push_back( 99 ); pfverts1.push_back( 99 );
  std::vector<int64_t> pfverts2;
  pfverts2.push_back( 99 ); pfverts2.push_back( 3 ); pfverts2.push_back( 99 ); pfverts2.push_back( 0 ); pfverts2.push_back( 99 ); pfverts2.push_back( 99 );
  std::vector<int64_t> pfverts3;
  pfverts3.push_back( 1 ); pfverts3.push_back( 99 ); pfverts3.push_back( 99 ); pfverts3.push_back( 99 ); pfverts3.push_back( 99 ); pfverts3.push_back( 99 );
  std::map<size_t,std::vector<int64_t> > pfverts;
  pfverts.insert( make_pair( 0, pfverts1 ) );
  pfverts.insert( make_pair( 1, pfverts2 ) );
  pfverts.insert( make_pair( 2, pfverts3 ) );
  vset4.AddPfverts( pfverts.begin(), pfverts.end() );
  _test( vset4.Pfvert( 0, 0 ) == 99 );  _test( vset4.Pfvert( 0, 1 ) == 1 );  _test( vset4.Pfvert( 0, 2 ) == 99 );
  _test( vset4.Pfvert( 0, 3 ) == 99 );  _test( vset4.Pfvert( 0, 4 ) == 99 );  _test( vset4.Pfvert( 0, 5 ) == 99 );
  _test( vset4.Pfvert( 1, 0 ) == 99 );  _test( vset4.Pfvert( 1, 1 ) == 3 );  _test( vset4.Pfvert( 1, 2 ) == 99 );
  _test( vset4.Pfvert( 1, 3 ) == 0 );  _test( vset4.Pfvert( 1, 4 ) == 99 );  _test( vset4.Pfvert( 1, 5 ) == 99 );
  _test( vset4.Pfvert( 2, 0 ) == 1 );  _test( vset4.Pfvert( 2, 1 ) == 99 );  _test( vset4.Pfvert( 2, 2 ) == 99 );
  _test( vset4.Pfvert( 2, 3 ) == 99 );  _test( vset4.Pfvert( 2, 4 ) == 99 );
  std::deque<std::vector<int64_t> > pfvertsDeque;
  pfvertsDeque.push_back( pfverts1 ); pfvertsDeque.push_back( pfverts2 ); pfvertsDeque.push_back( pfverts3 );
  vset4.AddPfverts( pfvertsDeque.begin(), pfvertsDeque.end() );
  _test( vset4.Pfvert( 0, 0 ) == 99 );  _test( vset4.Pfvert( 0, 1 ) == 1 );  _test( vset4.Pfvert( 0, 2 ) == 99 );
  _test( vset4.Pfvert( 0, 3 ) == 99 );  _test( vset4.Pfvert( 0, 4 ) == 99 );  _test( vset4.Pfvert( 0, 5 ) == 99 );
  _test( vset4.Pfvert( 1, 0 ) == 99 );  _test( vset4.Pfvert( 1, 1 ) == 3 );  _test( vset4.Pfvert( 1, 2 ) == 99 );
  _test( vset4.Pfvert( 1, 3 ) == 0 );  _test( vset4.Pfvert( 1, 4 ) == 99 );  _test( vset4.Pfvert( 1, 5 ) == 99 );
  _test( vset4.Pfvert( 2, 0 ) == 1 );  _test( vset4.Pfvert( 2, 1 ) == 99 );  _test( vset4.Pfvert( 2, 2 ) == 99 );
  _test( vset4.Pfvert( 2, 3 ) == 99 );  _test( vset4.Pfvert( 2, 4 ) == 99 );

  // BOUNDARY FLAGS
  vector<std::int8_t> bfmap( 13, INTERNAL );
  vset4.AddBFlags( bfmap.begin(), bfmap.end() );
  for( auto it = vset4.BFlagsBegin(); it != vset4.BFlagsEnd(); ++it )
    _test( (*it) == INTERNAL );


  // .)IO
  vset4.OutputTo( "vsetBIN", 0. );
  VSet<3U> vset5( vset4 ); double time( 0 );
  vset5.Erase();
  _test( vset5.Vertices() == 0 );
  vset5.InputFrom( "vsetBIN", time );
  _test( vset5 == vset4 );
// TODO: SKM: Parallel output needs to be refactored to work with PropertyData
  vset4.ParallelOutputTo( "vsetBINp", 0., 0 );
  vset5.Erase();
  _test( vset5.Elements() == 0 );
  size_t halo( 0 );
  vset5.ParallelInputFrom( "vsetBINp", time, halo );
  _test( vset5 == vset4 );

  // .)REDUCETO
  std::map<size_t,size_t> reduceMap;
  reduceMap.insert( std::make_pair( 1, 0 ) );
  reduceMap.insert( std::make_pair( 2, 1 ) );
  vset4.ReduceTo( reduceMap );
  _test( vset4.Elements() == 2 );
  _test( vset4.Vertices() == 9 );
  _test( vset4.PlistSize( 0 ) == 8 );
  _test( vset4.PlistSize( 1 ) == 5 );
  _test( vset4.ElementType( 0 ) == ISOPARAMETRIC_LINEAR_HEXAHEDRON );
  _test( vset4.ElementType( 1 ) == ISOPARAMETRIC_LINEAR_PYRAMID );

  // .)REMOVEDATA
  vset4.RemoveBflags();
  vset4.RemoveData( "scalar data" );	
  vset4.RemoveData( "vector data" );	
  vset4.RemoveData( "tensor data" );	
  _test( vset4.DataEmpty() == true );
  _test( vset4.PropertyValuesBegin() == vset4.PropertyValuesEnd() );

// NEW TESTS BY SKM
// comparitor test
  VSet<3U> vset6;
  test_Create_Prism_Hexa_VSet( vset6, false /* bSkewed */ );
  _test( vset6.MeshDimension() == 3 );
  vset6.OutputTo( "vsetBIN", 0. );
  VSet<3U> vset7;
  vset7.InputFrom( "vsetBIN", time );
  _test( vset7 == vset6 );

} // run

} // csmp


#endif // VSET_TEST_H
