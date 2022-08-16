#include "FiniteElement_Test.h"
#include "CSMP_highLevelUtilities.h"
#include "PL_Utilities.h"
#include "GlobalVerbose.h"
#include "compareFloats.h"

#include <fstream>

using namespace std;

namespace csmp{

FiniteElement_Test::FiniteElement_Test( FiniteElement* testee,
                                        const char* results_file,
                                        bool verbose )
  : fileName_( results_file ),
    femPtr_( testee ),
    femData_(verbose),
    verbose_(verbose)
{
  // debugging
  assert( testee != NULL );
}

FiniteElement_Test::~FiniteElement_Test()
{
    if( femPtr_!=NULL)
        delete femPtr_;
}




/// tests fem operations against provided data in text file
void FiniteElement_Test::run()
{
  bool global_verbose( GlobalVerbose::Instance().globalVerbose );

  // .) READING FEM TEST DATA & IO
  ifstream  femDataFile( fileName_.c_str(), ifstream::in );
  femDataFile >> femData_;
  string testName = "csmp::"; testName += femData_.FemName(); testName += "_Test";
  replaceWhiteSpaceBy( testName, '_' );
  setName( testName );
  vector<double> var_v;
  if ( verbose_ ) cout << "STARTING CORRESPONDING FEM TEST:\n";

  // .) ID
  if ( verbose_ ) cout << "Testing ID ops...\n";
  const size_t feID( 123 );
  femPtr_->CurrentID( feID );
  _test( femPtr_->CurrentID() == feID );

  // .) DIM
  if ( verbose_ ) cout << "Testing DIM ops...\n";
  _test( femData_.Dim() == femPtr_->Dim() );

  // .) NODES
  if ( verbose_ ) cout << "Testing Node ops...\n";
  _test( femData_.NodeCount() == femPtr_->Nodes() );

  // .) SEGMENTS
  if ( verbose_ ) cout << "Testing Segement ops...\n";
  _test( femData_.SegmentCount() == femPtr_->Segments() );

  // .) FACES
  if ( verbose_ ) cout << "Testing Face ops...\n";
  _test( femData_.FaceCount() == femPtr_->Faces() );

  // .) NEIGHBOURS
  if ( verbose_ ) cout << "Testing Neighbor ops...\n";
  _test( femData_.NeighborCount() == femPtr_->Neighbors() );

  // .) NODES PER FACE
  if ( verbose_ ) cout << "Testing nodes per faces...\n";
  for( auto i = 0; i < femPtr_->Faces(); ++ i )
    _test( femData_.NodesPerFace( i ) == femPtr_->NodesPerFace( i ) );

  // .) INTEGRATION POINTS
  if ( verbose_ ) cout << "Testing Integration points...\n";
  _test( femData_.IntegrationPointCount() == femPtr_->IntegrationPoints() );

  // .) INTERPOLATION NUMBER
  if ( verbose_ ) cout << "Testing interpolation...\n";
  _test( femData_.InterpolationNumber() == femPtr_->Interpolation() );

  // .) ISOPARAMETRIC
  if ( verbose_ ) cout << "Testing isoparametricity...\n";
  _test( femData_.Isoparametric() == femPtr_->Isoparametric() );

  // .) LOCAL/GLOBAL COORDS
  if ( verbose_ ) cout << "Testing local vs global coordinates...\n";
  _test( femData_.UsesLocalCoordinates() == femPtr_->UsesLocalCoordinates() );

  // .) ORDER OF SHAPE FUNCTIONS
  if ( verbose_ ) cout << "Testing shape function order...\n";
  _test( femData_.OrderOfShapeFunctions() == femPtr_->OrderOfShapeFunctions() );

  // .) LINE ELEMENT
  if ( verbose_ ) cout << "Testing element type...\n";
  _test( femData_.LineElement() == femPtr_->IsLine() );

  // .) SURFACE ELEMENT
  _test( femData_.SurfaceElement() == femPtr_->IsSurface() );

  // .) VOLUME ELEMENT
  _test( femData_.VolumeElement() == femPtr_->IsVolume() );

  // .) CSMP FEM TYPE
  _test( femData_.ElementType() == femPtr_->ElementType() );

  //    INITIALIZING NODES/COORDINATE MATRIX
  if ( verbose_ ) cout << "Initializing nodes...\n";
  femPtr_->XY.Resize( femData_.NodeCount(), 3 );
  for ( auto i = 0; i < femData_.NodeCount(); ++i )
    femPtr_->XY.AssignRow( i, femData_.NodePtr( i )->Coordinate() );

  // .) VOLUME
  if ( verbose_ ) cout << "Testing volume...\n";
  _equal( femData_.Volume(), femPtr_->Volume(), femData_.Tolerance() );

  // .) ASPECT RATIO
  if ( verbose_ ) cout << "Testing aspect ratio...\n";
  _equal( femData_.AspectRatio(), femPtr_->AspectRatio(), femData_.Tolerance() );

  // .) INNER RADIUS
  if ( verbose_ ) cout << "Testing inner radius...\n";
  _equal( femData_.InnerRadius(), femPtr_->InnerRadius(), femData_.Tolerance() );

  // .) EDGE LENGTHS
  if ( verbose_ ) cout << "Testing edges...\n";
  vector<double> edgeLengths;
  femPtr_->EdgeLengths( edgeLengths );
  _test( edgeLengths.size() == femData_.EdgeCount() );
  for( auto i = 0; i < edgeLengths.size(); ++i )
    _equal( edgeLengths.at( i ),  femData_.EdgeLength( i ), femData_.Tolerance() );

  // .) SEGMENT NODES
  if ( verbose_ ) cout << "Testing segment nodes...\n";
  vector<uint32_t> segmentNodes;
  for( auto i = 0; i < femPtr_->Segments(); ++i )
  {
    femPtr_->NodesOfSegment( i, segmentNodes );
    if( global_verbose )
      if ( verbose_ ) cout << "  Segment " << i << endl;
    for( auto j = 0; j < segmentNodes.size(); ++ j )
    {
      if( global_verbose )
        if ( verbose_ ) cout << "    Node " << j << ": " << segmentNodes.at( j ) << " versus " << femData_.NodeOfSegment( i, j ) << endl;
      _test( segmentNodes.at( j ) == femData_.NodeOfSegment( i, j ) );
    }
  }

  // .) FACE NODES
  if ( verbose_ ) cout << "Testing face nodes...\n";
  vector<uint32_t> faceNodes;
  for( auto i = 0; i < femPtr_->Faces(); ++i )
  {
    faceNodes = femPtr_->NodesOfFace( i );
    if( global_verbose )
      if ( verbose_ ) cout << "  Face " << i << endl;
    for( auto j = 0; j < faceNodes.size(); ++j )
    {
      if(global_verbose )
        if ( verbose_ ) cout << "    Node " << j << ": " << faceNodes.at( j ) << " versus " << femData_.NodeOfFace( i, j ) << endl;
      _test( faceNodes.at( j ) == femData_.NodeOfFace( i, j ) );
    }
  }

  // .) CORNER NODES
  if ( verbose_ ) cout << "Testing corner nodes...\n";
  vector<uint32_t> cornerNodes;
  femPtr_->CornerNodes( cornerNodes );
  for( auto i = 0; i < cornerNodes.size(); ++i )
  {
    _test( cornerNodes.at( i ) == femData_.NodeAtCorner( i ) );
    if( global_verbose )
      if ( verbose_ ) cout << "  Corner " << i << ": Node " << cornerNodes.at( i ) << " versus " << femData_.NodeAtCorner( i ) << endl;
  }
  _test( cornerNodes.size() == femPtr_->CornerNodes() );
  _test( femPtr_->CornerNodes() == femData_.CornerNodeCount() );

  // .) MIDSIDE NODES
  if ( verbose_ ) cout << "Testing midside nodes...\n";
  if( femData_.MidsideNodeCount() != 0 )
  {
    vector<uint32_t> midsideNodes;
    femPtr_->MidSideNodes( midsideNodes );
    for( auto i = 0; i < midsideNodes.size(); ++i )
    {
      _test( midsideNodes.at( i ) == femData_.NodeAtMidside( i ) );
      if( global_verbose )
        if ( verbose_ ) cout << "  Midside " << i << ": Node " << midsideNodes.at( i ) << " versus " << femData_.NodeAtMidside( i ) << endl;
    }
    _test( midsideNodes.size() == femPtr_->MidSideNodes() );
    _test( femPtr_->MidSideNodes() == femData_.MidsideNodeCount() );
  }

  // .) INTERIOR NODES
  if ( verbose_ ) cout << "Testing interior nodes...\n";
  if( femData_.InteriorNodeCount() != 0 )
  {
    vector<uint32_t> interiorNodes;
    femPtr_->InteriorNodes( interiorNodes );
    for( auto i = 0; i < interiorNodes.size(); ++i )
    {
      _test( interiorNodes.at( i ) == femData_.NodeAtInterior( i ) );
      if( global_verbose )
        if ( verbose_ ) cout << "  Interior " << i << ": Node " << interiorNodes.at( i ) << " versus " << femData_.NodeAtInterior( i ) << endl;
    }
    _test( interiorNodes.size() == femPtr_->InteriorNodes() );
    _test( femPtr_->InteriorNodes() == femData_.InteriorNodeCount() );
  }

  // .) FACE ELEMENT TYPES
  if ( verbose_ ) cout << "Testing face element types...\n";
  for( auto i = 0; i < femPtr_->Faces(); ++i )
    _test( femPtr_->ElementTypeOfFace( i ) == femData_.FaceElementType( i ) );

  // .) UNIT NORMAL
  if ( verbose_ ) cout << "Testing unit normal...\n";
  if( femPtr_->IsSurface() )
  {
    std::vector<double> unitNormal;
    std::vector<double> unitNormalTest( femData_.UnitNormal() );
    unitNormal = femPtr_->UnitNormal();
    for( auto i = 0; i < femData_.Dim(); ++i )
      _equal( unitNormal.at( i ), unitNormalTest.at( i ) ,femData_.Tolerance());
  }

  // .) COORDINATE MATRIX
  if ( verbose_ ) cout << "Testing coordinate matrix...\n";
  for( auto row = 0; row < femPtr_->Nodes(); ++row )
  {
    for( uint32_t dim = 0; dim < femPtr_->Dim(); ++dim )
    {
      _test( femPtr_->XYZ( row, dim ) == femData_.NodePtr( row )->operator []( dim ) );
    }
  }
  double cache( femPtr_->XYZ( 1,0 ) );
  femPtr_->XYZ( 1,0, 999.999 );
  _test( femPtr_->XYZ( 1,0 ) == 999.999 );
  _test( femPtr_->XYZ( 0,1 ) != 999.999 );
  femPtr_->XYZ( 1,0, cache );


  // .) EXTRAPOLATE INTEGRATION POINT TO NODE
  if ( verbose_ ) cout << "Testing gp to node extrapolation...\n";
  if( femData_.ExtrapolationVariableCount() != 0 )
  {
    vector<double> ivars( femPtr_->IntegrationPoints() * femData_.ExtrapolationVariableCount() );
    vector<double> nvars( femPtr_->Nodes() * femData_.ExtrapolationVariableCount() );
    int index( 0 );
    for( auto ip = 0; ip < femPtr_->IntegrationPoints(); ++ip )
    {
      for( auto v = 0; v < femData_.ExtrapolationVariableCount(); ++v )
        ivars.at( index++ ) = femData_.IntegrationPointVariable( ip, v );
    }
    femPtr_->ExtrapolateIntegrationPointVariableToNodes( femData_.ExtrapolationVariableCount(),
                                                         ivars, nvars );
    index = 0;
    for( auto np = 0; np < femPtr_->Nodes(); ++np )
    {
      for( auto v = 0; v < femData_.ExtrapolationVariableCount(); ++v )
        _equal( nvars.at( index++ ), femData_.NodePointVariable( np, v ), femData_.Tolerance() );
    }
  }

  // .) INTEGRATION POINT LOCAL TO GLOBAL
  if ( verbose_ ) cout << "Testing local-global transformation(integration points)...\n";
  if( femPtr_->IntegrationPoints() != 0 )
  {
    vector<double> ipGlobal;
    for( auto ip = 0; ip < femPtr_->IntegrationPoints(); ++ip )
    {
      femPtr_->IntegrationPoint( ip, ipGlobal );
      for( auto xyz = 0; xyz < femData_.Dim(); ++xyz )
        _equal( ipGlobal.at( xyz ), femData_.IpGlobal( ip, xyz ), femData_.Tolerance() );
    }
  }

  // .) INTEGRATION POINT WEIGHTS
  cout << "Testing integration point weights...\n";
  for( auto ip = 0; ip < femPtr_->IntegrationPoints(); ++ip )
    _equal( femPtr_->WeightAtIntegrationPoint( ip ), femData_.IntegrationPointWeight( ip ), femData_.Tolerance() );

  // .) SHAPE FUNCTION: N AT POINT
  if ( verbose_ ) cout << "Testing shape function N at point...\n";
  std::vector<double> shapeFunctionN;
  std::vector<double> shapeFunctionXYZ;
  for( auto xyz = 0; xyz < femPtr_->Dim(); ++xyz )
    shapeFunctionXYZ.push_back( femData_.ShapeFunctionXYZ( xyz ) );
  femPtr_->N( shapeFunctionN, shapeFunctionXYZ );
  for( auto n = 0; n < femPtr_->Nodes(); ++n )
    _equal( shapeFunctionN.at( n ), femData_.ShapeFunctionN( n ), femData_.Tolerance() );

  // .) SHAPE FUNCTION: N AT INTEGRATION POINT
  if ( verbose_ ) cout << "Testing shape function N at integration point...\n";
  for( auto ip = 0; ip < femPtr_->IntegrationPoints(); ++ip )
  {
    shapeFunctionN.clear();
    femPtr_->N_AtIntegrationPoint( ip, shapeFunctionN );
    for( auto n = 0; n < femPtr_->Nodes(); ++n )
      _equal( shapeFunctionN.at( n ), femData_.ShapeFunctionNatIP( ip, n ), femData_.Tolerance() );
  }

  // .) SHAPE FUNCTION: N AT BARY CENTER
  if ( verbose_ ) cout << "Testing shape function N at barycenter...\n";
  shapeFunctionN.clear();
  femPtr_->N_AtBaryCenter( shapeFunctionN );
  for( auto n = 0; n < femPtr_->Nodes(); ++n )
    _equal( shapeFunctionN.at( n ), femData_.ShapeFunctionNatBarycenter( n ), femData_.Tolerance() );

  /*
  // .) dN
  cout << "Testing shape function derivative...\n";
  DenseMatrix<DM_MIN> denseMatrix;
  femPtr_->dN( denseMatrix );
  for( size_t n = 0; n < femPtr_->Nodes(); ++n )
  {
    for( size_t xyz = 0; xyz < femPtr_->Dim(); ++xyz )
      _equal( denseMatrix( xyz, n ), femData_.ShapeFunctionDN( n, xyz ), femData_.Tolerance() );
  }
  */
  // .) dN AT XYZ
  if ( verbose_ ) cout << "Testing shape function derivative at point...\n";
  DenseMatrix<DM_MIN> denseMatrix;
  if( femData_.ShapeFunctionDNatXYZ() )
  {
    femPtr_->dN_At( denseMatrix, shapeFunctionXYZ );
    for( auto n = 0; n < femPtr_->Nodes(); ++n )
      {
        for( auto xyz = 0; xyz < femPtr_->Dim(); ++xyz )
          _equal( denseMatrix( xyz, n ), femData_.ShapeFunctionDNatXYZ( n, xyz ), femData_.Tolerance() );
      }
  }

  // .) dN AT IP
  if ( verbose_ ) cout << "Testing shape function derivative at integration point...\n";
  if( femData_.ShapeFunctionDNatIP() )
  {
    for( auto ip = 0; ip < femPtr_->IntegrationPoints(); ++ip )
    {
      femPtr_->dN_AtIntegrationPoint( denseMatrix, ip );
      for( auto n = 0; n < femPtr_->Nodes(); ++n )
        {
          for( auto xyz = 0; xyz < femPtr_->Dim(); ++xyz )
            _equal( denseMatrix( xyz, n ), femData_.ShapeFunctionDNatIP( ip, n, xyz ), femData_.Tolerance() );
        }
    }
  }

  // .) dN AT NODE
  if ( verbose_ ) cout << "Testing shape function derivative at nodes...\n";
  if( femData_.ShapeFunctionDNatNode() )
  {
    for( auto node = 0; node < femPtr_->Nodes(); ++node )
    {
      femPtr_->dN_AtNode( denseMatrix, node );
      for( auto n = 0; n < femPtr_->Nodes(); ++n )
      {
          for( auto xyz = 0; xyz < femPtr_->Dim(); ++xyz )
            _equal( denseMatrix( xyz, n ), femData_.ShapeFunctionDNatNode( node, n, xyz ), femData_.Tolerance() );
        }
    }
  }

  // .) dN AT BARY CENTER
  if ( verbose_ ) cout << "Testing shape function derivative at bary center...\n";
  if( femData_.ShapeFunctionDNatBaryCenter() )
  {
    femPtr_->dN_AtBarycenter( denseMatrix );
    for( auto n = 0; n < femPtr_->Nodes(); ++n )
    {
      for( auto xyz = 0; xyz < femPtr_->Dim(); ++xyz )
        _equal( denseMatrix( xyz, n ), femData_.ShapeFunctionDNatBaryCenter( n, xyz ), femData_.Tolerance() );
    }
  }

  // .) JACOBIAN AT IP
  if ( verbose_ ) cout << "Testing jacobian at integration points...\n";
  uint32_t dim_volume=1.0;

  if (femPtr_->IsSurface())
         dim_volume=2.0;
  if (femPtr_->IsVolume())
         dim_volume=3.0;

  if( femData_.JACOBIANatIP() )
  {
    for( auto ip = 0; ip < femPtr_->IntegrationPoints(); ++ip )
    {
      femPtr_->JacobianAtIntegrationPoint( ip );
      for( auto row = 0; row < dim_volume; ++row )
      {
        for( auto column = 0; column < femPtr_->Dim(); ++column )
          _equal( femPtr_->JAC( row, column ), femData_.JACOBIANatIP( ip, row, column ), femData_.Tolerance() );
      }
    }
  }


  // .) JACOBIAN AT RST
  if ( verbose_ ) cout << "Testing jacobian at rst point...\n";
  if( femData_.JACOBIANatRST() )
  {
    var_v.clear();
    for( auto rst = 0; rst < femPtr_->Dim(); ++rst )
      var_v.push_back( femData_.RST( rst ) );

    femPtr_->JacobianAt( var_v );

    for( auto column = 0; column < femPtr_->Dim(); ++column )
    {
      for( auto row = 0; row < dim_volume; ++row )
        _equal( femPtr_->JAC( row, column ), femData_.JACOBIANatRST( row, column ), femData_.Tolerance() );
    }
  }

  // .) TEST INTERPOLATION DERIVATIVES OVER A LINEAR SIMPLEX
  if ( femPtr_->IsSimplex() && femPtr_->Interpolation() == 1 ) {
    const auto iNrIps = femPtr_->IntegrationPoints();
    std::cerr << "Shape matrices for FE " << fileName_ << '\n';
    femPtr_->dN_AtBarycenter( denseMatrix );
    DenseMatrix<DM_MIN> denseMatrixIP;
    for (auto ip = 0; ip < iNrIps; ++ip) {
      femPtr_->dN_AtIntegrationPoint( denseMatrixIP, ip );
      _test( denseMatrix.Rows() == denseMatrixIP.Rows() );
      _test( denseMatrix.Cols() == denseMatrixIP.Cols() );
      for (auto i = 0; i < denseMatrix.Cols(); ++i) {
        for (auto j = 0; j < denseMatrix.Rows(); ++j) {
          _test( approximatelyEqual( denseMatrix(j,i), denseMatrixIP(j,i) ));
        }
      }
    }
  }

} // run

} // csmp
