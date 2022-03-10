#ifndef FINITE_ELEMENT_TESTDATA_H
#define FINITE_ELEMENT_TESTDATA_H

#include "FiniteElement.h"
#include "Node.h"

namespace csmp {

/// Container class for FEM Test/Validation Data
class FiniteElement_TestData {
   friend std::istream& operator >> ( std::istream&, FiniteElement_TestData& );

public:
  FiniteElement_TestData( bool verbose ) : verbose_(verbose) {}
  ~FiniteElement_TestData()
    {
      for( std::vector<Node<3U>*>::iterator it = nodes_.begin(); it != nodes_.end(); ++it )
        delete (*it);
        nodes_.clear();
    }

  std::string FemName() const { return femName_; }
  void        FemName( std::string femName ) { femName_ = femName; }
  int         Dim() const { return dim_; }
  void        Dim( int dim ) {
                   dim_ = dim; shapeFunctionXYZ_.assign( dim, 999. );
                   jacobianAtRST_.assign( dim, std::vector<double>( Dim(), 999. ) );
                   RST_.assign( dim, 999. );
                 }
  size_t NodeCount() const { return nodes_.size(); }
  void   NodeCount( int nodeCount ) {
                                        DeleteNodes();
                                        nodes_.assign( nodeCount, nullptr );
                                        counterNodes_.assign( nodeCount, 999. );
                                        shapeFunctionN_.assign( nodeCount, 999. );
                                        shapeFunctionNatBarycenter_.assign( nodeCount, 999. );
                                        shapeFunctionDN_.assign( nodeCount, std::vector<double>( Dim(), 999. ) );
                                        shapeFunctionDNatXZY_.assign( nodeCount, std::vector<double>( Dim(), 999. ) );
                                        shapeFunctionDNatNode_.assign( nodeCount, std::vector<std::vector<double> >( nodeCount, std::vector<double>( Dim(), 999. ) ) );
                                        shapeFunctionDNatBaryCenter_.assign( nodeCount, std::vector<double>( Dim(), 999. ) );
                                     }
  int SegmentCount() const { return segmentCount_; }
  void   SegmentCount( int segmentCount ) { segmentCount_ = segmentCount; }
  size_t FaceCount() const { return nodesPerFace_.size(); }
  void   FaceCount( int faceCount ) {
                                       nodesPerFace_.assign( faceCount, 0 );
                                       nodesOfFace_.assign( faceCount, std::vector<int>() );
                                       faceElementTypes_.assign( faceCount, UNKNOWN );
                                     }
  int NeighborCount() const { return neighborCount_; }
  void   NeighborCount( int neighborCount ) { neighborCount_ = neighborCount; }
  int NodesPerFace( int i ) const { return nodesPerFace_.at( i ); }
  void   NodesPerFace( int i, int nodesPerFace ) { nodesPerFace_.at( i ) = nodesPerFace; }
  int NodeOfFace( int face, int faceNode ) const { return nodesOfFace_.at( face ).at( faceNode ); }
  void   NodeOfFace( int face, int faceNode , int femNode ) { nodesOfFace_.at( face ).push_back( femNode ); }
  int IntegrationPointNeighborCount() const { return constraintPointNeighborCount_; }
  void   IntegrationPointNeighborCount( int constraintPointNeighborCount ) { constraintPointNeighborCount_ = constraintPointNeighborCount; }
  int IntegrationPointCount() const { return integrationPointCount_; }
  void   IntegrationPointCount( int count ) {
                                                 integrationPointCount_ = count;
                                                 ipGlobal_.assign( count, std::vector<double>( Dim() ) );
                                                 integrationPointWeights_.assign( count, 999. );
                                                 shapeFunctionNatIP_.assign( count, std::vector<double>( NodeCount(), 999. ) );
                                                 shapeFunctionDNatIP_.assign( count, std::vector<std::vector<double> >( NodeCount(), std::vector<double>( Dim(), 999. ) ) );
                                                 jacobianAtIP_.assign( count, std::vector<std::vector<double> >( Dim(), std::vector<double>( Dim(), 999. ) ) );
                                               }
  int InterpolationNumber() const { return interpolationNumber_; }
  void   InterpolationNumber( int interpolationNumber ) { interpolationNumber_ = interpolationNumber; }
  bool   Isoparametric() const { return isoparametric_; }
  void   Isoparametric( bool isoparametric ) { isoparametric_ = isoparametric; }
  bool   UsesLocalCoordinates() const { return localCoordinates_; }
  void   UsesLocalCoordinates( bool localCoordinates ) { localCoordinates_ = localCoordinates; }
  int OrderOfShapeFunctions() const { return orderOfShapeFunctions_; }
  void   OrderOfShapeFunctions( int orderOfShapeFunctions ) { orderOfShapeFunctions_ = orderOfShapeFunctions; }
  bool   LineElement() const { return lineElement_; }
  void   LineElement( bool lineElement ) { lineElement_ = lineElement; }
  bool   SurfaceElement() const { return surfaceElement_; }
  void   SurfaceElement( bool surfaceElement ) { surfaceElement_ = surfaceElement; }
  bool   VolumeElement() const { return volumeElement_; }
  void   VolumeElement( bool volumeElement ) { volumeElement_ = volumeElement; }
  CSMP_FEM_TYPE   ElementType() const { return elementType_; }
  void            ElementType( CSMP_FEM_TYPE elementType ) { elementType_ = elementType; }
  CSMP_FEM_TYPE   FaceElementType( int face ) const { return faceElementTypes_.at( face ); }
  void            FaceElementType( int face, CSMP_FEM_TYPE faceElementType ) { faceElementTypes_.at( face ) = faceElementType; }
  Node<3U>*       NodePtr( int i ) const { return nodes_.at( i ); }
  void            NodePtr( int i, Node<3U>* node ) { nodes_.at( i ) = node; }
  double Volume() const { return volume_; }
  void   Volume( double volume ) { volume_ = volume; }
  double Tolerance() const { return tolerance_; }
  void   Tolerance( double tolerance ) { tolerance_ = tolerance; }
  double AspectRatio() const { return aspectRatio_; }
  void   AspectRatio( double aspectRatio ) { aspectRatio_ = aspectRatio; }
  double InnerRadius() const { return innerRadius_; }
  void   InnerRadius( double innerRadius ) { innerRadius_ = innerRadius; }
  size_t EdgeCount() const { return edgeLengths_.size(); }
  void   EdgeCount( int edgeCount ) { edgeLengths_.clear(); edgeLengths_.assign( edgeCount, 0 ); }
  double EdgeLength( int i ) const { return edgeLengths_.at( i ); }
  void   EdgeLength( int i, double length ) { edgeLengths_.at( i ) = length; }
  size_t NodesPerSegment() const { return nodesOfSegment_.at( 0 ).size(); }
  void   NodesPerSegment( int nodesPerSegment ) { nodesOfSegment_.assign( SegmentCount(), std::vector<int>( nodesPerSegment, 999 ) ); }
  int NodeOfSegment( int segment, int segmentNode ) const { return nodesOfSegment_.at( segment ).at( segmentNode ); }
  void   NodeOfSegment( int segment, int segmentNode , int femNode ) { nodesOfSegment_.at( segment ).at( segmentNode ) = femNode; }
  size_t CornerNodeCount() const { return cornerNodes_.size(); }
  void   CornerNodeCount( int count ) { cornerNodes_.assign( count, 999 ); }
  int NodeAtCorner( int corner ) const { return cornerNodes_.at( corner ); }
  void   NodeAtCorner( int corner, int femNode ) { cornerNodes_.at( corner ) = femNode; }
  size_t MidsideNodeCount() const { return midsideNodes_.size(); }
  void   MidsideNodeCount( int count ) { midsideNodes_.assign( count, 999 ); }
  int NodeAtMidside( int midside ) const { return midsideNodes_.at( midside ); }
  void   NodeAtMidside( int midside, int femNode ) { midsideNodes_.at( midside ) = femNode; }
  size_t InteriorNodeCount() const { return interiorNodes_.size(); }
  void   InteriorNodeCount( int count ) { interiorNodes_.assign( count, 999 ); }
  int NodeAtInterior( int interior ) const { return interiorNodes_.at( interior ); }
  void   NodeAtInterior( int interior, int femNode ) { interiorNodes_.at( interior ) = femNode; }
  void   UnitNormal( std::vector<double> unitNormal ) { unitNormal_ = unitNormal; }
  std::vector<double> UnitNormal() const { return unitNormal_; }
  void   ExtrapolationVariableCount( int count ) { ivars_.assign( IntegrationPointCount(), std::vector<double>( count )  ); nvars_.assign( NodeCount(), std::vector<double>( count ) ); }
  size_t ExtrapolationVariableCount() const { return ivars_.at( 0 ).size(); }
  double IntegrationPointVariable( int ip, int var ) const { return ivars_.at( ip ).at( var ); }
  void   IntegrationPointVariable( int ip, int var, double value )  {  ivars_.at( ip ).at( var ) = value; }
  double NodePointVariable( int np, int var ) const { return nvars_.at( np ).at( var ); }
  void   NodePointVariable( int np, int var, double value )  {  nvars_.at( np ).at( var ) = value; }
  double IpGlobal( int ip, int xyz ) const { return ipGlobal_.at( ip ).at( xyz ); }
  void   IpGlobal( int ip, int xyz, double value )  {  ipGlobal_.at( ip ).at( xyz ) = value; }
  double IntegrationPointWeight( int ip ) const { return integrationPointWeights_.at( ip ); }
  void   IntegrationPointWeight( int ip, double value ) { integrationPointWeights_.at( ip ) = value; }
  int CounterClockWiseNode( int n ) const { return counterNodes_.at( n ); }
  void   CounterClockWiseNode( int n, int node ) { counterNodes_.at( n ) = node; }
  double ShapeFunctionXYZ( int xyz ) const { return shapeFunctionXYZ_.at( xyz ); }
  void   ShapeFunctionXYZ( int xyz, double value ) { shapeFunctionXYZ_.at( xyz ) = value; }
  double ShapeFunctionN( int n ) const { return shapeFunctionN_.at( n ); }
  void   ShapeFunctionN( int n, double value ) { shapeFunctionN_.at( n ) = value; }
  double ShapeFunctionNatIP( int ip, int n ) const { return shapeFunctionNatIP_.at( ip ).at( n ); }
  void   ShapeFunctionNatIP( int ip,  int n, double value ) { shapeFunctionNatIP_.at( ip ).at( n ) = value; }
  double ShapeFunctionNatBarycenter( int n ) const { return shapeFunctionNatBarycenter_.at( n ); }
  void   ShapeFunctionNatBarycenter( int n, double value ) { shapeFunctionNatBarycenter_.at( n ) = value; }
  double ShapeFunctionDN( int n, int xyz ) const { return shapeFunctionDN_.at( n ).at( xyz ); }
  void   ShapeFunctionDN( int n,  int xyz, double value ) { shapeFunctionDN_.at( n ).at( xyz ) = value; }
  double ShapeFunctionDNatXYZ( int n, int xyz ) const { return shapeFunctionDNatXZY_.at( n ).at( xyz ); }
  void   ShapeFunctionDNatXYZ( int n, int xyz, double value ) { shapeFunctionDNatXZY_.at( n ).at( xyz ) = value; }
  void   ShapeFunctionDNatXYZ( bool val ) { shapeFunctionDNatXZYbool_ = val; }
  bool   ShapeFunctionDNatXYZ() const { return shapeFunctionDNatXZYbool_; }
  double ShapeFunctionDNatIP( int ip, int n, int xyz ) const { return shapeFunctionDNatIP_.at( ip ).at( n ).at( xyz ); }
  void   ShapeFunctionDNatIP( int ip, int n,  int xyz, double value ) { shapeFunctionDNatIP_.at( ip ).at( n ).at( xyz ) = value; }
  void   ShapeFunctionDNatIP( bool val ) { shapeFunctionDNatIPbool_ = val; }
  bool   ShapeFunctionDNatIP() const { return shapeFunctionDNatIPbool_; }
  double ShapeFunctionDNatNode( int node, int n, int xyz ) const { return shapeFunctionDNatNode_.at( node ).at( n ).at( xyz ); }
  void   ShapeFunctionDNatNode( int node, int n,  int xyz, double value ) { shapeFunctionDNatNode_.at( node ).at( n ).at( xyz ) = value; }
  void   ShapeFunctionDNatNode( bool val ) { shapeFunctionDNatNodeBool_ = val; }
  bool   ShapeFunctionDNatNode() const { return shapeFunctionDNatNodeBool_; }
  double ShapeFunctionDNatBaryCenter( int n, int xyz ) const { return shapeFunctionDNatBaryCenter_.at( n ).at( xyz ); }
  void   ShapeFunctionDNatBaryCenter( int n,  int xyz, double value ) { shapeFunctionDNatBaryCenter_.at( n ).at( xyz ) = value; }
  void   ShapeFunctionDNatBaryCenter( bool val ) { shapeFunctionDNatBaryCenterBool_ = val; }
  bool   ShapeFunctionDNatBaryCenter() const { return shapeFunctionDNatBaryCenterBool_; }
  double JACOBIANatIP( int ip, int row, int column ) const { return jacobianAtIP_.at( ip ).at( row ).at( column ); }
  void   JACOBIANatIP( int ip, int row, int column, double value ) { jacobianAtIP_.at( ip ).at( row ).at( column ) = value; }
  void   JACOBIANatIP( bool val ) { JACOBIANatIPbool_ = val; }
  bool   JACOBIANatIP() const { return JACOBIANatIPbool_; }
  void   JACOBIANatRST( bool val ) { JACOBIANatRSTbool_ = val; }
  bool   JACOBIANatRST() const { return JACOBIANatRSTbool_; }
  void   JACOBIANatRST( int row, int column, double value ) { jacobianAtRST_.at( row ).at( column ) = value; }
  double JACOBIANatRST( int row, int column ) const { return jacobianAtRST_.at( row ).at( column ); }
  void   RST( int rst, double value ) { RST_.at( rst ) = value; }
  double RST( int rst ) const { return RST_.at( rst ); }



private:
  std::string femName_;
  int dim_;
  int segmentCount_;
  int neighborCount_;
  int constraintPointNeighborCount_;
  int integrationPointCount_;
  int interpolationNumber_;
  int orderOfShapeFunctions_;
  bool isoparametric_;
  bool localCoordinates_;
  bool lineElement_;
  bool surfaceElement_;
  bool volumeElement_;
  bool shapeFunctionDNatXZYbool_;
  bool shapeFunctionDNatIPbool_;
  bool shapeFunctionDNatNodeBool_;
  bool shapeFunctionDNatBaryCenterBool_;
  bool JACOBIANatIPbool_;
  bool JACOBIANatRSTbool_;
  double tolerance_;
  double volume_;
  double aspectRatio_;
  double innerRadius_;
  CSMP_FEM_TYPE elementType_;
  std::vector<int> nodesPerFace_;
  std::vector<int> cornerNodes_;
  std::vector<int> midsideNodes_;
  std::vector<int> interiorNodes_;
  std::vector<int> counterNodes_;
  std::vector<Node<3U>*> nodes_;
  std::vector<double> edgeLengths_;
  std::vector<double> unitNormal_;
  std::vector<double> integrationPointWeights_;
  std::vector<double> shapeFunctionXYZ_;
  std::vector<double> shapeFunctionN_;
  std::vector<double> shapeFunctionNatBarycenter_;
  std::vector<double> RST_;
  std::vector<std::vector<double> > ivars_;
  std::vector<std::vector<double> > nvars_;
  std::vector<std::vector<double> > ipGlobal_;
  std::vector<std::vector<double> > shapeFunctionNatIP_;
  std::vector<std::vector<double> > shapeFunctionDN_;
  std::vector<std::vector<double> > shapeFunctionDNatXZY_;
  std::vector<std::vector<double> > shapeFunctionDNatBaryCenter_;
  std::vector<std::vector<double> > jacobianAtRST_;
  std::vector<std::vector<std::vector<double> > > shapeFunctionDNatIP_;
  std::vector<std::vector<std::vector<double> > > shapeFunctionDNatNode_;
  std::vector<std::vector<std::vector<double> > > jacobianAtIP_;
  std::vector<std::vector<int> > nodesOfSegment_;
  std::vector<std::vector<int> > nodesOfFace_;
  std::vector<CSMP_FEM_TYPE> faceElementTypes_;
  
  const bool verbose_;

  void   DeleteNodes() {
             for( auto it = nodes_.begin(); it != nodes_.end(); ++it )
               delete (*it);
             nodes_.clear();
           } // DeleteNodes
};


/**
@class FiniteElement_TestData
@author Philipp Lang
@date Dec 2010

@todo (3) Verbose for read in feedback

Designed to hold test data of a given Finite Element type for later use in the generic FiniteElement_Test.
Overloaded istream operator to facilitate input from text file.

@section dataFile Data File Structure

Delimeters between variable names and corresponding values are tabs,
names may cointain white spaces. The given sequence is crucial. The
following sample test file should be selx explanatory, please scroll
further down for detailed information on the individual blocks.

@subsection Sample Test File Content and Structure

@code

entire test file here

@endcode

@subsection Extrapolation from Integration Point to Node

Testing extrapolation from integration(gauss) points to nodes. We choose
to test a fix number of two variables per IP. The first data line porvides
the input values in a form such that [ variable1 @ IP1, variable 2 @ IP1,
variable1 @ IP2, variable2 @ IP2 ... etc. for each IP ]
Similarly, the resulting extrapolated node values are given such that
[ variable1 @ N1, variable 2 @ N1, variable1 @ N2, variable2 @ N2 ... etc. for each Node ]

@code
Integration Points To Node Extrapolation(2 variables)
Integration Point Values(2x integration points)			1. 2. 1. 2. 1. 2. 1. 2.
Extrapolated Node Values(2x nodes)                    1. 2. 1. 2. 1. 2. 1. 2.
@endcode

*/

} // csmp

#endif // FINITEELEMENT_TESTDATA_H
