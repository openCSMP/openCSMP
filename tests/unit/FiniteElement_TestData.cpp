#include "FiniteElement_TestData.h"
#include "PL_Utilities.h"
#include "GlobalVerbose.h"
#include "Node.h"

#include <iomanip>
#include <sstream>
#include <limits>

namespace csmp {

/// PL 2010: overloaded istream operator for the FiniteElement_TestData class
std::istream& operator >> ( std::istream& iStream, FiniteElement_TestData& femTestData )
{
  const bool verbose(false);
  // dummies and delimeters, io feedback
  std::string cache;
  const std::streamsize dS( std::numeric_limits<std::streamsize>::max() );
  const char dcT( '\t' );
  bool var_b;
  using namespace std;
  bool global_verbose( GlobalVerbose::Instance().globalVerbose );
  if ( verbose ) cout << "\nReading FEM Test Data:\n";


  // reading fem name
  getline( iStream, femTestData.femName_, '\n' );
  getline( iStream, cache, '\n' );
  if ( verbose ) cout << "(" << femTestData.femName_ << ")" << endl;
  for( auto i = 0; i < femTestData.femName_.size()+2; ++ i )
    if ( verbose ) cout << "=";
  if ( verbose ) cout << endl;

  // reading verbose
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  global_verbose = var_b;

  // reading test tolerance
  iStream.ignore( dS, dcT );
  double tolerance;
  iStream >> tolerance;
  femTestData.Tolerance( tolerance );
  if ( verbose ) cout << "FEM Test Tolerance: " << tolerance << "(" << femTestData.Tolerance() << ")" << endl;

  // reading fem dim
  iStream.ignore( dS, dcT );
  uint32_t dim;
  iStream >> dim;
  femTestData.Dim( dim );
  if ( verbose ) cout << "FEM Spatial Dimension: " << dim << "(" << femTestData.Dim() << ")" << endl;

  // reading node count
  iStream.ignore( dS, dcT );
  size_t nodeCount;
  iStream >> nodeCount;
  femTestData.NodeCount( nodeCount );
  if ( verbose ) cout << "Node Count: " << nodeCount << "(" << femTestData.NodeCount() << ")"<< endl;

  // reading segments
  iStream.ignore( dS, dcT );
  size_t segmentCount;
  iStream >> segmentCount;
  femTestData.SegmentCount( segmentCount );
  if ( verbose ) cout << "Segment Count: " << segmentCount << "(" << femTestData.SegmentCount() << ")"<< endl;

  // reading faces
  iStream.ignore( dS, dcT );
  size_t faceCount;
  iStream >> faceCount;
  femTestData.FaceCount( faceCount );
  if ( verbose ) cout << "Face Count: " << faceCount << "(" << femTestData.FaceCount() << ")"<< endl;

  // reading neighbors
  iStream.ignore( dS, dcT );
  size_t neighborCount;
  iStream >> neighborCount;
  femTestData.NeighborCount( neighborCount );
  if ( verbose ) cout << "Neighbor Count: " << neighborCount << "(" << femTestData.NeighborCount() << ")"<< endl;

  // reading nodes per face
  size_t nodesPerFace;
  for( auto i = 0; i < femTestData.FaceCount(); ++i )
  {
    iStream.ignore( dS, dcT );
    iStream >> nodesPerFace;
    femTestData.NodesPerFace( i, nodesPerFace );
  }
  for( auto i = 0; i < femTestData.FaceCount(); ++i )
    if ( verbose ) cout << "Nodes Face " << i << ":" << femTestData.NodesPerFace( i ) << endl;

  // reading constraint point neighbors
  iStream.ignore( dS, dcT );
  size_t constraintPointNeighborCount;
  iStream >> constraintPointNeighborCount;
  femTestData.IntegrationPointNeighborCount( constraintPointNeighborCount );
  if ( verbose ) cout << "Constraint Point Neighbor Count: " << constraintPointNeighborCount << "(" << femTestData.IntegrationPointNeighborCount() << ")"<< endl;


  // reading integration points
  iStream.ignore( dS, dcT );
  size_t integrationPointCount;
  iStream >> integrationPointCount;
  femTestData.IntegrationPointCount( integrationPointCount );
  if ( verbose ) cout << "Integration Point Neighbor Count: " << integrationPointCount << "(" << femTestData.IntegrationPointCount() << ")"<< endl;

  // reading interpolation number
  iStream.ignore( dS, dcT );
  size_t interpolationNumber;
  iStream >> interpolationNumber;
  femTestData.InterpolationNumber( interpolationNumber );
  if ( verbose ) cout << "Interpolation Number: " << interpolationNumber << "(" << femTestData.InterpolationNumber() << ")"<< endl;

  // reading isoparametric
  iStream.ignore( dS, dcT );
  bool isoparametric;
  iStream >> isoparametric;
  femTestData.Isoparametric( isoparametric );
  if ( verbose ) cout << "Isoparametric: " << isoparametric << "(" << femTestData.Isoparametric() << ")"<< endl;

  // reading local/global coords
  iStream.ignore( dS, dcT );
  bool localCoordinates;
  iStream >> localCoordinates;
  femTestData.UsesLocalCoordinates( localCoordinates );
  if ( verbose ) cout << "Local Coordinates: " << localCoordinates << "(" << femTestData.UsesLocalCoordinates() << ")"<< endl;

  // reading order of shape functions
  iStream.ignore( dS, dcT );
  size_t orderOfShapeFunctions;
  iStream >> orderOfShapeFunctions;
  femTestData.OrderOfShapeFunctions( orderOfShapeFunctions );
  if ( verbose ) cout << "Order of shape functions: " << orderOfShapeFunctions << "(" << femTestData.OrderOfShapeFunctions() << ")"<< endl;

  // reading line element
  iStream.ignore( dS, dcT );
  bool lineElement;
  iStream >> lineElement;
  femTestData.LineElement( lineElement );
  if ( verbose ) cout << "Line Element: " << lineElement << "(" << femTestData.LineElement() << ")"<< endl;

  // reading surface element
  iStream.ignore( dS, dcT );
  bool surfaceElement;
  iStream >> surfaceElement;
  femTestData.SurfaceElement( surfaceElement );
  if ( verbose ) cout << "Surface Element: " << surfaceElement << "(" << femTestData.SurfaceElement() << ")"<< endl;

  // reading volume element
  iStream.ignore( dS, dcT );
  bool volumeElement;
  iStream >> volumeElement;
  femTestData.VolumeElement( volumeElement );
  if ( verbose ) cout << "Volume Element: " << volumeElement << "(" << femTestData.VolumeElement() << ")"<< endl;

  // reading element type
  iStream.ignore( dS, dcT );
  string elementType;
  iStream >> elementType;
  femTestData.ElementType( parseFiniteElementType( elementType ) );
  if ( verbose ) cout << "Element Type: " << elementType << "(" << femTestData.ElementType() << ")"<< endl;

  // reading nodes
  iStream.ignore( dS, dcT );
  //iStream >> cache >> cache;
  double x, y, z;
  for( auto i = 0; i < femTestData.NodeCount(); ++i )
  {
    iStream >> x; iStream >> y; iStream >> z;
    Node<3U> newNode;
    newNode.x( x ); newNode.y( y ); newNode.z( z );
    femTestData.NodePtr( i, &newNode );
    if ( verbose ) {
        cout << "Node(" << i << ") read in with: " << femTestData.NodePtr( i )->x() << "(" << x << ")" << " / ";
        cout                                       << femTestData.NodePtr( i )->y() << "(" << y << ")" << " / ";
        cout                                       << femTestData.NodePtr( i )->z() << "(" << z << ")" << endl;
      }
  }

  // reading fem volume
  iStream.ignore( dS, dcT );
  double femVolume;
  iStream >> femVolume;
  femTestData.Volume( femVolume );
  if ( verbose ) cout << "FEM Volume: " << femVolume << "(" << femTestData.Volume() << ")"<< endl;

  // reading aspect ratio
  iStream.ignore( dS, dcT );
  double aspectRatio;
  iStream >> aspectRatio;
  femTestData.AspectRatio( aspectRatio );
  if ( verbose ) cout << "FEM Aspect Ratio: " << aspectRatio << "(" << femTestData.AspectRatio() << ")"<< endl;

  // reading inner radius
  iStream.ignore( dS, dcT );
  double innerRadius;
  iStream >> innerRadius;
  femTestData.InnerRadius( innerRadius );
  if ( verbose ) cout << "FEM Inner Radius: " << innerRadius << "(" << femTestData.InnerRadius() << ")"<< endl;

  // reading edge count
  iStream.ignore( dS, dcT );
  size_t edgeCount;
  iStream >> edgeCount;
  femTestData.EdgeCount( edgeCount );
  if ( verbose ) cout << "FEM Edge Count: " << edgeCount << "(" << femTestData.EdgeCount() << ")"<< endl;

  // reading edge lengths
  //iStream >> cache >> cache;
  iStream.ignore( dS, dcT );
  double edgeLength;
  for( auto i = 0; i < femTestData.EdgeCount(); ++i )
  {
    iStream >> edgeLength;
    femTestData.EdgeLength( i, edgeLength );
    if ( verbose ) cout << "FEM Edge(" << i << ") Length: " << edgeLength << "(" << femTestData.EdgeLength( i ) << ")"<< endl;
  }

  // reading nodes per segment
  iStream.ignore( dS, dcT );
  size_t nodesPerSegment;
  iStream >> nodesPerSegment;
  femTestData.NodesPerSegment( nodesPerSegment );
  if ( verbose ) cout << "FEM Nodes per Segment: " << nodesPerSegment << "(" << femTestData.NodesPerSegment() << ")"<< endl;

  // reading segment nodes
  //iStream >> cache >> cache >> cache;
  iStream.ignore( dS, dcT );
  size_t femNodeOfSegment;
  for( auto i = 0; i < femTestData.SegmentCount(); ++i )
  {
    if ( verbose ) cout << "FEM Segment(" << i << ")" << endl;
    for( size_t j = 0; j < femTestData.NodesPerSegment(); ++j )
    {
      iStream >> femNodeOfSegment;
      femTestData.NodeOfSegment( i, j, femNodeOfSegment );
      if ( verbose ) cout << "Node(" << j << "): " << femNodeOfSegment << "(" << femTestData.NodeOfSegment( i, j ) << ") ";
    }
    if ( verbose ) cout << endl;
  }

  // reading face nodes
  //iStream >> cache >> cache >> cache;
  iStream.ignore( dS, dcT );
  size_t femNodeOfFace;
  for( auto i = 0; i < femTestData.FaceCount(); ++i )
  {
    if ( verbose ) cout << "FEM Face(" << i << ")" << endl;
    for( size_t j = 0; j < femTestData.NodesPerFace( i ); ++j )
    {
      iStream >> femNodeOfFace;
      femTestData.NodeOfFace( i, j, femNodeOfFace );
      if ( verbose ) cout << "Node(" << j << "): " << femNodeOfFace << "(" << femTestData.NodeOfFace( i, j ) << ") ";
    }
    if ( verbose ) cout << endl;
  }

  // reading corner nodes
  iStream.ignore( dS, dcT );
  size_t cornerNodeCount;
  iStream >> cornerNodeCount;
  femTestData.CornerNodeCount( cornerNodeCount );
  if ( verbose ) cout << "FEM Nodes at Corners: " << cornerNodeCount << "(" << femTestData.CornerNodeCount() << ")"<< endl;
  iStream >> cache >> cache >> cache;
  size_t cornerNode;
  for( auto i = 0; i < femTestData.CornerNodeCount(); ++i )
  {
    iStream >> cornerNode;
    femTestData.NodeAtCorner( i, cornerNode );
    if ( verbose ) cout << "Corner " << i << " Node " << cornerNode << "(" << femTestData.NodeAtCorner( i ) << ")" << endl;
  }

  // reading midside nodes
  iStream.ignore( dS, dcT );
  size_t midsideNodeCount;
  iStream >> midsideNodeCount;
  femTestData.MidsideNodeCount( midsideNodeCount );
  if ( verbose ) cout << "FEM Nodes at Interior: " << midsideNodeCount << "(" << femTestData.MidsideNodeCount() << ")"<< endl;
  iStream >> cache >> cache >> cache;
  size_t midsideNode;
  for( auto i = 0; i < femTestData.MidsideNodeCount(); ++i )
  {
    iStream >> midsideNode;
    femTestData.NodeAtMidside( i, midsideNode );
    if ( verbose ) cout << "Midside " << i << " Node " << midsideNode << "(" << femTestData.NodeAtMidside( i ) << ")" << endl;
  }

  // reading interior nodes
  iStream.ignore( dS, dcT );
  size_t interiorNodeCount;
  iStream >> interiorNodeCount;
  femTestData.InteriorNodeCount( interiorNodeCount );
  if ( verbose ) cout << "FEM Nodes at Interior: " << interiorNodeCount << "(" << femTestData.InteriorNodeCount() << ")"<< endl;
  iStream >> cache >> cache >> cache;
  size_t interiorNode;
  for( auto i = 0; i < femTestData.InteriorNodeCount(); ++i )
  {
    iStream >> interiorNode;
    femTestData.NodeAtInterior( i, interiorNode );
    if ( verbose ) cout << "Interior " << i << " Node " << interiorNode << "(" << femTestData.NodeAtInterior( i ) << ")" << endl;
  }

  // reading face element types
  //iStream >> cache >> cache >> cache;
  iStream.ignore( dS, dcT );
  string faceElementType;
  for( auto i = 0; i < femTestData.FaceCount(); ++i )
  {
    iStream >> faceElementType;
    femTestData.FaceElementType( i, parseFiniteElementType( faceElementType ) );
    if ( verbose ) cout << "Element Type of Face " << i << ": " << faceElementType << "(" << femTestData.FaceElementType( i ) << ")"<< endl;
  }

  // reading unit normal
  iStream >> cache >> cache;
  if( femTestData.SurfaceElement() )
  {
    if ( verbose ) cout << "Unit Normal: ";
    double comp;
    std::vector<double> unitNormal;
    for( auto i = 0; i < femTestData.Dim(); i++)
    {
      iStream >> comp;
      unitNormal.push_back( comp );
      if ( verbose ) cout << comp << "; ";
    }
    if ( verbose ) cout << endl;
    femTestData.UnitNormal( unitNormal );
  }

  // reading extrapolation validation data
  iStream.ignore( dS, dcT );
  size_t var_e;
  iStream >> var_e;
  femTestData.ExtrapolationVariableCount( var_e );
  if ( verbose ) cout<<"Number of variables to extrapolate: "<<var_e<<endl;
  //getline( iStream, cache );
  iStream.ignore( dS, dcT );
  double var;
  for( size_t ip = 0; ip < femTestData.IntegrationPointCount(); ++ip )
  {
    for( size_t v = 0; v < var_e; ++v )
    {
      iStream >> var;
      femTestData.IntegrationPointVariable( ip, v, var );
      if ( verbose ) cout << "Integration point " << ip << " variable " << v << ": " << var
                          << "(" << femTestData.IntegrationPointVariable( ip, v ) << ")\n";
    }
  }
  iStream.ignore( dS, dcT );
  for( size_t np = 0; np < femTestData.NodeCount(); ++np )
  {
    for( size_t v = 0; v < var_e; ++v )
    {
      iStream >> var;
      femTestData.NodePointVariable( np, v, var );
      if ( verbose ) cout << "Node point " << np << " variable " << v << ": " << var
                          << "(" << femTestData.NodePointVariable( np, v ) << ")\n";
    }
  }

  // reading global ip coordinates
  iStream.ignore( dS, dcT );
  for( size_t ip = 0; ip < femTestData.IntegrationPointCount(); ++ip )
  {
    if ( verbose ) cout << "Integration point " << ip << " global coordinates: ";
    for( size_t xyz = 0; xyz < femTestData.Dim(); ++xyz )
    {
      iStream >> var;
      femTestData.IpGlobal( ip, xyz, var );
      if ( verbose ) cout << var << "(" << femTestData.IpGlobal( ip, xyz ) << ")  ";
    }
    if ( verbose ) cout << endl;
  }

  // reading ip weights
  iStream.ignore( dS, dcT );
  for( size_t ip = 0; ip < femTestData.IntegrationPointCount(); ++ip )
  {
    iStream >> var;
    femTestData.IntegrationPointWeight( ip, var );
    if ( verbose ) cout << "Integration point " << ip << " weight: " << var << "("
                        << femTestData.IntegrationPointWeight( ip ) <<  ")\n";
  }

  // reading counter clock wise nodes
  iStream.ignore( dS, dcT );
  size_t var_t;
  cout << "Counter clock wise nodes: ";
  for( size_t n = 0; n < femTestData.NodeCount(); ++n )
  {
    iStream >> var_t;
    femTestData.CounterClockWiseNode( n, var_t );
    if ( verbose ) cout << var_t << "(" << femTestData.CounterClockWiseNode( n ) << ") ";
  }
  if ( verbose ) cout << endl;

  // reading shape function point
  //iStream >> cache >> cache >> cache >> cache >> cache >> cache;
  iStream.ignore( dS, dcT );
  for( size_t xyz = 0; xyz < femTestData.Dim(); ++ xyz )
  {
    iStream >> var;
    femTestData.ShapeFunctionXYZ( xyz, var );
    if ( verbose ) cout << "Shape function eval point " << xyz << ": " << var << "("
                        << femTestData.ShapeFunctionXYZ( xyz ) <<  ")\n";
  }

  // reading shape function N
  iStream.ignore( dS, dcT );
  for( size_t n = 0; n < femTestData.NodeCount(); ++ n )
  {
    iStream >> var;
    femTestData.ShapeFunctionN( n, var );
    if ( verbose ) cout << "Shape function N(" << n << "): " << var << "("
                        << femTestData.ShapeFunctionN( n ) <<  ")\n";
  }

  // reading shape functions at IP's
  for( size_t ip = 0; ip < femTestData.IntegrationPointCount(); ++ip )
  {
    iStream.ignore( dS, dcT );
    for( size_t n = 0; n < femTestData.NodeCount(); ++n )
    {
      iStream >> var;
      femTestData.ShapeFunctionNatIP( ip, n, var );
      if ( verbose ) cout << "Shape function N at(" << ip << "," << n << "): " << var << "("
                          << femTestData.ShapeFunctionNatIP( ip, n ) <<  ")\n";
    }
  }

  // reading shape function N at Barycenter
  iStream.ignore( dS, dcT );
  for( size_t n = 0; n < femTestData.NodeCount(); ++ n )
  {
    iStream >> var;
    femTestData.ShapeFunctionNatBarycenter( n, var );
    if ( verbose ) cout << "Shape function N(" << n << ") at Barycenter: " << var << "("
                        << femTestData.ShapeFunctionNatBarycenter( n ) <<  ")\n";
  }
  /*
  // reading shape functions derivative
  iStream.ignore( dS, dcT );
  for( size_t xyz = 0; xyz < femTestData.Dim(); ++xyz )
  {
    for( size_t n = 0; n < femTestData.NodeCount(); ++n )
    {
      iStream >> var;
      femTestData.ShapeFunctionDN( n, xyz, var );
      cout << "Shape function dN at(" << n << "," << xyz << "): " << var << "("
           << femTestData.ShapeFunctionDN( n, xyz ) <<  ")\n";
    }
  }
  */

  // reading shape functions derivative at XYZ
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  femTestData.ShapeFunctionDNatXYZ( var_b );
  if( femTestData.ShapeFunctionDNatXYZ() )
  {
    for( size_t xyz = 0; xyz < femTestData.Dim(); ++xyz )
    {
      for( size_t n = 0; n < femTestData.NodeCount(); ++n )
      {
        iStream >> var;
        femTestData.ShapeFunctionDNatXYZ( n, xyz, var );
        if ( verbose ) cout << "Shape function dN at Point(" << n << "," << xyz << "): " << var << "("
                            << femTestData.ShapeFunctionDNatXYZ( n, xyz ) <<  ")\n";
      }
    }
  }

  // reading shape functions derivative at IP
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  femTestData.ShapeFunctionDNatIP( var_b );
  if( femTestData.ShapeFunctionDNatIP() )
  {
    for( size_t ip = 0; ip < femTestData.IntegrationPointCount(); ++ip )
    {
      iStream.ignore( dS, dcT );
      for( size_t xyz = 0; xyz < femTestData.Dim(); ++xyz )
      {
        for( size_t n = 0; n < femTestData.NodeCount(); ++n )
        {
          iStream >> var;
          femTestData.ShapeFunctionDNatIP( ip, n, xyz, var );
          if ( verbose ) cout << "Shape function dN at Integration Point " << ip << "(" << n << "," << xyz << "): " << var << "("
                              << femTestData.ShapeFunctionDNatIP( ip, n, xyz ) <<  ")\n";
        }
      }
    }
  }

  // reading shape functions derivative at Nodes
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  femTestData.ShapeFunctionDNatNode( var_b );
  if( femTestData.ShapeFunctionDNatNode() )
  {
    for( size_t node = 0; node < femTestData.NodeCount(); ++node )
    {
      iStream.ignore( dS, dcT );
      for( size_t xyz = 0; xyz < femTestData.Dim(); ++xyz )
      {
        for( size_t n = 0; n < femTestData.NodeCount(); ++n )
        {
          iStream >> var;
          femTestData.ShapeFunctionDNatNode( node, n, xyz, var );
          if ( verbose ) cout << "Shape function dN at Node " << node << "(" << n << "," << xyz << "): " << var << "("
                              << femTestData.ShapeFunctionDNatNode( node, n, xyz ) <<  ")\n";
        }
      }
    }
  }

  // reading shape functions derivative at bary center
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  femTestData.ShapeFunctionDNatBaryCenter( var_b );
  if( femTestData.ShapeFunctionDNatBaryCenter() )
  {
      iStream.ignore( dS, dcT );
    for( size_t xyz = 0; xyz < femTestData.Dim(); ++xyz )
    {
      for( size_t n = 0; n < femTestData.NodeCount(); ++n )
      {
        iStream >> var;
        femTestData.ShapeFunctionDNatBaryCenter( n, xyz, var );
        if ( verbose ) cout << "Shape function dN at bary center(" << n << "," << xyz << "): " << var << "("
                            << femTestData.ShapeFunctionDNatBaryCenter( n, xyz ) <<  ")\n";
      }
    }
  }

  // reading jacobians at integration points
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  femTestData.JACOBIANatIP( var_b );

  uint32_t dim_volume=1;

  if( femTestData.SurfaceElement() )
      dim_volume=2;
  if( femTestData.VolumeElement() )
      dim_volume=3;

  if( femTestData.JACOBIANatIP() )
  {
    for( size_t ip = 0; ip < femTestData.IntegrationPointCount(); ++ip )
    {
      iStream.ignore( dS, dcT );
      for( size_t row = 0; row < dim_volume; ++row )
      {
        for( size_t column = 0; column < femTestData.Dim(); ++column )
        {
          iStream >> var;
          femTestData.JACOBIANatIP( ip, row, column, var );
          if ( verbose ) cout << "Jacobian at integration point " << ip << "(" << column << "," << row << "): " << var << "("
                              << femTestData.JACOBIANatIP( ip, row, column ) <<  ")\n";
        } // columns
      } // rows
    } // integration points
  } // if

  // reading jacobian at RST
  iStream.ignore( dS, dcT );
  iStream >> var_b;
  femTestData.JACOBIANatRST( var_b );
  if( femTestData.JACOBIANatRST() )
  {iStream.ignore( dS, dcT );
      // reading rst coordinates
      for( size_t rst = 0; rst < femTestData.Dim(); ++rst )
      {

          iStream >> var;
          femTestData.RST( rst, var );
          if ( verbose ) cout << "RST (" << rst << "): " << var << " (" << femTestData.RST( rst ) << ")\n";
      }


      // reading jacobian
      for( size_t row = 0; row < dim_volume; ++row )
      {
          for( size_t column = 0; column < femTestData.Dim(); ++column )
          {
              iStream >> var;
              femTestData.JACOBIANatRST( row, column, var );
             if ( verbose ) cout << "Jacobian at RST(" << row << "," << column << "): " << var << "("
                                 << femTestData.JACOBIANatRST( row, column ) <<  ")\n";
          }
      }
}



  // io feedback done
  for( auto i = 0; i < femTestData.femName_.size()+2; ++ i ) cout << "=";
    if ( verbose ) cout << endl << endl;

  return iStream;

} // std::istream& operator >>

} // csmp
