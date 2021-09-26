#include "PL_Utilities.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<>
double64 maximumDifference<3U,Node>( Model<3U>& model, const char* prop1Name, const char* prop2Name, const char* region )
{
  // keys to database for props
  const Index prop1Key( model.Database().StorageKey( prop1Name ) );
  const Index prop2Key( model.Database().StorageKey( prop2Name ) );

  // max difference
  double64 maxDifference( 0. ), difference( 0. );
  double64 prop1, prop2;

  // looping over all nodes of the region of interest
  Region<3>& rref( model.Region( region ) );
  const std::vector<Node<3U>*>::const_iterator nodesEnd( rref.NodesEnd() );
  for( std::vector<Node<3U>*>::iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
  {
    prop1 = (*it)->Read( prop1Key );
    prop2 = (*it)->Read( prop2Key );
    // calculating difference at current node
    difference = fabs( prop1 - prop2 );
    // and assign to maximum difference if larger
    maxDifference = maxDifference <  difference ? difference : maxDifference;

  } // nodes of region

  std::cout << "\n\nMaximumDifference: " << maxDifference << std::endl;
  cout.flush();
  return maxDifference;
}


template<>
double64 averageDifference<3U,Node>( Model<3U>& model, const char* prop1Name, const char* prop2Name, const char* region )
{
  // keys to database for props
  const Index prop1Key( model.Database().StorageKey( prop1Name ) );
  const Index prop2Key( model.Database().StorageKey( prop2Name ) );

  // max difference
  double64 cumulativeDifference( 0. ), averageDifference;
  double64 prop1, prop2;

  // looping over all nodes of the region of interest
  Region<3>& rref( model.Region( region ) );
  const std::vector<Node<3U>*>::const_iterator nodesEnd( rref.NodesEnd() );
  for( std::vector<Node<3U>*>::iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
  {
    prop1 = (*it)->Read( prop1Key );
    prop2 = (*it)->Read( prop2Key );
    // adding difference from current node
    cumulativeDifference += fabs( prop1 - prop2 );
  } // nodes of region

  averageDifference = cumulativeDifference / rref.Nodes();
  std::cout << "\n\nAverageDifference: " << averageDifference << std::endl;
  cout.flush();
  return averageDifference;
}


template<>
double64 maximumOfProperty<3U,Node>( Model<3U>& model, const char* propName, const char* region )
{
  // keys to database for props
  const Index propKey( model.Database().StorageKey( propName ) );

  // max
  double64 maxOfProp( 0. ), prop;

  // looping over all nodes of the region of interest
  Region<3>& rref( model.Region( region ) );
  const std::vector<Node<3U>*>::const_iterator nodesEnd( rref.NodesEnd() );
  for( std::vector<Node<3U>*>::iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
  {
    prop = (*it)->Read( propKey );
    maxOfProp = maxOfProp < prop ? prop : maxOfProp;
  } // nodes of region

  std::cout << "\nMax of " <<  model.Database().Name(propKey) << ": " << maxOfProp << std::endl;
  cout.flush();
  return maxOfProp;
}

/// discards given number of istream objects
void ignoreIstream( std::istream& iStream, size_t noOfTerms )
{
  static std::string cache;
  for( size_t i = 0; i < noOfTerms; ++i )
    iStream >> cache;
} // ignoreIstream


bool copyFile( string fileSource, string fileCopy )
{
   // file paths
   ifstream inputFile( fileSource.c_str(), ios::in|ios::binary );
   ofstream outputFile( fileCopy.c_str(), ios::out|ios::binary );

   // buffer
   inputFile.seekg( 0, ios::end );
   size_t fileSize = inputFile.tellg();

   // copy files
   if( inputFile.is_open() && outputFile.is_open() )
   {
      short* buffer = new short[fileSize/2];
      inputFile.seekg( 0, ios::beg );
      inputFile.read( reinterpret_cast<char*>( buffer ), fileSize );
      outputFile.write( reinterpret_cast<char*>( buffer ), fileSize );
      delete[] buffer;
   }
   else
   {
     throw csmp::Exception( ERROR,
                            "copyFile",
                            "Unable to open file",
                            "unknown" );
     return false;
   }

   inputFile.close();
   outputFile.close();

   return true;
}

// needs debugging for last/first nodes
void bubbleSortNodes1D( vector<Node<1U>*>::iterator NODES_END,
                        vector<Node<1U>*>::iterator NODES_BEGIN )
{  
  size_t cacheIdx( 9999 );
  const vector<Node<1U>*>::iterator LAST_NODE( NODES_END );

  for( vector<Node<1U>*>::iterator nodesUp = NODES_BEGIN; nodesUp < LAST_NODE; ++nodesUp )
  {
    for( vector<Node<1U>*>::iterator nodesDown = LAST_NODE; nodesDown > nodesUp; --nodesDown )
    {
      if( (*(nodesDown-1))->x() > (*nodesDown)->x() )
      {
        cacheIdx = (*nodesDown)->Idx();
        (*nodesDown)->Idx( (*(nodesDown-1))->Idx() );
        (*(nodesDown-1))->Idx( cacheIdx );
      } // swap
    } // inner, descending loop
  } // outer, ascending loop
} // bubbleSortNodes1D


void scaleModelByFactor( Model<3U>& model, double64 xFactor, double yFactor, double zFactor )
{
  Region<3>& mref( model.Region( "Model" ) );
  const vector<Node<3U>*>::const_iterator nodesEnd( mref.NodesEnd() );

  for( vector<Node<3U>*>::iterator it = mref.NodesBegin(); it != nodesEnd; ++it )
  {
    (*it)->x( (*it)->x() * xFactor );
    (*it)->y( (*it)->y() * yFactor );
    (*it)->z( (*it)->z() * zFactor );
  }

} // scaleModelByFactor


void minMaxXYZ( Point<3U>& min, Point<3U>& max, const Model<3U>& model, const char* regionName )
{
  double64 xMin( 0. ), xMax( 0. ), yMin( 0. ), yMax( 0. ), zMin( 0. ), zMax( 0. );
  const Region<3>& region( model.Region( regionName ) );

  const std::vector<csmp::Node<3U>*>::const_iterator regionNodesEnd( region.NodesEnd() );
  for( std::vector<csmp::Node<3U>*>::const_iterator it = region.NodesBegin(); it != regionNodesEnd; ++it )
  {
    xMax = (*it)->x() > xMax ? (*it)->x() : xMax;
    yMax = (*it)->y() > yMax ? (*it)->y() : yMax;
    zMax = (*it)->z() > zMax ? (*it)->z() : zMax;

    xMin = (*it)->x() < xMin ? (*it)->x() : xMin;
    yMin = (*it)->y() < yMin ? (*it)->y() : yMin;
    zMin = (*it)->z() < zMin ? (*it)->z() : zMin;
  }


  // storing values in point references
  min[0] = xMin; min[1] = yMin; min[2] = zMin;
  max[0] = xMax; max[1] = yMax; max[2] = zMax;
}


Point<3U> modelMidpoint( const Model<3U>& model, const char* regionName )
{
  Point<3U> midpoint, min, max;
  minMaxXYZ( min, max, model, regionName );
  midpoint[0] = ( min[0] + max[0] ) / 2;
  midpoint[1] = ( min[1] + max[1] ) / 2;
  midpoint[2] = ( min[2] + max[3] ) / 2;

  return midpoint;
}



} // csmp
