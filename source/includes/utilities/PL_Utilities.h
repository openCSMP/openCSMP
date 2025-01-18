#ifndef PL_UTILITIES_H
#define PL_UTILITIES_H

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>

#include "Model.h"
#include "PropertyHandle.h"

namespace csmp{

/**
@file PL_Utilities.h
@author P. Lang
*/

/// P. Lang 2010
template<bool> struct CompileTimeError;
template<> struct CompileTimeError<true> {};


#define STATIC_ASSERT(expr) \
   ( CompileTimeError<(expr) != 0>() )

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// returns whether values are within tolerance
template<typename T>
inline bool withinTolerance( T value1, T value2, T tolerance )
{
  return( std::abs( value1 - value2 ) <= tolerance );
}

/// returns whether values are within tolerance
template<typename T>
inline bool smallerThanWithinTolerance( T valueInQuestion, T referenceValue, T tolerance )
{
  return( valueInQuestion < (referenceValue-tolerance) );
}

/// returns the maximum difference found on any node/element between two properties property
template<uint32_t dim, template<uint32_t> class NodeOrElement>
double maximumDifference( Model<dim>& model, const char* nodeProp1, const char* nodeProp2, const char* region = "Model" );

/// returns the averaged difference found on all nodes/elements between two properties property
template<uint32_t dim, template<uint32_t> class NodeOrElement>
double averageDifference( Model<dim>& model, const char* nodeProp1, const char* nodeProp2, const char* region = "Model" );


/// returns the maximum value on any node/element of a property
template<uint32_t dim, template<uint32_t> class NodeOrElement>
double maximumOfProperty( Model<dim>& model, const char* prop, const char* region = "Model" );

/// returns arithmetic average of parent elements prop to node
template<uint32_t dim>
double elementToNodeProperty( Node<dim>* node, Index key )
  {
    double cacheDouble = 0.;
    for( size_t parent = 0; parent < node->Parents(); ++parent )
      cacheDouble +=  node->Parent( parent )->Read( key );
    return cacheDouble / node->Parents();
  }


/// prints the content of a vector to std::cout, optional title std::string 'message'
template<typename T>
inline void vectorOut( const std::vector<T>& v, std::string message = "" )
 {
  std::cout << "\nvectorOut(" << message << "):\n";
  for( auto i = 0; i < v.size(); ++i )
    std::cout << v.at( i ) << std::endl;
 }

/// discards given number of istream objects
void ignoreIstream( std::istream& iStream, size_t noOfTerms );


/// copies any given file
bool copyFile( std::string fileSource, std::string fileCopy );

/// orders node indices so that the follow the x position(distance from 0)
void bubbleSortNodes1D( std::vector<Node<1U>*>::iterator NODES_BEGIN,
                        std::vector<Node<1U>*>::iterator NODES_END );


/// adjusts the size of the model by the provided factor
void scaleModelByFactor( Model<3U>& model, double xFactor, double yFactor, double zFactor );

/// returns points which CONTAIN the lowest and largest x,y,z coordinates. not equal to MinMaxCoordinates() !!!
void minMaxXYZ( Point<3U>& min, Point<3U>& max, const Model<3U>& model, const char* regionName = "Model" );

/// returns the true center of a csmp model
Point<3U> modelMidpoint( const Model<3U>& model, const char* regionName = "Model" );


/// NCFVT bug workaround to calculate inflow across a boundary with PP conditions. Requires add. 'nodal volume flux'
template<uint32_t dim>
double boundaryInflow( Model<dim>& model, csmp::Boundary<dim>& boundary )
  {
    model.ExtrapolateCellToNodeProperty("volume flux", "nodal volume flux" );
    return  boundary.SurfaceIntegral( model.Database(), "nodal volume flux" );
  }


template<uint32_t dim>
double regionInflow( Model<dim>& model, csmp::Region<dim>& region )
  {
  model.ExtrapolateCellToNodeProperty("volume flux", "nodal volume flux" );
  return  region.VolumeIntegral( "nodal volume flux" ,false);
  }

/// establishes 'velocity' and 'volume flux' for a single phase system ('conductivity'), beware of surface elements
template<uint32_t dim>
void singlePhaseVelocity( Model<dim>& model, const std::string& regionName, 
                         const std::string velocityName = "velocity", 
                         const std::string pressureName = "fluid pressure", 
                         const std::string conductivityName = "conductivity", 
                         const std::string volumeFluxName = "volume flux" )
  {
    VectorVariable<dim>   velo;
    ScalarVariable        flux;
    DenseMatrix<DM_MIN>   DERIV(dim,dim);
    Index velocityKey( model.Database().StorageKey( velocityName.c_str() ) );
    Index fluidPresssureKey( model.Database().StorageKey( pressureName.c_str() ) );
    Index conductivityKey( model.Database().StorageKey( conductivityName.c_str() ) );
    Index volumeFluxKey( model.Database().StorageKey( volumeFluxName.c_str() ) );

    Region<dim>&  rref(  model.Region( regionName.data() ) );

    const auto elementsEnd = rref.CellsEnd();
    for ( auto it = rref.CellsBegin(); it != elementsEnd; ++it )
      {
        // vt = -k (lt grad p)
        const double conductivity = (*it)->Read( conductivityKey );
        velo = 0.;
        (*(*it)).dN_AtBaryCenter( DERIV, 1U );
        for ( auto i = 0; i < (*it)->Nodes(); ++i )
          {
            double pf = (*it)->N(i)->Read( fluidPresssureKey );
            for( auto xyz = 0; xyz < dim; ++xyz )
              velo( xyz ) += pf * -DERIV( xyz, i ) * conductivity;
          }
        // storing the computed velocity
        (*it)->Store( velocityKey, velo );
        (*it)->Store(  volumeFluxKey , makeScalar( PLAIN, velo.Length() ) );
      } // elements region
  }


/// looks up value for given key and passes it to parameter; returns true if found, false otherwise
template<typename KeyType,typename ValueType>
bool mapFind( const KeyType& KEY, ValueType& VALUE, const std::map<KeyType,ValueType>& mapToWriteTo )
  {
    const typename std::map<KeyType,ValueType>::const_iterator IT( mapToWriteTo.find( KEY ) );
    if( IT == mapToWriteTo.end() )
      return false;
    VALUE = IT->second;
    return true;
  }


/// returns if value is to be found
template<typename KeyType,typename ValueType>
bool mapFind( const KeyType& KEY, const std::map<KeyType,ValueType>& mapToWriteTo )
  {
    const typename std::map<KeyType,ValueType>::const_iterator IT( mapToWriteTo.find( KEY ) );
    if( IT == mapToWriteTo.end() )
      return false;
    return true;
  }


/// Writes map to tab-seperated ascii file, (given types have '<<' overload) and returns success
template<typename KeyType,typename ValueType>
bool mapOut( const std::map<KeyType,ValueType>& mapToOutput, const char* filename )
{
  std::ofstream outputFile;
  outputFile.open(filename);
  if( !outputFile.is_open() )
    return false;
  for( typename std::map<KeyType,ValueType>::const_iterator it( mapToOutput.begin() ); it != mapToOutput.end(); ++it )
    outputFile << it->first << "\t" << it->second << std::endl;
  return true;
}


template<uint32_t dim,template<uint32_t> class Domain>
bool connected( const Domain<dim>& d1, const Domain<dim>& d2 )
{
  const typename std::vector<typename Domain<dim>::Simplex*>::const_iterator simplicesEnd = d1.CellsEnd();
  for ( typename std::vector<typename Domain<dim>::Simplex*>::const_iterator it = d1.CellsBegin(); it != simplicesEnd; ++it )
  {
    const typename std::vector<typename Domain<dim>::Simplex*>::const_iterator simplicesEndInner = d2.CellsEnd();
    for ( typename std::vector<typename Domain<dim>::Simplex*>::const_iterator iit = d2.CellsBegin(); iit != simplicesEndInner; ++iit )
      for( size_t i(0); i < (*iit)->Neighbors(); ++i )
        if( (*iit)->Neighbor(i) == (*it) )
          return true;
  } // simplices domain 1
  return false;
}


template<uint32_t dim,template<uint32_t> class Domain>
double averageDistanceBetween( const Domain<dim>& d1, const Domain<dim>& d2 )
  {
    double averageDistance( 0. );
    size_t distanceCount( 0 );

    const typename std::vector<Node<dim>*>::const_iterator nodesEnd = d1.NodesEnd();
    for ( typename std::vector<Node<dim>*>::const_iterator it = d1.NodesBegin(); it != nodesEnd; ++it )
      {
      const typename std::vector<Node<dim>*>::const_iterator nodesEndInner = d2.NodesEnd();
      for ( typename std::vector<Node<dim>*>::const_iterator iit = d2.NodesBegin(); iit != nodesEndInner; ++iit )
        {
          ++distanceCount;
          averageDistance += (*it)->Coordinate().DistanceTo( (*iit)->Coordinate() );
        } // nodes domain 2
      } // nodes domain 1

    return static_cast<double>( averageDistance / distanceCount );
  }


/// @todo (2-F) (2-C) Should return a vector instead
template<uint32_t dim,template<uint32_t> class Domain>
double minimumDistanceBetween( const Domain<dim>& d1, const Domain<dim>& d2 )
{
  double minumumDistance( std::numeric_limits<double>::max() );

  const typename std::vector<Node<dim>*>::const_iterator nodesEnd = d1.NodesEnd();
  for ( typename std::vector<Node<dim>*>::const_iterator it = d1.NodesBegin(); it != nodesEnd; ++it )
  {
    const typename std::vector<Node<dim>*>::const_iterator nodesEndInner = d2.NodesEnd();
    for ( typename std::vector<Node<dim>*>::const_iterator iit = d2.NodesBegin(); iit != nodesEndInner; ++iit )
      {
        double const currentDistance = (*it)->Coordinate().DistanceTo( (*iit)->Coordinate() );
        minumumDistance = currentDistance < minumumDistance ? currentDistance : minumumDistance;
      }
  } // nodes domain 1

  return minumumDistance;
}


 /**
 @}
 */

} // csmp

#endif // PL_UTILITIES_H
