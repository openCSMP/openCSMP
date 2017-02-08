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
template<size_t dim, template<size_t> class NodeOrElement>
double64 maximumDifference( Model<dim>& model, const char* nodeProp1, const char* nodeProp2, const char* region = "Model" );

/// returns the averaged difference found on all nodes/elements between two properties property
template<size_t dim, template<size_t> class NodeOrElement>
double64 averageDifference( Model<dim>& model, const char* nodeProp1, const char* nodeProp2, const char* region = "Model" );


/// returns the maximum value on any node/element of a property
template<size_t dim, template<size_t> class NodeOrElement>
double64 maximumOfProperty( Model<dim>& model, const char* prop, const char* region = "Model" );

/// conversion from float number to string
template<typename T>
void floatnumberToString( T number, std::string& recipient )
{
  std::ostringstream stringStream;
#ifdef __APPLE__
  if ( fabs(number) < 1.0e-30 ) number=0.0;
  else if ( number >1.0e+30 ) number=1.0e+30;
#endif
  stringStream << number;
  recipient = stringStream.str();
}

/// conversion from float number to string
template<typename T>
std::string floatnumberToString( T number )
{
  std::string cacheString;
  floatnumberToString<T>( number, cacheString );
  return cacheString;
}

/// conversion from number to string, no formatting specified
template<typename T>
void numberToString( T number, std::string& recipient )
{
  std::ostringstream stringStream;
  // SKM FIX
  stringStream.setf(std::ios::scientific);
  stringStream << number;
  recipient = stringStream.str();
}

// Specialisations for avoiding numbers with an exponent that
// has 3 digits (for PARAVIEW on Mac)
#ifdef __APPLE__
template<>
void numberToString<float>( float number, std::string& recipient )
{
  std::ostringstream stringStream;
  // SKM FIX
  stringStream.setf(std::ios::scientific);
//  stringStream.precision(std::numeric_limits<float>::digits10);
  if ( fabsf(number) < 1.0e-30f ) number=0.f;
  else if ( fabsf(number) > 1.0e+30f ) number=1.0e+30f;
  stringStream << number;
  recipient = stringStream.str();
}

template<>
void numberToString<double>( double number, std::string& recipient )
{
  std::ostringstream stringStream;
  // SKM FIX
  stringStream.setf(std::ios::scientific);
//  stringStream.precision(std::numeric_limits<double>::digits10);
  if ( fabs(number) < 1.0e-30 ) number=0.;
  else if ( fabs(number) > 1.0e+30 ) number=1.0e+30;
  stringStream << number;
  recipient = stringStream.str();
}

#endif

/// conversion from number to string, no formatting specified
template<typename T>
std::string numberToString( T number )
{
  std::string cacheString;
  numberToString<T>( number, cacheString );
  return cacheString;
}


template <typename T>
T stringToNumber ( const std::string &Text )
  {                               
    std::stringstream ss(Text);
    // SKM FIX
    ss.setf(std::ios::scientific);
//    ss.precision(std::numeric_limits<double>::digits10);

    T result;
    return ss >> result ? result : 0;
  }

/// returns arithmetic average of parent elements prop to node
template<size_t dim>
double64 elementToNodeProperty( Node<dim>* node, Index key )
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
  for( size_t i = 0; i < v.size(); ++i )
    std::cout << v.at( i ) << std::endl;
 }

/// extrapolation for volume relative element variables to nodal variable
template<size_t dim>
void extrapolateElementToNodalVariable( Model<dim>& mref,
                                        const char* region,
                                        const char* eprop,
                                        const char* nprop )
{
  Region<dim>& rref = mref.Region( region );
  PropertyHandle<dim> targetProp( mref, nprop );
  targetProp = 0.;
  const Index eKey( mref.Database().StorageKey( eprop ) );
  Index nKey( mref.Database().StorageKey( nprop ) );

  ScalarVariable elementContribution( PLAIN, 0. ), nodeValue( PLAIN, 0. ) ;
  const typename std::vector<Element<dim>* >::const_iterator elementsEnd( rref.ElementsEnd() );
  for( typename std::vector<Element<dim>* >::iterator it = rref.ElementsBegin(); it != elementsEnd; ++it )
  {
    elementContribution = (*it)->Read( eKey ) * (*it)->Volume() / (*it)->Nodes();

    for( typename std::vector<Node<dim>* > ::iterator iit = (*it)->NodesBegin(); iit != (*it)->NodesEnd(); ++iit )
    {
      (*iit)->Read( nKey, nodeValue );
      nodeValue += elementContribution;
      (*iit)->Store( nKey, nodeValue );
    }
  }

  std::cout << "extrapolateElementToNodalVariable<" << dim << "> extrapolated "
            << eprop << " to " << nprop << " in " << region << std::endl;

} // extrapolateElementToNodalVariable


/// discards given number of istream objects
void ignoreIstream( std::istream& iStream, size_t noOfTerms );


/// copies any given file
bool copyFile( std::string fileSource, std::string fileCopy );

/// orders node indices so that the follow the x position(distance from 0)
void bubbleSortNodes1D( std::vector<Node<1U>*>::iterator NODES_BEGIN,
                        std::vector<Node<1U>*>::iterator NODES_END );


/// adjusts the size of the model by the provided factor
void scaleModelByFactor( Model<3U>& model, double64 xFactor, double yFactor, double zFactor );

/// returns points which CONTAIN the lowest and largest x,y,z coordinates. not equal to MinMaxCoordinates() !!!
void minMaxXYZ( Point<3U>& min, Point<3U>& max, const Model<3U>& model, const char* regionName = "Model" );

/// returns the true center of a csmp model
Point<3U> modelMidpoint( const Model<3U>& model, const char* regionName = "Model" );


/// NCFVT bug workaround to calculate inflow across a boundary with PP conditions. Requires add. 'nodal volume flux'
template<size_t dim>
double boundaryInflow( Model<dim>& model, csmp::Boundary<dim>& boundary )
  {
    model.ExtrapolateElementToNodeProperty("volume flux", "nodal volume flux" );
    return  boundary.SurfaceIntegral( model.Database(), "nodal volume flux" );
  }


template<size_t dim>
double regionInflow( Model<dim>& model, csmp::Region<dim>& region )
  {
  model.ExtrapolateElementToNodeProperty("volume flux", "nodal volume flux" );
  return  region.VolumeIntegral( "nodal volume flux" ,false);
  }

/// establishes 'velocity' and 'volume flux' for a single phase system ('conductivity'), beware of surface elements
template<size_t dim>
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

    const typename std::vector<Element<dim>*>::const_iterator elementsEnd = rref.ElementsEnd();
    for ( typename std::vector<Element<dim>*>::iterator it = rref.ElementsBegin(); it != elementsEnd; ++it )
      {
        // vt = -k (lt grad p)
        const double64 conductivity = (*it)->Read( conductivityKey );
        velo = 0.;
        (*(*it)).dN_AtBaryCenter( DERIV, 1U );
        for ( size_t i = 0; i < (*it)->Nodes(); ++i )
          {
            double64 pf = (*it)->N(i)->Read( fluidPresssureKey );
            for( size_t xyz = 0; xyz < dim; ++xyz )
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


template<size_t dim,template<size_t> class Domain>
bool connected( const Domain<dim>& d1, const Domain<dim>& d2 )
{
  const typename std::vector<typename Domain<dim>::Simplex*>::const_iterator simplicesEnd = d1.ElementsEnd();
  for ( typename std::vector<typename Domain<dim>::Simplex*>::const_iterator it = d1.ElementsBegin(); it != simplicesEnd; ++it )
  {
    const typename std::vector<typename Domain<dim>::Simplex*>::const_iterator simplicesEndInner = d2.ElementsEnd();
    for ( typename std::vector<typename Domain<dim>::Simplex*>::const_iterator iit = d2.ElementsBegin(); iit != simplicesEndInner; ++iit )
      for( size_t i(0); i < (*iit)->Neighbors(); ++i )
        if( (*iit)->Neighbor(i) == (*it) )
          return true;
  } // simplices domain 1
  return false;
}


template<size_t dim,template<size_t> class Domain>
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
template<size_t dim,template<size_t> class Domain>
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
