#include "IndexTracker.h"
#include "Index.h"
#include "Exception.h"

using namespace std;

namespace csmp{

  IndexTracker::IndexTracker()
    : trackedIndices_()
    {
    }


  IndexTracker::~IndexTracker()
    {
    DetachFromAll();
    }


  void IndexTracker::DetachFromAll()
    {
    for( map<csmp::Index*,string>::const_iterator it( IndicesBegin() ); it != IndicesEnd(); ++it )
      if(it->first)
        it->first->Attach(NULL);
    trackedIndices_.clear();
    }


  void IndexTracker::Attach( csmp::Index* newIndex, std::string parameterName )
    {
    if(newIndex)
      trackedIndices_[newIndex] = parameterName;
    }


  void IndexTracker::Attach( csmp::Index* newIndex, const csmp::Index* existingIndex )
    {
    map<csmp::Index*,string>::const_iterator it = trackedIndices_.find( const_cast<csmp::Index*>(existingIndex) );
    if( it == trackedIndices_.end() ) {
         cerr <<"\n\texisting Index: ";
         if ( existingIndex != NULL ) cerr << (*existingIndex) << endl;
         cerr <<"\n\tnew Index: ";
         if ( newIndex != NULL ) cerr << (*newIndex) << endl;
         throw csmp::Exception( CSMP_ERROR, "IndexTracker::Attach:", "Existing Index not registered." );
      }
    string parameterName = it->second;
    Attach( newIndex, parameterName );
    }


  void IndexTracker::Detach( const csmp::Index* existingIndex )
    {
    map<csmp::Index*,string>::iterator it = trackedIndices_.find( const_cast<csmp::Index*>(existingIndex) );
    if( it == trackedIndices_.end() )
      throw csmp::Exception( CSMP_ERROR, "IndexTracker::Detach", "Index not registered" );
    trackedIndices_.erase(it);
    }


  void IndexTracker::Detach( string parameterName )
    {
    map<csmp::Index*,string>::iterator it( trackedIndices_.begin() );
    while( it != trackedIndices_.end() )
      {
      if( it->second == parameterName )
        {
        it->first->Attach(NULL);
        trackedIndices_.erase(it);
        it = trackedIndices_.begin();
        }
      else
        ++it;
      }
    }




  } // csmp
