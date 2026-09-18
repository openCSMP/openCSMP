// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IndexTracker.h"
#include "Index.h"
#include "Exception.h"
#include "ErrorHandler.h"

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



std::map<csmp::Index*,std::string>::const_iterator IndexTracker::IndicesBegin() const
  {
    return trackedIndices_.begin();
  }

    
std::map<csmp::Index*,std::string>::const_iterator IndexTracker::IndicesEnd() const
  {
    return trackedIndices_.end();
  }


void IndexTracker::DetachFromAll()
  {
  for( map<csmp::Index*,string>::const_iterator it( IndicesBegin() ); it != IndicesEnd(); ++it )
    if(it->first)
      it->first->Attach(nullptr);
  trackedIndices_.clear();
  }


void IndexTracker::Attach( csmp::Index* newIndex, std::string parameterName )
  {
    if ( newIndex != nullptr )
      trackedIndices_[newIndex] = parameterName;
  }


/**
    Attaches an existing index object to a new index.
    Reports if the existing index object is orphan, i.e. has no name because 
    it has not been tracked before.
    
    In this case, the parameter is called 'unspecified'.
    This will be picked up by the property database at the next update
    attempt.
*/
void IndexTracker::Attach( csmp::Index* newIndex, const csmp::Index* existingIndex )
  {
      string parameterName("unspecified");
    
      map<csmp::Index*,string>::const_iterator it = trackedIndices_.find( const_cast<csmp::Index*>(existingIndex) );
      if ( it == trackedIndices_.end() ) {
           cerr <<"\nIndexTracker::Attach (detach from existing instance to new Index instance):";
           cerr <<"\n\tpre-existing Index object is unregistered: ";
           if ( existingIndex != nullptr ) existingIndex->Out();
           cerr <<"\n\tnew Index: ";
           if ( newIndex != nullptr ) newIndex->Out();
           // error reporting
           ErrorHandler&  csmp_error( ErrorHandler::Instance() );
           cerr <<"\nregistered index objects:";
           Out();
           csmp_error.Note( ERROR, "IndexTracker::Attach:", "Existing Index object could not be re-attached as it was not registered." );
        }
      else parameterName = it->second;
      Attach( newIndex, parameterName );
  }



/**
     Detach index object as identified by a (unique) pointer to it.
*/
void IndexTracker::Detach( const csmp::Index* existingIndex )
  {
#ifdef DEBUG
     const bool verbose(true);
#else
     const bool verbose(false);
#endif
     map<csmp::Index*,string>::iterator it = trackedIndices_.find( const_cast<csmp::Index*>(existingIndex) );
    
     if ( it == trackedIndices_.end() ) {
          ErrorHandler&  csmp_error( ErrorHandler::Instance() );
          if ( verbose ) {
               cerr <<"\nIndexTracker::Detach:";
               cerr <<"\nunregistered index:\n";
               if ( existingIndex != nullptr ) existingIndex->Out();
               cerr <<"\nregistered index objects: "<< trackedIndices_.size();
               Out();
            }
          csmp_error.Note( WARNING, "IndexTracker::Detach:", "target Index was not registered." );
       }
     else trackedIndices_.erase(it);
  }



/**
    Detaches all index objects associated with a particular variable
    that is identified by its (unique) name.
*/
void IndexTracker::Detach( string parameterName )
  {
    map<csmp::Index*,string>::iterator it( trackedIndices_.begin() );
    while( it != trackedIndices_.end() )
      {
      if( it->second == parameterName )
        {
          it->first->Attach(nullptr);
          trackedIndices_.erase(it);
          it = trackedIndices_.begin();
        }
      else
        ++it;
      }
  }
  
  
  
  
void IndexTracker::Out() const
 {
     cout <<"\nIndexTracker::Out: currently tracked index objects and their values:";
     int counter(1);
     for ( auto it=trackedIndices_.begin(); it!=trackedIndices_.end(); ++it ) {
          cout <<"\n"<< counter++ <<": "<< *(*it).first <<" of variable '"<< (*it).second <<"'";
       }
     cout << endl;
 }
 
 } // end namespace csmp







