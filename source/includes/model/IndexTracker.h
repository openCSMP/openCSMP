#ifndef INDEX_TRACKER_H
#define INDEX_TRACKER_H

#include <string>
#include <map>

namespace csmp {

  struct Index;

  /**
   @class IndexTracker
  
   @brief Keeps track of all Index objects once they have been registered. Serves PropertyDatabase to update all existing Index objects if necessary.
  
   @author  P. Lang
   @date  9/25/2012
   */
  class IndexTracker
    {
    public:
      IndexTracker();
      ~IndexTracker();

      std::map<csmp::Index*,std::string>::const_iterator IndicesBegin() const;
      std::map<csmp::Index*,std::string>::const_iterator IndicesEnd() const;

      /// Called when a new Index is issued
      void Attach( csmp::Index* newIndex, std::string parameterName );
      /// Called when a new Index is copied using an existing Index
      void Attach( csmp::Index* newIndex, const csmp::Index* existingIndex );
      /// Called when an existing Index is destroyed
      void Detach( const csmp::Index* existingIndex );
      /// Called when an existing Index is destroyed
      void Detach( std::string parameterName );
      /// Detaches itself from all indices, called upon destruction
      void DetachFromAll();
      
    protected:
      // We do not want to issue a copy of this. Instead, create a new instance and assign variables.
      IndexTracker( const IndexTracker& );

    private:
     std::map<csmp::Index*,std::string> trackedIndices_; ///< Pointers to all indices and corresponding parameter name
    };


  inline std::map<csmp::Index*,std::string>::const_iterator IndexTracker::IndicesBegin() const
    {
      return trackedIndices_.begin();
    }

      
  inline std::map<csmp::Index*,std::string>::const_iterator IndexTracker::IndicesEnd() const
    {
      return trackedIndices_.end();
    }


} // csmp

#endif