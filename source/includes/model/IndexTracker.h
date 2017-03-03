#ifndef INDEX_TRACKER_H
#define INDEX_TRACKER_H

#include <string>
#include <map>

namespace csmp {

  struct Index;

  /**
   @class IndexTracker
  
   @brief Keeps track of all Index objects once they have been registered. 
   
   Serves PropertyDatabase to keep all existing Index objects 
   objects up to date, which is necessary when new variables are created 
   at runtime.
   
   @note index tracking is necessary because some classes have indices as members.
   These objects can get invalidated when new variables of the same kind and place.
   are created. If so, for instance, the 5th scalar on the node might become
   the 10th. 
   
   By tracking all Index objects that exist in a model, these can be updated
   by the PropertyDatabase everytime a variable is created or deleted.
  
   @author  P. Lang
   @date  9/25/2012
   
   TODO: @todo SKM 2/2017 - The IndexTracker -> update methodology seems to fail
   if ArrayVariable objects of different length are created at runtime; fix
   or rethink approach.
   
*/
class IndexTracker {
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
      
      /// Detaches itself from all Index objects; this method is called upon destruction
      void DetachFromAll();
  
      /// print all the currently registered Index objects
      void Out() const;
      
    protected:
      // We do not want to issue a copy of this. Instead, create a new instance and assign variables.
      IndexTracker( const IndexTracker& );

    private:
     std::map<csmp::Index*,std::string> trackedIndices_; ///< Pointers to all indices and corresponding parameter name
 };

} // csmp

#endif
