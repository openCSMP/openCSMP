#ifndef SPLIT_BOUNDARY_INTERFACE_H
#define SPLIT_BOUNDARY_INTERFACE_H

#include "SplitBoundary.h"

namespace csmp {

/**
    Creation, management and deletion of SplitBoundary objects.
 
    Split boundaries can be created for:
      1. internal model boundaries
      2. node-matched disconnected boundaries of the mesh
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
class SplitBoundaryInterface {
  public:
    SplitBoundaryInterface();
    ~SplitBoundaryInterface();

    // -----------------------------------------------
    // Access of SplitBoundaries objects
    // -----------------------------------------------
    csmp::SplitBoundary<dim>&         SplitBoundary( const std::string& spbName );
    const csmp::SplitBoundary<dim>&   SplitBoundary( const std::string& spbName ) const;
    bool                              ContainsSplitBoundary( const std::string& bname ) const;

    typedef typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator         splitBoundaryIterator;
    typedef typename std::map<std::string,csmp::SplitBoundary<dim> >::const_iterator   splitBoundaryConstIterator;

    splitBoundaryIterator        SplitBoundariesBegin();
    splitBoundaryIterator        SplitBoundariesEnd();
    splitBoundaryConstIterator   SplitBoundariesBegin() const;
    splitBoundaryConstIterator   SplitBoundariesEnd() const;
    size_t                       SplitBoundaries() const;

    // -----------------------------------------------
    // SplitBoundary creation and deletion
    // -----------------------------------------------
    std::string  CreateSplitBoundaryName( const std::pair<std::string,std::string>& juxtaposed_regions ) const;
        
    /// Creates SplitBoundary detecting and connecting node-matched disconnected perimeter element faces in mesh; these are grouped and named for regions
    bool DetectAndCreateSplitBoundaries();

//  TODO: a separate method is needed for the one-to-one conversion of a Boundary into a SplitBoundary
    bool CreateFrom( const Boundary<dim>& );

    /// Creation of SplitBoundary around region
    bool InsertSplitBoundary( const std::string& region );

    /// Creation of SplitBoundary between regions
    bool InsertSplitBoundary( const std::string& region1, const std::string& region2, bool createRegionBetween = false );

    /// Creates lower-dimensional regions that lie in between the mesh patches that are separated by the split boundaries; creates unique names indicating region juxtaposition
    bool RegionsFromSplitBoundaries( std::set<std::string>& newly_created_regions );

    /// Removes splitboundary
    void RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary, bool deleteElements = false );

    /// prints current split boundaries
    void SplitBoundariesOut() const;
  
    // -----------------------------------------------------------
    // Binary input/output
    // The properties can also be read from a subset of variables
    // -----------------------------------------------------------

    bool OutputSplitBoundariesToBinary( const char* fileName ) const;
    bool InputSplitBoundariesFromBinary( const char* fileName, const std::set<std::string>* subset_variables = nullptr );


  protected:
    std::map<std::string,csmp::SplitBoundary<dim> >  splitBoundaryMap_; ///< boundary name & boundary container of key-value pairs

  private:
    SplitBoundaryInterface( const SplitBoundaryInterface& bd );
};

} // csmp

#endif
