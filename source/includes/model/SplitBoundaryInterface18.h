#ifndef SPLIT_BOUNDARY_INTERFACE_18_H
#define SPLIT_BOUNDARY_INTERFACE_18_H

#include "SplitBoundary.h"

namespace csmp {

/**
    Creation, management and deletion of SplitBoundary objects.
 
    Split boundaries can be created for:
      1. internal model boundaries
      2. node-matched disconnected boundaries of the mesh
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
class SplitBoundaryInterface18 {
  public:
    SplitBoundaryInterface18();
    ~SplitBoundaryInterface18();

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
    // SplitBoundaries creation and deletion
    // -----------------------------------------------
    std::string  CreateSplitBoundaryName( const std::pair<std::string,std::string>& juxtaposed_regions ) const;

    /// checks whether name string contains SPLIT_BOUNDARY
    bool IsEligibleSplitBoundaryRegionName( const std::string& regionName ) const;
  
    /// Creates SplitBoundary detecting and connecting node-matched disconnected perimeter element faces in mesh; these are grouped and named for regions
    bool DetectAndCreateSplitBoundaries();

    /// Creation of SplitBoundary around region
    bool InsertSplitBoundary( const std::string& region );

    /// Creation of SplitBoundary between regions
    bool InsertSplitBoundary( const std::string& region1, const std::string& region2, bool createRegionBetween = false );

    /// Removes splitboundary
    void RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary, bool deleteElements = false );

    /// prints current split boundaries
    void SplitBoundariesOut() const;
  
    // -----------------------------------------------
    // Binary input/output
    // -----------------------------------------------

    bool OutputSplitBoundariesToBinary( const char* fileName ) const;
    bool InputSplitBoundariesFromBinary( const char* fileName );


  protected:
    std::map<std::string,csmp::SplitBoundary<dim> >  splitBoundaryMap_; ///< boundary name & boundary container of key-value pairs

  private:
    SplitBoundaryInterface18( const SplitBoundaryInterface18& bd );
};

} // csmp

#endif
