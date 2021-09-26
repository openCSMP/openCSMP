#ifndef SPLIT_BOUNDARY_INTERFACE_H
#define SPLIT_BOUNDARY_INTERFACE_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Boundary;
template<size_t> class SplitBoundary;

/**
    Creation, management and deletion of SplitBoundary objects.
 
    Split boundaries can be created for:
      1. internal model boundaries
      2. node-matched disconnected boundaries of the mesh
*/
template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
class SplitBoundaryInterface {
  public:
    SplitBoundaryInterface() {}
    SplitBoundaryInterface( const SplitBoundaryInterface& bd ) = delete;
    ~SplitBoundaryInterface() {}

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
    /// Creates SplitBoundaries detecting and connecting node-matched disconnected perimeter element faces in mesh (already created in ANSYS or other); these are grouped and named for regions
    std::pair<std::set<std::string>,bool>  DetectAndCreateSplitBoundaries();
    
    /// creation of one or multiple SplitBoundaries from a lower dimensional region 
    std::pair<std::set<std::string>,bool>  CreateSplitBoundaryFrom( const char* dim_1_region );

    /// one-to-one conversion of a model Boundary into a SplitBoundary, non-constant because Boundary gets removed
    std::pair<std::string,bool>  CreateSplitBoundaryFrom( Boundary<dim>& );

    /// Creation of SplitBoundary between regions via boundary that gets deleted afterwards
    std::pair<std::string,bool>  InsertSplitBoundary( const char* region1, const char* region2, bool createRegionBetween = false );

    /// inserts a lower-dimensional Region inside of the SplitBoundary, assigning its elements to the InterveningElement() pointers of its interfaces; the name will be that of the SplitBoundary followed by _REGION
    std::pair<std::string,bool>  InsertRegionIntoSplitBoundary( const char* split_boundary, int32 material_id_for_new_elements );

    /// Creates isolated lower-dimensional mesh regions between split boundaries with unique names matching those of the SplitBoundary objects; set will be empty if none created
    std::set<std::string>  InsertLowerDimensionalRegionsIntoSplitBoundaries( int32 material_id_for_new_elements );

    /// Creates a single lower-dimensional mesh region taking into account all split boundaries objects; returning its name and whether this operation was successful
    bool  SingleRegionFromAllSplitBoundaries( const char* name_of_new_region );

    /// Removes splitboundary including interfaces, but does not fuse the mesh back together again
    void RemoveSplitBoundary( const char* split_boundary );

    /// Removes splitboundary including interfaces, but does not fuse the mesh back together again
    void RemoveSplitBoundary( csmp::SplitBoundary<dim>& );

    // -----------------------------------------------------------
    // Input/output
    // The properties can also be read from a subset of variables
    // -----------------------------------------------------------
    
    /// prints current split boundaries to screen
    void SplitBoundariesOut() const;
    
    bool OutputSplitBoundariesToBinary( const char* fileName ) const;
    bool InputSplitBoundariesFromBinary( const char* fileName, const std::set<std::string>& subset_variables );

  protected:

    /// creating name for the case when the SplitBoundary was already present in the input mesh
    std::string CreateSplitBoundaryName( const std::pair<std::string, std::string>& juxtaposed_regions ) const;
  
    /// JCK: multiplicates the nodes of the boundaryin order to create split boundary TODO: method does not seem to be used!
    void SplitNodes( const Boundary<dim>& boundary, csmp::SplitBoundary<dim>& splitboundary );
    
    /// Since the construction of a new SplitBoundary may have affected existing ones, this method updates the connectivity of all SplitBoundary objects;
    void UpdateSplitBoundaryComplex();
 
  protected:
    std::map<std::string,csmp::SplitBoundary<dim> >  splitBoundaryMap_; ///< boundary name & boundary container of key-value pairs
};

} // csmp

#endif
