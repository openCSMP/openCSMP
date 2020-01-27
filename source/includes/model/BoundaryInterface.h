#ifndef CSMP_BOUNDARY_INTERFACE_H
#define CSMP_BOUNDARY_INTERFACE_H

#include "CSMP_definitions.h"
#include "Box.h"

namespace csmp {

template<size_t> class Face;
template<size_t> class Element;
template<size_t> class Boundary;
template<size_t> class Region;
template<size_t> class VSet;
class FaceConstructionData;

/// finds the neighbors of dim-1 element, and their faces that connect to it; index records neighbor materials
template<size_t dim>
FaceConstructionData  higherDimensionalNeighbors( const csmp::Element<dim>&, const csmp::Index& );

/// finds inside neighbor of dim-1 element, and the face that connects to it; index records neighbor materials
template<size_t dim>
const csmp::Element<dim>* const  higherDimensionalNeighbor( const csmp::Element<dim>&, const csmp::Index&,
                                                            size_t& local_face_number_of_e, double64& material_ID  );
template<size_t dim>
bool  higherDimensionalNeighbors( const Element<dim>& , std::vector<Element<dim>*>& );

/**

@brief Policy of Model class for the management of internal and external
boundaries which are of a lower spatial dimensional than the model and consist of Face objects.

The BoundaryInterface manages the boundaries of the model, including 
their construction and removal.

@section implementation Implementation

By contrast with the RegionInterface whose interfaces for Region creation
rely on the existance of a master region that contains all elements of the model, no such 
region exists for the boundaries, but boundaries are typically constructed from lower-dimensional
regions tagged via the name string 'BOUNDARY' as future boundaries in the geomodelling process.

@note the BoundaryInterface is not responsible for the management of Face objects.
These are handled by the MeshManager. 
The only thing that the BoundaryInterface does is to allow the creation / modification of Boundaries
and to give access to them. 

@author S.K. Matthai
@author P. Lang
@date 2010, 2017

*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
class BoundaryInterface {
  public:
    BoundaryInterface();
    BoundaryInterface( const BoundaryInterface& ) = delete;
    virtual ~BoundaryInterface();
    
    friend class Boundary_Test;
    friend class BoundaryInterface_Test; ///< to gain access to protected member functions for testing

    csmp::Boundary<dim>&        Boundary( const std::string& bName );
    const csmp::Boundary<dim>&  Boundary( const std::string& bName ) const;
    bool                        ContainsBoundary( const std::string& bName ) const;
      
    typedef typename std::map<std::string,csmp::Boundary<dim> >::iterator         boundaryIterator;
    typedef typename std::map<std::string,csmp::Boundary<dim> >::const_iterator   boundaryConstIterator;

    boundaryIterator       BoundariesBegin();
    boundaryIterator       BoundariesEnd();
    boundaryConstIterator  BoundariesBegin() const;
    boundaryConstIterator  BoundariesEnd() const;
    boundaryIterator       Boundary( const csmp::Boundary<dim>& );
    size_t                 Boundaries() const;

    /// combines strings (including 'BOUNDARY') into a unique name of the boundary
    std::string CreateBoundaryNameFrom( const FaceConstructionData&, const std::vector<std::string>& region_names ) const;

    /// searches for a boundary that intersects the supplied higher-dimensional region(s); returns null string ('\0') if not found
    std::string  FindBoundaryName( const std::set<std::string>& intersected_regions ) const;

    /// searches for the boundaries whose name contains the supplied strings (region names etc.)
    size_t  FindBoundaryNames( const std::set<std::string>& intersected_regions,
                               std::set<std::string>& region_patches_found ) const;

    /// checks whether name contains the strings BOUNDARY or any of the predefined boundary names
    bool IsBoundaryName( const std::string& regionName ) const;
    
    /// checks whether the model contains the BOX_BOUNDARY=equivalent boundaries TOP, BOTTOM etc.
    bool BoxShaped() const;


    // -----------------------------------------------
    // Boundary creation, modification & removal
    // -----------------------------------------------
    
    /// creates uniquely named boundary patches, returning their names if successful; the patches are created from meshed surface inside of model which will be removed by default
    std::pair<std::set<std::string>,bool>  CreateInternalBoundaryFrom( const char* dimension_minus1_region, 
                                                                       bool remove_dim_minus1_region=true );

    /// insert Boundary<Face> between two equidimensional unique regions, first on inside by convention returns name
    std::pair<std::string,bool>  InsertBoundary( const char* region1, const char* region2, bool createRegionBetween = false );

    /// converts lower dimensional element regions surrounding the target region and containing strings like BOUNDARY in their name into a Boundary<Face> object
    bool InsertBoundary( BOX_BOUNDARY boxBoundary, const char* region = "Model" );

    /// tries to create Box Boundary objects surrounding 'Model' into TOP, BOTTOM, IRREGULAR if possible; updates BOX_BOUNDARY flags
    bool EstablishBoxBoundariesFromOrientation();
    
    /// using the assigned box boundary flags, tries to create corresponding Boundary objects
    bool EstablishBoxBoundariesFromFlags();
    
    /// Removes boundary with optional deletion of its faces by the MeshManager
    void RemoveBoundary( csmp::Boundary<dim>& boundary, bool deleteElements = false );

    // -----------------------------------------------
    //  input/output
    // -----------------------------------------------
    
    /// prints boundary names and other stats to screen 
    void BoundariesOut() const;

    /// method used in the storage of a model to binary file
    bool OutputAllBoundariesToBinary( const char* file_name ) const;
  
    /// reads boundaries stored in CSMP native binary file written by OutputAllBoundariesToBinary; it can also read only a subset of variables
    void InputAllBoundariesFromBinary( const char* file_name, const std::set<std::string>* subset_variables = nullptr );


 protected:
 
    /// creates boundary from already interconnected faces that also know their parent elements
    bool AddBoundary( const char* name,
                      typename std::vector<Face<dim>*>::iterator facesBegin,
                      typename std::vector<Face<dim>*>::iterator facesEnd,
                      BOX_BOUNDARY );

    /// creates new boundary from existing faces
    bool InsertBoundary( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                         const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                         const std::string& bName );

    /// inserts csmp::Boundary for box boundaries
    bool EstablishBoxBoundaries();
  
    /// creates edge Boundary objects for box-shaped model from side boundaries
    bool EstablishEdgeBoundariesOfBoxShapedModel();

      /// @todo DEPRECATE use boundary information from VSet to create lower-dimensional regions of edge elements
    bool EstablishEdgeRegionsOfBoxShapedModel( const VSet<dim>&  );

    /// inserts irregular csmp::Boundary for all eligible regions in the model
    bool EstablishBoundariesFromRegions( bool remove_original_lower_dimensional_regions );

    /// inserts irregular csmp::Boundary for all eligible regions in the discontiguous model
    bool EstablishBoundariesFromDiscontiguousModel(bool remove_original_lower_dimensional_regions);

    /// construct Boundary<Face> objects around a region
    bool AddFaces(const char* region);

    /// Splits boundary based on regions into new boundaries
    bool DivideBoundary( typename std::map<std::string, csmp::Boundary<dim> >::iterator boundary,
                         const typename std::map<std::string, Region<dim> >::const_iterator subRegion );

    /// Splits boundary into remainder and new name with parameter name
    bool DivideBoundary(  typename std::map<std::string, csmp::Boundary<dim> >::iterator boundary,
                          const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                          const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                          const std::string& bName );
 protected:
   std::map<std::string,csmp::Boundary<dim> >   faceBoundaryMap_; ///< storage of the boundaries
};

} // csmp


#endif
