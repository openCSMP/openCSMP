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

@brief Policy of the Model class for the management of its internal and external
boundaries which are of a lower spatial dimensional and consist of Face objects.

@author P. Lang
@author S.K. Matthai
@date 2010

The BoundaryInterface manages the boundaries of the model, including 
their construction and removal.

@section motivation Motivation

The BoundaryInterface has been designed as a policy of the model providing
access and the technology to work with boundaries in CSMP computations.

@section consequences Consequences

@section implementation Implementation

By contrast with the RegionInterface whose creational operations
rely on a Master Region that encompasses all elements, no such 
region exists for the boundaries, but they are constructed from 
pre-existing regions whose existance will always precede that of 
the boundaries.

@note the BoundaryInterface is not responsible for the management of the Face objects
which are handled by the MeshManager. 
By contrast, the BoundaryInterface only manages boundaries=vectors of pointers to the 
faces.

@section applications Applications

@todo (3) storage in binary models
@todo (3) Include Edges into the boundaries

*/
template<size_t dim, template<size_t> class BOUNDARY_COMPLEX>
class BoundaryInterface {
  public:
    BoundaryInterface();
    virtual ~BoundaryInterface();

    csmp::Boundary<dim>&        Boundary( const std::string& bName );
    const csmp::Boundary<dim>&  Boundary( const std::string& bName ) const;
    bool                        ContainsBoundary( const std::string& bName ) const;
  
    /// searches for a boundary that intersects the supplied higher-dimensional region(s); returns null string ('\0') if not found
    std::string  FindBoundaryName( const std::set<std::string>& intersected_regions ) const;

    /// searches for the boundaries whose name contains the supplied strings (region names etc.)
    size_t  FindBoundaryNames( const std::set<std::string>& intersected_regions,
                               std::set<std::string>& region_patches_found ) const;

    /// checks whether name contains the strings BOUNDARY or any of the predefined boundary names
    bool IsBoundaryName( const std::string& regionName ) const;
    
    /// checks whether the model contain the BOX_BOUNDARY=equivalent boundaries TOP, BOTTOM etc.
    bool BoxShaped() const;

    // ----------------------------------------------------------------
    // Boundaries access and manipulations with the list of Boundaries
    // ----------------------------------------------------------------

    typedef typename std::map<std::string,csmp::Boundary<dim> >::iterator         boundaryIterator;
    typedef typename std::map<std::string,csmp::Boundary<dim> >::const_iterator   boundaryConstIterator;

    boundaryIterator       BoundariesBegin();
    boundaryIterator       BoundariesEnd();
    boundaryConstIterator  BoundariesBegin() const;
    boundaryConstIterator  BoundariesEnd() const;
    boundaryIterator       Boundary( const csmp::Boundary<dim>& );
    size_t                 Boundaries() const;

    /// Removes boundary
    void RemoveBoundary( csmp::Boundary<dim>& boundary, bool deleteElements = false );

    // -----------------------------------------------
    // Boundary creation
    // -----------------------------------------------
    /// create boundary from already interconnected faces that also know their parent elements
    bool AddBoundary( const char* name,
                      typename std::vector<Face<dim>*>::iterator facesBegin,
                      typename std::vector<Face<dim>*>::iterator facesEnd,
                      BOX_BOUNDARY );
    
    /// returns number of created boundary segments = to region juxtapositions encountered
    size_t  CreateInternalBoundaryFrom( const char* dimension_minus1_region, bool remove_dim_minus1_region=true );

    /// construct Boundary<Face> from lower dimensional boundary flagged elements of given region
    bool InsertBoundary( BOX_BOUNDARY boxBoundary, const char* region = "Model" );

    /// insert Boundary<Face> between two equidimensional unique regions, first on the inside by convention
    bool InsertBoundary( const char* region1, const char* region2, bool createRegionBetween = false );

    /// Creates new boundary from existing faces
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

    /// partitions encompassing boundary 'Model' into TOP, BOTTOM, IRREGULAR if possible; updates BOX_BOUNDARY flags
    bool EstablishRegularities();
    bool EstablishRegularitiesForEclipse();

    /// construct Boundary<Face> objects around a region
    bool AddFaces(const char* region);

    /// Splits boundary based on regions into new boundaries
    bool DivideBoundary( typename std::map<std::string, csmp::Boundary<dim> >::iterator boundary,
                         const typename std::map<std::string, Region<dim> >::const_iterator subRegion );

    /// Splits boundary into remainder and new name with parameter name
    bool DivideBoundary(typename std::map<std::string, csmp::Boundary<dim> >::iterator boundary,
		const typename std::vector<Face<dim>*>::const_iterator facesBegin,
		const typename std::vector<Face<dim>*>::const_iterator facesEnd,
		const std::string& bName);

    // -----------------------------------------------
    // Binary input/output
    // -----------------------------------------------

    /// method used in the storage of a model to binary file
    bool OutputAllBoundariesToBinary( const char* file_name ) const;
  
    /// reads boundaries stored in CSMP native binary file written by OutputAllBoundariesToBinary; it can also read only a subset of variables
    void InputAllBoundariesFromBinary( const char* file_name, const std::set<std::string>* subset_variables = nullptr );

  protected:
    std::map<std::string,csmp::Boundary<dim> >   faceBoundaryMap_; ///< storage of the boundaries

  private:
    BoundaryInterface( const BoundaryInterface& );

 };

} // csmp


#endif
