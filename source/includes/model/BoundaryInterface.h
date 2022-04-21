#ifndef CSMP_BOUNDARY_INTERFACE_H
#define CSMP_BOUNDARY_INTERFACE_H

#include "CSMP_definitions.h"
#include "Box.h"

namespace csmp {

class ModelTopology;
template<uint32_t> class Face;
template<uint32_t> class Element;
template<uint32_t> class Boundary;
template<uint32_t> class Region;
template<uint32_t> class VSet;
class FaceConstructionData;

/// finds the neighbors of dim-1 element, and their faces that connect to it; index records neighbor materials
template<uint32_t dim>
FaceConstructionData  higherDimensionalNeighbors( const csmp::Element<dim>&, const csmp::Index& );

/// finds inside neighbor of dim-1 element, and the face that connects to it; index records neighbor materials
template<uint32_t dim>
const csmp::Element<dim>* const  higherDimensionalNeighbor( const csmp::Element<dim>&, const csmp::Index&,
                                                            size_t& local_face_number_of_e, double& material_ID  );
template<uint32_t dim>
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
template<uint32_t dim, template<uint32_t> class BOUNDARY_COMPLEX>
class BoundaryInterface {
  public:
    BoundaryInterface() {}   
    BoundaryInterface( const BoundaryInterface& ) = delete;
    ~BoundaryInterface() {}
    
    friend class Boundary_Test;
    friend class BoundaryInterface_Test; ///< to gain access to protected member functions for testing

    csmp::Boundary<dim>&        Boundary( const std::string& bName );
    const csmp::Boundary<dim>&  Boundary( const std::string& bName ) const;
    bool                        ContainsBoundary( const std::string& bName ) const;
      
    typedef typename std::map<std::string,csmp::Boundary<dim> >::iterator  boundaryIterator;
    typedef typename std::map<std::string,csmp::Boundary<dim> >::const_iterator  boundaryConstIterator;

    boundaryIterator       BoundariesBegin();
    boundaryIterator       BoundariesEnd();
    boundaryConstIterator  BoundariesBegin() const;
    boundaryConstIterator  BoundariesEnd() const;
    boundaryIterator       Boundary( const csmp::Boundary<dim>& );
    size_t                 Boundaries() const;

    /// creates a name for a new boundary that shall be created between two touching unique (non-overlapping) regions
    std::string CreateBoundaryName( const std::string& inside_region, const std::string& outside_region );

    /// searches for a boundary that intersects the supplied higher-dimensional region(s); returns null string ('\0') if not found
    std::string  FindBoundaryByName( const std::set<std::string>& intersected_regions ) const;

    /// searches for the boundaries whose name contains the supplied strings (region names etc.)
    size_t  FindBoundaryByNames( const std::set<std::string>& intersected_regions,
                                 std::set<std::string>& region_patches_found ) const;

    /// checks whether name contains the strings BOUNDARY or any of the predefined boundary names
    bool IsBoundaryName( const std::string& regionName ) const;
    
    /// checks whether the model contains the BOX_BOUNDARY=equivalent boundaries TOP, BOTTOM etc.
    bool BoxShaped() const;


    // -----------------------------------------------
    // Boundary creation, modification & removal
    // -----------------------------------------------
    
    /// uses the Face ids stored in the model topology object to form boundaries with corresponding names; returns number of boundaries formed
    size_t FormBoundariesFrom( const ModelTopology& );
    
    /// creates Faces and uniquely named boundary patches, returning their names if successful; the patches are created from meshed surface inside of model which will be removed by default
    std::pair<std::set<std::string>,bool>  CreateInternalBoundaryFrom( const char* dimension_minus1_region, 
                                                                       bool remove_dim_minus1_region=true );
                                                                       
    /// creates boundary on the outside of the model; no partitioning based on contacting regions will occur
    bool CreateExternalBoundaryFrom( const char* dimension_minus1_region, BOX_BOUNDARY boxBoundary );

    /// creates INTERNAL boundary, ignoring already existing boundaries or split-boundaries as well as lower-dimensional regions, region will be on inside
    bool CreateBoundaryAround( const char* region );

    /// insert lower-dimensional Region between two equidimensional unique regions, and then converts it into Boundary; returns boundary name
    std::pair<std::string,bool>  CreateBoundaryBetween( const char* region1, const char* region2 );
    
    /// Removes boundary with  deletion of its faces in the MeshManager
    void RemoveBoundary( const char* boundary, bool eerase_faces );
    
    /// Remove boundary, also deleting associated face objects if so requested
    void RemoveBoundary( csmp::Boundary<dim>&, bool erase_faces );


    // -----------------------------------------------
    //  input/output
    // -----------------------------------------------
    
    /// prints boundary names and other stats to screen 
    void BoundariesOut() const;

    /// method used in the storage of a model to binary file
    bool OutputBoundariesToBinary( const char* file_name ) const;
  
    /// reads boundaries stored in CSMP native binary file written by OutputBoundariesToBinary; if set non-empty will only read subset of variables
    void InputBoundariesFromBinary( const char* file_name, const std::set<std::string>& subset_variables );


 protected:
 
    /// combines strings (including 'BOUNDARY') into a unique name of the boundary
    std::string CreateBoundaryNameFrom( const FaceConstructionData&, const std::vector<std::string>& region_names ) const;

    /// creates boundary from already interconnected faces that also know their parent elements
    bool AddBoundary( const char* name,
                      typename std::vector<Face<dim>*>::iterator facesBegin,
                      typename std::vector<Face<dim>*>::iterator facesEnd,
                      BOX_BOUNDARY );

    /// creates csmp::Boundary objects replacing lower-dimensional BOX_BOUNDARY named regions with boundaries with the same names
    bool EstablishBoxBoundaries();
  
    /// inserts  box boundary or irregular csmp::Boundary objects for all eligible regions in the model; returns the names of the created boudaries
    std::set<std::string>  EstablishBoundariesFromRegions();

    /// tries to create Box Boundary objects surrounding 'Model' into TOP, BOTTOM, IRREGULAR if possible; updates BOX_BOUNDARY flags
    bool EstablishBoxBoundariesFromOrientation();
    
    /// creates BOX boundaries using the node flags to identify sides, edges, and corners; use for simple models where corresponding lines or surfaces are missing
    void EstablishBoxBoundariesFromNodeFlags( bool recreate_box_boundary_flags_before );

    // EDGES
    
    /// creates edge Boundary objects for box-shaped model from side boundaries
    bool EstablishEdgeBoundariesOfBoxShapedModel();

 protected:
   std::map<std::string,csmp::Boundary<dim> >  boundaryMap_; ///< storage of the boundaries
};

} // csmp


#endif
