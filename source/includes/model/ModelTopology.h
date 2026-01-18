#ifndef CSMP_MODEL_TOPOLOGY_H
#define CSMP_MODEL_TOPOLOGY_H

#include "CSMP_definitions.h"
#include "ANSYS_ElementSpecifications.h"
#include "CSMP_ElementSpecifications.h"

namespace csmp {

template<uint32_t> class VSet;

/** 
     Model topology  associates elements with regions and boundaries,
     recording their names, element types, and their dimensionality.
     This is important in model construction process.
     The ModelTopology thus stores this "topologic" information, which
     cannot be stored in the VSet.
     
     In CSMP regions are model subdomains, just as boundaries and split boundaries. All of these can be stored in the ModelTopology.
     They are distinguished only by their names. Thus, a Boundary either contains the string BOUNDARY in its name or has one of the standard names
     for boundaries and or edges of box-shaped models (BOTTOM, TOP, LEFT, RIGHT, BACK, FRONT, IRREGULAR etc., EDGE1..EDGE12.
     Likewise, SplitBoundary objects contain the string SPLITBOUNDARY in their name, typically separated by underscores from other more specific
     parts of the name.
    
     @attention It follows that the user has to assign appropriate names to the geometrical entities in their model that they want to become boundaries.
  
     Model topology can also be used to perform a consistency check on the cells stored in the VSet.
     Thus, it can check the cell numbering and the neighbor connectivity.
     To some degree it can also fix such potential problems with the VSet.
    
     When the model is supposed box-shaped, this is tested and potentially missing information is restored in collaboration with the VSet / VData and Box.
     In this case it can deduce the correct BOX_BOUNDARY flags and assign them to the nodes.
     
     When a model is created by ANSYS and output using the CSMP interface of ANSYS
     the information needed to initialise the ModelToplogy class is contained in the '.asc' file.

    @section motivation Motivation

    ModelTopology was created to convey topological information from
    an ANSYS model to CSMP so that corresponding named groups of elements 
    can be created and assigned specific material properties.  

    For this purpose, ModelTopology stores named maps of the element
    IDs and types that make up specific model regions. By default, the 
    element names are the ANSYS element names as captured inside
    the class ANSYS_ElementSpecifications. Functionality is there
    in the public interface to convert these into CSMP names.

    Since one does not always want to use all topological information and
    or finite elements in a simulation, ModelTopology allows the user
    to reduce the initial topology as, for instance, obtained from ANSYS's
    meshing tools, to a few target regions - or element types, such as
    only the volume elements. Several interfaces are provided for this
    purpose.  

    @section structure Structure
     
    ModelTopology assumes a complementary role to the VSet which 
    stores the connectivity between elements and their nodes. Thus, there
    must always be a supporting VSet in order to build a Model.
    When the ModelTopology is reduced to a subset of the original 
    model, corresponding operations must be performed on the VSet.
     
    @author S.K. Matthaei
    @date 2001
    
 */
class ModelTopology {
  public:
    explicit ModelTopology( bool isoparametric_element_mesh=false );
    ModelTopology( const char* model_name, bool isoparametric_element_mesh );
 
    /// general model info
    void        ModelName( const char* name );
    std::string ModelName() const;
    size_t      ModelDomains() const;
    
    /// model shape and spatial dimension
    bool        LineModel() const;
    bool        SurfaceModel() const; // true if only surface and/or line elements
    bool        SolidModel() const;
    size_t      MinimumSpatialDimensionOfDomain( const char* region ) const;

    /// type of elements
    bool        IsoparametricFiniteElements() const;
    void        UseIsoparametricFiniteElementTypes();

    /// FEM interpolation function properties
    size_t      InterpolationOrder() const;
    bool        LinearElementMesh() const;
    bool        QuadraticElementMesh() const;
    
    /// return total number and CSMP names of FE-types in the model
    size_t      FiniteElementTypes( std::set<std::string>& etypes ) const;
    size_t      FiniteElementTypes( std::set<int32_t>& etypes ) const;
    
    // Since Regions consist of Elements, Boundaries of Faces, and SplitBoundaries of InterFaces, all these are just referred to as Cells
    
    /// changes the finite - element type associated with a particular region; see FiniteElement.h for available types
    void        ChangeCellType( const std::string& old_cell_type, std::string new_cell_type_same_as_enum_names );
    void        EliminateCellType( const char* etype );
    void        EliminateCellTypes( const std::list<std::string>& etypes );
    
    /// remove all curves
    void        EliminateLineCells();
    
    /// remove all surfaces
    void        EliminateSurfaceCells();
    
    /// remove all volumetric elements from a model
    void        EliminateVolumeCells();

    /// manipulation with elements within region
    size_t      Cells() const;
    void        Cells( std::set<size_t>& elmt_ids ) const;
    size_t      CellsWithinDomain( const char* region ) const;
    std::vector<size_t>::const_iterator  CellsOfDomainBegin( const char* region ) const;
    std::vector<size_t>::const_iterator  CellsOfDomainEnd( const char* region ) const;
    void        CellTypesOfDomain( const char* region, std::set<std::string>& etypes ) const;
    bool        IsWithinDomain( const char* region, size_t elmt_id ) const;

    /// operations on regions
    /// output subset as another model topology
    void        ExportSelectionTo( const std::list<std::string>& regions,
                                   ModelTopology& mt ) const;
       
    /// lists all the contained subdomains (regions, boundaries and splitboundaries
    void        OutputAll( std::list<std::string>& regions ) const;
    void        OutputRegions( std::list<std::string>& regions ) const;
    void        OutputBoundaries( std::list<std::string>& regions ) const;
    void        OutputSplitBoundaries( std::list<std::string>& regions ) const;

    /// check whether region is alreday included
    bool        Contains( const char* region ) const;
    /// remove certain regions
    void        Erase(); ///< all regions
    void        RemoveDomains( const std::set<std::string>& regions );
    void        RemoveDomain( const char* name );
    void        RemoveEmptyDomains();
    /// eliminate all model regions other than the ones specified in '*-regions.txt' file
    void        ReduceToDomains( const char* regions_file );
    void        ReduceToDomains( const std::set<std::string>& regions );
    /// add new regions
    bool        AddDomain( const char* rname,
                           const std::set<std::string>& fem_types,
                           const std::vector<size_t>& elms );
                           
    /// if the region only consists of a single element type, choose this method to add it (enums and strings are defined in FiniteElement.h)
    void        AddDomainCellType( const char* rname, const std::string& fem_types );
    
    /// add element types for regions that consist of multiple ones (enums and strings are defined in FiniteElement.h)
    void        AddDomainCellTypes( const char* rname, const std::set<std::string>& fem_types );
    
    /// element IDs range from 0..n-1
    void        AddDomainCellId( const char* rname, size_t elm );
    
    void        AddDomainCellIds( const char* rname, const std::vector<size_t>& elms );
                                     
    /// adds unique (non-overlapping) regions to topology, using the name, element types and element indices supplied
    bool        AddDomains( const std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >& unique_regions );
    
    /// adds unique (non-overlapping) regions to topology, using the region specifications (name and element types) and element index lists supplied
    bool        AddDomains( const std::multimap<std::string,std::string>& object_specs,
                            const std::multimap<std::string,std::vector<size_t> >& object_elements );
     
    /// adds regions to topology, eliminating lower-dimensional regions with the same name as the ones with the same dimension as the model
    bool        AddDomainsWithEquidimensionalCheck( const std::multimap<std::string,std::string>& object_specs,
                                                    const std::multimap<std::string,std::vector<size_t> >& object_elements );
    template<uint32_t dim>
    void        RemoveLowDimCellsFromDomains( csmp::VSet<dim>& vset );

	  /// region, boundary or split boundary names
	  void		    DomainNames( std::vector<std::string>& ) const;
   
    /// changes the name of a domain in the most efficient way; reports whether operation was successful
    bool        ChangeDomainName( const std::string& old_name, const std::string& new_name );
    
    /// iterate and modify regions
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator DomainsBegin() { return model_domains_.begin(); }
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator DomainsEnd() { return model_domains_.end(); }

    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator DomainsBegin() const { return model_domains_.begin(); }
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator DomainsEnd() const { return model_domains_.end(); }

    /// properties of regions
    void        PropertiesOfDomains( const char* regions_file,
                                     std::list<std::string>& properties,
                                     std::map<std::string,std::list<double> >& props ) const;


    // ---------------------------------------------------------
    // consistency checks and restoration of missing information
    // ---------------------------------------------------------

    /// initialises topology object and fixes pontentially wrong surface element orientations, non-consecutive numbering, orphan nodes, neighbor connectivity etc.
    template<uint32_t dim>
    bool        EstablishTopology( VSet<dim>& vset,
                                   const std::multimap<std::string,std::string>& object_specs,
                                   const std::multimap<std::string,std::vector<size_t> >& object_elements,
                                   bool require_unique_names_for_vol_surf_lines  = true,
                                   bool interactive_property_assignment = false,
                                   bool correct_orientation_of_surface_elements = false,
                                   bool reassign_boundary_flags = true );
  
    /// calls EstablishTopology with a reduced set of options
    template<uint32_t dim>
    bool        EstablishTopology( VSet<dim>& vset,
                                   bool require_unique_names_for_vol_surf_lines = true,
                                   bool correct_orientation_of_surface_elements = false,
                                   bool reassign_boundary_flags = true );
  
    /// checks that 2D model contains the boundaries LEFT, RIGHT, BOTTOM, TOP
    bool RectangleShapedModel() const;
    
    /// checks that 3D model contains the boundaries LEFT, RIGHT, BOTTOM, FRONT, BACK; TOP omitted because it may be IRREGULAR
    bool BoxShapedModel() const;
    
    /// recreates the BOX_BOUNDARY node flags if a problem was detected
    template<uint32_t dim> 
    bool  AssignBoxShapedModelFlags( VSet<dim>& ) const;

    template<uint32_t dim>
    void  AssignMaterialProperties( VSet<dim>&,const std::multimap<std::string,std::vector<size_t> >& object_elements );

    bool  CheckCellNumbering() const;

    /// numbering / repair
    void  CreateNewCellNumbers( std::map<size_t,size_t>& old_to_new_mapping, bool check_output=true );

    /// output info
    void  Out() const;
    void  Out( const char* output_file, const std::string header = std::string() ) const;

    /// initialises ModelTopology from text file as written by Out()
    bool InputFromTextFile( const char* file_dot_asc );

  private:  
    template<uint32_t dim>
    void  RenumberCells( VSet<dim>& vset, bool check_whether_already_correct );
    void  RenumberCells( const std::map<size_t /* old */,size_t /* new */>& eid_mapping );

    /// tests that the corner elements are indeed present
    bool Infer_BOX_BOUNDARY_EdgeAndCornerFlagsFromSideFlags( const VSet<2U>& ) const;
    
    /// deduces node boundary flags from BOX_BOUNDARY and other regions the name of which contains 'BOUNDARY'
    bool  FlagNodesUsingBoundaryDomains( VSet<2U>& vset ) const;
    bool  FlagNodesUsingBoundaryDomains( VSet<3U>& vset ) const;

    /// Boundary Flag nodes on lowerdimensional elements to INTERNAL - ATTENTION - THIS OVERWRITES, SO IT SHOULD BE DONE BEFORE BOXBOUNDARY FLAGS ARE SET
    //E.P Experimental 08.2022 - This is not used or tested, but may be needed to set TOPO Flags correctly, since they rely on INTERNAL flag being set.
    bool  FlagNodesOnLowerDimensionalElementsAsINTERNAL( VSet<2U>& vset ) const;
    bool  FlagNodesOnLowerDimensionalElementsAsINTERNAL( VSet<3U>& vset ) const;


  private:
    //       region name           etypes-of-region      ids of elements in region
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >  model_domains_;
    // public information on csmp element types
    typedef CSMP_ElementSpecifications  fem_specs;
    std::string  model_name_;
    bool  isoparametric_mesh_; // default is false
};


/// read regions from file
bool doesDomainsFileExist( const char* regions_file );
void readDesiredDomains( const char* regions_file, std::set<std::string>& desired_regions );

/// enlist the boundary names of box-shaped models
void BoundariesOfBoxShapedModel( std::set<std::string>& bs );

/// enlist the boundary names of rectangular models
void BoundariesOfRectangleShapedModel( std::set<std::string>& );

} // csmp


#endif 

