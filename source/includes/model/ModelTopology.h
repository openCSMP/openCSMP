#ifndef CSMP_MODEL_TOPOLOGY_H
#define CSMP_MODEL_TOPOLOGY_H

#include "CSMP_definitions.h"
#include "ANSYS_ElementSpecifications.h"
#include "CSMP_ElementSpecifications.h"

namespace csmp {

template<size_t> class VSet;

/** 
     Model topology  associates elements with regions and (future) boundaries,
     recording the region and boundary names, element types, and their dimensionality.
     This is important in model construction process.
     The ModelTopology stores this "topologic" information, which
     cannot be stored in the VSet.
     
     Model topology is used also to perform a consistency check on the VSet. It checks the neighbor connectivity and gets VSet to fix it if there is a problem.
    
     When the model is supposed box-shaped, this is tested and potentially missing information is restored in collaboration with the VSet / VDataand Box.
     
     When the model was created by ANSYS and output using its CSP interface,
     the information needed to initialise the ModelToplogy class is contained in the '.asc' file.

    @author S.K. Matthaei
    @date 2001
     

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
     
     
    @section participants Participants
     
    The current implementation depends on the ANSYS_ElementSpecifications 
    object for the interpretation of finite element names.  
     
     
    @section examples Application Examples
     
    The topology class is used inside of ANSYS_Model3D  and  2D.
 */
class ModelTopology {
  public:
    explicit ModelTopology( bool isoparametric_element_mesh=false );
    ModelTopology( const char* model_name, bool isoparametric_element_mesh );
    ModelTopology( const ModelTopology& mt );
    ModelTopology& operator=( const ModelTopology& mt );
    ~ModelTopology();

    /// general model info
    void        ModelName( const char* name );
    std::string ModelName() const;
    size_t      ModelRegions() const;
    
    /// model shape and spatial dimension
    bool        LineModel() const;
    bool        SurfaceModel() const; // true if only surface and/or line elements
    bool        SolidModel() const;
    size_t      MinimumSpatialDimensionOfRegion( const char* region ) const;

    /// type of elements
    bool        IsoparametricElements() const;
    void        TreatElementsAsIsoparametric();

    /// FEM interpolation function properties
    size_t      InterpolationOrder() const;
    bool        LinearElementMesh() const;
    bool        QuadraticElementMesh() const;
    
    /// return total number and CSMP names of FE-types in the model
    size_t      FiniteElementTypes( std::set<std::string>& etypes ) const;
    size_t      FiniteElementTypes( std::set<int32>& etypes ) const;
    void        ChangeElementType( const std::string& old_element_type, std::string new_element_type );
    void        EliminateElementType( const char* etype );
    void        EliminateElementTypes( const std::list<std::string>& etypes );
    void        EliminateLineElements();
    void        EliminateSurfaceElements();
    void        EliminateVolumeElements();

    /// manipulation with elements within region
    size_t      Elements() const;
    void        Elements( std::set<size_t>& elmt_ids ) const;
    size_t      ElementsOfRegion( const char* region ) const;
    std::vector<size_t>::const_iterator  ElementsOfRegionBegin( const char* region ) const;
    std::vector<size_t>::const_iterator  ElementsOfRegionEnd( const char* region ) const;
    void        ElementTypesOfRegion( const char* region, std::set<std::string>& etypes ) const;
    bool        IsWithinRegion( const char* region, size_t elmt_id ) const;

    /// operations on regions
    /// output subset as another model topology
    void        ExportSelectionTo( const std::list<std::string>& regions,
                                   ModelTopology& mt ) const;
    void        Out( std::list<std::string>& regions ) const;
    /// check whether region is alreday included
    bool        Contains( const char* region ) const;
    /// remove certain regions
    void        Erase(); ///< all regions
    void        RemoveRegions( const std::set<std::string>& regions );
    void        RemoveRegion( const char* name );
    /// eliminate all model regions other than the ones specified in '*-regions.txt' file
    void        ReduceToRegions( const char* regions_file );
    void        ReduceToRegions( const std::set<std::string>& regions );
    /// add new regions
    bool        AddRegion( const char* rname,
                           const std::set<std::string>& fem_types,
                           const std::vector<size_t>& elms );
    void        AddRegionElementType( const char* rname,
                                      const std::string& fem_types );
    void        AddRegionElementTypes( const char* rname,
                                      const std::set<std::string>& fem_types );
    void        AddRegionElementId( const char* rname,
                                    size_t elm );
    void        AddRegionElementIds( const char* rname,
                                     const std::vector<size_t>& elms );
    bool        AddRegions( const std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >& unique_regions );
    
    bool        AddRegionsWithoutEquidimensionalCheck( const std::multimap<std::string,std::string>& object_specs,
                                                       const std::multimap<std::string,std::vector<size_t> >& object_elements );
                                                       
    bool        AddRegionsWithEquidimensionalCheck( const std::multimap<std::string,std::string>& object_specs,
                                                    const std::multimap<std::string,std::vector<size_t> >& object_elements );
    template<size_t dim>
    void        RemoveLowDimElementsFromRegions( csmp::VSet<dim>& vset );

	  /// region names
	  void		    RegionNames( std::vector<std::string>& ) const;

    /// properties of regions
    void        PropertiesOfRegions( const char* regions_file,
                                     std::list<std::string>& properties,
                                     std::map<std::string,std::list<double64> >& props ) const;


    // ---------------------------------------------------------
    // consistency checks and restoration of missing information
    // ---------------------------------------------------------

    /// checks and fixes pontentially wrong surface element orientations, non-consecutive numbering, orphan nodes, neighbor connectivity etc.
    template<size_t dim>
    bool        CheckTopology( VSet<dim>& vset,
                               const std::multimap<std::string,std::string>& object_specs,
                               const std::multimap<std::string,std::vector<size_t> >& object_elements,
                               bool require_unique_names_for_vol_surf_lines  = true,
                               bool interactive_property_assignment = false,
                               bool correct_orientation_of_surface_elements = false,
                               bool non_box_boundary = true );
  
    /// calls CheckTopology with a reduced set of options
    template<size_t dim>
    bool        CheckTopology( VSet<dim>& vset,
                               bool require_unique_names_for_vol_surf_lines = true,
                               bool correct_orientation_of_surface_elements = false,
                               bool non_box_boundary = true );
  
    /// checks that 2D model contains the boundaries LEFT, RIGHT, BOTTOM, TOP
    bool RectangleShapedModel() const;
    
    /// checks that 3D model contains the boundaries LEFT, RIGHT, BOTTOM, FRONT, BACK; TOP omitted because it may be IRREGULAR
    bool BoxShapedModel() const;
    
    /// recreates the BOX_BOUNDARY node flags if a problem was detected
    template<size_t dim> 
    bool  AssignBoxShapedModelFlags( VSet<dim>& ) const;

    template<size_t dim>
    void  AssignMaterialProperties( VSet<dim>&,const std::multimap<std::string,std::vector<size_t> >& object_elements);

    bool  CheckElementNumbering() const;

    /// numbering / repair
    void  CreateNewElementNumbers( std::map<size_t,size_t>& old_to_new_mapping, bool check_output=true );

    /// output info
    void  Out() const;
    void  Out( const char* output_file ) const;


  private:
  
    template<size_t dim>
    void  RenumberElements( VSet<dim>& vset, bool check_whether_already_correct );
    void  RenumberElements( const std::map<size_t /* old */,size_t /* new */>& eid_mapping );

    /// tests that the corner elements are indeed present
    bool Infer_BOX_BOUNDARY_EdgeAndCornerFlagsFromSideFlags( const VSet<2U>& ) const;
    
    /// deduces node boundary flags from BOX_BOUNDARY and other regions the name of which contains 'BOUNDARY'
    bool  FlagNodesUsingBoundaryRegions( VSet<2U>& vset ) const;
    bool  FlagNodesUsingBoundaryRegions( VSet<3U>& vset ) const;


  private:
    //       region name          etypes-of-region       ids of elements in region
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >  model_regions;
    // public information on csmp element types
    typedef CSMP_ElementSpecifications  fem_specs;
    std::string  model_name;
    bool  isoparametric_mesh; // default is false
};


/// read regions from file
bool doesRegionsFileExist( const char* regions_file );
void readDesiredRegions( const char* regions_file, std::set<std::string>& desired_regions );

/// box-shaped models
void BoundariesOfBoxShapedModel( std::set<std::string>& bs );
void BoundariesOfRectangleShapedModel( std::set<std::string>& );

/// surface elements orientation ( works soo far only for linear elements )
void CorrectSurfaceElementOrientations( VSet<3U>& vset );
void CorrectSurfaceElementOrientations( VSet<2U>& vset );
void CorrectSurfaceElementOrientations( VSet<1U>& vset );

} // csmp


#endif 

