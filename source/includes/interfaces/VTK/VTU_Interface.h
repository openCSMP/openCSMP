#ifndef VTU_INTERFACE_H
#define VTU_INTERFACE_H

#include "CSMP_definitions.h"
#include "VTK_Type.h"

#include "Exception.h"
#include "CSMP_highLevelUtilities.h"
#include "PL_Utilities.h"
#include "XML_Document.h"

#include "Model.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

namespace csmp {

// forwards
struct Index;
template<uint32_t> class Model;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Node;
template<uint32_t> class VTU_Interface;
template<uint32_t,template<uint32_t> class CELL> class ModelSubDomain;

// hack forwards
double zCoordinate( Point<3> const& );
double zCoordinate( Point<2> const& );
double zCoordinate( Point<1> const& );
double yCoordinate( Point<3> const& );
double yCoordinate( Point<2> const& );
double yCoordinate( Point<1> const& );

/// writes model boundaries to a binary file
void outputBoundariesToVTU( const char* modelBinFIleName );

/// interactive, reads csmp binary model file and writes user-specified variable to VTU file.
void csmpBinaryToVTU( const char* modelBinFIleName );

// we need this since VTK_Interface does not support 1D; definitions in VTU_Interface.hpp
namespace vtuInterfaceDispatch
{
  template<uint32_t dim,class T>
  void dispatchCompileTimeToRuntimeVTK_Output( const csmp::Model<dim>& model, const std::string& subDomainName, const std::string& propertyName, T timeStep );
}

/**
    @brief Interface for the creation of XML files for the visualisation toolkit (VTK).

    @author P. Lang
    @date 2010

    @attention still issues with nan values(when attempting to visualized uninitialized vars in paraview)

    Interface for serial unstructured VTK XML files (.vtu).
    The VTU file format is considered to be of greater flexibility
    and to offer more options than the legacy format(VTK). The main
    advantage of this interface lies in the fact that the connectivity
    part is only build once for each region, and is reused in consecutive
    outputs. This results in major speed increase during repetitive
    outputs.

    It is possible to output multiple variables to one file, the interface
    automatically separates element placed vectors and tensors to another
    file, since a different connectivity is required (point cloud output).

    @section features Features
    The key features of this interface are:
    -Direct mapping from CSMP to VTK elements / nodes
    -Linear and Quadratic elements
    -Node and Element properties
    -Caching of connectivity data for speed up in consecutive outputs
    -Writing several property data into one file
    -Output of 2D and 3D models

    @section application Application Example

    The methods of this interface resemble those of the VTK_Interface.
    To instantiata a VTU_Interface, dimension and a reference to the
    model has to be provided

    In the simplest case, all you need is
    @code
    VTU_Interface<SPACE> vtu( model );
    vtu.OutputDataToVTU( "BOX", "fluid pressure" );
    vtu.OutputDataToVTU( "BOX", "fluid pressure", "FRACTURES" );
    vtu.OutputDataToVTU( "BOX", "fluid pressure", "FRACTURES", model_time );
    @endcode

    For further details:

    @code
    // dimension, csmp::Model and an optional session/simulation title
    VTU_Interface<3U> vtu( model, "Debug Session" );
    @endcode

    If you elect to output one variable only:

    @code
    // 'permeability' of entire model at time 0
    // file will be named 'ModelPermeability0.vtu'
    vtu.OutputDataToVTU( "ModelPermeability", "permeability" );

    // as above, but for region 'FRACTURES'
    // file will be named 'FRACTURES_Permeability1333.vtu' for model_time = 1333
    vtu.OutputDataToVTU( "Permeability", "permeability", "FRACTURES",  model_time );
    @endcode

    If you wish to output multiple properties:

    @code
    // create a list of the variables at the beginning
    list<string> vtuVariables;
    props.push_back( "fluid pressure" );
    props.push_back( "total mobility" );
    props.push_back( "saturation oil" );

    //...

    // to output all properties to single VTU file for entire model
    // file will be named 'Simulation0.vtu'
    vtu.OutputDataToVTU( "Simulation", vtuVariables );

    // as above but for specific region
    // file will be named 'FRACTURES_Simulation0.vtu' etc
    vtu.OutputDataToVTU( "Simulation", vtuVariables, "FRACTURES" );

    // as above but for specific region and model time
    // file will be named 'FRACTURES_Simulation1333.vtu' etc
    vtu.OutputDataToVTU( "Simulation", vtuVariables, "FRACTURES",  model_time );
    @endcode

    @section elements Featured CSMP Elements

    csmp Elements are directly mapped to their vtk counterpart(no simplification
    of geometry/transformation of nodes/element types)

    @code
    ISOPARAMETRIC_LINEAR_BAR
    LINEAR_BAR
    LINEAR_TRIANGLE3D
    ISOPARAMETRIC_LINEAR_TRIANGLE
    LINEAR_TRIANGLE
    ISOPARAMETRIC_LINEAR_TETRAHEDRON
    LINEAR_TETRAHEDRON
    ISOPARAMETRIC_LINEAR_PYRAMID
    ISOPARAMETRIC_LINEAR_PRISM
    ISOPARAMETRIC_LINEAR_QUADRILATERAL
    ISOPARAMETRIC_LINEAR_HEXAHEDRON
    ISOPARAMETRIC_QUADRATIC_BAR
    QUADRATIC_BAR
    ISOPARAMETRIC_QUADRATIC_TRIANGLE
    QUADRATIC_TRIANGLE
    ISOPARAMETRIC_QUADRATIC_QUADRILATERAL
    ISOPARAMETRIC_QUADRATIC_TETRAHEDRON
    QUADRATIC_TETRAHEDRON
    ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20
    ISOPARAMETRIC_QUADRATIC_PRISM15
    ISOPARAMETRIC_QUADRATIC_PYRAMID13
    @endcode

    @section implementation Implementation

    Requires 'VTK_Interface.h' for VTK element enumeration.

    For the sake of user readability, binary output was postponed to a later
    point in development, as was offset/appended data arrays.

    Furthermore, Point locations are described by a single data array, holding their
    x, y and z position(in this order). This approach also should increase readability
    compared to the alternative of having 3 seperate DataArrays, each holding the x, y and
    z coordinates for all points.

    Instantiating creates the connectivity file that may be reused consequently.
    Calling OutputToVTU opens a new file, writes Headers, Mesh and Data and closes
    the file, leaving vtuFile_ == 0.

    The node numbering convention of VTK is consistent with CSMP's for linear FEs, but
    not so for quadratic ones.

    Uses member function pointers for visitor like node/element loops when writing
    data arrays.

    Convention is that all outputs end with a new line each.

    @todo (3) Fix tensor output( if dim == 2 )
    @todo (3) Exceptions
    @todo (3) Implement true/false checks in bool members
    @todo (3) In insertConnectivityFile<dim> check whether insertion sucessfull
    @todo (3) Update doc
    @todo (3) think about bringing in mem alloc class to manage vtu files
    @todo (3) Additional info in xml header(see legacy vtu interface)
    @todo (3) Replace list by set as property container 
 
    @note minor refactoring by SKM.
    
    @todo break class into smaller functional units and remove code bloat by duplication.
    @todo simplify interfaces.
 
*/
template<uint32_t dim>
class VTU_Interface {
  public:

    VTU_Interface( const Model<dim>&, const std::string& problemTitle = "CSMP_Simulation", bool use_problem_title_as_output_folder_name = false );
    VTU_Interface( const Model<dim>&, const std::string& problemTitle, const std::string& subFolderTitle, bool use_problem_title_as_output_folder_name = false );

    ~VTU_Interface();

    // CSMP DATA OUTPUT
    // ----------------

    /// output a single property to VTU for given region
    template<class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::string& propertyName,
                          const std::string& regionName = "Model",
                          T timestep = static_cast<T>(0) );
    
    /// output a list of properties to VTU for given region
    template<class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::list<std::string>& propertyNames,
                          const std::string& regionName = "Model",
                          T timestep = static_cast<T>(0) );
    
    /// output a vector of properties to VTU for given region
    template<class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::vector<std::string>& propertyNames,
                          const std::string& regionName = "Model",
                          T timestep = static_cast<T>(0) );
    
    /// output a set of properties to VTU for given region
    template<class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::set<std::string>& propertyNames,
                          const std::string& regionName = "Model",
                          T timestep = static_cast<T>(0) );

    /// output a single property to VTU for given model sub domain (Region, Boundary or SplitBoundary)
    template<template <uint32_t> class CELL,class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::string& propertyName,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );

    /// output a list of properties to VTU for given model sub domain (Region, Boundary or SplitBoundary)
    template<template <uint32_t> class CELL,class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::list<std::string>& propertyNames,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );
    
    /// output a vector of properties to VTU for given model sub domain (Region, Boundary or SplitBoundary)
    template<template <uint32_t> class CELL,class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::vector<std::string>& propertyNames,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );
    
    /// KEY METHOD: output a set of properties to VTU for given model sub domain (Region, Boundary or SplitBoundary)
    template<template <uint32_t> class CELL,class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::set<std::string>& propertyNames,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );

    // SINGLE VERTEX OUTPUT
    // --------------------

    /// output vectors in origin of coordinate system to VTU file
    static bool OutputVectorsToVTU( const std::string& fileName, const std::string& propertyCaption,
                                    const std::vector<std::vector<double> >& vectors );
    
    /// output principle axes with given length in origin of coordinate system to VTU file
    static bool OutputPrincipalVectorsToVTU( const std::string& fileName, const std::string& propertyCaption,
                                             const std::vector<double>& xyzLengths );

    /// output tensor placed in origin of coordinate system to VTU file
    static bool OutputTensorToVTU( const std::string& fileName, const std::string& propertyCaption,
                                   const std::vector<std::vector<double> >& tensor );


    // SETUP PARTICULAR OUTPUT OPTIONS
    // -------------------------------

    /// provide aliases for csmp variables that will be shown in vtu instead
    bool VariableNameAliases( const std::map<std::string,std::string>& variableNameAliases );

    /**
        allows the creation of an index-to-name correspondance.  If the resulting container is
        empty or the array variable name does not exist upon the output call, the behaviour is the default one (output of [#] suffix).
        If the component name exists,it gets prefixed to the variable name for the vtu file. Works for flagged arrays, too.
        Julian, July 2014
    */
    void CreateArrayComponentPrefixNames(std::string array_var_name, std::vector<std::string>& index_to_name );

    /// Allows adding text after the timestep number to the vtu filename. By default, this string is empty.
    void SetSuffixText(std::string text);
    const std::string& GetProblemTitle() const;

    /// control whether '0' is appended
    void OmitZeroInFileName( bool omitZeroInFileName );
    bool OmitZeroInFileName() const;

    /// setup
    void OutputElementVectorAndTensorDataAtCellCenters ( bool element_data_at_cell_centers);
    void OutputRegionVectorAndTensorDataAtRegionCenters( bool region_data_at_cell_centers );

    /// allows to force recreation of connectivity map (FOR THE UPCOMING OUTPUT ONLY!!!)
    void DeleteConnectivity();

  protected:

    VTU_Interface();

    std::string FindVariableOutputAlias( const std::string& csmpVariableName ) const;
    std::string OutputFileNamePrefix( const std::string& );
    std::string DomainSpecificOutputFileNamePrefix( const std::string& fileName, const std::string& domainName );
    template<class T>
    std::string FullOutputFileName( const std::string& fileName, T timestep );

    template<template <uint32_t> class CELL>
    void OutputMultiBlockVTU( const std::string& fileName,
                          const std::vector<std::string>& fileNames,
                          const ModelSubDomain<dim,CELL>& subDomain );
                          
    template<template <uint32_t> class CELL>
    bool OutputFieldNodesAndElementDataToVTU( const std::string& fileName,
                          const std::list<csmp::Index>& fieldDataIndices,
                          const std::list<csmp::Index>& nodeIndices,
                          const std::list<csmp::Index>& elementIndices,
                          const ModelSubDomain<dim,CELL>& subDomain );
                          
    template<template <uint32_t> class CELL>
    bool OutputElementBarycentricDataToVTU( const std::string& fileName,
                          const std::list<csmp::Index>& elementMatrixIndices,
                          const ModelSubDomain<dim,CELL>& subDomain );
                          
    template<template <uint32_t> class CELL>
    bool OutputRegionDataToVTU( const std::string& fileName,
                          const std::list<csmp::Index>& regionIndices,
                          const ModelSubDomain<dim,CELL>& subDomain );
                          
    template<template <uint32_t> class CELL>
    bool OutputFiniteElementIntegrationPointsDataToVTU( const std::string& fileName,
                          const std::list<csmp::Index>& feipIndices,
                          const ModelSubDomain<dim,CELL>& subDomain );
                          
    template<template <uint32_t> class CELL>
    bool OutputFiniteVolumeSectorIntegrationPointsDataToVTU( const std::string& fileName,
                          const std::list<csmp::Index>& fvsipIndices,
                          const ModelSubDomain<dim,CELL>& subDomain );
                          
    template<template <uint32_t> class CELL>
    bool OutputFiniteVolumeFacetIntegrationPointsDataToVTU( const std::string& fileName,
                          const std::list<csmp::Index>& fvfipIndices,
                          const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    void OutputFieldDataToVTU( XML_Document& outputFile,
                               const ModelSubDomain<dim,CELL>& subDomain,
                               const std::list<Index>& indices );

    template<template <uint32_t> class CELL>
    void OutputPointDataToVTU( XML_Document& outputFile,
                               const ModelSubDomain<dim,CELL>& subDomain,
                               const std::list<Index>& indices );

    template<template <uint32_t> class CELL>
    void OutputCellDataToVTU( XML_Document& outputFile,
                              const ModelSubDomain<dim,CELL>& subDomain,
                              const std::list<Index>& indices );
    /// write variables
    void WriteScalar( XML_Document& vtu,const size_t& MAX_ENTRIES_PER_LINE,
                      double variable, size_t& entriesOfLine, bool& newLine ) const;
    void WriteVector( XML_Document& vtu,const size_t& MAX_ENTRIES_PER_LINE,
                      const VectorVariable<dim>& variable, size_t& entriesOfLine, bool& newLine ) const;
    void WriteTensor( XML_Document& vtu,const size_t& MAX_ENTRIES_PER_LINE,
                      const TensorVariable<dim>& variable, size_t& entriesOfLine, bool& newLine ) const;

    /// write field data
    template<class Var>
    void WriteFieldDataArray( const Index& key, XML_Document& vtu) const;

    /// write point data
    template<template <uint32_t> class CELL>
    void WritePointDataArrayScalar( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayVector( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayTensor( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayScalarArray( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayScalarFlaggedArray( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// write element data
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayScalar( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayVector( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayTensor( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayScalarArray( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayScalarFlaggedArray( const Index&, XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Connnectivity File's @todo how about some typedefs?
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMap( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMapBCS( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMapRCS( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMapFEIPS( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMapFVSIPS( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMapFVFIPS( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>& GetConnectivityMapMultiBlock( const ModelSubDomain<dim,Element>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMap( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMapBCS( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMapRCS( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMapFEIPS( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMapFVSIPS( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMapFVFIPS( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>& GetConnectivityMapMultiBlock( const ModelSubDomain<dim,Face>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMap( const ModelSubDomain<dim,InterFace>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMapBCS( const ModelSubDomain<dim,InterFace>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMapRCS( const ModelSubDomain<dim,InterFace>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMapFEIPS( const ModelSubDomain<dim,InterFace>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMapFVSIPS( const ModelSubDomain<dim,InterFace>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMapFVFIPS( const ModelSubDomain<dim,InterFace>& subDomain );
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& GetConnectivityMapMultiBlock( const ModelSubDomain<dim,InterFace>& subDomain );

    template<template <uint32_t> class CELL>
    XML_Document* findConnectivityFile( std::map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap, const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool insertConnectivityFile( std::map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap, XML_Document* newConnectivityFile, const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    XML_Document* ConnectivityFile( std::map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap, const ModelSubDomain<dim,CELL>& subDomain );

    bool CloseFile( const std::string& fileName, const std::string& extension, XML_Document& );
    void EstablishConnectivityFileHeader( XML_Document& ) const;

    /// Elements and Nodes
    template<template <uint32_t> class CELL>
    void EstablishConnectivityFile( XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Cell centers
    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileBCPC( XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Region's Point Clouds
    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileRPC( XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Finite Element Integration Points
    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileFEIP( XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Finite Volume Sector Integration Points
    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileFVSIP( XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Finite Volume Facet Integration Points
    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileFVFIP( XML_Document&, const ModelSubDomain<dim,CELL>& subDomain ) const;

    /// Quadratic Elements connectivity
    template<template <uint32_t> class CELL>
    void QuadraticWedgeConnectivity( const CELL<dim>* const, std::vector<size_t>& data ) const;
    template<template <uint32_t> class CELL>
    void QuadraticHexahedronConnectivity( const CELL<dim>* const, std::vector<size_t>& data ) const;

    template<template <uint32_t> class CELL>
    VTK_TYPE ElementType( const CELL<dim>* const elmt ) const;

    std::string   DomainName( const ModelSubDomain<dim,Element>& subDomain ) const;
    std::string   DomainName( const ModelSubDomain<dim,Face>& subDomain ) const;
    std::string   DomainName( const ModelSubDomain<dim,InterFace>& subDomain ) const;

  private:

    const Model<dim>&                 model_;                     ///< reference to csmp model of interest
    bool                              toFolder_;                  ///< output to folder
    bool                              toSubFolder_;               ///< output to subfolder
    std::string                       suffix_text_;
    std::string                       problemTitle_;              ///< title that goes into xml header
    std::string                       subFolderName_;             ///< name of subfolder
    bool                              omitZeroInFileName_;        ///< if false(default), '0' is appended to filename if no time privided (i.e. timestep=0)
    bool                              elementVecAndTensDataAtCellCenters_;    /// by default: false, if true then a separate file will be created
    bool                              regionVecAndTensDataAtCellCenters_;     /// by default: false( assigned as field data ), if true then separate file will be created
    std::map<std::string,std::string> variableNameAliases_;                   ///< if not NULL(default) this list contains aliases for csmp vars as brought into vtu
    std::map<std::string, std::vector<std::string> > array_index_to_name_;    /// if empty for a particular variable name, default behaviour,
                                                                              /// otherwise, it should contain the corresponding name for each array component
                                                                              /// (array variable name is the key of the outer map). Should work for FlaggedArray, too.

    /// files for field, node and element data
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_;                   ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_;                 ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_;            ///< xml connectivity files for all split boundaries of interest

    /// files for element barycentric data
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_bcd_;               ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_bcd_;             ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_bcd_;        ///< xml connectivity files for all split boundaries of interest

    /// files for regions data
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_rcd_;               ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_rcd_;             ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_rcd_;        ///< xml connectivity files for all split boundaries of interest

    /// files for finite element integration points data
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_feipsd_;            ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_feipsd_;          ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_feipsd_;     ///< xml connectivity files for all split boundaries of interest

    /// files for finite volume sector integration points data
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_fvsipsd_;           ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_fvsipsd_;         ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_fvsipsd_;    ///< xml connectivity files for all split boundaries of interest

    /// files for finite volume faect integration points data
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_fvfipsd_;           ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_fvfipsd_;         ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_fvfipsd_;    ///< xml connectivity files for all split boundaries of interest

    /// collection of all "*.vtu" files
    std::map<const ModelSubDomain<dim,Element>*,XML_Document*>        regionConnectivityFiles_multiblock_;        ///< xml connectivity files for all regions of interest
    std::map<const ModelSubDomain<dim,Face>*,XML_Document*>           boundaryConnectivityFiles_multiblock_;      ///< xml connectivity files for all boundaries of interest
    std::map<const ModelSubDomain<dim,InterFace>*,XML_Document*>      splitBoundaryConnectivityFiles_multiblock_; ///< xml connectivity files for all split boundaries of interest
};


/**
*/

}// csmp


#endif // VTU_INTERFACE_H

