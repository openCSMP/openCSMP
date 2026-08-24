#ifndef VTU_INTERFACE_H
#define VTU_INTERFACE_H

#include "CSMP_definitions.h"
#include "VTK_Type.h"
#include "XML_Document.h"

namespace csmp {

// forwards
struct Index;
template<uint32_t> class Model;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Node;
template<uint32_t> class VTU_Interface;
template<uint32_t, template<uint32_t> class CELL> class ModelSubDomain;

// hack forwards
double zCoordinate( Point<3> const& );
double zCoordinate( Point<2> const& );
double zCoordinate( Point<1> const& );
double yCoordinate( Point<3> const& );
double yCoordinate( Point<2> const& );
double yCoordinate( Point<1> const& );

/// writes model boundaries to a binary file
void outputBoundariesToVTU( const char* modelBinFileName );

/// interactive, reads csmp binary model file and writes user-specified variable to VTU file
void csmpBinaryToVTU( const char* modelBinFileName );

// we need this since VTK_Interface does not support 1D
namespace vtuInterfaceDispatch
{
    template<uint32_t dim, class T>
    void dispatchCompileTimeToRuntimeVTK_Output(
        const csmp::Model<dim>& model,
        const std::string& subDomainName,
        const std::string& propertyName,
        T timeStep );
}

/**
    @brief Interface for the creation of XML files for the visualisation toolkit (VTK).

    @author P. Lang
    @date 2010

    Interface for serial unstructured VTK XML files (.vtu).

    @see Full documentation in original header.
*/
template<uint32_t dim>
class VTU_Interface {
public:

    // ------------------------------------------------------------------
    // Convenience alias templates
    //
    // ConnectivityMap<CELL>  replaces:
    //   std::map<const ModelSubDomain<dim,CELL>*, std::unique_ptr<XML_Document>>
    //
    // Usage inside the class:  ConnectivityMap<Element>
    // ------------------------------------------------------------------
    template<template <uint32_t> class CELL>
    using ConnectivityMap = std::map<const ModelSubDomain<dim, CELL>*, std::unique_ptr<XML_Document>>;

    // ------------------------------------------------------------------
    // Construction / destruction
    // ------------------------------------------------------------------

    VTU_Interface( const Model<dim>&,
                   const std::string& problemTitle = "CSMP_Simulation",
                   bool use_problem_title_as_output_folder_name = false );

    VTU_Interface( const Model<dim>&,
                   const std::string& problemTitle,
                   const std::string& subFolderTitle,
                   bool use_problem_title_as_output_folder_name = false );

    ~VTU_Interface();

    // ------------------------------------------------------------------
    // CSMP DATA OUTPUT
    // ------------------------------------------------------------------

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

    /// output a single property to VTU for given model sub domain
    template<template <uint32_t> class CELL, class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::string& propertyName,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );

    /// output a list of properties to VTU for given model sub domain
    template<template <uint32_t> class CELL, class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::list<std::string>& propertyNames,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );

    /// output a vector of properties to VTU for given model sub domain
    template<template <uint32_t> class CELL, class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::vector<std::string>& propertyNames,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );

    /// KEY METHOD: output a set of properties to VTU for given model sub domain
    template<template <uint32_t> class CELL, class T>
    bool OutputDataToVTU( const std::string& fileName,
                          const std::set<std::string>& propertyNames,
                          const ModelSubDomain<dim,CELL>& subDomain,
                          T timestep = static_cast<T>(0) );

    // ------------------------------------------------------------------
    // SINGLE VERTEX OUTPUT
    // ------------------------------------------------------------------

    static bool OutputVectorsToVTU( const std::string& fileName,
                                    const std::string& propertyCaption,
                                    const std::vector<std::vector<double>>& vectors );

    static bool OutputPrincipalVectorsToVTU( const std::string& fileName,
                                             const std::string& propertyCaption,
                                             const std::vector<double>& xyzLengths );

    static bool OutputTensorToVTU( const std::string& fileName,
                                   const std::string& propertyCaption,
                                   const std::vector<std::vector<double>>& tensor );

    // ------------------------------------------------------------------
    // SETUP / OPTIONS
    // ------------------------------------------------------------------

    bool VariableNameAliases( const std::map<std::string,std::string>& variableNameAliases );

    void CreateArrayComponentPrefixNames( std::string array_var_name,
                                          std::vector<std::string>& index_to_name );

    void SetSuffixText( std::string text );
    const std::string& GetProblemTitle() const;

    void OmitZeroInFileName( bool omitZeroInFileName );
    bool OmitZeroInFileName() const;

    void OutputElementVectorAndTensorDataAtCellCenters ( bool element_data_at_cell_centers );
    void OutputRegionVectorAndTensorDataAtRegionCenters( bool region_data_at_cell_centers );

    /// forces recreation of connectivity map (FOR THE UPCOMING OUTPUT ONLY)
    void DeleteConnectivity() noexcept;

protected:

    VTU_Interface();

    std::string FindVariableOutputAlias( const std::string& csmpVariableName ) const;
    std::string OutputFileNamePrefix( const std::string& );
    std::string DomainSpecificOutputFileNamePrefix( const std::string& fileName,
                                                    const std::string& domainName );
    template<class T>
    std::string FullOutputFileName( const std::string& fileName, T timestep );

    template<template <uint32_t> class CELL>
    void OutputMultiBlockVTU( const std::string& fileName,
                              const std::vector<std::string>& fileNames,
                              const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool OutputFieldNodesAndElementDataToVTU(
                              const std::string& fileName,
                              const std::list<csmp::Index>& fieldDataIndices,
                              const std::list<csmp::Index>& nodeIndices,
                              const std::list<csmp::Index>& elementIndices,
                              const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool OutputElementBarycentricDataToVTU(
                              const std::string& fileName,
                              const std::list<csmp::Index>& elementMatrixIndices,
                              const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool OutputRegionDataToVTU( const std::string& fileName,
                                const std::list<csmp::Index>& regionIndices,
                                const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool OutputFiniteElementIntegrationPointsDataToVTU(
                              const std::string& fileName,
                              const std::list<csmp::Index>& feipIndices,
                              const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool OutputFiniteVolumeSectorIntegrationPointsDataToVTU(
                              const std::string& fileName,
                              const std::list<csmp::Index>& fvsipIndices,
                              const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool OutputFiniteVolumeFacetIntegrationPointsDataToVTU(
                              const std::string& fileName,
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

    // write variables
    void WriteScalar( XML_Document& vtu, const size_t& MAX_ENTRIES_PER_LINE,
                      double variable, size_t& entriesOfLine, bool& newLine ) const;
    void WriteVector( XML_Document& vtu, const size_t& MAX_ENTRIES_PER_LINE,
                      const VectorVariable<dim>& variable,
                      size_t& entriesOfLine, bool& newLine ) const;
    void WriteTensor( XML_Document& vtu, const size_t& MAX_ENTRIES_PER_LINE,
                      const TensorVariable<dim>& variable,
                      size_t& entriesOfLine, bool& newLine ) const;

    // write field data
    template<class Var>
    void WriteFieldDataArray( const Index& key, XML_Document& vtu ) const;

    // write point data
    template<template <uint32_t> class CELL>
    void WritePointDataArrayScalar( const Index&, XML_Document&,
                                    const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayVector( const Index&, XML_Document&,
                                    const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayTensor( const Index&, XML_Document&,
                                    const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayScalarArray( const Index&, XML_Document&,
                                         const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WritePointDataArrayScalarFlaggedArray( const Index&, XML_Document&,
                                                const ModelSubDomain<dim,CELL>& ) const;

    // write element data
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayScalar( const Index&, XML_Document&,
                                      const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayVector( const Index&, XML_Document&,
                                      const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayTensor( const Index&, XML_Document&,
                                      const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayScalarArray( const Index&, XML_Document&,
                                           const ModelSubDomain<dim,CELL>& ) const;
    template<template <uint32_t> class CELL>
    void WriteElementDataArrayScalarFlaggedArray( const Index&, XML_Document&,
                                                  const ModelSubDomain<dim,CELL>& ) const;

    // ------------------------------------------------------------------
    // Connectivity map accessors
    // Each group of three overloads dispatches on Element / Face / InterFace
    // ------------------------------------------------------------------
    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMap( const ModelSubDomain<dim,CELL>& );

    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMapBCS( const ModelSubDomain<dim,CELL>& );

    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMapRCS( const ModelSubDomain<dim,CELL>& );

    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMapFEIPS( const ModelSubDomain<dim,CELL>& );

    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMapFVSIPS( const ModelSubDomain<dim,CELL>& );

    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMapFVFIPS( const ModelSubDomain<dim,CELL>& );

    template<template <uint32_t> class CELL>
    ConnectivityMap<CELL>& GetConnectivityMapMultiBlock( const ModelSubDomain<dim,CELL>& );

    // ------------------------------------------------------------------
    // Connectivity file management
    // ------------------------------------------------------------------
    template<template <uint32_t> class CELL>
    XML_Document* findConnectivityFile( ConnectivityMap<CELL>& connectivityMap,
                                        const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    bool insertConnectivityFile( ConnectivityMap<CELL>& connectivityMap,
                                 std::unique_ptr<XML_Document> newConnectivityFile,
                                 const ModelSubDomain<dim,CELL>& subDomain );

    template<template <uint32_t> class CELL>
    XML_Document* ConnectivityFile( ConnectivityMap<CELL>& connectivityMap,
                                    const ModelSubDomain<dim,CELL>& subDomain );

    bool CloseFile( const std::string& fileName,
                    const std::string& extension,
                    XML_Document& );

    void EstablishConnectivityFileHeader( XML_Document& ) const;

    template<template <uint32_t> class CELL>
    void EstablishConnectivityFile( XML_Document&,
                                    const ModelSubDomain<dim,CELL>& ) const;

    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileBCPC( XML_Document&,
                                        const ModelSubDomain<dim,CELL>& ) const;

    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileRPC( XML_Document&,
                                       const ModelSubDomain<dim,CELL>& ) const;

    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileFEIP( XML_Document&,
                                        const ModelSubDomain<dim,CELL>& ) const;

    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileFVSIP( XML_Document&,
                                         const ModelSubDomain<dim,CELL>& ) const;

    template<template <uint32_t> class CELL>
    void EstablishConnectivityFileFVFIP( XML_Document&,
                                         const ModelSubDomain<dim,CELL>& ) const;

    template<template <uint32_t> class CELL>
    void QuadraticWedgeConnectivity( const CELL<dim>* const,
                                     std::vector<size_t>& data ) const;

    template<template <uint32_t> class CELL>
    void QuadraticHexahedronConnectivity( const CELL<dim>* const,
                                          std::vector<size_t>& data ) const;

    template<template <uint32_t> class CELL>
    VTK_TYPE ElementType( const CELL<dim>* const elmt ) const;

    std::string DomainName( const ModelSubDomain<dim,Element>&   ) const;
    std::string DomainName( const ModelSubDomain<dim,Face>&      ) const;
    std::string DomainName( const ModelSubDomain<dim,InterFace>& ) const;

private:

    // ------------------------------------------------------------------
    // Helper: clears one connectivity map, destroying all owned documents
    // ------------------------------------------------------------------
    template<template <uint32_t> class CELL>
    static void clearConnectivityMap( ConnectivityMap<CELL>& m ) noexcept
    {
        m.clear();   // unique_ptr destructors run here; no manual delete needed
    }

    // ------------------------------------------------------------------
    // Model reference and scalar options
    // ------------------------------------------------------------------
    const Model<dim>&                 model_;
    bool                              toFolder_;
    bool                              toSubFolder_;
    std::string                       suffix_text_;
    std::string                       problemTitle_;
    std::string                       subFolderName_;
    bool                              omitZeroInFileName_;
    bool                              elementVecAndTensDataAtCellCenters_;
    bool                              regionVecAndTensDataAtCellCenters_;
    std::map<std::string,std::string> variableNameAliases_;
    std::map<std::string, std::vector<std::string>> array_index_to_name_;

    // ------------------------------------------------------------------
    // Connectivity caches — 7 purposes x 3 cell types = 21 maps
    //
    // Each map owns its XML_Document objects via unique_ptr.
    // Destruction is automatic; DeleteConnectivity() provides early release.
    // ------------------------------------------------------------------

    /// Field, node and element data
    ConnectivityMap<Element>   regionConnectivityFiles_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_;

    /// Element barycentre data
    ConnectivityMap<Element>   regionConnectivityFiles_bcd_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_bcd_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_bcd_;

    /// Region data
    ConnectivityMap<Element>   regionConnectivityFiles_rcd_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_rcd_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_rcd_;

    /// Finite Element integration point data
    ConnectivityMap<Element>   regionConnectivityFiles_feipsd_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_feipsd_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_feipsd_;

    /// Finite Volume sector integration point data
    ConnectivityMap<Element>   regionConnectivityFiles_fvsipsd_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_fvsipsd_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_fvsipsd_;

    /// Finite Volume facet integration point data
    ConnectivityMap<Element>   regionConnectivityFiles_fvfipsd_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_fvfipsd_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_fvfipsd_;

    /// Multiblock dataset index files
    ConnectivityMap<Element>   regionConnectivityFiles_multiblock_;
    ConnectivityMap<Face>      boundaryConnectivityFiles_multiblock_;
    ConnectivityMap<InterFace> splitBoundaryConnectivityFiles_multiblock_;
};

} // namespace csmp

#endif // VTU_INTERFACE_H
