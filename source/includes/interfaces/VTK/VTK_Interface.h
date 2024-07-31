#ifndef CSMP_VTK_INTERFACE_H
#define CSMP_VTK_INTERFACE_H

#include "FiniteElement.h"
#include "VTK_Type.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Region;
template<uint32_t> class Element;
template<uint32_t> class VTK_Interface;

/// for the given CSMP element type, finds the matching VTK geometric primitive
VTK_TYPE parseElementType( CSMP_FEM_TYPE );

/// prints the supplied node variable data DATA(dim x nodes) to file
template<uint32_t dim>
void outputNodeDataToVTK( const Element<dim>&, const std::string& file_name,
                          const std::string& var_name, DenseMatrix<DM_MIN>& DATA );

/// prints node variable data to file
template<uint32_t dim>
void outputNodeDataToVTK( const Element<dim>&, const csmp::Index&,
                          const char* file_name, const char* var_name );

/// prints element or node data to file
void outputDataToVTK( const char* file_name, const char* var_name,
                      VTK_TYPE vtk_type,
                      size_t cell_idx,
                      const DenseMatrix<DM_MIN>& XY,
                      const DenseMatrix<DM_MIN>& DATA );
                      
/// prints the face normals scaled by element size to file
template<uint32_t dim>
void outputFaceNormalsToVTK( const Element<dim>&, const char* file );

/// prints the variable values at the integration points to file
template<uint32_t dim>
void outputIntegrationPointDataToVTK( const Element<dim>&, const csmp::Index&, const char* var_name, const char* file );

/// prints integration points to file, scalar data are the integration point numbers
template<uint32_t dim>
void outputIntegrationPointsToVTK( const Element<dim>&, const char* file );

/// output the perimeter line of a region to VTK file
void outputRegionBoundaryToVTK( const Model<2U>&, const char* region, const char* file );

/// output the perimeter surface of a region to VTK file
void outputRegionBoundaryToVTK( const Model<3U>&, const char* region, const char* file );

/// all properties discretised on the model; region by region
template<uint32_t dim>
void  outputPropertiesOfRegionToVTK( const Model<dim>&,
                                     VTK_Interface<dim>&,
                                     const char* region, long output_time );

/// interactive function that allows extraction of desired variable from csmp binary file into a VTK file
void csmpBinaryToVTK( const char* modelBinFIleName );


/** Outputs csmp data to VTK Ascii text format

Known problems:

@todo (3) Outputs an inconsistent number of elements and number of points for a quadratic element in case of vector output (A)
@todo (3) Fails to output constraint point properties correctly (A)
@todo (3) Display quadratic FEM with the corresponding VTK primitives (A)
@todo (3) Display 2D regions using VTK polygons (A)
@todo (3) Save data as VTK binaries (A)
@todo (3) Create and store connectivity only once and then just add new data blocks to files (A)
*/

template<uint32_t dim>
class VTK_Interface {

  public:

    VTK_Interface( const std::string& problemTitle = "CSMP_Simulation", bool use_problem_title_as_output_folder_name = false );
    VTK_Interface( const std::string& problemTitle, const std::string& subFolderTitle, bool use_problem_title_as_output_folder_name = false );

    ~VTK_Interface();

    /// for all node property data associated with Model object
    template<class T>
    void OutputNodeDataToVTK( const Model<dim>&,
                              const std::string& file_name,
                              T timestep );

    /// for all node property data associated with a specific region
    template<class T>
    void OutputNodeDataToVTK( const Model<dim>&,
                              const std::string& domain_name,
                              const std::string& file_name,
                              T timestep );

    /// for entire models represented by Model objects
    template<class T>
    void OutputDataToVTK( const Model<dim>&,
                          const std::string& file_name,
                          const std::string& var_name,
                          T timestep,
                          bool update_topology=true );

    /// for a specific subregion of a model
    template<class T>
    void OutputDataToVTK( const Model<dim>&,
                          const std::string& domain_name,
                          const std::string& file_name,
                          const std::string& var_name,
                          T timestep,
                          bool update_topology=true );
                          
    /// all subregion of model, with file_name preface and var-name attached
    template<class T>
    void OutputRegionByRegionToVTK( const Model<dim>&,
                                    const std::string& file_name,
                                    const std::string& var_name,
                                    T timestep,
                                    bool update_topology=true );

    bool NodeOutputOfElementData() const;
    void NodeOutputOfElementData( bool yes_or_no );

 private:
    struct PointDescriptor {
        size_t element_or_node_;
        size_t sector_or_facet_;
    };
   bool toFolder_;                                           ///< output to folder
   bool toSubFolder_;                                        ///< output to subfolder
   std::string problemTitle_;                                ///< title that goes into xml header
   std::string subFolderName_;                               ///< name of subfolder
   std::deque<VTK_TYPE>                     geometric_primitives_VTK_;
   std::map<size_t,std::vector<size_t> >    plist_;
   std::deque<std::vector<size_t> >         transformed_plist_;
   std::map<size_t,PointDescriptor>         node_mapping_;
   std::map<size_t,std::vector<double> >    pxyz_data_;
   PLACEMENT                                last_visualized_;
   bool                                     node_output_of_element_data_;

   std::string OutputFileAndSubFolderName( const std::string& );

   void NodeBasedTopology( const Region<dim>&,
                           const std::vector<size_t>& elmt_ids,
                           std::map<size_t,std::vector<size_t> >& plist,
                           std::map<size_t,PointDescriptor>&  node_nums );

   void ElmtIntegrationPointBasedTopology( const Region<dim>&,
                                      const std::vector<size_t>& elmt_ids,
                                      std::map<size_t,std::vector<size_t> >& clist,
                                      std::map<size_t,PointDescriptor>&  cpoint_nums );

   void FacetIntegrationPointBasedTopology( const Region<dim>&,
                                      const std::vector<size_t>& elmt_ids,
                                      std::map<size_t,std::vector<size_t> >& clist,
                                      std::map<size_t,PointDescriptor>&  cpoint_nums );

   void ElementBasedTopology( const Region<dim>&,
                              const std::vector<size_t>& elmt_ids,
                              std::map<size_t,std::vector<size_t> >&  plist,
                              std::map<size_t,size_t>&  node_nums );

   void PointBasedTopology( const Region<dim>&,
                            const std::vector<size_t>& elmt_ids,
                            std::map<size_t,std::vector<size_t> >&   plist,
                            std::map<size_t,PointDescriptor>&  node_nums );
                                
   void TransformPlist( const Region<dim>&,
                        const std::map<size_t,std::vector<size_t> >& plist,
                        std::deque<std::vector<size_t> >& tplist );

   void RetrieveData( const Region<dim>&,
                      const csmp::Index&,
                      const std::map<size_t,PointDescriptor>& obj_nums,
                      std::map<size_t,std::vector<double> >& sgdata );

   void NodeCoordinates( const Region<dim>&,
                         const std::map<size_t,PointDescriptor>& node_nums,
                         std::map<size_t,std::vector<double> >& pxyz_data );

   void IntegrationPointData( const Region<dim>&,
                             const csmp::Index&,
                             std::map<size_t,std::vector<double> >& pxyz_data );   

   void ElmtIntegrationPointData( const Region<dim>&,
                             const csmp::Index&,
                             std::map<size_t,std::vector<double> >& pxyz_data );   

   void FacetIntegrationPointData( const Region<dim>&,
                             const csmp::Index&,
                             std::map<size_t,std::vector<double> >& pxyz_data );   

   void NodeData( const Region<dim>&,
                  const csmp::Index&,
                  const std::map<size_t,PointDescriptor>& node_nums,
                  std::map<size_t,std::vector<double> >& pxyz_data );

   void ElementData( const Region<dim>&,
                     const csmp::Index&,
                     const std::map<size_t,size_t>& node_nums,
                     std::map<size_t,std::vector<double> >& pxyz_data );

   void ElementPointData( const Region<dim>&,
                          const csmp::Index&,
                          const std::map<size_t,PointDescriptor>& node_nums,
                          std::map<size_t,std::vector<double> >& pxyz_data );

   void CellData( const Region<dim>&,
                  const csmp::Index&,
                  const std::map<size_t,std::vector<size_t> >& plist,
                  std::map<size_t,std::vector<double> >& pxyz_data );
};

} // csmp

#endif




