#ifndef CGNS_INTERFACE_H
#define CGNS_INTERFACE_H

#include "CSMP_definitions.h"
#include "Model.h"
#include "VSet.h"
#include "ModelTopology.h"
#include "CGNS_ElementSpecifications.h"
#include "CSMP_ElementSpecifications.h"

namespace csmp {

/**
 
@brief CFD General Notation System (CGNS cfd data standard, NASA) interface to CSMP

cgns.org.

@author R. Manasipov
@date 2015

*/
class CGNS_Interface {

public:
    explicit CGNS_Interface( bool isoparametric_mesh = true );

    ~CGNS_Interface();
    
    template<uint32_t dim>
    int  Write_CGNS_Mesh( const std::string&    filename,
                          const csmp::Model<dim>& model );

    template<uint32_t dim>
    int  Read_CGNS_Mesh( const std::string&     filename,
                         csmp::VSet<dim>&       vset,
                         csmp::ModelTopology&   mesh_topology );

    int  SetFileType( int );
    int  GetFileType( ) const;
    void AssignMaterialPropertiesInteractively( bool assign = true );

private:

    /// writing mesh file
    template<uint32_t dim>
    int  WriteUnstructMesh( const std::string& filename,int cgfile,const csmp::Model<dim>& model );
    template<uint32_t dim>
    void WriteBase( const std::string& filename, int cgfile, int& cgbase );
    template<uint32_t dim>
    void WriteZone( const std::string& filename, int cgfile, int cgbase, int& cgzone, const csmp::Model<dim>& model );
    template<uint32_t dim>
    void WriteCoords( int cgfile, int cgbase, int cgzone, const csmp::Model<dim>& model );
    template<uint32_t dim>
    void WriteSubDomains( int cgfile, int cgbase, int cgzone, const csmp::Model<dim>& model );
    template<uint32_t dim,template <uint32_t> class SIMPLEX>
    void WriteElements( int cgfile, int cgbase, int cgzone, const csmp::ModelSubDomain<dim,SIMPLEX>& subDomain );

    /// reading mesh file
    template<uint32_t dim>
    int  ReadUnstructMesh( int cgfile, const std::string& filename, csmp::VSet<dim>& vset, csmp::ModelTopology&   mesh_topology );
    template<uint32_t dim>
    void ReadCoords( int cgfile, int cgbase, int cgzone, int total_num_coords, csmp::VSet<dim>& vset, csmp::ModelTopology&   mesh_topology );
    template<uint32_t dim>
    void ReadElements( int cgfile, int cgbase, int cgzone, int total_num_elements, csmp::VSet<dim>& vset, csmp::ModelTopology&   mesh_topology );

    void Clear();

  private:

    /// one of avaliable formats: CG_FILE_HDF5, CG_FILE_ADF2, CG_FILE_ADF
    int file_type_;

    /// property assignment
    bool interactive_property_assignment_;

    /// element specs
    bool  isoparametric_;
    CGNS_ElementSpecifications         elmt_specs_;
    using CSMP_ElementSpecifications   csmp_elmt_specs_;  ///< this is a static object
    size_t global_eid_;
    //             element nodes       element type  element id
    std::map<std::vector<uint32_t>,std::pair<size_t,size_t> > element_ids_;
};

/// CGNS - to - CSMP FE's
template<uint32_t dim>
void Convert_CGNS_To_CSMP_FiniteElementTypes( csmp::VSet<dim>&,bool);

class CGNS_ModelSettings
{
public:

    CGNS_ModelSettings( const std::string& mesh_file_prefix );
    ~CGNS_ModelSettings();

    CGNS_ModelSettings( const CGNS_ModelSettings& );
    CGNS_ModelSettings& operator=( const CGNS_ModelSettings& );

    void MeshSetup( const std::string& regions_file_prefix );
    void MeshSetup( const std::set<std::string>& regions );

public:

    std::string mesh_file_prefix_;    ///< name of model or main mesh file
    std::set<std::string> regions_;   ///< name of regions to be used
};

} // end namespace csmp

#endif
