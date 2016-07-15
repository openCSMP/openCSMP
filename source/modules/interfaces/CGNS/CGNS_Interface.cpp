#include "CGNS_Interface.h"
#include "VSet.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "ModelTopology.h"
#include "Model.h"
#include "ErrorHandler.h"

/* General IO include section */

#if defined(_WIN32) && !defined(__NUTC__)
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif

/* CGNS include section */

/* cgnslib.h file must be located in directory specified by -I during compile: */
#include "cgnslib.h"

namespace csmp {

CGNS_Interface::CGNS_Interface( bool isoparametric_mesh )
    :file_type_(CG_FILE_ADF2),
     isoparametric_( isoparametric_mesh ),
     interactive_property_assignment_(false)
 {
 }
 
CGNS_Interface::~CGNS_Interface()
 {
 }

/** Avaliable file formats: CG_FILE_HDF5, CG_FILE_ADF2, CG_FILE_ADF
*/
int CGNS_Interface::SetFileType( int file_type )
{
    file_type_ = file_type;

    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    //std::string errmsg;
    if (cg_set_file_type( file_type_ ) )
    {
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::SetFileFormat():","cg_set_file_type");
        return 1;
    }

    return 0;
}

int CGNS_Interface::GetFileType( ) const
{
    return file_type_;
}

void CGNS_Interface::AssignMaterialPropertiesInteractively( bool ass )
{
    interactive_property_assignment_ = ass;
}

void CGNS_Interface::Clear()
{
    if ( !element_ids_.empty() )
      element_ids_.erase( element_ids_.begin(), element_ids_.end() );
}

/**
Reads simple 3-D unstructured grid from a CGNS file
*/
template<size_t dim>
int CGNS_Interface::Read_CGNS_Mesh( const std::string&   filename,
                                    csmp::VSet<dim>&     vset,
                                    csmp::ModelTopology& mesh_topology )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    std::string infile = filename;
    infile += ".cgns";

    /// open CGNS file for read-only
    int cgfile;
    if ( cg_open(infile.c_str(),CG_MODE_READ,&cgfile) )
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::Read_CGNS_Mesh():","cg_open");

    /// read mesh
    ReadUnstructMesh<dim>(cgfile,filename,vset,mesh_topology);

    /// close CGNS file
    cg_close(cgfile);

    return 0;
}


template int CGNS_Interface::Read_CGNS_Mesh( const std::string&,csmp::VSet<1U>&,csmp::ModelTopology& );
template int CGNS_Interface::Read_CGNS_Mesh( const std::string&,csmp::VSet<2U>&,csmp::ModelTopology& );
template int CGNS_Interface::Read_CGNS_Mesh( const std::string&,csmp::VSet<3U>&,csmp::ModelTopology& );


template<size_t dim>
int CGNS_Interface::ReadUnstructMesh( int cgfile, const std::string& filename, csmp::VSet<dim>& vset, csmp::ModelTopology&   mesh_topology )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    Clear();

    int cgbase;
    int cgzone;
    int nbase;
    int nzone;
    char zonename[33];

    /// check that there is only one base
    cg_nbases(cgfile,&nbase);
    if( nbase!=1 ){
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::ReadUnctructMesh():","Unexpected number of bases, works only for 1 base.");
        return 1;
    }
    cgbase=1;
    /// since cgns format does not provide explicit element numbering
    /// at least the way that the reading is implemented here requres
    /// explicit element numbering
    global_eid_ = 0;
    /// assign name to ModelTopology
    mesh_topology.ModelName( filename.c_str() );
    /// assign hybrid mesh type to VSet
    vset.HybridElementTypeMesh( true );

    /// read number of zones ( regions )
    cg_nzones(cgfile,cgbase,&nzone);

    /// do loop over the zones
    for( cgzone = 1; cgzone <= 1 /*nzone*/; ++cgzone ){
        cgsize_t size[3];
        /// get zone size (and name - although not needed here)
        cg_zone_read(cgfile,cgbase,cgzone,zonename,size);
        /// read grid coordinates
        ReadCoords<dim>( cgfile, cgbase, cgzone, size[0], vset, mesh_topology );
        /// read elements
        ReadElements<dim>( cgfile, cgbase, cgzone, size[1], vset, mesh_topology );
    }

    /// construct ModelTopology
    bool require_unique_names_of_volumes_surfaces_and_lines = true;
    bool correct_orientation_of_surface_elements = false;
    bool non_box_boundary = true;
    mesh_topology.CheckTopology( vset,
                                 require_unique_names_of_volumes_surfaces_and_lines,
                                 correct_orientation_of_surface_elements,
                                 non_box_boundary );

    Clear();

    if( error_handler.Verbose() )
        std::cout<<"\nCGNS_Interface::ReadUnstructMesh():Successfully read unstructured grid from file"<<filename<<"\n";

    return 0;
}

template int CGNS_Interface::ReadUnstructMesh<1U>( int,const std::string&,csmp::VSet<1U>&,csmp::ModelTopology& );
template int CGNS_Interface::ReadUnstructMesh<2U>( int,const std::string&,csmp::VSet<2U>&,csmp::ModelTopology& );
template int CGNS_Interface::ReadUnstructMesh<3U>( int,const std::string&,csmp::VSet<3U>&,csmp::ModelTopology& );


template<size_t dim>
void CGNS_Interface::ReadCoords( int cgfile, int cgbase, int cgzone, int total_num_coords, csmp::VSet<dim>& vset, csmp::ModelTopology&   mesh_topology  )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    float* xcoord = new float[ total_num_coords ];
    float* ycoord = new float[ total_num_coords ];
    float* zcoord = new float[ total_num_coords ];

    /// lower range index
    cgsize_t irmin = 1;
    /// upper range index of vertices
    cgsize_t irmax = total_num_coords;
    /// read grid coordinates
    cg_coord_read(cgfile,cgbase,cgzone,"CoordinateX",
                  RealSingle,&irmin,&irmax,xcoord);
    cg_coord_read(cgfile,cgbase,cgzone,"CoordinateY",
                  RealSingle,&irmin,&irmax,ycoord);
    cg_coord_read(cgfile,cgbase,cgzone,"CoordinateZ",
                  RealSingle,&irmin,&irmax,zcoord);

    /// adding xcoord,ycoord,zcoord to VSet
    vset.ResizeNodes( total_num_coords );
    for ( int i=0U; i<total_num_coords; ++i ) {
         vset.Px( i, xcoord[i] );
         vset.Py( i, ycoord[i] );
         vset.Pz( i, zcoord[i] );
      }
    delete[] xcoord;
    delete[] ycoord;
    delete[] zcoord;
}

template void CGNS_Interface::ReadCoords<1U>( int,int,int,int,csmp::VSet<1U>&,csmp::ModelTopology& );
template void CGNS_Interface::ReadCoords<2U>( int,int,int,int,csmp::VSet<2U>&,csmp::ModelTopology& );
template void CGNS_Interface::ReadCoords<3U>( int,int,int,int,csmp::VSet<3U>&,csmp::ModelTopology& );

template<size_t dim>
void CGNS_Interface::ReadElements( int cgfile, int cgbase, int cgzone, int total_num_elements, csmp::VSet<dim>& vset, csmp::ModelTopology&   mesh_topology  )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// reading sections
    int cgsect;

    /// find out how many sections
    int  nsections;
    cg_nsections(cgfile,cgbase,cgzone,&nsections);
    if( error_handler.Verbose() )
        std::cout<<"\nCGNS_Interface::ReadCoordsAndElements(): number of sections=" << nsections << std::endl;

    /// read element connectivity
    cgsize_t* elements( NULL );
    cgsize_t* parent_data( NULL );
    cgsize_t  istart,iend;
    cgsize_t  element_data_size;
    int  num_elements;
    int  nperelmt;
    int  iparent_flag;
    int  nbndry;
    ElementType_t cgns_fem_type;
    csmp::CSMP_FEM_TYPE csmp_fem_type;
    std::string csmp_fem_type_name;
    std::string section_name;
    std::pair<typename std::map<std::vector<size_t>,std::pair<size_t,size_t> >::iterator,bool> eit;
    for ( cgsect=1; cgsect <= nsections; ++cgsect )
    {
        char sname[33];
        cg_section_read(cgfile,cgbase,cgzone,cgsect,sname,
                      &cgns_fem_type,&istart,&iend,&nbndry,&iparent_flag);
        section_name = sname;
        //cg_ElementDataSize(cgfile,cgbase,cgzone,cgsect,&element_data_size);
        num_elements = iend - istart + 1;
        std::vector<size_t> section_element_ids;
        if( error_handler.Verbose() )
        {
            std::cout<<"\nCGNS_Interface::ReadCoordsAndElements():Reading section data...\n";
            std::cout<<"   section name=" << section_name << std::endl;
            std::cout<<"   section type=" << ElementTypeName[cgns_fem_type] << std::endl;
            std::cout<<"   istart,iend=" << (int)istart << ", "<< (int)iend << std::endl;
        }
        if( cgns_fem_type == CGNS_ENUMV(MIXED) )
        {
            cg_ElementDataSize(cgfile,cgbase,cgzone,cgsect,&element_data_size);
            elements = new cgsize_t[ element_data_size ];
            if( error_handler.Verbose() )
                std::cout<<"   reading element data for this element\n";
            /// filling fem types and element nodes
            cg_elements_read(cgfile,cgbase,cgzone,cgsect,elements,parent_data);
            std::set<std::string>    section_element_types;
            size_t i = 0;
            while( i < element_data_size )
            {
                csmp_fem_type = elmt_specs_.CSMP_TypeFrom_CGNS_Type( elements[i++], isoparametric_, dim );
                csmp_fem_type_name = csmp_elmt_specs_.CSMP_TypeName( csmp_fem_type );
                nperelmt = csmp_elmt_specs_.NodesPerElementOfType( csmp_fem_type );
                std::vector<size_t> nids( nperelmt );
                size_t iend = i + nperelmt;
                size_t j = 0;
                while( i < iend )
                    nids[j++] = elements[i++] - 1;
                eit = element_ids_.insert( std::make_pair( nids, std::make_pair(static_cast<size_t>(csmp_fem_type), 0 ) ) );
                /// if element was not yet considered
                if ( eit.second )
                {
                    (*eit.first).second.second = global_eid_++;
                    /// filling VSet
                    vset.ResizeElementTypes( global_eid_ );
                    vset.ResizePlist( global_eid_, nperelmt );
                    vset.ElementType( (*eit.first).second.second, csmp_fem_type );
                    for( size_t nid = 0 ; nid < nperelmt; ++nid )
                        vset.Plist( (*eit.first).second.second, nid, nids[nid] );
                }
                /// filling MeshTopology
                section_element_types.insert( csmp_fem_type_name );
                section_element_ids.push_back( (*eit.first).second.second );
            }
            delete[] elements;
            /// adding section information to MeshTopology
            mesh_topology.AddRegionElementTypes( section_name.c_str(), section_element_types );
            mesh_topology.AddRegionElementIds( section_name.c_str(), section_element_ids );
        }
        else
        {
            csmp_fem_type = elmt_specs_.CSMP_TypeFrom_CGNS_Type( cgns_fem_type, isoparametric_, dim );
            csmp_fem_type_name = csmp_elmt_specs_.CSMP_TypeName( csmp_fem_type );
            if ( csmp_fem_type != csmp::UNKNOWN )
            {
                //cg_ElementDataSize(cgfile,cgbase,cgzone,cgsect,&element_data_size);
                nperelmt = csmp_elmt_specs_.NodesPerElementOfType( csmp_fem_type );
                element_data_size = num_elements * nperelmt;
                elements = new cgsize_t[ element_data_size ];
                if( error_handler.Verbose() )
                    std::cout<<"   reading element data for this element\n";
                cg_elements_read(cgfile,cgbase,cgzone,cgsect,elements,parent_data);
                /// filling fem types and element nodes
                std::vector<size_t> nids( nperelmt );
                size_t i = 0;
                while( i < element_data_size )
                {
                    size_t iend = i + nperelmt;
                    size_t j = 0;
                    while( i < iend )
                        nids[j++] = elements[i++] - 1;
                    eit = element_ids_.insert( std::make_pair( nids, std::make_pair(static_cast<size_t>(csmp_fem_type), 0 ) ) );
                    if ( eit.second )
                    {
                        (*eit.first).second.second = global_eid_++;
                        /// filling VSet
                        vset.ResizeElementTypes( global_eid_ );
                        vset.ResizePlist( global_eid_, nperelmt );
                        vset.ElementType( (*eit.first).second.second, csmp_fem_type );
                        for( size_t nid = 0 ; nid < nperelmt; ++nid )
                            vset.Plist( (*eit.first).second.second, nid, nids[nid] );
                    }
                    /// filling MeshTopology
                    section_element_ids.push_back( (*eit.first).second.second );
                }
                delete[] elements;
                /// adding section information to MeshTopology
                mesh_topology.AddRegionElementType( section_name.c_str(), csmp_fem_type_name );
                mesh_topology.AddRegionElementIds( section_name.c_str(), section_element_ids );
            }
            else
            {
                if( error_handler.Verbose() )
                    std::cout<<"   not reading element data for this element\n";
            }
        }
    }
}

template void CGNS_Interface::ReadElements<1U>( int,int,int,int,csmp::VSet<1U>&,csmp::ModelTopology& );
template void CGNS_Interface::ReadElements<2U>( int,int,int,int,csmp::VSet<2U>&,csmp::ModelTopology& );
template void CGNS_Interface::ReadElements<3U>( int,int,int,int,csmp::VSet<3U>&,csmp::ModelTopology& );


/**

Replaces the element identifier information in the supplied VSet into the
corresponding CSMP finite element types taking into account whether
an isoparametric model shall be created or not.

The method converts the argument VSet.
*/
template<size_t dim>
void Convert_CGNS_To_CSMP_FiniteElementTypes( csmp::VSet<dim>& vset, bool isoparametric )
 {
    CGNS_ElementSpecifications  elmt_specs;

    if ( !vset.HybridElementTypeMesh() )
      vset.ElementType( 0U, elmt_specs.CSMP_TypeFrom_CGNS_Type( vset.ElementType(0U), isoparametric, dim ) );
    else
      for ( size_t i=0U; i<vset.ElementTypes(); i++ )
        vset.ElementType( i, elmt_specs.CSMP_TypeFrom_CGNS_Type( vset.ElementType(i), isoparametric, dim ) );

 } // end

template void Convert_CGNS_To_CSMP_FiniteElementTypes( csmp::VSet<1U>&,bool );
template void Convert_CGNS_To_CSMP_FiniteElementTypes( csmp::VSet<2U>&,bool );
template void Convert_CGNS_To_CSMP_FiniteElementTypes( csmp::VSet<3U>&,bool );

/**
Writess simple 3-D unstructured grid to CGNS file
*/

template<size_t dim>
int CGNS_Interface::Write_CGNS_Mesh( const std::string&      filename,
                                     const csmp::Model<dim>& model )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    /// prepare file
    int cgfile;
    std::string outfile = filename;
    outfile += ".cgns";

    /// set file type
    if (cg_set_file_type( GetFileType() ) )
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::Write_CGNS_Mesh():","cg_set_file_type");

    /// open CGNS file for write
    unlink(outfile.c_str());
    if (cg_open(outfile.c_str(), CG_MODE_WRITE, &cgfile))
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::Write_CGNS_Mesh():","cg_open");

    /// ----------------------------------------------------------

    if( WriteUnstructMesh<dim>( filename, cgfile, model ) )
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::Write_CGNS_Mesh():","Error while writing mesh.");

    /// ----------------------------------------------------------
    /// close CGNS file

    if ( cg_close(cgfile) )
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::Write_CGNS_Mesh():","cg_close");

    return 0;

}

template int CGNS_Interface::Write_CGNS_Mesh( const std::string&,const csmp::Model<1U>& );
template int CGNS_Interface::Write_CGNS_Mesh( const std::string&,const csmp::Model<2U>& );
template int CGNS_Interface::Write_CGNS_Mesh( const std::string&,const csmp::Model<3U>& );

template<size_t dim>
int CGNS_Interface::WriteUnstructMesh( const std::string& filename, int cgfile, const csmp::Model<dim>& model )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// write base
    int cgbase;
    WriteBase<dim>(filename,cgfile,cgbase);

    /// write zone
    int cgzone;
    WriteZone<dim>(filename,cgfile,cgbase,cgzone,model);

    /// write coordinates
    WriteCoords<dim>(cgfile,cgbase,cgzone,model);

    /// write subdomains
    WriteSubDomains<dim>(cgfile,cgbase,cgzone,model);

    return 0;
}

template int CGNS_Interface::WriteUnstructMesh<1U>( const std::string&,int,const csmp::Model<1U>& );
template int CGNS_Interface::WriteUnstructMesh<2U>( const std::string&,int,const csmp::Model<2U>& );
template int CGNS_Interface::WriteUnstructMesh<3U>( const std::string&,int,const csmp::Model<3U>& );

template<size_t dim>
void CGNS_Interface::WriteBase( const std::string& filename, int cgfile, int& cgbase )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// create base (user can give any name)
    if( error_handler.Verbose() )
        std::cout<< "CGNS_Interface::Write_CGNS_Mesh(): writing unstructured base\n";
    fflush (stdout);

    /// mesh info
    int CellDim = dim;
    int PhysDim = dim;
    std::string basename = "Base";
    //std::string basename( filename.c_str() );
    //basename += "-Base";
    if (cg_base_write(cgfile, basename.c_str(), CellDim, PhysDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor", "Multi-block Unstructured Grid"))
            error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::Write_CGNS_Mesh():","unstructured base");

    /// additional data
    //        cg_dataclass_write(CGNS_ENUMV(NormalizedByDimensional)
    //        cg_units_write(CGNS_ENUMV(Kilogram),        // MassUnits_t
    //                       CGNS_ENUMV(Meter),           // LengthUnits_t
    //                       CGNS_ENUMV(Second),          // TimeUnits_t
    //                       CGNS_ENUMV(Kelvin),          // TemperatureUnits_t
    //                       CGNS_ENUMV(Radian)))         // AngleUnits_t
    //        cg_unitsfull_write(CGNS_ENUMV(Kilogram),    // MassUnits_t
    //                       CGNS_ENUMV(Meter),           // LengthUnits_t
    //                       CGNS_ENUMV(Second),          // TimeUnits_t
    //                       CGNS_ENUMV(Kelvin),          // TemperatureUnits_t
    //                       CGNS_ENUMV(Radian),          // AngleUnits_t
    //                       CGNS_ENUMT(Ampere),          // ElectricCurrentUnits_t
    //                       CGNS_ENUMT(Mole),            // SubstanceAmountUnits_t
    //                       CGNS_ENUMT(Candela)));       // LuminousIntensityUnits_t

}
template void CGNS_Interface::WriteBase<1U>( const std::string&,int,int& );
template void CGNS_Interface::WriteBase<2U>( const std::string&,int,int& );
template void CGNS_Interface::WriteBase<3U>( const std::string&,int,int& );


template<size_t dim>
void CGNS_Interface::WriteZone( const std::string& filename, int cgfile, int cgbase, int& cgzone, const csmp::Model<dim>& model )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// write zone
    cgsize_t size[3];
    for (size_t n = 0; n < 3; ++n )
        size[n] = 0;
    size[0] = model.Region("Model").Nodes();
    size[1] = model.Region("Model").Elements();
    if( model.Boundaries() > 0 )
    {
        typename csmp::Model<dim>::boundaryConstIterator bit    = model.BoundariesBegin();
        typename csmp::Model<dim>::boundaryConstIterator bitEnd = model.BoundariesEnd();
        for ( ; bit!=bitEnd; ++bit )
            size[1] += (*bit).second.Elements();
    }

    //std::string zonename = "Zone";
    std::string zonename( filename.c_str() );
    //zonename += "-Zone";
    if (cg_zone_write( cgfile, cgbase, zonename.c_str(), size, CGNS_ENUMV(Unstructured), &cgzone ) )
    {
        errmsg  = "Cannot write Zone('";
        errmsg += filename;
        errmsg += "')";
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::WriteZone():",errmsg.c_str() );
    }
}
template void CGNS_Interface::WriteZone<1U>( const std::string&,int,int,int&,const csmp::Model<1U>& );
template void CGNS_Interface::WriteZone<2U>( const std::string&,int,int,int&,const csmp::Model<2U>& );
template void CGNS_Interface::WriteZone<3U>( const std::string&,int,int,int&,const csmp::Model<3U>& );


template<size_t dim>
void CGNS_Interface::WriteCoords( int cgfile, int cgbase, int cgzone, const csmp::Model<dim>& model )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// write coordinates

    int cgcoord;

    const csmp::Region<dim>& subDomain( model.Region("Model") );
    typename csmp::Region<dim>::vertexConstIterator nit    = subDomain.NodesBegin();
    typename csmp::Region<dim>::vertexConstIterator nitEnd = subDomain.NodesEnd();
    int num_coord = subDomain.Nodes();
    float* xcoord = new float[num_coord];
    float* ycoord = new float[num_coord];
    float* zcoord = new float[num_coord];
    int nid = 0;
    if( dim == 3U )
        for ( ; nit!=nitEnd; ++nit ){
            xcoord[ nid ] = static_cast<float>( (*nit)->operator[](0U) );
            ycoord[ nid ] = static_cast<float>( (*nit)->operator[](1U) );
            zcoord[ nid ] = static_cast<float>( (*nit)->operator[](2U) );
            nid++;
        }
    else if( dim == 2U )
        for ( ; nit!=nitEnd; ++nit ){
            xcoord[ nid ] = static_cast<float>( (*nit)->operator[](0U) );
            ycoord[ nid ] = static_cast<float>( (*nit)->operator[](1U) );
            zcoord[ nid ] = 0.0;
            nid++;
        }
    else
        for ( ; nit!=nitEnd; ++nit ){
            xcoord[ nid ] = static_cast<float>( (*nit)->operator[](0U) );
            ycoord[ nid ] = 0.0;
            zcoord[ nid ] = 0.0;
            nid++;
        }

    if (cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateX", xcoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone, "GridCoordinates_t", 1,
            "CoordinateX", 0, NULL) ||
        cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateY", ycoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone, "GridCoordinates_t", 1,
            "CoordinateY", 0, NULL) ||
        cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateZ", zcoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone, "GridCoordinates_t", 1,
            "CoordinateZ", 0, NULL)) {
        errmsg += "Region('";
        errmsg += subDomain.Name();
        errmsg += "') coordinates couldn't be written.";
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::WriteCoords():",errmsg.c_str());
    }
    delete[] xcoord;
    delete[] ycoord;
    delete[] zcoord;
}
template void CGNS_Interface::WriteCoords<1U>( int,int,int,const csmp::Model<1U>& );
template void CGNS_Interface::WriteCoords<2U>( int,int,int,const csmp::Model<2U>& );
template void CGNS_Interface::WriteCoords<3U>( int,int,int,const csmp::Model<3U>& );

template<size_t dim>
void CGNS_Interface::WriteSubDomains( int cgfile, int cgbase, int cgzone, const csmp::Model<dim>& model )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// write regions
    typename csmp::Model<dim>::regionConstIterator rit    = model.UniqueRegionsBegin();
    typename csmp::Model<dim>::regionConstIterator ritEnd = model.UniqueRegionsEnd();
    std::string regionname;
    for ( ; rit!=ritEnd; ++rit ) {
        regionname = (*rit).first;
        if( regionname != "Model" )
            WriteElements<dim>(cgfile,cgbase,cgzone,(*rit).second);
    }

    /// write boundaries
    typename csmp::Model<dim>::boundaryConstIterator bit    = model.BoundariesBegin();
    typename csmp::Model<dim>::boundaryConstIterator bitEnd = model.BoundariesEnd();
    for ( ; bit!=bitEnd; ++bit )
        WriteElements<dim>(cgfile,cgbase,cgzone,(*bit).second);
}
template void CGNS_Interface::WriteSubDomains<1U>( int,int,int,const csmp::Model<1U>& );
template void CGNS_Interface::WriteSubDomains<2U>( int,int,int,const csmp::Model<2U>& );
template void CGNS_Interface::WriteSubDomains<3U>( int,int,int,const csmp::Model<3U>& );

template<size_t dim,template <size_t> class SIMPLEX>
void CGNS_Interface::WriteElements( int cgfile, int cgbase, int cgzone, const csmp::ModelSubDomain<dim,SIMPLEX>& subDomain )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );
    std::string errmsg;

    /// write elements
    int cgsect;

    typename csmp::ModelSubDomain<dim,SIMPLEX>::simplexConstIterator eit;
    typename csmp::ModelSubDomain<dim,SIMPLEX>::simplexConstIterator eitEnd = subDomain.ElementsEnd();
    int num_element = 0;
    int eid  = 0;
    int enodes;
    for ( eit = subDomain.ElementsBegin(); eit!=eitEnd; ++eit )
        num_element += ( (*eit)->Nodes() + 1 ); // nodes + element type
    cgsize_t* elements = new cgsize_t[ num_element ];
    num_element = subDomain.Elements();
    for ( eit = subDomain.ElementsBegin(); eit!=eitEnd; ++eit ){
        /// type of element
        elements[ eid++ ] = static_cast<int>(elmt_specs_.CGNS_TypeFrom_CSMP_Type( (*eit)->FE_Type() ) );
        /// element nodes
        enodes = (*eit)->Nodes();
        for ( size_t enid = 0; enid < enodes; ++enid )
            elements[ eid++ ] = (*eit)->N( enid )->Idx() + 1;
    }

    if (cg_section_write(cgfile, cgbase, cgzone, subDomain.Name().c_str(), CGNS_ENUMV(MIXED),
                         1, num_element, 0/* last boundary element id, 0 if unsorted elements*/, elements, &cgsect) ) {
        errmsg += "Region('";
        errmsg += subDomain.Name();
        errmsg += "') elements couldn't be written.";
        error_handler.notice(csmp::FATAL_ERROR, "CGNS_Interface::WriteCoordsAndElements():",errmsg.c_str() );
    }
    delete[] elements;
}

template void CGNS_Interface::WriteElements<1U>( int,int,int,const csmp::ModelSubDomain<1U,csmp::Element>& );
template void CGNS_Interface::WriteElements<2U>( int,int,int,const csmp::ModelSubDomain<2U,csmp::Element>& );
template void CGNS_Interface::WriteElements<3U>( int,int,int,const csmp::ModelSubDomain<3U,csmp::Element>& );

template void CGNS_Interface::WriteElements<1U>( int,int,int,const csmp::ModelSubDomain<1U,csmp::Face>& );
template void CGNS_Interface::WriteElements<2U>( int,int,int,const csmp::ModelSubDomain<2U,csmp::Face>& );
template void CGNS_Interface::WriteElements<3U>( int,int,int,const csmp::ModelSubDomain<3U,csmp::Face>& );



// CGNS MODEL INTERFACE

CGNS_ModelSettings::CGNS_ModelSettings( const std::string& mesh_file_prefix )
    : mesh_file_prefix_     ( mesh_file_prefix )
{
    if( csmp::isRegionsFileExist( mesh_file_prefix.c_str() ) ){
        regions_.clear();
        csmp::readDesiredRegions( mesh_file_prefix.c_str(), regions_ );
    }
}

CGNS_ModelSettings::CGNS_ModelSettings( const CGNS_ModelSettings& s )
    : mesh_file_prefix_     ( s.mesh_file_prefix_ ),
      regions_              ( s.regions_ )
{
}

CGNS_ModelSettings& CGNS_ModelSettings::operator=( const CGNS_ModelSettings& s )
{
    if( &s != this )
    {
        mesh_file_prefix_     = s.mesh_file_prefix_;
        regions_              = s.regions_;
    }
    return *this;
}

CGNS_ModelSettings::~CGNS_ModelSettings()
{
}

void CGNS_ModelSettings
::MeshSetup( const std::string& regions_file_prefix )
{
    if( csmp::isRegionsFileExist( regions_file_prefix.c_str() ) )
    {
        regions_.clear();
        csmp::readDesiredRegions( regions_file_prefix.c_str(), regions_ );
    }
}

void CGNS_ModelSettings
::MeshSetup( const std::set<std::string>& regions )
{
    regions_ = regions;
}


} // end namespace csmp
