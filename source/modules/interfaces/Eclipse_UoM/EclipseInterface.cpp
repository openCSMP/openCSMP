#include "EclipseInterface_UoM.h"

#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "ErrorHandler.h"
#include "TextFileInterface.h"
#include "vectorOperations.h"

namespace csmp {

namespace eclipse {

EclipseInterface::EclipseInterface()
: NX_(0),
  NY_(0),
  NZ_(0)
{
}


EclipseInterface::~EclipseInterface()
{
}

#if 0
void EclipseInterface::ClearBefore()
{
    /// grid specs
    NX_ = 0;
    NY_ = 0;
    NZ_ = 0;

    /// topology data
    regions_.clear();
    faults_.clear();
    wells_.clear();

    ClearAfter();
}


void EclipseInterface::ClearAfter()
{
    /// grid
    grid_.Clear();
    block_grid_.Clear();
    box_.clear();

    /// regions
    satnum_.clear();
    pvtnum_.clear();
    rocknum_.clear();
    eqlnum_.clear();
    fipnum_.clear();

    /// cell properties
    ntg_.clear();
    poro_.clear();
    permxyz_.clear();
    pressure_.clear();
    swat_.clear();
    soil_.clear();
    sgas_.clear();

    /// face properties
//    multflt_.clear();
//    multxyz_.clear();
//    tranxyz_.clear();
}


void EclipseInterface
::AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids )
{
    grid_.AddWellFacePath( well_name, cell_ids );
}


void EclipseInterface
::AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids, const std::vector<std::pair<size_t,size_t> >& face_ids )
{
    grid_.AddWellFacePath( well_name, cell_ids, face_ids );
}


void EclipseInterface
::AddWellEdgePath( const std::string& well_name, const std::vector<size_t>& cell_ids, const std::vector<std::pair<size_t,size_t> >& edge_ids )
{
    grid_.AddWellEdgePath( well_name, cell_ids, edge_ids );
}


#endif





/**
    MASTER METHOD of EclipseInterface which does everything:
    
    1. Reads sequence of files which all have the .grdecl extension
*/

bool EclipseInterface::ReadFile( csmp::VSet<3U>& vset,
                                 csmp::ModelTopology& model_topology,
                                 const std::string& fname,
                                 bool exclude_inactive_cells, bool tetra_mesh )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( csmp_error.Verbose() )
        std::cout <<"\nEclipseInterface<::ReadFile: Reading ECLIPSE ASCII file '"<<fname<<"'..."<< std::endl;

    std::ifstream  ifs;
    size_t    line_length(256);
    std::vector<std::string> file_extensions;
    file_extensions.push_back( ".GRDECL" );
    file_extensions.push_back( ".grdecl" );

    /// 0. Open fname.grdecl file
    csmp::openFile( ifs, fname, file_extensions );

    /// 1. Assign object members and reader settings
    model_name_     = fname;
    vset_           = &vset;
    model_topology_ = &model_topology;
    model_topology_->ModelName( fname.c_str() );

    // 2. Read file
#if 0
    ClearBefore();
#endif

    if( !ReadFile( ifs, line_length ) ) {
        csmp_error.notice( csmp::FATAL_ERROR, "EclipseInterface::","File could not be properly read!!!");
        return false;
    }

    if( csmp_error.Verbose() )
        std::cout <<"\nEclipseInterface<::ReadFile: file '"<< fname <<"' read successfully!"<< std::endl;

    // 3. convert grid from cell-centered to corner point format
    //    by creating new pillars in the center of the cells
#if 0
    if ( grid_.GetPillars().empty() &&
        !block_grid_.GetCellDepths().empty() &&
        !block_grid_.GetCellSizes(0).empty() &&
        !block_grid_.GetCellSizes(1).empty() &&
        !block_grid_.GetCellSizes(2).empty() )
      {
          block_grid_.AssignCellCoordinatesToPillars( grid_.GetPillars() );
          block_grid_.Clear();
      }
#endif

// DEBUGGING
#if 0
std::cerr <<"\nEclipseInterface:: printing current pillars\n";
for ( auto it=grid_.GetPillars().begin(); it!=grid_.GetPillars().end(); ++it )
  for ( auto rit=(*it).begin(); rit!=(*it).end(); ++rit )
    (*rit).Out();
#endif
  
    // KEY METHOD
    grid_.CreateModel( model_name_,*vset_, *model_topology_, zcorn_, regions_, faults_, wells_,
                       tetra_mesh, exclude_inactive_cells );

   /// 4. Write properties to VSet
    WritePropertiesToVSet();
#if 0
    ClearAfter();
#endif

    return true;
}


/**
     ECLIPSE_ROCKNUM is not recognised?
*/

void EclipseInterface::WritePropertiesToVSet()
{
    /// 1. write rock types
    if( !satnum_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   satnum_, properties_[ECLIPSE_SATNUM].name, properties_[ECLIPSE_SATNUM].key.place );

    if( !pvtnum_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   pvtnum_, properties_[ECLIPSE_PVTNUM].name, properties_[ECLIPSE_PVTNUM].key.place );

    if( !rocknum_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   rocknum_, properties_[ECLIPSE_ROCKNUM].name, properties_[ECLIPSE_ROCKNUM].key.place );

    if( !eqlnum_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   eqlnum_, properties_[ECLIPSE_EQLNUM].name, properties_[ECLIPSE_EQLNUM].key.place );

    if( !fipnum_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   fipnum_, properties_[ECLIPSE_FIPNUM].name, properties_[ECLIPSE_FIPNUM].key.place );

    /// 2. write ntg
    if( !ntg_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   ntg_, properties_[ECLIPSE_NTG].name, properties_[ECLIPSE_NTG].key.place );

    /// 3. write porosity
    if( !poro_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   poro_, properties_[ECLIPSE_PORO].name, properties_[ECLIPSE_PORO].key.place );

    /// 4. write permeability
    if( !permxyz_.empty() )
    {
        /// converting to SI system
        double conversion_factor( 1.0 );
        if( properties_[ECLIPSE_PERM].unit == "mD")
            conversion_factor = 1.0e-15;
        else if( properties_[ECLIPSE_PERM].unit == "D")
            conversion_factor = 1.0e-12;
        if( conversion_factor != 1.0 )
        {
            const size_t num_cells( permxyz_.size() );
            for( size_t i=0; i<num_cells; ++i)
                permxyz_[i] *= conversion_factor;
        }
        if( properties_[ECLIPSE_PERM].key.type != csmp::TENSOR )
            WriteScalarPropertyToVSet( *vset_, grid_,
                                       permxyz_, properties_[ECLIPSE_PERM].name, properties_[ECLIPSE_PERM].key.place );
        else
            WriteTensorPropertyToVSet( *vset_, grid_,
                                       permxyz_, properties_[ECLIPSE_PERM].name, properties_[ECLIPSE_PERM].key.place );
    }

    /// 5. write pressure
    if( !pressure_.empty() )
    {
        /// converting to SI system
        double conversion_factor( 1.0 );
        if( properties_[ECLIPSE_PRESSURE].unit == "kPa")
            conversion_factor = 1.0e3;
        else if( properties_[ECLIPSE_PRESSURE].unit == "MPa")
            conversion_factor = 1.0e6;
        else if( properties_[ECLIPSE_PRESSURE].unit == "GPa")
            conversion_factor = 1.0e9;
        else if( properties_[ECLIPSE_PRESSURE].unit == "bar")
            conversion_factor = 1.0e5;
        else if( properties_[ECLIPSE_PRESSURE].unit == "psi")
            conversion_factor = 6.8948e3;
        if( conversion_factor != 1.0 )
        {
            const size_t num_cells( pressure_.size() );
            for( size_t i=0; i<num_cells; ++i)
                pressure_[i] *= conversion_factor;
        }
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   pressure_, properties_[ECLIPSE_PRESSURE].name, properties_[ECLIPSE_PRESSURE].key.place );
    }

    /// 6. write water saturation
    if( !swat_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   swat_, properties_[ECLIPSE_SWAT].name, properties_[ECLIPSE_SWAT].key.place );

    /// 7. write oil saturation
    if( !soil_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   soil_, properties_[ECLIPSE_SOIL].name, properties_[ECLIPSE_SOIL].key.place );

    /// 8. write gas saturation
    if( !sgas_.empty() )
        WriteScalarPropertyToVSet( *vset_, grid_,
                                   sgas_, properties_[ECLIPSE_SGAS].name, properties_[ECLIPSE_SGAS].key.place );

}


void EclipseInterface::WriteScalarPropertyToVSet( csmp::VSet<3U>&  vset,
                                                   const CornerPointGrid_UoM& grid,
                                                   const std::vector<csmp::ScalarVariable>& scalar_data,
                                                   const std::string& property_name,
                                                   const csmp::PLACEMENT& place )
 {
    // Add scalar data to vset
    grid.WritePropertyToVSet(vset,scalar_data,property_name,place);
    return;
 }



void EclipseInterface::WriteScalarPropertyToVSet( csmp::VSet<3U>&  vset,
                                                   const CornerPointGrid_UoM& grid,
                                                   const std::vector<csmp::VectorVariable<3U> >& vector_data,
                                                   const std::string& property_name,
                                                   const csmp::PLACEMENT& place )
 {
    // Convert vector data to scalar data
    const size_t data_size( vector_data.size() );
    std::vector<csmp::ScalarVariable> scalar_data( data_size, csmp::ScalarVariable(csmp::PLAIN,0.0) );
    for( size_t i = 0; i < data_size; i++ )
        scalar_data[ i ] = valueAverage( vector_data[i] );
    // Add scalar data to vset
    grid.WritePropertyToVSet(vset,scalar_data,property_name,place);
    return;
 }



void EclipseInterface::WriteScalarPropertyToVSet( csmp::VSet<3U>&  vset,
                                                   const CornerPointGrid_UoM& grid,
                                                   const std::vector<csmp::TensorVariable<3U> >& tensor_data,
                                                   const std::string& property_name,
                                                   const csmp::PLACEMENT& place )
 {
    // Convert tensor data to scalar data
    const size_t data_size( tensor_data.size() );
    std::vector<csmp::ScalarVariable> scalar_data( data_size, csmp::ScalarVariable(csmp::PLAIN,0.0) );
    for( size_t i = 0; i < data_size; i++ )
        scalar_data[ i ] = tensor_data[ i ].Trace()/3.0;
    // Add scalar data to vset
    grid.WritePropertyToVSet(vset,scalar_data,property_name,place);
    return;
 }




void EclipseInterface::WriteVectorPropertyToVSet( csmp::VSet<3U>&  vset,
                                                 const CornerPointGrid_UoM& grid,
                                                 const std::vector<csmp::VectorVariable<3U> >& vector_data,
                                                 const std::string& property_name,
                                                 const csmp::PLACEMENT& place )
 {
    // Add vector data to vset
    grid.WritePropertyToVSet(vset,vector_data,property_name,place);
    return;
 }




void EclipseInterface::WriteTensorPropertyToVSet( csmp::VSet<3U>&  vset,
                                                 const CornerPointGrid_UoM& grid,
                                                 const std::vector<csmp::TensorVariable<3U> >& tensor_data,
                                                 const std::string& property_name,
                                                 const csmp::PLACEMENT& place )
 {
    // Add tensor data to vset
    grid.template WritePropertyToVSet<csmp::TensorVariable<3U> >(vset,tensor_data,property_name,place);
    return;
 }

// READING PROCESS

/**
    Configures the interface for the case where everything is supplied as separated files
*/

bool EclipseInterface::ReadFile( std::ifstream&  ifs, size_t line_length )
{
    TextFileInterface<EclipseInterface > reader( line_length );

    // keyword INCLUDE allows to include other file to be read
    reader.AddKeyword("INCLUDE");

    /// section
    reader.AddKeyword("GRID");
    reader.AddKeyword("RUNSPEC");
    reader.AddKeyword("PROPS");
    reader.AddKeyword("REGIONS");

    /// grid
    reader.AddKeyword("DIMENS");
    reader.AddKeyword("DX");
    reader.AddKeyword("DY");
    reader.AddKeyword("DZ");
    reader.AddKeyword("SPECGRID");
    reader.AddKeyword("TOPS");
    reader.AddKeyword("COORD");
    reader.AddKeyword("ZCORN");
    reader.AddKeyword("FAULTS");
//    reader.AddKeyword("BOX");
//    reader.AddKeyword("ENDBOX");

    /// regions
    reader.AddKeyword("ACTNUM");
    reader.AddKeyword("SATNUM");
    reader.AddKeyword("PVTNUM");
    reader.AddKeyword("ROCKNUM");
    reader.AddKeyword("EQLNUM");
    reader.AddKeyword("FIPNUM");

    /// cell properties
    reader.AddKeyword("NTG");
    reader.AddKeyword("PORO");
    reader.AddKeyword("PERMX");
    reader.AddKeyword("PERMY");
    reader.AddKeyword("PERMZ");
    reader.AddKeyword("PRESSURE");
    reader.AddKeyword("SWAT");
    reader.AddKeyword("SOIL");
    reader.AddKeyword("SGAS");

    /// face properties
//    reader.AddKeyword("MULTFLT");
//    reader.AddKeyword("MULTX");
//    reader.AddKeyword("MULTY");
//    reader.AddKeyword("MULTZ");
//    reader.AddKeyword("MULTX-");
//    reader.AddKeyword("MULTY-");
//    reader.AddKeyword("MULTZ-");
//    reader.AddKeyword("TRANX");
//    reader.AddKeyword("TRANY");
//    reader.AddKeyword("TRANZ");

    /// well properties
    reader.AddKeyword("WELSPECS");
    reader.AddKeyword("COMPDAT");
    reader.AddKeyword("COMPDATEF");  /// advanced well data assignment ( explicit path through faces )
    reader.AddKeyword("COMPDATEN");  /// advanced well data assignment ( explicit path through nodes )

    reader.ReadFile( ifs, this, (&isEclipseCommentLine), (&EclipseInterface::ReadBlock) );

    return true;

} // end ReadFile




bool EclipseInterface::ReadIncludeFile( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    // 0. Read file name
    std::string fname;
    if( readEclipseIncludeFileName( ifs, text_line, line_length, fname ) == 0 )
        return false;

    // 1. Open the file
    std::ifstream ifs_include;
    csmp::openFile( ifs_include, fname );

    // 2. Read file
    if( csmp_error.Verbose() )
        std::cout <<"\nEclipseInterface::ReadIncludeFile: Reading ECLIPSE ASCII file '"<<fname<<"'..."<< std::endl;
    bool result( ReadFile( ifs_include, line_length ) );
    if( csmp_error.Verbose() )
        std::cout <<"\nEclipseInterface::ReadIncludeFile: file '"<< fname <<"'read successfully!"<< std::endl;
    return result;
}




bool EclipseInterface
::ReadUnknownBlock( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( csmp_error.Verbose() )
        std::cout <<"\nEclipseInterface::ReadUnknownBlock: reading "<<keyword_<<" block of data...";

    if( skipEclipseBlock( ifs, text_line, line_length, csmp_error.Verbose() ) == 0)
        return false;

    if( csmp_error.Verbose() )
        std::cout <<"\nEclipseInterface::ReadUnknownBlock: block "<<keyword_<<" has been skipped!\n";

    return true;
}



/**

Generic function which reads blocks of data identified by keywords

@author Roman, 2014
@todo add properly BOX and ENDBOX keywords
@todo add face properties: MULTXYZ, MULTFLT, TRANXYZ

@attention SKM:  KEY FUNCTION ALWAYS EXECUTED BY INTERFACE

*/


bool EclipseInterface::ReadBlock( std::ifstream& ifs, char* text_line, size_t line_length )
{

    // read keyword and it's parameters
    if(! readEclipseKeyword( ifs, text_line, line_length,
                             keyword_, keyword_parameters_ ) )
        return false;

    // PROCESSING KEYWORDS

    /// GRID SECTION
    if( keyword_ == "DIMENS" )
    {
        return ReadDimensions( ifs, text_line, line_length );
    }
    else if( keyword_ == "DX" || keyword_ == "DY" || keyword_ == "DZ" )
    {
        return ReadCellSizes( ifs, text_line, line_length );
    }
    else if ( keyword_ == "TOPS" )
    {
        return ReadCellDepths( ifs, text_line, line_length );
    }
    if( keyword_ == "SPECGRID" )
    {
        return ReadGridSpecs( ifs, text_line, line_length );
    }
    else if( keyword_ == "COORD" )
    {
// SKM fix:        return ReadPillarCoordinates( ifs, text_line, line_length );
        return Read_COORD( ifs, text_line, line_length );
    }
    else if( keyword_ == "ZCORN" )
    {
// SKM fix        return ReadCornerDepths( ifs, text_line, line_length );
        return Read_ZCORN( ifs, text_line, line_length );
    }
    else if( keyword_ == "FAULTS" )
    {
        return ReadFaultsData( ifs, text_line,line_length );
    }

    /// REGIONS SECTION
    else if( keyword_ == "ACTNUM" )
    {
        return ReadActiveCells( ifs, text_line,line_length );
    }
    else if( keyword_ == "SATNUM" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, satnum_ );
    }
    else if( keyword_ == "PVTNUM" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, pvtnum_ );
    }
    else if( keyword_ == "ROCKNUM" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, rocknum_ );
    }
    else if( keyword_ == "EQLNUM" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, eqlnum_ );
    }
    else if( keyword_ == "FIPNUM" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, fipnum_ );
    }

    /// CELL PROPERTIES
    else if( keyword_ == "NTG" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, ntg_ );
    }
    else if( keyword_ == "PORO" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, poro_ );
    }
    else if( keyword_ == "PERMX" || keyword_ == "PERMY" || keyword_ == "PERMZ" )
    {
        return ReadTensorProperty( ifs, text_line,line_length, permxyz_, "PERMX", "PERMY", "PERMZ" );
    }
    else if( keyword_ == "PRESSURE" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, pressure_ );
    }
    else if( keyword_ == "SWAT" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, swat_ );
    }
    else if( keyword_ == "SOIL" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, soil_ );
    }
    else if( keyword_ == "SGAS" )
    {
        return ReadScalarProperty( ifs, text_line,line_length, sgas_ );
    }

    /// FACE PROPERTIES
//    else if( keyword_ == "MULTFLT" )
//    {
//        return ReadFaultTransmissibilityMultipliers( ifs, text_line,line_length );
//    }
//    else if(   keyword_ == "MULTX" || keyword_ == "MULTY" || keyword_ == "MULTZ"
//            || keyword_ == "MULTX-" || keyword_ == "MULTY-" || keyword_ == "MULTZ-" )
//    {
//    }
//    else if(   keyword_ == "TRANX" || keyword_ == "TRANY" || keyword_ == "TRANZ" )
//    {
//    }

    /// WELL PROPERTIES
    else if( keyword_ == "WELSPECS" )
    {
        return ReadWellSpecs( ifs, text_line,line_length );
    }
    else if( keyword_ == "COMPDAT" )
    {
        return ReadWellCompletionsData( ifs, text_line,line_length );
    }
    else if( keyword_ == "COMPDATEF" )
    {
        return ReadExplicitFaceWellCompletionsData( ifs, text_line,line_length );
    }
    else if( keyword_ == "COMPDATEN" )
    {
        return ReadExplicitNodeWellCompletionsData( ifs, text_line,line_length );
    }

    /// INCLUDE EXTERNAL FILE
    else if( keyword_ == "INCLUDE" )
    {
        return ReadIncludeFile( ifs, text_line,line_length );
    }
    else if( ( keyword_ == "GRID" ) || ( keyword_ == "RUNSPEC" ) || ( keyword_ == "PROPS" ) || ( keyword_ == "REGIONS" ) )
    {
        return true;
    }
    else
    {
        return ReadUnknownBlock( ifs, text_line, line_length );
    }

    return true;
}






bool EclipseInterface::ReadDimensions( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( readEclipseDimensions( NX_,NY_,NZ_,
                               ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
        return false;

    AssignGridDimensions();

    return true;
}





bool EclipseInterface::ReadGridSpecs( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( readEclipseGridSpecs( NX_,NY_,NZ_,
                              ifs, text_line, line_length, csmp_error.Verbose()) == 0 )
        return false;

    AssignGridDimensions();

    return true;
}


bool EclipseInterface::ReadCellSizes( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t direction( ( keyword_ == "DX") ? 0U : ( keyword_ == "DY") ? 1U : 2U );
    if( readEclipseCellData( block_grid_.GetCellSizes( direction ),
                             ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
        return false;
    return true;
}


bool EclipseInterface::ReadCellDepths( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( readEclipseCellData( block_grid_.GetCellDepths(),
                             ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
        return false;
    return true;
}


bool EclipseInterface::ReadPillarCoordinates( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if ( readEclipsePillarCoordinates( NX_,NY_,
                                       grid_,
                                       ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
        return false;
    AssignGridDimensions();
    return true;
}



bool EclipseInterface::ReadCornerDepths( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if ( readEclipseCornerDepths( NX_,NY_,NZ_,
                                  grid_,
                                  ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    AssignGridDimensions();
    return true;
}




bool EclipseInterface::ReadActiveCells( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseActiveCells( grid_.GetCellActivity(),
                                 ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}




bool EclipseInterface::ReadWellSpecs( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if ( readEclipseWellSpecs( well_data_,
                               ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}



bool EclipseInterface::ReadWellCompletionsData( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseWellCompletionsData( NX_,NY_,NZ_,
                                         well_data_,
                                         well_face_path_,
                                         ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}



bool EclipseInterface::ReadExplicitFaceWellCompletionsData( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseWellCompletionsData( NX_,NY_,NZ_,
                                         well_data_,
                                         well_face_path_,
                                         ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}



bool EclipseInterface::ReadExplicitNodeWellCompletionsData( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseWellCompletionsData( NX_,NY_,NZ_,
                                         well_data_,
                                         well_face_path_,
                                         ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}



bool EclipseInterface::ReadFaultsData( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseFaultData( NX_,NY_,NZ_,
                               faults_data_,
                               ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}




bool EclipseInterface::ReadFaultTransmissibilityMultipliers( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseFaultTransmissibilityMultipliers( NX_,NY_,NZ_,
                                                      multflt_,
                                                      ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}


bool EclipseInterface
::ReadBoxData( std::ifstream& ifs, char* text_line, size_t line_length )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
  std::vector<size_t> box_data;
    if ( readEclipseBoxData( box_data,
                             ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
       return false;
    return true;
}


bool EclipseInterface::ReadScalarProperty( std::ifstream& ifs, char* text_line, size_t line_length,
                      std::vector<csmp::ScalarVariable>& prop_data )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( readEclipseCellData( prop_data,
                              ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
        return false;
    return true;
}



bool EclipseInterface::ReadTensorProperty( std::ifstream& ifs, char* text_line, size_t line_length,
                                          std::vector<csmp::TensorVariable<3U> >& prop_data,
                                          const std::string& compx_name, const std::string& compy_name, const std::string& compz_name )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    std::vector<csmp::ScalarVariable> prop_x_or_y_or_z;
    if( readEclipseCellData( prop_x_or_y_or_z,
                             ifs, text_line, line_length, csmp_error.Verbose() ) == 0 )
        return false;
    SaveTensorProperty( ( ( keyword_ == compx_name) ? 0U : ( keyword_ == compy_name) ? 1U : 2U ),
                            prop_x_or_y_or_z,
                            prop_data );
    prop_x_or_y_or_z.clear();
    return true;
}



void EclipseInterface::SaveVectorProperty( size_t component,
                                          const std::vector<csmp::ScalarVariable>& scalar_data ,
                                          std::vector<csmp::VectorVariable<3U> >&  vector_data )
 {
    const size_t data_size( scalar_data.size() );
    if( component == 0 )
        vector_data.resize( data_size, csmp::VectorVariable<3U>( csmp::PLAIN, 0.0 ) );
    for( size_t i = 0; i < data_size; i++ )
        for( size_t j = component; j < 3U; j++ )
            vector_data[ i ]( j ) = scalar_data[ i ]();
 }



void EclipseInterface::SaveTensorProperty( size_t component,
                                          const std::vector<csmp::ScalarVariable>& scalar_data ,
                                          std::vector<csmp::TensorVariable<3U> >&  tensor_data )
 {
    const size_t data_size( scalar_data.size() );
    if( component == 0 )
        tensor_data.resize( data_size, csmp::TensorVariable<3U>( csmp::PLAIN, 0.0 ) );
    for( size_t i = 0; i < data_size; i++ )
        for( size_t j = component; j < 3U; j++ )
            tensor_data[ i ]( j, j ) = scalar_data[ i ]();
 }

void EclipseInterface::AssignGridDimensions()
{
    grid_.AssignDimensions( NX_, NY_, NZ_ );
#if 0
    block_grid_.AssignDimensionX( NX_ );
    block_grid_.AssignDimensionY( NY_ );
    block_grid_.AssignDimensionZ( NZ_ );
#endif
}



int readEclipseDimensions( size_t& NX, size_t& NY, size_t& NZ,
                           std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                         )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseDimensions: reading dimensions NX,NY,NZ...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    bool        endOfblock(false);
    size_t      counter( 0U );
    size_t      num;
    size_t      value;
    bool        default_value;
    std::vector<size_t> values;

    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            token = strtok( text_line, delims );
            do{
                if ( !readEclipseValue<size_t>(std::string(token),default_value,num,value) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseDimension:",
                                      "value cannot be read!");
                    return 0;
                }
                else
                {
                    for( size_t i = 0; i < num; i++ )
                    {
                        values.push_back( value );
                        counter++;
                    }
                }
                token = strtok( NULL, delims );
                endOfblock = isEclipseEndOfBlock(token);
            }
            while( ( counter < 3U ) && (token != NULL ) && !endOfblock );
        }
        if( counter < 3U ){
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
            if( !endOfblock && csmp::isBlankLine(text_line) )
            {
                ifs.getline( text_line, line_length );
                endOfblock = isEclipseEndOfBlock(text_line);
            }
        }
    }
    while ( counter< 3U && !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    /// assign dimensions
    if( values.size() > 3U ){
        csmp_error.notice(csmp::ERROR,
                          "readEclipseDimension:",
                          "found more than 3 dimensions!");
        return 0;
    }

    NX = values[0U];
    NY = values[1U];
    NZ = values[2U];

    if( verbose )
        std::cout <<"\nreadEclipseDimensions: dimensions have been read successfully.\n";

    return 1;
}

int readEclipseGridSpecs( size_t& NX, size_t& NY, size_t& NZ,
                          std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                        )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseGridSpecs: reading dimensions NX,NY,NZ...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    bool        endOfblock(false);
    size_t      counter( 0U );
    size_t      num;
    size_t      value;
    bool        default_value;
    std::vector<size_t> values;

    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            token = strtok( text_line, delims );
            do{
                if ( !readEclipseValue<size_t>(std::string(token),default_value,num,value) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseGridSpecs:",
                                      "value cannot be read!");
                    return 0;
                }
                else
                {
                    for( size_t i = 0; i < num; i++ )
                    {
                        values.push_back( value );
                        counter++;
                    }
                }
                token = strtok( NULL, delims );
                endOfblock = isEclipseEndOfBlock(token);
            }
            while( ( counter < 3U ) && ( token != NULL ) && !endOfblock );
        }
        if( counter < 3U ){
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
            if( !endOfblock && csmp::isBlankLine(text_line) )
            {
                ifs.getline( text_line, line_length );
                endOfblock = isEclipseEndOfBlock(text_line);
            }
        }
    }
    while ( counter< 3U && !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    /// assign dimensions
    if( values.size() > 3U ){
        csmp_error.notice(csmp::ERROR,
                              "readEclipseGridSpecs:",
                              "found more than 3 dimensions!");
        return 0;
    }

    NX = values[0U];
    NY = values[1U];
    NZ = values[2U];

    if ( verbose ) {
         std::cout <<"\nreadEclipseGridSpecs: dimensions have been read successfully.";
         std::cout <<"(NX="<< NX <<", NY="<< NY <<", NZ="<< NZ <<").\n";
      }

    return 1;
}


/**
    Reading COORD block.
    Reading top and bottom point coordinates of pillars ( coordinate lines ).
    Format:
    COORD
    x(1,1)top   y(1,1)top   z(1,1)top   x(1,1)btm   y(1,1)btm   z(1,1)btm
    x(2,1)top   y(2,1)top   z(2,1)top   x(2,1)btm   y(2,1)btm   z(2,1)btm
    x(3,1)top   y(3,1)top   z(3,1)top   x(3,1)btm   y(3,1)btm   z(3,1)btm
    .
    .
    .
    x(NX+1,NY+1)top   y(NX+1,NY+1)top   z(NX+1,NY+1)top
    x(NX+1,NY+1)btm   y(NX+1,NY+1)btm   z(NX+1,NY+1)btm\
*/
int readEclipsePillarCoordinates( size_t& NX, size_t& NY,
                                  CornerPointGrid_UoM& grid, /* matrix of Pillars=cells? */
                                  std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose ) std::cout <<"\nReadEclipsePillarCoordinates: reading COORD...";

    char* token(0);
    const char*  delims =" ,:,\t,\n,\r";
    size_t idx(0);
    size_t idy(0);
    size_t cycle(0);
    size_t position(0);
    size_t direction(0);
  csmp::Point<3u> pillar_coord[2];

    std::string message;
    bool same_line( false );
    bool endOfblock(false);
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    grid.Resize(NX, NY);
  
    do {
        same_line = false;
        do{
            if( !endOfblock && !isEclipseCommentLine( text_line ) )
            {
                if( !same_line )
                    token = strtok( text_line, delims );

                position  = (cycle/3);
                direction = cycle%3;

                /// initialize coordinate point
                if( direction == 0 )
                    pillar_coord[position] = 0.;

                if( token == NULL )
                {
                    message  = "Cannot read ";
                    message += ( direction == 0 ? "x" : ( direction == 1 ? "y" : "z" ) );
                    message += "(";
                    message  = ( position == 0 ? "top" : "bottom" );
                    message += ") ";
                    message += "value!";
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipsePillarCoordinates:",
                                      message.c_str() );
                    return 0;
                }
                else if( !isRealNumber(token) )
                {
                    message  = ( direction == 0 ? "x" : ( direction == 1 ? "y" : "z" ) );
                    message += "(";
                    message  = ( position == 0 ? "top" : "bottom" );
                    message += ") ";
                    message += "value is not a digit!";
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipsePillarCoordinates:",
                                      message.c_str() );
                    return 0;
                }
                else
                {
                    assert( position < 2 );
                    assert( direction < 3 ); // spatial dimensions
                    pillar_coord[position][direction] = atof(token);
                }

                /// add new pillar
                if( position == 1 && direction == 2 )
                {
                  Pillar pillar(pillar_coord[0], pillar_coord[1]);
                  grid(idx, idy) = std::move(pillar);
                    assert(idy < NY);
                    ++idx;
                    cycle = 0;
                }
                else
                    ++cycle;

                token = strtok( NULL, delims );
                endOfblock = isEclipseEndOfBlock(token);
                if( !endOfblock )
                {
                    if( token != NULL )
                        same_line = true;
                    else
                    {
                        same_line = false;
                        ifs.getline( text_line, line_length );
                        endOfblock = isEclipseEndOfBlock(text_line);
                    }
                }
            }
            else if( !endOfblock  )
            {
                same_line = false;
                ifs.getline( text_line, line_length );
            }
        }
        while ( !endOfblock && !csmp::isBlankLine(text_line) );

        if( !endOfblock )
        {
            ++idy;
            idx = 0;
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
        }
    }
//    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );
    while ( !endOfblock && !ifs.eof() );

    /// set grid specs
#if 0
    NY = pillars.size() - 1;
    NX = 0;
    for( size_t i = 0; i < NY; i++ )
        if( pillars[i].size() > NX )
            NX = pillars[i].size();
    NX -= 1;
#endif
    if( verbose ) std::cout <<"\nreadEclipsePillarCoordinates: COORD block has been read successfully.\n";

    return 1;
}




/**
    Reading ZCORN block (indices of the pillars).
 
    @attention Recall that X points to the east, Y points to the south, Z point downwards.
    The coordinates of cell corners are written layer by layer and
    can be denoted by NW,NE,SW,SE.

                   ______________X ( EAST )
                  |\
                  | \    NW --- NE
                  |  \    \      \
                  |   \    SW --- SE
                  |    \
    ( DOWNWARD )Z |     \Y ( SOUTH )

    Format:
    ZCONR
    z(1,1,1)top,NW    z(1,1,1)top,NE    z(2,1,1)top,NW    z(2,1,1)top,NE   ... z(NX,1,1)top,NW    z(NX,1,1)top,NE "= cycle 0 with idy =0"
    z(1,1,1)top,SW    z(1,1,1)top,SE    z(2,1,1)top,SW    z(2,1,1)top,SE   ... z(NX,1,1)top,SW    z(NX,1,1)top,SE "= cycle 1 with idx =0"
    z(1,2,1)top,NW    z(1,2,1)top,NE    z(2,2,1)top,NW    z(2,2,1)top,NE   ... z(NX,2,1)top,NW    z(NX,2,1)top,NE "= cycle 0 with idy =1"
    z(1,2,1)top,SW    z(1,2,1)top,SE    z(2,2,1)top,SW    z(2,2,1)top,SE   ... z(NX,2,1)top,SW    z(NX,2,1)top,SE "= cycle 1 with idx =1"
    .
    .
    z(1,NY,1)top,SW   z(1,NY,1)top,SE   z(2,NY,1)top,SW   z(2,NY,1)top,SE  ... z(NX,NY,1)top,SW   z(NX,NY,1)top,SE
    z(1,1,1)btm,NW    z(1,1,1)btm,NE    z(2,1,1)btm,NW    z(2,1,1)btm,NE   ... z(NX,1,1)btm,NW    z(NX,1,1)btm,NE
    .
    .
    z(1,NY,1)btm,SW   z(1,NY,1)btm,SE   z(2,NY,1)btm,SW   z(2,NY,1)btm,SE  ... z(NX,NY,1)btm,SW   z(NX,NY,1)btm,SE
    z(1,1,2)top,NW    z(1,1,2)top,NE    z(2,1,2)top,NW    z(2,1,2)top,NE   ... z(NX,1,2)top,NW    z(NX,1,2)top,NE
    .
    .
    z(1,NY,NZ)btm,SW  z(1,NY,NZ)btm,SE  z(2,NY,NZ)btm,SW  z(2,NY,NZ)btm,SE ... z(NX,NY,NZ)btm,SW  z(NX,NY,NZ)btm,SE
*/
int readEclipseCornerDepths( size_t NX, size_t NY, size_t& NZ,
                             CornerPointGrid_UoM& grid,
                             std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    std::string message;

    if( verbose ) std::cout <<"\nreadEclipseCornerDepths: reading ZCORN...";

    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";


    bool endOfblock(false);
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    /// 4 cycles:
    /// cycle 0: top NW-NE
    /// cycle 1: top SW-SE
    /// cycle 2: btm NW-NE
    /// cycle 3: btm SW-SE
    size_t cycle( 0 );
    /// 8 points of hexahedron ( TNW, TNE, TSW, TSE, BNW, BNE, BSE, BSW )
    csmp::Point<3U> pt(0.0);
    size_t idx(0), idy(0), idz(0);
    size_t pidx(0), pidy(0);
    size_t east_west_position(0);
    bool same_line( false );
    do {
        same_line = false;
        do {
            if( csmp::isBlankLine(text_line) )
            {
                ifs.getline( text_line, line_length );
 //std::cerr << text_line <<"\n";
                endOfblock = isEclipseEndOfBlock(text_line);
            }
            if( !endOfblock && !isEclipseCommentLine( text_line ) )
            {
                pidy  = idy;
                pidy += ( cycle%2 );
                pidx  = idx;
                pidx += east_west_position;

                if( !same_line )
                    token = strtok( text_line, delims );
                if( token == NULL )
                {
                    message  = "Cannot read ";
                    message += ( cycle/2 == 0 ? "top" : "bottom" );
                    message += " ";
                    message += ( cycle%2 == 0 ? "north" : "south" );
                    message += "-";
                    message += ( east_west_position == 0 ? "east" : "west" );
                    message += " ";
                    message += "value!";
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseCornerDepths:",
                                      message.c_str() );
                    return 0;
                }
                else if( !isRealNumber(token) )
                {
                    message = ( cycle/2 == 0 ? "top" : "bottom" );
                    message += " ";
                    message += ( cycle%2 == 0 ? "north" : "south" );
                    message += "-";
                    message += ( east_west_position == 0 ? "east" : "west" );
                    message += " ";
                    message += "value is not a digit!";
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseCornerDepths:",
                                      message.c_str() );
                    return 0;
                }
                else
                {
                    double64 z = atof(token);
                    assert( pidy < NY );
                    assert( pidx < NX );
                    grid(pidx,pidy).AddZCoord( z );
                }

                /// renew east-west position
                if( east_west_position == 1 )
                {
                    ++idx;
                    east_west_position = 0;
                }
                else
                    east_west_position = 1;

                token = strtok( NULL, delims );
                endOfblock = isEclipseEndOfBlock(token);
                if( !endOfblock )
                {
                    if( token != NULL )
                        same_line = true;
                    else
                    {
                        same_line = false;
                        ifs.getline( text_line, line_length );
                        endOfblock = isEclipseEndOfBlock(text_line);
                    }
                }
            }
            else if( !endOfblock  )
            {
                same_line = false;
                ifs.getline( text_line, line_length );
//std::cerr << text_line <<"\n";
             }
        }
        while ( idx < NX && !endOfblock );

        idx = 0;
        east_west_position = 0;

        /// northern points
        if( cycle%2 == 0 )
        {
            ++cycle;
        }
        /// southern points
        else
        {
            ++idy;
            if( idy < NY && !endOfblock )
            {
                --cycle;
            }
            else
            {
                idy = 0;
                cycle = (++cycle)%4;

                if( cycle == 0 && !endOfblock )
                {
                    ++idz;
                    if( csmp::isBlankLine(text_line) )
                    {
                        ifs.getline( text_line, line_length );
                        endOfblock = isEclipseEndOfBlock(text_line);
//std::cerr << text_line <<"\n";
                     }
                }
                else if( cycle == 0 && endOfblock )
                {
                    ++idz;
                }
                else if( endOfblock )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseCornerDepths:",
                                      "Unexpected end!!! Not all the blocks were read!" );
                    return 0;
                }
            }
        }
    }
    while ( !endOfblock && !ifs.eof() );



// DEBUGGING
//std::cerr <<"\nEclipseInterface::readEclipseCornerDepth: printing current pillars\n";
//for ( auto it=pillars.begin(); it!=pillars.end(); ++it )
//  for ( auto rit=(*it).begin(); rit!=(*it).end(); ++rit )
//    (*rit).Out();

    /// set grid specs
    NZ = idz;
    const size_t num_points_in_z_direction_internal_pillars( 8*NZ );
    const size_t num_points_in_z_direction_boundary_pillars( 4*NZ );
    const size_t num_points_in_z_direction_corner_pillars( 2*NZ );
    for( size_t i=1; i<NX; ++i )
    {
        for( size_t j=1; j<NY; ++j )
        {
            if( grid(i,j).GetNumPoints() != num_points_in_z_direction_internal_pillars )
            {
                std::stringstream msg;
                msg << "ZCONR were not read correctly for internal pillars!!! ";
                msg << grid(i,j).GetNumPoints();
                msg << " point instead of expected ";
                msg << num_points_in_z_direction_internal_pillars;
                msg << " points!!!";
                csmp_error.notice( csmp::ERROR,
                                   "readEclipseCornerDepths:",
                                   msg.str() );
                return 0;
            }
        }
    }
    for( size_t i=1; i<NX; ++i )
    {
        if( grid(i,0).GetNumPoints()  != num_points_in_z_direction_boundary_pillars ||
            grid(i,NY).GetNumPoints() != num_points_in_z_direction_boundary_pillars    )
        {
            std::stringstream msg;
            msg << "ZCONR were not readed correctly for boundary pillars!!! ";
            msg << grid(i,0).GetNumPoints() << " , " << grid(i,NY).GetNumPoints();
            msg << " point instead of expected ";
            msg << num_points_in_z_direction_boundary_pillars;
            msg << " points!!!";
            csmp_error.notice( csmp::ERROR,
                               "readEclipseCornerDepths:",
                               msg.str() );
            return 0;
        }
    }
    for( size_t j=1; j<NY; ++j )
    {
        if( grid(0,j).GetNumPoints()  != num_points_in_z_direction_boundary_pillars ||
            grid(NX,j).GetNumPoints() != num_points_in_z_direction_boundary_pillars    )
        {
            std::stringstream msg;
            msg << "ZCONR were not readed correctly for boundary pillars!!! ";
            msg << grid(0,j).GetNumPoints() << " , " << grid(NX,j).GetNumPoints();
            msg << " point instead of expected ";
            msg << num_points_in_z_direction_boundary_pillars;
            msg << " points!!!";
            csmp_error.notice( csmp::ERROR,
                               "readEclipseCornerDepths:",
                               msg.str() );
            return 0;
        }
    }
    if( grid(0,0).GetNumPoints()   != num_points_in_z_direction_corner_pillars ||
        grid(NX,0).GetNumPoints()  != num_points_in_z_direction_corner_pillars ||
        grid(0,NY).GetNumPoints()  != num_points_in_z_direction_corner_pillars ||
        grid(NX,NY).GetNumPoints() != num_points_in_z_direction_corner_pillars    )
    {
        std::stringstream msg;
        msg << "ZCONR were not readed correctly for corner pillars!!! ";
        msg << grid(0,0).GetNumPoints() << " , " << grid(NX,0).GetNumPoints() << grid(0,NY).GetNumPoints() << " , " << grid(NX,NY).GetNumPoints();
        msg << " point instead of expected ";
        msg << num_points_in_z_direction_corner_pillars;
        msg << " points!!!";
        csmp_error.notice( csmp::ERROR,
                           "readEclipseCornerDepths:",
                           msg.str() );
        return 0;
    }

    if( verbose )
        std::cout <<"\nreadEclipseCornerDepths: ZCORN block has been read successfully.\n";

    return 1;
}


/**
      Reading ACTNUM block.
      1 - active cell, 0 -inactive cell
      Format:
      1 0 1 0 0 0 0 0 1 1 1 1 1 ...\
 */
int readEclipseActiveCells( std::vector<uint8_t>& cell_activity, std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if ( verbose ) std::cout <<"\nreadEclipseActiveCells: reading ACTNUM...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    bool        endOfblock(false);
    size_t      active_num( 0U );
    bool        default_value;
    size_t      num;

    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    cell_activity.clear();
    size_t active_cells(0U);
  
    do {
        if ( !isEclipseCommentLine( text_line ) )
        {
            token = strtok( text_line, delims );
            do{
                if ( !readEclipseValue<size_t>(std::string(token),default_value,num,active_num) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                          "readEclipseActiveCells:",
                                          "index value cannot be read!");
                    return 0;
                }
                else
                {
                    if( active_num == 0 || active_num == 1 )
                    {
                        if ( active_num == 1 ) active_cells++;
                        for( size_t i = 0; i < num; i++ )
                            cell_activity.push_back( active_num );
                    }
                    else
                    {
                        std::cout <<"\n"<< token << std::endl;
                        csmp_error.notice(csmp::ERROR,
                                          "readEclipseActiveCells:",
                                          "index value is not valid! Should be 0 or 1.");
                        return 0;
                    }
                }
                token = strtok( NULL, delims );
                endOfblock = isEclipseEndOfBlock(token);
            }
            while( ( token != NULL ) && !endOfblock );
        }
        if( !endOfblock )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
            if( !endOfblock && csmp::isBlankLine(text_line) )
            {
                ifs.getline( text_line, line_length );
                endOfblock = isEclipseEndOfBlock(text_line);
            }
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );
  
    if ( verbose )
      std::cout <<"\nreadEclipseActiveCells: ACTNUM block with "<< cell_activity.size() <<" cells ("<< active_cells <<"=active) has been read successfully.\n";

    return 1;
}


/**
    Reads scalar property values from Eclipse file.
*/
int readEclipseCellData( std::vector<csmp::ScalarVariable>& values,
                         std::ifstream& ifs, char* text_line,
                         size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseCellData: reading cell data...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    bool        endOfblock(false);
    size_t      counter   ( 0U );
    size_t      num;
    double      value(-999.);
    bool        default_value;

    values.clear();
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            token = strtok( text_line, delims );
            do{
                if ( !readEclipseValue<double>(std::string(token),default_value,num,value) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseCellData:",
                                      "value cannot be read!");
                    return 0;
                }
                else
                {
                    for( size_t i = 0; i < num; i++ )
                    {
                        values.push_back( csmp::ScalarVariable(csmp::PLAIN,value) );
                        counter++;
                    }
                }
                token = strtok( NULL, delims );
                endOfblock = isEclipseEndOfBlock(token);
            }
            while( ( token != NULL ) && !endOfblock );
        }
        if( !endOfblock )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
            if( !endOfblock && csmp::isBlankLine(text_line) )
            {
                ifs.getline( text_line, line_length );
                endOfblock = isEclipseEndOfBlock(text_line);
            }
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        std::cout <<"\nreadEclipseCellData: cell data has been read successfully.\n";

    return 1;
}




/**
    Reads the definitions of the wells from file.
    
    Well head, bottom hole pressure etc.
*/
int readEclipseWellSpecs( std::map<std::string,EclipseWell>& well_data,
                          std::ifstream& ifs, char* text_line,
                          size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseWellSpecs: reading well specs...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    const char* symbolsToremove ="'";
    bool        endOfblock(false);

    /// temp data
    std::vector<size_t> indices( 2, 0 );
    std::string  name;
    double value;

    /// phase map
    std::map<std::string,size_t> phase;
    phase.insert( std::make_pair( "OIL", ECLIPSE_OIL ) );
    phase.insert( std::make_pair( "WATER", ECLIPSE_WATER ) );
    phase.insert( std::make_pair( "GAS", ECLIPSE_GAS ) );
    phase.insert( std::make_pair( "LIQ", ECLIPSE_LIQ ) );

    /// read first line
    well_data.clear();
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            EclipseWell well;

            /// 1. read well name
            token = strtok( text_line, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseWellSpecs:",
                                  "cannot read well name!");
                return 0;
            }
            else
            {
                name = token;
                removeSymbolsFromString( name, symbolsToremove );
                well.well_name_ = name;
            }

            /// 2. read group name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseWellSpecs:",
                                  "cannot read group name!");
                return 0;
            }
            else
            {
                name = token;
                removeSymbolsFromString( name, symbolsToremove );
                well.group_name_ = name;
            }

            /// 2. read i,j indices
            for( size_t i = 0; i<2; ++i )
            {
                token = strtok( NULL, delims );
                if ( token == NULL )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseWellSpecs:",
                                      "cannot read index!");
                    return 0;
                }
                else if( !isIntegerNumber(token) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseWellSpecs:",
                                      "index is not a digit!");
                    return 0;
                }
                else
                {
                    indices[ i ] = atoi(token);
                    /// start indexing form 0
                    indices[ i ] -= 1;
                }
            }
            well.i_start_ = indices[0];
            well.j_start_ = indices[1];

            /// 3. read bottom hole reference pressure
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cout <<"\n"<< token << std::endl;
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseWellSpecs:",
                                  "cannot read reference pressure value!");
                return 0;
            }
            else if( !isRealNumber(token) )
            {
                std::cout <<"\n"<< token << std::endl;
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseWellSpecs:",
                                  "reference pressure is not a digit!");
                return 0;
            }
            else
            {
                value = atof(token);
                well.z_bhp_ = value;
            }

            /// 4. read well phase
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cout <<"\n"<< token << std::endl;
                csmp_error.notice( csmp::ERROR,
                                   "readEclipseWellSpecs:",
                                   "cannot read phase name!");
                return 0;
            }
            else
            {
                name = token;
                removeSymbolsFromString( name, symbolsToremove );
                well.phase_ = static_cast<ECLIPSE_PHASE>( phase[ name ] );
            }

            /// 5. assign well data
            well_data[ well.well_name_ ] = well;
        }
        ifs.getline( text_line, line_length );
        endOfblock = isEclipseEndOfBlock(text_line);
        if( !endOfblock && csmp::isBlankLine(text_line) )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        std::cout <<"\nreadEclipseWellSpecs: well specifications have been read successfully.\n";

    return 1;
}



/**
     Following the reading of the well specifications, this function reads 
     information about the completion.
     
     Top and bottom of completion; open or closed etc.
*/
int readEclipseWellCompletionsData( size_t NX, size_t NY, size_t NZ,
                                    std::map<std::string,EclipseWell>& well_data,
                                    std::map<std::string,EclipseWellPath>& well_path,
                                    std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseWellCompletionsData: reading well completions data...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    const char* symbolsToremove ="'";
    bool        endOfblock(false);

    /// temp data
    const size_t NX_x_NY( NX * NY );
    std::vector<size_t> indices( 4, 0 );
    std::map<int,CORNER_POINT_CELL_FACE_INDEX> face_map;
    face_map.emplace( -1,       CORNER_POINT_CELL_FACE_Xminus );
    face_map.emplace( +1,       CORNER_POINT_CELL_FACE_Xplus  );
    face_map.emplace( -NX,      CORNER_POINT_CELL_FACE_Yminus );
    face_map.emplace( +NX,      CORNER_POINT_CELL_FACE_Yplus  );
    face_map.emplace( -NX_x_NY, CORNER_POINT_CELL_FACE_Zminus );
    face_map.emplace( +NX_x_NY, CORNER_POINT_CELL_FACE_Zplus  );

    std::string well_name;

    /// read first line
    well_data.clear();
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            /// 1. read well name
            token = strtok( text_line, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseWellCompletionsData:",
                                  "cannot read well name!");
                return 0;
            }
            else
            {
                well_name = token;
                removeSymbolsFromString( well_name, symbolsToremove );
            }

            EclipseWell& well( well_data[ well_name ] );
            EclipseWellCompletion wellcomp;

            /// 2. read i,j,k_top,k_bot indices
            for( size_t i = 0; i<4; ++i )
            {
                token = strtok( NULL, delims );
                if ( token == NULL )
                {
                     csmp_error.notice(csmp::ERROR,
                                      "readEclipseWellCompletionsData:",
                                      "cannot read index!");
                    return 0;
                }
                else if( !isIntegerNumber(token) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseWellCompletionsData:",
                                      "index is not a digit!");
                    return 0;
                }
                else
                {
                    indices[ i ] = atoi(token);
                    /// start indexing form 0
                    indices[ i ] -= 1;
                }
            }
            wellcomp.ic_    = indices[0];
            wellcomp.jc_    = indices[1];
            wellcomp.k_top_ = indices[2];
            wellcomp.k_bot_ = indices[3];

            /// 3. assign well path
            std::vector<ijk> cell_ids;
            EclipseWellPath wpath;
            for( size_t k = wellcomp.k_top_; k <= wellcomp.k_bot_; ++k )
            {
                cell_ids.emplace_back( wellcomp.ic_, wellcomp.jc_, k );
                well.well_data_.push_back( wellcomp );
            }
            addWellPath( NX, NY, NZ, well_name, cell_ids, well_path );

            /// 4. correct well path
            size_t num_cells( wellcomp.k_bot_ - wellcomp.k_top_ + 1 );
            size_t num_wells( wpath.path.size() );
            size_t wid( num_wells - num_cells );
            for( ; wid < num_wells; ++wid )
            {
                well.well_data_[wid].w_start_ = wpath.path[wid].from;
                well.well_data_[wid].w_end_   = wpath.path[wid].to;
            }
            well_path.emplace( well_name, wpath );
        }
        ifs.getline( text_line, line_length );
        endOfblock = isEclipseEndOfBlock(text_line);
        if( !endOfblock && csmp::isBlankLine(text_line) )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        std::cout <<"\nreadEclipseWellCompletionsData: well completions data has been read successfully.\n";

    return 1;
}


/**
*/
int readEclipseFaultData( size_t NX, size_t NY, size_t NZ,
                          std::map<std::string,EclipseFault>& fault_data,
                          std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseFaultData: reading fault data...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    const char* symbolsToremove ="'";
    bool        endOfblock(false);

    /// fault data
    std::vector<std::pair<size_t,size_t> > data;
    std::vector<size_t> index_range( 6, 0 );
    std::string fault_name;
    std::string face;

    /// face map
    std::unordered_map<std::string,CORNER_POINT_CELL_FACE_INDEX> facemap;
    facemap.emplace( "X-", CORNER_POINT_CELL_FACE_Xminus );
    facemap.emplace( "Y-", CORNER_POINT_CELL_FACE_Yminus );
    facemap.emplace( "Z-", CORNER_POINT_CELL_FACE_Zminus );
    facemap.emplace( "X+", CORNER_POINT_CELL_FACE_Xplus );
    facemap.emplace( "Y+", CORNER_POINT_CELL_FACE_Yplus );
    facemap.emplace( "Z+", CORNER_POINT_CELL_FACE_Zplus );
    facemap.emplace( "X", CORNER_POINT_CELL_FACE_Xplus );
    facemap.emplace( "Y", CORNER_POINT_CELL_FACE_Yplus );
    facemap.emplace( "Z", CORNER_POINT_CELL_FACE_Zplus );

    /// read first line
    fault_data.clear();
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            /// 1. read fault name
            token = strtok( text_line, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseFaultData:",
                                  "cannot read fault name!");
                return 0;
            }
            else
            {
                fault_name = token;
                removeSymbolsFromString( fault_name, symbolsToremove );
            }

            /// 2. read index ranges
            for( size_t i = 0; i<6; ++i )
            {
                token = strtok( NULL, delims );
                if ( token == NULL )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseFaultData:",
                                      "cannot read index!");
                    return 0;
                }
                else if( !isIntegerNumber(token) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseFaultData:",
                                      "index is not a digit!");
                    return 0;
                }
                else
                {
                    index_range[ i ] = atoi(token);
                    /// start indexing form 0
                    index_range[ i ] -= 1;
                }

            }

            /// 3. read face
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cout <<"\n"<< token << std::endl;
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseFaultData:",
                                  "cannot read face value!");
                return 0;
            }
            else
            {
                face = token;
                removeSymbolsFromString( face, symbolsToremove );
            }

            /// 4. assign fault data
            EclipseFault fault;
            fault.fault.reserve((index_range[1] - index_range[0] + 1)
                                * (index_range[3] - index_range[2] + 1)
                                * (index_range[5] - index_range[4] + 1));
            for( size_t i = index_range[0]; i <= index_range[1]; ++i )
                for( size_t j = index_range[2]; j <= index_range[3]; ++j )
                    for( size_t k = index_range[4]; k <= index_range[5]; ++k )
                    {
                      fault.Add( i, j, k, facemap[face] );
                    }
            fault_data.emplace(fault_name, fault);
        }
        ifs.getline( text_line, line_length );
        endOfblock = isEclipseEndOfBlock(text_line);
        if( !endOfblock && csmp::isBlankLine(text_line) )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        std::cout <<"\nreadEclipseFaultData: fault data has been read successfully.\n";

    return 1;
}



/**
    Reads property multipliers that allow users to treat the fault thickness.
*/
int readEclipseFaultTransmissibilityMultipliers( size_t NX, size_t NY, size_t NZ,
                                                 std::map<std::string,double>& multflt,
                                                 std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                               )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseFaultTransmissibilityMultipliers: reading fault transmissibility multipliers ...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    const char* symbolsToRemove ="'";
    bool        endOfblock(false);

    /// fault data
    std::string fault_name;
    double      mult;

    /// read first line
    multflt.clear();
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            /// 1. read fault name
            token = strtok( text_line, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseFaultTransmissibilityMultipliers:",
                                  "cannot read fault name!");
                return 0;
            }
            else
            {
                fault_name = token;
                removeSymbolsFromString( fault_name, symbolsToRemove );
            }

            /// 2. read fault transmissibility multiplier
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseFaultTransmissibilityMultipliers:",
                                  "cannot read index!");
                return 0;
            }
            else if( !isRealNumber(token) )
            {
                std::cout <<"\n"<< token << std::endl;
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseFaultTransmissibilityMultipliers:",
                                  "index is not a digit!");
                return 0;
            }
            else
                mult = atof(token);

            /// 3. assign multiplier
            multflt[ fault_name ] = mult;

        }
        ifs.getline( text_line, line_length );
        endOfblock = isEclipseEndOfBlock(text_line);
        if( !endOfblock && csmp::isBlankLine(text_line) )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        std::cout <<"\nreadEclipseFaultTransmissibilityMultipliers: fault transmissibility multipliers has been read successfully.\n";

    return 1;
}


/**
    To retrieve boundaries of the model?
*/
int readEclipseBoxData( std::vector<size_t>& box_data,
                        std::ifstream& ifs, char* text_line,
                        size_t line_length, bool verbose )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if( verbose )
        std::cout <<"\nreadEclipseBoxData: reading box data...";

    char*       token(0);
    const char* delims =" ,:,\t,\n,\r";
    bool        endOfblock(false);

    /// box data
    box_data.clear();
    box_data.resize(6,0);

    /// read first line
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    do{
        if ( !isEclipseCommentLine( text_line ) )
        {
            /// 1. read first box index
            token = strtok( text_line, delims );
            if ( token == NULL )
            {
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseBoxData:",
                                  "cannot read index!");
                return 0;
            }
            else if( !isIntegerNumber(token) )
            {
                std::cout <<"\n"<< token << std::endl;
                csmp_error.notice(csmp::ERROR,
                                  "readEclipseBoxData:",
                                  "index is not a digit!");
                return 0;
            }
            else
                box_data[ 0 ] = atoi(token);

            /// 1. read index ranges
            for( size_t i = 1; i<6; ++i )
            {
                token = strtok( NULL, delims );
                if ( token == NULL )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseBoxData:",
                                      "cannot read index!");
                    return 0;
                }
                else if( !isIntegerNumber(token) )
                {
                    std::cout <<"\n"<< token << std::endl;
                    csmp_error.notice(csmp::ERROR,
                                      "readEclipseBoxData:",
                                      "index is not a digit!");
                    return 0;
                }
                else
                {
                    box_data[ i ] = atoi(token);
                    /// start indexing form 0
                    box_data[ i ] -= 1;
                }
            }
        }
        ifs.getline( text_line, line_length );
        endOfblock = isEclipseEndOfBlock(text_line);
        if( !endOfblock && csmp::isBlankLine(text_line) )
        {
            ifs.getline( text_line, line_length );
            endOfblock = isEclipseEndOfBlock(text_line);
        }
    }
    while ( !endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        std::cout <<"\nreadEclipseBoxData: box data has been read successfully.\n";

    return 1;
}


bool isEclipseCommentLine( char* str )
{
    if ( str    == NULL )
        return false;
    if ( str[0] == '#' || str[0] == '%' )
        return true;
    if ( ( str[0] == '-' ) && ( str[1] == '-' ) )
        return true;

    const size_t strlength( strlen(str) );
    for ( size_t i = 0U; i<strlength; i++ )
      if ( str[i] == '#' || str[i] == '%' ){
           str[i] = '\0';
           break;
      }
      else if ( str[i] == '-' ) {
         if( i != strlength-1 )
             if ( str[i+1] == '-' ){
                 str[i] = '\0';
                 break;
             }
      }
    return false;
}

/// Eclipse ends its definitions with a forward slash
bool isEclipseEndOfBlock( char* str )
{
    if ( str == NULL ) return false;
    if ( str[0] == '/' )
        return true;
    const size_t strlength( strlen(str) );
    for ( size_t i = 0U; i<strlength; i++ )
      if ( str[i] != ' ' && str[i] != '/' )
          return false;
      else if( str[i] == '/' )
      {
          if( i+1 < strlength )
              str[i+1] = '\0';
          return true;
      }
    return false;
}

bool isEclipseLineWithEndOfBlock( char* str )
{
    if ( str == NULL ) return false;

    if ( str[0] == '/' )
        return true;

    const size_t strlength( strlen(str) );
    for ( size_t i = 0U; i<strlength; i++ )
      if ( str[i] == '/' ){
          if( i+1 < strlength )
              str[i+1] = '\0';
          return true;
      }

    return false;
}

bool readEclipseNonBlankLine( std::ifstream& ifs, char* text_line, size_t line_length, bool& endOfblock )
{
    ifs.getline( text_line, line_length );
    endOfblock = isEclipseEndOfBlock( text_line );
    while ( csmp::isBlankLine(text_line) && !ifs.eof() ) {
        ifs.getline( text_line, line_length );
        endOfblock = isEclipseEndOfBlock( text_line );
    }
    return true;
}

int readEclipseFirstLineInBlock( std::ifstream& ifs, char* text_line, size_t line_length, bool& endOfblock )
{
    // read first line in block
    if ( !readEclipseNonBlankLine( ifs, text_line, line_length, endOfblock ) )
        return 0;

    if ( ifs.eof() )
    {
        ifs.close();
        std::cout <<"\nreadFirstLineInBlock: Reading completed!" << std::endl;
        return 2;
    }

    return 1;
}

bool readEclipseKeyword( std::ifstream& ifs, char* text_line, size_t line_length,
                         std::string& keyword, std::vector<std::string>&  keyword_parameters )
{
    // Read keyword and it's parameters
    char*       token( 0 );
    //SKM fix: adding space to keyword delimiter
    //    const char* keyword_delims       = ":,\t,\n,\r";
    //    const char* keyword_param_delims = " :,\t,\n,\r";
    const char* keyword_delims       = " :,\t,\n,\r";
    const char* keyword_param_delims = " :,\t,\n,\r";

    // read keyword and make it upper case
    // in order to avoid case sensitive mistakes
    keyword = strtok( text_line, keyword_delims);
    std::transform( keyword.begin(), keyword.end(), keyword.begin(), ::toupper );

    bool endOfblock = false;
    // read keyword parameters
    keyword_parameters.clear();
    do
    {
        token = strtok( NULL, keyword_param_delims );
        endOfblock = isEclipseEndOfBlock( token );
        if( token != NULL && !endOfblock )
            keyword_parameters.push_back( token );
    }
    while( token != NULL && !endOfblock );

    return true;
}

int readEclipseIncludeFileName( std::ifstream& ifs, char* text_line, size_t line_length,
                                std::string& fname )
{
    // Read keyword and it's parameters
    const char* filename_delims       = ":,\t,\n,\r'`/";

    bool endOfblock = false;
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    // read file name
    fname = strtok( text_line, filename_delims);

    return 1;
}


/**
    Reads property values associated with the Eclipse grid.
*/
template<class varType>
bool readEclipseValue( const std::string& str, bool& default_value, size_t& num, varType& value )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    if ( !isRealNumber(str) )
    {
        std::cout <<"\n"<< str << std::endl;
        csmp_error.notice(csmp::ERROR,
                          "readEclipseValue:",
                          "value is not a digit!");
        return false;
    }
    else
    {
        const size_t star_position( str.find('*') );
        const bool multiple_values(star_position!=std::string::npos);
        if( !multiple_values )
        {
            num = 1;
            std::istringstream(str) >> value;
            default_value = false;
        }
        else
        {
            std::string number;

            // read number of assignments
            number = str.substr(0,star_position);
            std::istringstream(number) >> num;

            // read value to assign
            const size_t str_len( str.length() );
            if( star_position + 1 < str_len )
            {
                number = str.substr(star_position+1);
                std::istringstream(number) >> value;
                default_value = false;
            }
            else
            {
                default_value = true;
            }
        }
    }
    return true;
}
bool readEclipseValue(const std::string&,bool&,size_t&,size_t&);
bool readEclipseValue(const std::string&,bool&,size_t&,int&);
bool readEclipseValue(const std::string&,bool&,size_t&,double&);

int skipEclipseBlock( std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    bool endOfblock(false);
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    while ( !isEclipseLineWithEndOfBlock(text_line) && !ifs.eof() )
    {
        ifs.getline( text_line, line_length );
    }

    return 1;
}

void removeSymbolsFromString( std::string &str, const char* symbolsToRemove )
{
   const size_t length( std::strlen(symbolsToRemove) );
   for ( unsigned int i = 0; i < length; ++i ) {
      str.erase( std::remove( str.begin(), str.end(), symbolsToRemove[i] ), str.end() );
   }
}



// ECLIPSE WELLS

EclipseWellCompletion::EclipseWellCompletion()
:ic_(0),
 jc_(0),
 k_top_(0),
 k_bot_(0),
 w_start_(0),
 w_end_(0),
 is_open_(false)
{
}

EclipseWellCompletion
::~EclipseWellCompletion()
{
}


EclipseWell
::EclipseWell():
    well_name_ ("unknown"),
    group_name_("FIELD"),
    i_start_(0),j_start_(0),
    z_bhp_(0.0),
    phase_(ECLIPSE_OIL),
    drainage_radius_(0.0),
    inflow_type_(0),
    stop_or_shut_in_(false),
    cross_flow_(true)
{
}

EclipseWell
::~EclipseWell()
{
}



// ECLIPSE MODEL SETTINGS

EclipseModelSettings::EclipseModelSettings( const std::string& mesh_file_prefix )
    : mesh_file_prefix_   ( mesh_file_prefix ),
      regions_file_prefix_( mesh_file_prefix ),
      create_boundaries_  ( false ),
      tetra_mesh_         ( false )
{
    if( csmp::isRegionsFileExist( mesh_file_prefix.c_str() ) ){
        regions_.clear();
        csmp::readDesiredRegions( mesh_file_prefix.c_str(), regions_ );
    }
}

EclipseModelSettings::~EclipseModelSettings()
{
}


void EclipseModelSettings
::MeshSetup( bool exclude_inactive_cells,
             bool tetra_mesh,
             bool create_boundaries )    /* true = creates boundaries around model, false = does not create boundaries */
{
    exclude_inactive_cells_ = exclude_inactive_cells;
    tetra_mesh_             = tetra_mesh;
    create_boundaries_      = create_boundaries;
}

void EclipseModelSettings
::MeshSetup( const std::string& regions_file_prefix,
             bool exclude_inactive_cells,
             bool tetra_mesh,
             bool create_boundaries )    /* true = creates boundaries around model, false = does not create boundaries */
{
    exclude_inactive_cells_ = exclude_inactive_cells;
    tetra_mesh_             = tetra_mesh;
    create_boundaries_      = create_boundaries;
    if( csmp::isRegionsFileExist( regions_file_prefix.c_str() ) )
    {
        regions_.clear();
        csmp::readDesiredRegions( regions_file_prefix.c_str(), regions_ );
    }
}

void EclipseModelSettings
::MeshSetup( const std::set<std::string>& regions,
             bool exclude_inactive_cells,
             bool tetra_mesh,
             bool create_boundaries )    /* true = creates boundaries around model, false = does not create boundaries */
{
    exclude_inactive_cells_ = exclude_inactive_cells;
    tetra_mesh_         = tetra_mesh;
    create_boundaries_  = create_boundaries;
    regions_            = regions;
}


// PROPERTIES SETUP

void EclipseModelSettings
::SatNumPropertySetup( const std::string& prop_name,
                       const csmp::VARIABLE_TYPE& prop_type,
                       const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_SATNUM,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "X";
    prop.min   = 0.0;
    prop.max   = 1.0e3;
}

void EclipseModelSettings
::PvtNumPropertySetup( const std::string& prop_name,
                       const csmp::VARIABLE_TYPE& prop_type,
                       const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_PVTNUM,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "X";
    prop.min   = 0.0;
    prop.max   = 1.0e3;
}

void EclipseModelSettings
::RockNumPropertySetup( const std::string& prop_name,
                        const csmp::VARIABLE_TYPE& prop_type,
                        const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_ROCKNUM,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "X";
    prop.min   = 0.0;
    prop.max   = 1.0e3;
}

void EclipseModelSettings
::EqlNumPropertySetup( const std::string& prop_name,
                       const csmp::VARIABLE_TYPE& prop_type,
                       const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_EQLNUM,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "X";
    prop.min   = 0.0;
    prop.max   = 1.0e3;
}

void EclipseModelSettings
::FipNumPropertySetup( const std::string& prop_name,
                       const csmp::VARIABLE_TYPE& prop_type,
                       const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_FIPNUM,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "X";
    prop.min   = 0.0;
    prop.max   = 1.0e3;
}

void EclipseModelSettings
::NtgPropertySetup( const std::string& prop_name,
                    const csmp::VARIABLE_TYPE& prop_type,
                    const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_NTG,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "m^3 m^-3";
    prop.min   = 0.0;
    prop.max   = 1.0;
}

void EclipseModelSettings
::PoroPropertySetup( const std::string& prop_name,
                         const csmp::VARIABLE_TYPE& prop_type,
                         const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_PORO,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "m^3 m^-3";
    prop.min   = 0.0;
    prop.max   = 1.0;
}

void EclipseModelSettings
::PermPropertySetup( const std::string& prop_name,
                     const csmp::VARIABLE_TYPE& prop_type,
                     const csmp::PLACEMENT& prop_place,
                     const std::string&     prop_unit )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_PERM,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = prop_unit;
    prop.min   = 1.0e-1;
    prop.max   = 1.0e-20;
}

void EclipseModelSettings
::PressurePropertySetup( const std::string& prop_name,
                         const csmp::VARIABLE_TYPE& prop_type,
                         const csmp::PLACEMENT& prop_place,
                         const std::string&     prop_unit )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_PRESSURE,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = prop_unit;
    prop.min   = 0.0;
    prop.max   = 1.0e9;
}

void EclipseModelSettings
::SwatPropertySetup( const std::string& prop_name,
                     const csmp::VARIABLE_TYPE& prop_type,
                     const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_SWAT,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "m^3 m^-3";
    prop.min   = 0.0;
    prop.max   = 1.0;
}

void EclipseModelSettings
::SoilPropertySetup( const std::string& prop_name,
                     const csmp::VARIABLE_TYPE& prop_type,
                     const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_SOIL,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "m^3 m^-3";
    prop.min   = 0.0;
    prop.max   = 1.0;
}

void EclipseModelSettings
::SgasPropertySetup( const std::string& prop_name,
                     const csmp::VARIABLE_TYPE& prop_type,
                     const csmp::PLACEMENT& prop_place )
{
    std::pair<std::map<int,csmp::Parameter>::iterator,bool>
            it = properties_.insert( std::make_pair( ECLIPSE_SGAS,csmp::Parameter() ) );
    csmp::Parameter& prop = (*it.first).second;
    prop.name  = prop_name;
    prop.key.type  = prop_type;
    prop.key.place = prop_place;
    prop.unit  = "m^3 m^-3";
    prop.min   = 0.0;
    prop.max   = 1.0;
}

// ************************************************************************************************

// STEPHAN & CARO's REPLACEMENT CODE

// ************************************************************************************************


/// provide continous stream of number tokens
char* const popToken( std::ifstream& ifs, char* text_line, size_t line_length )
 {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const char*  delims =" ,:,\t,\n,\r";
    static bool first_call(true);
    char*  token = (first_call == true ) ? strtok( text_line, delims ) : strtok( NULL, delims );
   
    while ( token == NULL and !ifs.eof() ) {
         ifs.getline( text_line, line_length );
         token = strtok( text_line, delims );
      }
    assert( token != NULL );

     // if something else than a number is read
     if ( ifs.eof() ) {
            csmp_error.notice(csmp::ERROR, "EclipseInterface::Read_COORD:", "reached end of file." );
            return NULL;
        }
     if ( !isRealNumber(token) ) {
            std::cout <<"\n"<< token << std::endl;
            csmp_error.notice(csmp::ERROR, "EclipseInterface::Read_COORD:", "value is not a digit!" );
            return NULL;
        }

    first_call = false;
    return token;
 }
 
 


/**
    ReadPillarCoordinates -> Read_COORD
     
    COORD = vertical grid top->bottom coordinate lines, where x,y,z... x pointing east, y south, and z down.
    These lines become the grid pillars, extending all the way from the top to the bottom of the grid.
    Subsequently, the pillars are subdivided by cell corners stored in ZCORN, for which x,y is found by
    interpolation.

    The ordering of the pillars and z-cordinates is in “book order”, line after line until a page has been read,
    then the next page until all data has been read. i.e.
    first the top layer (K=1) is read, line by line 
    (begin with J=1, read I=1,...,NX, then J=2, etc.) 
    Then repeat for K=2,...,NZ.
 
     NX number of cells in x-direction (E)
     NY number of cells in y-direction (S)
     
     pillar grid is y-first (rows), then x (columns) then z (3rd dim= downwards)
     
     called by ReadBlock(s) after ReadGridSpecs()  so dimensions of grid are already known
     
     @todo SKM write Out() const function for CornerPointGrid testing
*/
bool EclipseInterface::Read_COORD( std::ifstream& ifs, char* text_line, size_t line_length )
 {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    if ( csmp_error.Verbose() ) std::cout <<"\nEclipseInterface::Read_COORD: reading COORD...\n";
    assert( NX_ > 0 );
    assert( NZ_ > 0 );
    assert( NY_ > 0 );


    bool endOfblock(false);
    int firstLine( readEclipseFirstLineInBlock( ifs, text_line, line_length, endOfblock) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;
   
    const size_t i_stride(NX_ + 1);
    const size_t j_stride(NY_ + 1);
 
    grid_.Resize( i_stride, j_stride );
  
    // keeping j constant i go over the i's
    for ( size_t j=0U; j<j_stride; j++ )
      {
         // reading the pillar coordinates
         for ( size_t i=0U; i<i_stride; i++ )
           {
             // read pillar top
             double64 pillar_top_x = atof( popToken( ifs, text_line, line_length ) );
             double64 pillar_top_y = atof( popToken( ifs, text_line, line_length ) );
             double64 pillar_top_z = atof( popToken( ifs, text_line, line_length ) );
             double64 pillar_bot_x = atof( popToken( ifs, text_line, line_length ) );
             double64 pillar_bot_y = atof( popToken( ifs, text_line, line_length ) );
             double64 pillar_bot_z = atof( popToken( ifs, text_line, line_length ) );
             
             // creating pillar
             Pillar pillar(Point<3u>(pillar_top_x, pillar_top_y, pillar_top_z),
                           Point<3u>(pillar_bot_x, pillar_bot_y, pillar_bot_z));
             grid_(i,j) = std::move(pillar);
           }
      }
   
// TODO: grid_.Out();

    if ( csmp_error.Verbose() ) std::cout <<"\nreadEclipsePillarCoordinates: COORD block has been read successfully.\n";
   
    return true;
   
 }  // end Read_COORD
  
  
  
  
  
  
/** ReadCornerDepths -> Read_ZCORN, i.e. puts the beads on the pillars, dividing them by the cell corners.

    The cell coorners are stored in grid, an i,j-array of beads. However, here the CornerPointGrid class is used for this purpose.

    Each row contains NX + 1 coordinate entries, but for adjacent cells the values are duplicated.
    Inside the Pillar grids the rows are duplicated as well to create a 1:1 mapping of cell corners
    to entries.
    
    For each pillar point only the Z coordinate is read the others are deduced from the 
    the direction of the Pillar. This assumes that the pillar is a straight line.

    Corner depths are defined by the keyword ZCORN. The intersection between a (non-horizontal) coordinate line 
    and a depth value is unique, such that from coordinate lines and corner depths, 
    all coordinates can be calculated.
    It would perhaps seem natural to define corners cell by cell, 
    but Eclipse sticks strictly to the book page format, so (with x pointing east,
    y south) first the top northern edge of all cells are read from west to east, 
    then the southern edge, then advancing to next row of cells. 
    When the top of a layer has been read, the bottom is read in a similar fashion, and then we are ready for the next layer. 
    Recall that corners associated with the same coordinate line may have different depths, but are equal in continuous areas, 
    such that all corners must be defined – nothing is “implicitly assumed”. When explaining the syntax we will use indices 
    NW, NE, SW, SE to denote the corners (assuming a standard orientation of the grid), and T, B for Top, Bottom.
 
    ZCORN
    z(1,1,1)T,NW z(1,1,1)T,NE z(2,1,1)T,NW z(2,1,1)T,NE ... z(NX,1,1)T,NE z(1,1,1)T,SW z(1,1,1)T,SE z(2,1,1)T,SW z(2,1,1)T,SE ... 
    z(NX,1,1)T,SE z(1,2,1)T,NW z(1,2,1)T,NE z(2,2,1)T,NW z(2,2,1)T,NE ... 
    z(NX,2,1)T,NE z(1,2,1)T,SW z(1,2,1)T,SE z(2,2,1)T,SW z(2,1,1)T,SE ... z(NX,2,1)T,SE .
    .
    .
    z(1,NY,1)T,SW z(1,NY,1)T,SE z(2,NY,1)T,SW z(2,NY,1)T,SE ... z(NX,NY,1)T,SE z(1,1,1)B,NW z(1,1,1)B,NE z(2,1,1)B,NW z(2,1,1)B,NE ... z(NX,1,1)B,NE
    .
    .
    .
    z(1,NY,1)B,SW z(1,NY,1)B,SE z(2,NY,1)B,SW z(2,NY,1)B,SE ... z(NX,NY,1)B,SE z(1,1,2)T,NW z(1,1,2)T,NE z(2,1,2)T,NW z(2,1,2)T,NE ... z(NX,1,2)T,NE
    .
    .
    .
    z(1,NY,NZ)T,SW z(1,NY,NZ)T,SE z(2,NY,NZ)T,SW z(2,NY,NZ)T,SE ... z(NX,NY,NZ)T,SE
*/

  bool EclipseInterface::Read_ZCORN( std::ifstream& ifs, char* text_line, size_t line_length )
  {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    
    if ( csmp_error.Verbose() ) std::cout <<"\nEclipseInterface::Read_ZCORN: reading ZCORN...\n";
    assert( NX_ > 0 );
    assert( NZ_ > 0 );
    assert( NY_ > 0 );
    
	zcorn_.clear();
    zcorn_.resize(NX_ * NY_ * NZ_ * 8, std::numeric_limits<double64>::quiet_NaN());
    
    // 1. Read z coordinates
    struct ZCoords {
      double64 z[8];
    };
    size_t NXxNY = NX_ * NY_;

	//JC: To handle with Eclipse's shorthand notation
	//    ex) n*val means that the value val shall be repeated n times.
	std::vector<double64> values;

	char*       token(0);
	const char* delims = " ,:,\t,\n,\r";
	bool        endOfblock(false);
	size_t      counter(0U);
	size_t      num;
	double      value(-999.);
	bool        default_value;

	int firstLine(readEclipseFirstLineInBlock(ifs, text_line, line_length, endOfblock));
	if (firstLine == 0 || firstLine == 2)
		return firstLine;

	do {
		if (!isEclipseCommentLine(text_line))
		{
			token = strtok(text_line, delims);
			do {
				if (!readEclipseValue<double64>(std::string(token), default_value, num, value))
				{
					std::cout << "\n" << token << std::endl;
					csmp_error.notice(csmp::ERROR,
						"readEclipseCellData:",
						"value cannot be read!");
					return 0;
				}
				else
				{
					for (size_t i = 0; i < num; i++)
					{
						values.push_back(value);
						counter++;
					}
				}
				token = strtok(NULL, delims);
				endOfblock = isEclipseEndOfBlock(token);
			} while ((token != NULL) && !endOfblock);
		}
		if (!endOfblock)
		{
			ifs.getline(text_line, line_length);
			endOfblock = isEclipseEndOfBlock(text_line);
			if (!endOfblock && csmp::isBlankLine(text_line))
			{
				ifs.getline(text_line, line_length);
				endOfblock = isEclipseEndOfBlock(text_line);
			}
		}
	} while (!endOfblock && !csmp::isBlankLine(text_line) && !ifs.eof());

	//JC: verify # of values in the block ZCORN
	size_t num_zcorn = (NX_ * 4) * NY_ * 2 * NZ_;
	if (num_zcorn != values.size())
	{
		std::cout << "\n" << token << std::endl;
		csmp_error.notice(csmp::ERROR,
			"Read_ZCORN:",
			"# of data is invalid!");
		return false;
	}

	size_t pos = 0;
    for ( size_t k=0U; k<NZ_; k++ )
    {
      for ( size_t j=0U; j<NY_; j++ )
      {
		for ( size_t i=0U; i<NX_; i++ )
        {			
		  double64 t_nw = values.at(pos++);
		  double64 t_ne = values.at(pos++);
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 0] = t_nw;		  
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 1] = t_ne;
        }
        for ( size_t i=0U; i<NX_; i++ )
        {			
		  double64 t_sw = values.at(pos++);
		  double64 t_se = values.at(pos++);
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 2] = t_sw;
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 3] = t_se;
        }
      }
      for ( size_t j=0U; j<NY_; j++ )
      {
        for ( size_t i=0U; i<NX_; i++ )
        {
		  double64 b_nw = values.at(pos++);			
		  double64 b_ne = values.at(pos++);
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 4] = b_nw;
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 5] = b_ne;
        }
        for ( size_t i=0U; i<NX_; i++ )
        {
		  double64 b_sw = values.at(pos++);
		  double64 b_se = values.at(pos++);
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 6] = b_sw;
          zcorn_[(i + j * NX_ + k * NXxNY)*8 + 7] = b_se;
        }
      }
    }
    
    if ( csmp_error.Verbose() ) std::cout <<"\nRead_ZCORN: ZCORN block has been read successfully.\n";
    
    return true;
    
  } // Read_ZCORN
  

void EclipseInterface::SetProperties(const std::map<int,csmp::Parameter>& p )
{
    properties_ = p;
}

template<class Container>
void EclipseInterface
::GetRegions( Container& data )
{
  Container newdata( regions_.begin(), regions_.end() );
  std::swap(data, newdata);
}

template void EclipseInterface::GetRegions( std::vector<std::string>& );
template void EclipseInterface::GetRegions( std::list<std::string>& );
template void EclipseInterface::GetRegions( std::set<std::string>& );





template<class Container>
void EclipseInterface
::GetFaults( Container& data )
{
  Container newdata( faults_.begin(), faults_.end() );
  std::swap(data, newdata);
}

template void EclipseInterface::GetFaults( std::vector<std::string>& );
template void EclipseInterface::GetFaults( std::list<std::string>& );
template void EclipseInterface::GetFaults( std::set<std::string>& );





template<class Container>
void EclipseInterface
::GetWells( Container& data )
{
  Container newdata( wells_.begin(), wells_.end() );
  std::swap(data, newdata);
}

template void EclipseInterface::GetWells( std::vector<std::string>& );
template void EclipseInterface::GetWells( std::list<std::string>& );
template void EclipseInterface::GetWells( std::set<std::string>& );
  
  
  
  std::multimap<ijk,size_t>&
  EclipseInterface::IJKMap()
  {
    return grid_.IJKMap();
  }

 } // eclipse

} // end namespace csmp








