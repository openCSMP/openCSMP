#include "ANSYS_Interface.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "ModelTopology.h"
#include "VSet.h"
#include "Model.h"

using namespace std;

namespace csmp {

/** Constructs the ANSYS_Interface interface object.

@section arguments Input Arguments 

The default argument of the constructor (false) means that it will
try to output finite elements using a global interpolation coordinate
system to the client objects. However there are currently only 
line, triangular and tetrahedral elements with linear basis functions
available in CSMP because the elements must be straight-sided for such
an approach to work. Instead isoparametric elements can be used.  

In order to create isoparametric element models use 'true'.  

At the moment, CSMP supports analytical element forms only for tetra-triangular
poly-element type meshes.   */
ANSYS_Interface::ANSYS_Interface( bool create_isoparametric_element_mesh )
   : file_header_("not initialized"),
     boundary_tag_( "BOUNDARY" ),
     isoparametric_(create_isoparametric_element_mesh)
 {
 }

 
ANSYS_Interface::~ANSYS_Interface()
 {
 }



void ANSYS_Interface::Clear()
{
    if ( !object_specs_.empty() )
      object_specs_.erase( object_specs_.begin(), object_specs_.end() );
    if ( !object_elements_.empty() )
      object_elements_.erase( object_elements_.begin(), object_elements_.end() );
}



/// Builds Model after it was constructed with the default constructor.
template<uint32_t dim>
void ANSYS_Interface::Read_ANSYS_Mesh( const string& filename,
                                       VSet<dim>&         vset,
                                       ModelTopology&     mesh_topology,
                                       bool               binary_input_file,
                                       bool               reassign_boundary_flags )
{
   /// read ansys mesh
   if ( binary_input_file )
       ReadMeshBinary( filename, vset, mesh_topology, reassign_boundary_flags );
   else
       ReadMeshASCII( filename, vset, mesh_topology, reassign_boundary_flags );
}


template void ANSYS_Interface::Read_ANSYS_Mesh( const string&,VSet<1U>&,ModelTopology&,bool,bool);
template void ANSYS_Interface::Read_ANSYS_Mesh( const string&,VSet<2U>&,ModelTopology&,bool,bool);
template void ANSYS_Interface::Read_ANSYS_Mesh( const string&,VSet<3U>&,ModelTopology&,bool,bool);
































/**
 
Reads ANSYS finite element meshes created by any of the tools accessible
through MED. If box shaped model geometry is assumed
the method will attempt to set flags that will later be used
to assign boundary counditions in CSMP (see method documentation for
Model::AssignBoundaryValues3D() ).  

@section arguments Input Arguments 

The method needs the filename (without extension) of the '*.asc',
'*.dat' file pair that is output by the Input utility in the
ANSYS mesher.  

@param vset The correctly read mesh and the family definitions are returned into
the user supplied vset and model topology objects, respectively.  

@param mesh_topology  the elements that make up the different model regions

@section implementation Implementation

Reads Ulrike's binary output from Tetra 30/7/2003 which is written as 
follows:  

@code
// 32 bit integers and unsigned integers
unsigned int nnodes;
unsigned int records;
unsigned int lplist, lpfverts; // overall length of plist, pfvert records
unsigned int lelem;

// 32 bit integers and unsigned integers arrays
int          *pbflags, *pfverts;
unsigned int *pelmt, *pmtrl, *plist;

// 64 bit double arrays
double *pnts, *px, *py, *pz;
double *pbvals;

records = (unsigned int) nnodes; // changed

// writing the nnodes-sized data blocks to file
// size = nnodes (Knotenpunkte)
fwrite(&records, sizeof(unsigned int), 1, fdat); 
fwrite(px, sizeof(double), nnodes, fdat);
fwrite(py, sizeof(double), nnodes, fdat);
fwrite(pz, sizeof(double), nnodes, fdat);
fwrite(pbflags, sizeof(int), nnodes, fdat);
fwrite(pbvals, sizeof(double), nnodes, fdat);
// size = nelements (Finite Elemente)
fwrite(&lelem, sizeof(unsigned int), 1, fdat);
fwrite(pelmt, sizeof(int), lelem, fdat);
// size = sum over nodes per element
fwrite(&lplist, sizeof(unsigned int), 1, fdat);
fwrite(plist, sizeof(unsigned int), lplist, fdat);
// size = sum over neighbors per element
fwrite(&lpfverts, sizeof(unsigned int), 1, fdat);
fwrite(pfverts, sizeof(int), lpfverts, fdat);
// size = nelements (Finite Elemente)
fwrite(&lelem, sizeof(unsigned int), 1, fdat);
fwrite(pmtrl, sizeof(int), lelem, fdat);
@endcode

NB: Element IDs in the ANSYS dataset are 0..n-1 by default. 

@section application Application

Creates an interface the ANSYS suite of meshing tools. This interface
can also be used for mixed (poly-element type meshes).  

@section messages Messages 

To protect the
user against the erratic reading of transferred binary files a 
detailed reporting of the reading progress is output and 
errors are raised is encountered.  

*/

template<uint32_t dim>
void ANSYS_Interface::ReadMeshBinary( const string&  meshfile,
                                      VSet<dim>&  vset,
                                      ModelTopology&  mesh_topology,
                                      bool reassign_boundary_flags )
 {
    ErrorHandler& csmp_error(ErrorHandler::Instance());

    Clear();

    // 1. Reading the '*.asc' file with the mesh topology description
    // ---------------------------------------------------------------
    string   asc_name  = meshfile; asc_name += ".asc";
    ifstream ifs_asc( asc_name.c_str() );
    if ( !ifs_asc.is_open() )
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                    "ASCII geometry input file with extension '.asc' could not be opened");
    
    if ( !ReadTitleASCII( ifs_asc, file_header_ ) )
      csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                 "File header not read correctly");

    // initialise object_specs_ that is used in the construction of the model topology further below
    if ( !ReadRegionsAndElementTypesASCII( ifs_asc ) )
      csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                 "Region and element type information not read correctly");
   
    convert_ANSYS_To_CSMP_FiniteElementTypes( object_specs_, isoparametric_, dim );
    ifs_asc.close(); // '*.asc' geometry file

    // 2.0 Reading the '*.dat' file with the mesh data
    // -----------------------------------------------
    FILE*  ifs_dat(0);
    string dat_name = meshfile; dat_name += ".dat";

    if ( (ifs_dat=fopen( dat_name.c_str(), "rb" )) == NULL )
      csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                     "DAT input file with extension '.dat' could not be opened");

    if( csmp_error.Verbose() )
    {
        cout <<"\n\nANSYS_Interface::ReadMeshBinary: Reading connectivity file: '"<< dat_name;
        cout <<"'"<< endl << endl;
    }

    if ( !ReadNodeCoordinatesBinary( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                       "Node coordinates not read correctly");
      }
    if ( !ReadBoundaryFlagsAndConditionsBinary( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                       "Boundary flags and conditions not read correctly");
      }
    if ( !ReadPelementBinary( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                       "Element types not read correctly");
      }
    convert_ANSYS_To_CSMP_FiniteElementTypes( vset, isoparametric_ );

    if ( !ReadPlistBinary( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                       "Nodes per element information not read correctly");
      }
    if ( !ReadPfvertsBinary( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                       "Element neighbor information not read correctly");
      }
    if ( !ReadPmaterialBinary( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshBinary",
                                       "Material property identifiers for elements not read correctly");
      }
    fclose( ifs_dat ); // '*.dat' pdata file

    if ( object_specs_.size() != object_elements_.size() )
    {
      csmp_error.Note( ERROR, "ANSYS_Interface::ReadMeshBinary():",
                                "Mismatch in object names and object specifiers");
    }
    else if( csmp_error.Verbose() )
    {
         cout <<"\n\nANSYS_Interface::ReadMeshBinary(): Input (asc & dat) files '";
         cout << meshfile <<"' read successfully."<< endl;
    }

    // build the model topology
    const bool require_unique_names_of_volumes_surfaces_and_lines = true;
    const bool correct_orientation_of_surface_elements = true;
    const bool interactive_property_assignment = false;
    mesh_topology.ModelName( meshfile.c_str() );
    
    mesh_topology.EstablishTopology( vset,
                                     object_specs_, object_elements_, // region names, element types and elements per region
                                     require_unique_names_of_volumes_surfaces_and_lines,
                                     interactive_property_assignment,
                                     correct_orientation_of_surface_elements,
                                     reassign_boundary_flags );
    Clear();

 } // ReadMeshBinary

template void ANSYS_Interface::ReadMeshBinary( const string&,VSet<1U>&,ModelTopology&,bool);
template void ANSYS_Interface::ReadMeshBinary( const string&,VSet<2U>&,ModelTopology&,bool);
template void ANSYS_Interface::ReadMeshBinary( const string&,VSet<3U>&,ModelTopology&,bool);






















/** Reads ANSYS '*.asc' and '*.dat' files as created by the Input module
in the MED interface. The dat files which store the mesh connectivity
must be in ASCII format.

This method reads models of an arbitrary shape and returns them into
the supplied VSet and mesh topology file. In order to assign boundary
conditions to irregular shaped models one can use the outer surfaces
if they are named appropriately. This must be done using the CAD
and ANSYS tools.

@section arguments Input Arguments

The method needs the model name which it will automatically expand into
the names of the '*.asc' and '*.dat' files.

@section application Application

To read the ASCII file '*dat' specifications that make sense for small
box-shaped models for which the mesh data are to be transferred between
different platforms so that binary files cannot be used.

@section messages Messages

The output is rather verbose, alerting the user to all potential errors
that may have occurred in the reading and processing phase of the
input data.
*/

template<uint32_t dim>
void ANSYS_Interface::ReadMeshASCII( const string& meshfile,
                                     VSet<dim>&  vset,
                                     ModelTopology&  mesh_topology,
                                     bool reassign_boundary_flags )
 {
    ErrorHandler& csmp_error(ErrorHandler::Instance());

    Clear();

    // 1. Reading the '*.asc' file with the mesh topology description
    // --------------------------------------------------------------
    string    asc_name = meshfile; asc_name += ".asc";
    ifstream  ifs_asc( asc_name.c_str() );
    if ( !ifs_asc.is_open() )
       csmp_error.Note( ERROR, "ANSYS_Interface::ReadMeshASCII",
                                 "ASCII geometry input file with extension '.asc' could not be opened");

    if ( !ReadTitleASCII( ifs_asc, file_header_ ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "File header not read correctly");
      }
    // initialises 'object_specs_' that is used in EstablishTopology, see below
    if ( !ReadRegionsAndElementTypesASCII( ifs_asc ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Region and element type information not read correctly");
      }
    convert_ANSYS_To_CSMP_FiniteElementTypes( object_specs_, isoparametric_ , dim );
    ifs_asc.close(); // '*.asc' geometry file


    // 2. Reading the '*.dat' file with the mesh data
    // -------------------------------------------------------------
    string    dat_name = meshfile; dat_name += ".dat";
    ifstream  ifs_dat( dat_name.c_str() );

    if ( !ifs_dat.is_open() )
         csmp_error.Note( ERROR, "ANSYS_Interface::ReadMeshASCII",
                                   "DAT pdata input file with extension '.dat' could not be opened");

    if( csmp_error.Verbose() )
    {
        cout <<"\n\nANSYS_Interface::ReadMeshASCII: Reading connectivity file: '"<< dat_name;
        cout <<"'"<< endl << endl;
    }

    if ( !ReadNodeCoordinatesASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Node coordinates not read correctly");
      }
    if ( !ReadBoundaryFlagsAndConditionsASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Boundary flags and conditions not read correctly");
      }
    if ( !ReadPelementASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Element types not read correctly");
      }
    convert_ANSYS_To_CSMP_FiniteElementTypes( vset, isoparametric_ );

    if ( !ReadPlistASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Nodes per element information not read correctly");
      }
    if ( !ReadPfvertsASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Element neighbor information not read correctly");
      }
    if ( !ReadPmaterialASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "ANSYS_Interface::ReadMeshASCII",
                                    "Material property identifiers for elements not read correctly");
      }
    ifs_dat.close(); // '*.data' pdata file

    if ( object_specs_.size() != object_elements_.size() )
    {
      csmp_error.Note( ERROR, "ANSYS_Interface::ReadMeshASCII():",
                                "Mismatch in object names and object specifiers");
    }
    else if( csmp_error.Verbose() )
    {
         cout <<"\n\nANSYS_Interface::ReadMeshASCII(): Input (asc & dat) files '";
         cout << meshfile <<"' read successfully."<< endl;
    }

    /// check topology
    const bool require_unique_names_of_volumes_surfaces_and_lines = true;
    const bool correct_orientation_of_surface_elements = true;
    const bool interactive_property_assignment = false;
    mesh_topology.ModelName( meshfile.c_str() );
    mesh_topology.EstablishTopology( vset,
                                     object_specs_, object_elements_,
                                     require_unique_names_of_volumes_surfaces_and_lines,
                                     interactive_property_assignment,
                                     correct_orientation_of_surface_elements,
                                     reassign_boundary_flags );
    Clear();

 } // ReadMeshASCII


template void ANSYS_Interface::ReadMeshASCII( const string&,VSet<1U>&,ModelTopology&,bool );
template void ANSYS_Interface::ReadMeshASCII( const string&,VSet<2U>&,ModelTopology&,bool );
template void ANSYS_Interface::ReadMeshASCII( const string&,VSet<3U>&,ModelTopology&,bool );










/**

Replaces the element identifier information in the supplied VSet into the
corresponding CSMP finite element types taking into account whether
an isoparametric model shall be created or not.

The method converts the argument VSet.
*/
template<uint32_t dim>
void convert_ANSYS_To_CSMP_FiniteElementTypes( VSet<dim>& vset, bool isoparametric )
 {
    if ( !vset.HybridElementTypeMesh() )
      vset.ElementType( 0U, ANSYS_ElementSpecifications::CSMP_TypeFrom_ANSYS_Type( vset.ElementType(0U), isoparametric, dim ) );
    else
      for ( size_t i{0U}; i<vset.TotalNumberOfCells(); i++ )
        vset.ElementType( i, ANSYS_ElementSpecifications::CSMP_TypeFrom_ANSYS_Type( vset.ElementType(i), isoparametric, dim ) );

 } // end

template void convert_ANSYS_To_CSMP_FiniteElementTypes( VSet<1U>&,bool );
template void convert_ANSYS_To_CSMP_FiniteElementTypes( VSet<2U>&,bool );
template void convert_ANSYS_To_CSMP_FiniteElementTypes( VSet<3U>&,bool );

void convert_ANSYS_To_CSMP_FiniteElementTypes( multimap<string,string>& object_specs,
                                               bool isoparametric, uint32_t dim )
{
    string elmt_type;
    multimap<string,string>::iterator itEnd = object_specs.end();
    for ( multimap<string,string>::iterator
       it=object_specs.begin(); it!=itEnd; it++ )
    {
        elmt_type = (*it).second;
        (*it).second = ANSYS_ElementSpecifications::CSMP_TypeNameFrom_ANSYS_TypeName( elmt_type, isoparametric, dim );
    }
}










// READING COMMENTS


/** If string commences with # sign true is returned, else
false is returned.
*/
bool ANSYS_Interface::IsCommentLine( char* str ) const
 {
    if ( str[0] == '#' )                  return true;
    if ( str[0] == ' ' && str[1] == '#' ) return true;

    return false;
 }



/** If a comment occurs in the current line everything after the comment
is skipped. The functions returns true if such a comment was found.
*/
bool ANSYS_Interface::SkipPotentialComment( ifstream& ifs ) const
 {
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    char    c('#'); ifs.get(c);
    bool    is_comment(false);
    string  comment;

    while ( !ifs.eof() && !isdigit(c) && c != '\n' && c != '\r' )
      {
         // checking for comments
         if ( c == '#' ) is_comment = true;
         if ( is_comment ) comment += c;
         ifs.get(c);
      }
    ifs.unget();

    if ( is_comment && csmp_error.Verbose() )
        cout <<"\n\tSkipped comment: "<< comment << endl;

    return is_comment;

 } // end SkipPotentialComment




/** Advances the file stream to behind the comment line.
  */
void ANSYS_Interface::AdvancePastCommentLine( ifstream& ifs ) const
 {
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    char  text_line[256];
    do {
          ifs.getline( text_line, 256 );
       }
    while ( (!IsCommentLine(text_line) && !ifs.eof()) );

    if( csmp_error.Verbose() )
        cout <<"\n\tSkipped comment line: "<< text_line << endl;

 } // end SkipPotentialCommentLines













/// ASCII REGION'S FILE


/** Reads the title line of the '*.asc' file and echoes it to screen.

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The file title is returned.  
*/
bool ANSYS_Interface::ReadTitleASCII( ifstream& ifs, string& title )
 {
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    char  text_line[256];
    
    // 1. Reading file header 
    ifs.getline( text_line, 256 ); // title line
    if( csmp_error.Verbose() )
        cout <<"\nANSYS_Interface::ReadTitleASCII: File header: "<< text_line << endl;
    title = text_line;
    
    ifs.getline( text_line, 256 ); // second explanatory line
    if( csmp_error.Verbose() )
        cout <<"\n\tFile specifications: "<< text_line << endl;

    if ( title.empty() ) return false;
    return true;
 
  } // ReadTitleASCII



/** Reads regions indiscriminately of their element types that they contain
and stores all the information in corresponding multimaps.

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return If the ascii file was read correctly, the method returns true, else 
false.  

@section implementation Implementation

This method is based on the assumption that object names never start 
with a digit. 

@section application Application

Used by the public interfaces of the class.  
*/
bool ANSYS_Interface::ReadRegionsAndElementTypesASCII( ifstream& ifs )
 {
    char                       text_line[INFO_STRING];
    const char* const          delims =" ,\t,:,\n,\r";
    string                object_name;
    string                elmt_specifier;
    int32_t                    n_regions(0), region(0);
    size_t                     n, n_elements;
    int32_t                    material;
    vector<size_t>        empty_list;
    set<string>      excluded_elmts;
    
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    // 0. excluding some line element types and polygons from the reading process
    excluded_elmts.insert("POLYGON");

    // 1. Finding number of regions (while no blank line is encountered)
    do {
          // checking for comment lines & lines which start with a number 
          if ( isdigit(ifs.get()) ) {
               ifs.unget(); 
               // parsing number of objects
               ifs.getline( text_line, INFO_STRING );
               n_regions = atoi( strtok( text_line, delims ) );
               if ( n_regions <= 0 )
                 throw csmp::Exception( ERROR, 
                                       "ANSYS_Interface::ReadRegionsAndElementTypesASCII", 
                                       "'*.asc' file seems to lack region descriptors"); 
               break;
            }
          ifs.unget(); 
          ifs.getline( text_line, NAME_STRING );
       }
    while ( strlen(text_line) > 0U && !IsCommentLine(text_line) );

    // 2. reading the descriptions of the model regions, if a region name is duplicated
    //    its element_type is appended to the name to make the region name unique
    while ( region < n_regions )
      {
          ifs.getline( text_line, INFO_STRING );

          if ( !IsCommentLine(text_line) ) {
               // 2.1 model region name (does not have to be a volumetric object) 
               object_name = strtok( text_line, delims );
               // 2.2 type of finite element used in the region 
               elmt_specifier = strtok(NULL,delims);
               // 2.3 material type identifier (material key)
               material = atoi(strtok(NULL,delims));
               // 2.4 number of elements which make up this region
               n_elements = atol( strtok(NULL,delims) );
               if ( n_elements == 0U ) 
                 throw csmp::Exception( ERROR, 
                                "ANSYS_Interface::ReadRegionsAndElementTypesASCII", 
                                "Region does not contain any finite elements:", 
                                 object_name.c_str() );
                                   
               else {
                    // 2.5 inserting the region description into 'object_specs'
                    object_specs_.insert( make_pair(object_name,elmt_specifier) );
                    region++;
                 }
            }
       } 
    
    // 3. Reading which individual elements belong to the subregions
    set<string>  processed_regions;
    region=0; 
    
    while ( region < n_regions )
      {
          ifs.getline( text_line, INFO_STRING );

          if ( !IsCommentLine(text_line) ) {
               // 3.1 reading region name & number of elements of region
               object_name = strtok(text_line,delims);
               
               // 3.2 reading element type specifier 
               elmt_specifier = strtok(NULL,delims);

               // 3.3 testing whether region of that name has already been processed 
               if ( (processed_regions.find(object_name)) != processed_regions.end() ) 
                  processed_regions.insert(object_name);

               // 3.4 checking name for its consistency with previous list
               if ( object_specs_.find(object_name) == object_specs_.end() )
                    throw csmp::Exception( ERROR, "ANSYS_Interface::ReadRegionsAndElementTypesASCII",
                                "Region name is not contained in earlier tabulation:",
                                 object_name.c_str() ); 
                 
               // 3.5 reading number of elements                   
               n_elements  = atol( strtok(NULL,delims) );
               if ( n_elements <= 0 ) 
                 throw csmp::Exception( ERROR, 
                                "ANSYS_Interface::ReadRegionsAndElementTypesASCII", 
                                "Region does not contain any finite elements:", 
                                 object_name.c_str() );
               
               // 3.6 adding region and filling the elements into a list 
               //     If the region does already exist a new name is created using the element                   
               else { 
                    auto rit = object_elements_.insert( make_pair(object_name,empty_list) );
                    (*rit).second.reserve( n_elements );
                    for ( size_t i{0U}; i<n_elements; i++ ) {
                         // expecting that elements are numbered 0...n-1
                         ifs >> n; 
                         (*rit).second.push_back( n );
                      }
                    ifs.getline( text_line, 256 ); // get pending line break
                 }
               region++;
            }
       } 

    // final checks   
    auto it2(object_elements_.begin());
    if ( !object_specs_.empty() ) {
        if( csmp_error.Verbose() )
        {
            cout <<"\nANSYS_Interface::ReadRegionsAndElementTypesASCII ";
            cout <<"Input file contains the geometric objects: "<< endl;
            for ( multimap<string,string>::const_iterator
               it=object_specs_.begin(); it!=object_specs_.end(); it++, it2++ )
            {
                cout <<"\n\t"<< (*it2).second.size() <<" "<< (*it).second <<" elements ";
                cout <<"\tforming region: '"<< (*it).first <<"'";
            }
            cout << endl << endl;
        }
        return true;
      }
      
    return false;
 
 } // ReadRegionsAndElementTypesASCII
 






 

/// ASCII DATA FILE


/** The method reads the node coordinate block from the supplied file
stream.  

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The node coordinates are stored in the px, py, and pz records of the
VSet.  
*/
template<uint32_t dim>
bool ANSYS_Interface::ReadNodeCoordinatesASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    size_t  n_nodes(0U);

    // getting the record size    
    ifs >> n_nodes;
    SkipPotentialComment( ifs );
    
    // Px recornd
    deque<double> X(n_nodes);
    if ( n_nodes == 0U ) {
         throw csmp::Exception( ERROR, "ANSYS_Interface::ReadNodeCoordinatesASCII", 
                                    "No node coordinates are specified in file" );
         return false;
      }
    for ( size_t i{0U}; i<n_nodes; i++ ) ifs >> X[i];

    // Py record
    deque<double> Y(n_nodes);
    for ( size_t i{0U}; i<n_nodes; i++ ) ifs >> Y[i];
    
    // Pz record
    deque<double> Z(n_nodes);
    for ( size_t i{0U}; i<n_nodes; i++ ) ifs >> Z[i];
     
    if( csmp_error.Verbose() )
    {
        cout <<"\nANSYS_Interface::ReadNodeCoordinatesASCII: ";
        cout <<"Read x,y,z coordinates of: "<< n_nodes <<" nodes."<< endl;
    }

    vset.AddXYZ( X, Y, Z );  
      
    if( csmp_error.Verbose() )
    {
        cout <<"\nANSYS_Interface::ReadNodeCoordinatesASCII: ";
        cout <<"Coordinates of "<< n_nodes <<" nodes read successfully."<< endl;
    }

    return true;
    
 } // ReadNodeCoordinatesASCII

template bool ANSYS_Interface::ReadNodeCoordinatesASCII( ifstream&,VSet<1U>& );
template bool ANSYS_Interface::ReadNodeCoordinatesASCII( ifstream&,VSet<2U>& );
template bool ANSYS_Interface::ReadNodeCoordinatesASCII( ifstream&,VSet<3U>& );


/**
 
The method reads the boundary condition and flag records from the 
input stream and stores them in the pbound and pbcond records
in the VSet.  

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The boundary flags and values at where nodes were flagged as boundary
nodes are stored in the second method argument.  
*/
template<uint32_t dim>
bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    int8_t        flag;
    double      bvalue;
    vector<bool>  bconds(vset.Vertices(),false);
    
    // Pbflags record
    if ( bconds.size() == 0 ) {
         throw csmp::Exception( WARNING, "ANSYS_Interface::ReadBoundaryFlagsAndConditionsASCII", 
                                  "vset is not initialized with coordinates yet" );
         return false;
      }

    AdvancePastCommentLine( ifs );

    for ( size_t i{0U}; i<bconds.size(); i++ ) {
         ifs >> flag;
         if ( flag != 0 ) {
              bconds[i] = true;
              // '0...n-1' since nodes are numbered this way
              vset.BFlag( i, flag );
           }
         else vset.BFlag( i, NOT );
      }
      
    AdvancePastCommentLine( ifs );

    // Pbvals record
    for ( size_t i{0U}; i<bconds.size(); i++ ) {
         ifs >> bvalue;
//         if ( bconds[i] ) vset.AddBValue( i+1, bvalue ); ignoring these values
      }

    if( csmp_error.Verbose() )
    {
        cout <<"\nANSYS_Interface::ReadBoundaryFlagsAndConditionsASCII: ";
        cout <<" Boundary flags and boundary values read successfully."<< endl;
    }

    return true;
    
 } // ReadBoundaryFlagsAndConditionsASCII


template bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsASCII( ifstream&,VSet<1U>& );
template bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsASCII( ifstream&,VSet<2U>& );
template bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsASCII( ifstream&,VSet<3U>& );


/**

Reads element-type integer-based identifiers for each element which are
later parsed to determine the number of nodes per element and neighbors
per element for the reading process.

@section arguments Input Arguments

A reference to the initialized ASCII input file stream.
 */
template<uint32_t dim>
bool ANSYS_Interface::ReadPelementASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    unsigned long  records(0U);
    int8_t         etype;

    // Getting record size
    ifs >> records;

    SkipPotentialComment( ifs );

    if ( records == 0U ) {
         csmp_error.Note( WARNING, "ANSYS_Interface::ReadPelementASCII",
                                  "Record of element types per element appears empty." );
         return false;
      }

    vector<int8_t> elmt_types;
    elmt_types.reserve(records);
    for ( size_t i{0U}; i<records; i++ ) {
        ifs >> etype;
        elmt_types.push_back( etype );
      }

    if ( elmt_types.size() < vset.Elements() ) {
         throw csmp::Exception( ERROR, "ANSYS_Interface::ReadPelementASCII",
                                    "Element type information could not be obtained for all elements" );
         return false;
      }

    vset.ElementTypes( elmt_types );

    return true;

 } // ReadPelementASCII

template bool ANSYS_Interface::ReadPelementASCII( ifstream&,VSet<1U>& );
template bool ANSYS_Interface::ReadPelementASCII( ifstream&,VSet<2U>& );
template bool ANSYS_Interface::ReadPelementASCII( ifstream&,VSet<3U>& );


/**
 
The method reads the 'plist' record from the supplied ASCII file. The
plist contains the nodel id's of the consecutive finite elements that 
make up the mesh. These are output to the VSet.  

@section arguments Input Arguments 

A reference to the initialized ASCII file stream that points to the 
beginning of the plist data block.  

@return The VSet into which the plist data shall be stored is supplied as the 
second method argument.  
*/
template<uint32_t dim>
bool ANSYS_Interface::ReadPlistASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    map<size_t,vector<int64_t> >  plist;
    vector<int64_t>          dummy;
    size_t                   total_items,
                             element(0), item(0),
                             id, nodes;
    const size_t             n_nodes(vset.Vertices());
    
    pair<map<size_t,vector<int64_t> >::iterator,bool>  it;
    pair<size_t,vector<int64_t> > data;
    
    // now the vset can be resized according to the new information
    assert( vset.HybridElementTypeMesh() );
    deque<uint32_t>  ndele(vset.TotalNumberOfCells());
                       
    for ( size_t i{0U}; i<vset.TotalNumberOfCells(); i++ )
      ndele[i] = csmp_elmt_specs::NodesPerElementOfType( vset.ElementType(i) );
    vset.ResizePlist( ndele );

    // reading number of data identifiers in the record
    ifs >> total_items;

    SkipPotentialComment( ifs );

    if ( total_items == 0 ) {
         csmp_error.Note( WARNING, "ANSYS_Interface::ReadPlistASCII", 
                                     "File record of nodes per element (plist) appears empty" );
         return false;
      }
    
    // reading plist
    while ( item < total_items )
      {
         // getting the number of nodes of the element to be read
         nodes = ndele[element];
         assert( nodes >= 2 && nodes <= 32 );
         
         data.first  = ++element;
         data.second =   dummy;
         it = plist.insert(data); // insertion of empty vector
         assert( it.second );
         (*it.first).second.reserve(nodes);
         
         // node ID's in file range 0...nodes-1
         for ( size_t i{0U}; i<nodes; i++ ) {
              ifs >> id;
              assert( id < n_nodes );
              (*it.first).second.push_back( id );
              item++;
           }
      }

    if ( item != total_items ) {
         throw csmp::Exception( ERROR, "ANSYS_Interface::ReadPlistASCII", 
                                        "File record of nodes per element (plist) was not correctly read" );
         return false;
      }

    vset.AddPlist( plist.begin(), plist.end() );
   
    if( csmp_error.Verbose() )
    {
        cout <<"\nANSYS_Interface::ReadPlistASCII: ";
        cout <<"Member node IDs read for: "<< element <<" elements."<< endl;
    }

    return true;
    
 } // ReadPlistASCII


template bool ANSYS_Interface::ReadPlistASCII( ifstream&,VSet<1U>& );
template bool ANSYS_Interface::ReadPlistASCII( ifstream&,VSet<2U>& );
template bool ANSYS_Interface::ReadPlistASCII( ifstream&,VSet<3U>& );


/**
 
The 'pfverts' record stores the IDs of the finite elements that sit 
opposite to the faces of the current finite element. If the face is 
located at the model boundary a negative boundary flag is stored in
stead of the neighbor ID.  

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The pfvert information is added to the second argument VSet. 
 
*/
template<uint32_t dim>
bool ANSYS_Interface::ReadPfvertsASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    // making an array of numbers of neighbors of each element
    assert( vset.HybridElementTypeMesh() );
    deque<uint32_t>  nbors( vset.TotalNumberOfCells() );
    for ( size_t i{0U}; i<vset.TotalNumberOfCells(); ++i )
      nbors[i] = csmp_elmt_specs::NeighborsPerElementOfType( vset.ElementType(i) );
    vset.ResizePfverts( nbors );

    // reading how many neighbor-element data identifiers are in the file record
    size_t total_items;
    ifs >> total_items;
    assert( total_items >= 1 );

    SkipPotentialComment( ifs );

    if ( total_items == 0 ) {
         csmp_error.Note( WARNING, "ANSYS_Interface::ReadPfvertsASCII", 
                                     "File record of neighbors per element (pfverts) appears empty" );
         return false;
      }
    
    // reading the pfvert file record
    size_t                       element(0), item(0);
    map<size_t,vector<int64_t> >  pfverts;
    vector<int64_t>               dummy;

    while ( item < total_items )
      {
         // getting the number of neighbors of the element to be read
         size_t neighbors = nbors[element];
         assert( neighbors >= 2  and  neighbors <= 32 );
         pair<size_t,vector<int64_t> >  data(element,dummy);
        
         // insertion of empty vector
         pair<map<size_t,vector<int64_t> >::iterator,bool>  it(pfverts.insert(data));
         assert( it.second );
         (*it.first).second.reserve(neighbors);

         // element ID's in file range 0...elements-1
         for ( size_t i{0U}; i<neighbors; ++i ) {
              int32_t  idx;
              ifs >> idx;
              if ( ifs.bad() ) {
                   cerr <<"\nread 'pfvert' record for element: "<< element <<", neighbor: "<< i <<", value: "<< idx;
                   throw csmp::Exception( ERROR, "ANSYS_Interface::ReadPfvertsASCII:",
                                         "file stream went bad, when reading 'pfverts' record; may be not enough entries." );
                }
              // element number must not be larger than the number of elements in the mesh
              assert( idx < vset.Elements() );
              // ascertaining that one of the possible options of boundary identifiers was used
              if ( idx < 0 ) assert( idx >= REGION_BOUNDARY );
              (*it.first).second.push_back( idx );
              item++;
           }
         element++;
      }

    if ( item != total_items ) {
         throw csmp::Exception( ERROR, "ANSYS_Interface::ReadPfvertsASCII:",
                                    "File record of neighbors per element (pfverts) was not correctly read." );
         return false;
      }

    // negative element numbers must now be identified as model boundaries
    // according to the families of the object neighbors which have these
    // element id's
    vset.AddPfverts( pfverts.begin(), pfverts.end() );
   
    if( csmp_error.Verbose() )
    {
        cout <<"\nANSYS_Interface::ReadPfvertsASCII: ";
        cout <<"Neighbor IDs read for: "<< element <<" elements."<< endl;
    }

    return true;

 } // ReadPfvertsASCII

template bool ANSYS_Interface::ReadPfvertsASCII( ifstream&,VSet<1U>& );
template bool ANSYS_Interface::ReadPfvertsASCII( ifstream&,VSet<2U>& );
template bool ANSYS_Interface::ReadPfvertsASCII( ifstream&,VSet<3U>& );

/**
 
Method reads the material ID record from the file stream and assigns 
the data to the argument VSet.  

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The material information is assigned to the supplied VSet.  
*/
template<uint32_t dim>
bool ANSYS_Interface::ReadPmaterialASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // Getting record size that should be equivalent to the number of elements in the model
    size_t records;
    ifs >> records;
    assert( records == vset.Elements() );

    SkipPotentialComment( ifs );

    if ( records == 0 ) {
         csmp_error.Note( WARNING, "ANSYS_Interface::ReadPmaterialASCII:",
                                  "Record of material types per element appears empty." );
         return false;
      }
    
    // only for the volumetric elements material data are read from file
    vector<int32_t> elmt_mtrls;
    elmt_mtrls.reserve(records);
    int32_t   emtrl;
    size_t  i(0);
    for ( i=0; i<records; ++i ) {
         ifs >> emtrl;
         elmt_mtrls.push_back(emtrl);
      }

    if ( elmt_mtrls.size() < vset.Elements() ) {
         cerr <<"\nFor the "<< vset.Elements() <<" elements, material parameter values were provided for only "<< records <<"\n";
         csmp_error.Note( ERROR, "ANSYS_Interface::ReadPmaterialASCII:",
                           "Element material information could not be obtained for all elements." );
         return false;
      }
    
    vset.AddPmtrl( elmt_mtrls.begin(), elmt_mtrls.end() );

    return true;

 } // ReadPmaterialASCII


template bool ANSYS_Interface::ReadPmaterialASCII( ifstream&,VSet<1U>& );
template bool ANSYS_Interface::ReadPmaterialASCII( ifstream&,VSet<2U>& );
template bool ANSYS_Interface::ReadPmaterialASCII( ifstream&,VSet<3U>& );









// BINARY data file

template<uint32_t dim>
bool ANSYS_Interface::ReadNodeCoordinatesBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  uibytes = sizeof( uint32_t);
    const size_t  dbytes  = sizeof(double);
    size_t        entries(0);

    // reading node coordinates
    // -------------------------------

    // reading number of nodes
    assert( uibytes <=sizeof(size_t) );
    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX ); // max size of unsigned 32 bit integer
    if( csmp_error.Verbose() )
    {
        cout <<"\n\treading "<< entries <<" node coordinates 'px,py,pz'..."<< endl;
        cout.flush();
    }
    // reading node coordinates 'px', 'py', 'pz' (double)
    double* px = new double[ entries ];
    double* py = new double[ entries ];
    double* pz = new double[ entries ];
    if ( fread( (void*) px, dbytes, entries, fp ) != entries )
     throw csmp::Exception( ERROR, "ANSYS_Interface::ReadNodeCoordinatesBinary",
                              "'px' array did not read correctly");
    if ( fread( (void*) py, dbytes, entries, fp ) != entries )
     throw csmp::Exception( ERROR, "ANSYS_Interface::ReadNodeCoordinatesBinary",
                              "'py' array did not read correctly");
    if ( fread( (void*) pz, dbytes, entries, fp ) != entries )
     throw csmp::Exception( ERROR, "ANSYS_Interface::ReadNodeCoordinatesBinary",
                              "'pz' array did not read correctly");

    // resizing the node array
    vset.ResizeNodes( entries );

    // assign node coordinates
    for ( size_t i{0U}; i<entries; i++ ) {
         vset.Px( i, px[i] );
         vset.Py( i, py[i] ); // ANSYS models will always have a three coordinate's
         vset.Pz( i, pz[i] ); // ANSYS models will always have a three coordinate's
      }
    delete[] px;
    delete[] py;
    delete[] pz;

    return true;

    /*
    // coordinate range checking
    double* xmin=min_element( px, px + entries );
    double* xmax=max_element( px, px + entries );
    double* ymin=min_element( py, py + entries );
    double* ymax=max_element( py, py + entries );
    double* zmin=min_element( pz, pz + entries );
    double* zmax=max_element( pz, pz + entries );

    if ( (*xmax) - (*xmin) < numeric_limits<double>::epsilon() )
      throw csmp::Exception( FATAL_ERROR, "ANSYS_Interface::ReadNodeCoordinatesBinary:",
                                          "all nodes have same X-coordinate. 2D model must lie in YZ plane.");

    if ( (*ymax) - (*ymin) < numeric_limits<double>::epsilon() )
      throw csmp::Exception( FATAL_ERROR, "ANSYS_Interface::ReadNodeCoordinatesBinary:",
                                          "all nodes have same Y-coordinate. 2D model must lie in XZ plane.");

    if ( (*zmax) - (*zmin) < numeric_limits<double>::epsilon() )
      throw csmp::Exception( FATAL_ERROR, "ANSYS_Interface::ReadNodeCoordinatesBinary:",
                                          "all nodes have same Z-coordinate. 2D model must lie in XY plane.");
    */
}

template bool ANSYS_Interface::ReadNodeCoordinatesBinary( FILE*,VSet<1U>&);
template bool ANSYS_Interface::ReadNodeCoordinatesBinary( FILE*,VSet<2U>&);
template bool ANSYS_Interface::ReadNodeCoordinatesBinary( FILE*,VSet<3U>&);


template<uint32_t dim>
bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  dbytes  = sizeof(double);

    // 1. reading boundary flags 'pbflags' (int)
    // ------------------------------------------

    if( csmp_error.Verbose() )
    {
        cout <<"\n\treading boundary flags 'pbflags'..."<< endl;
        cout.flush();
    }
    int32_t        ival;
    const int32_t  min29(MULTIPLE_BOUNDARIES);
    const size_t   nodes(vset.Vertices());
    size_t         flag_value_errors{0U};
    for ( size_t i{0U}; i<nodes; i++ ){
         fread( (void*) &ival, ibytes, 1U, fp );
         if ( ival < min29 ) {
              ival = IRREGULAR;
              flag_value_errors++;
           }
         vset.BFlag( i, static_cast<int8_t>(ival) );
     }

    if ( flag_value_errors > 0 ) {
         cerr <<"\n\n\t"<<"encountered "<< flag_value_errors <<" BOX_BOUNDARY flag values that were out of the range defined for this enum [-29,0].\n";
         csmp_error.Note( INFO, "ANSYS_Interface::ReadBoundaryFlagsAndConditionsBinary","'pbflag' value out of range, setting to IRREGULAR.");
      }


    // 2. reading boundary condition values 'pbounds' (double) and ignoring them
    // --------------------------------------------------------------------------
    // CURRENTLY IGNORED because users never assign any in ANSYS
    if ( csmp_error.Verbose() ) {
        cout <<"\n\treading boundary values 'pbounds'..."<< endl;
        cout.flush();
      }
    double dval;
    for ( size_t i{0U}; i<nodes; i++ )
      fread( (void*) &dval, dbytes, 1U, fp );

    return true;
}

template bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsBinary( FILE*,VSet<1U>&);
template bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsBinary( FILE*,VSet<2U>&);
template bool ANSYS_Interface::ReadBoundaryFlagsAndConditionsBinary( FILE*,VSet<3U>&);





template<uint32_t dim>
bool ANSYS_Interface::ReadPelementBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  uibytes = sizeof( uint32_t);
    size_t        entries(0);

    // reading element-type information record 'pelement' (unsigned int)
    // ---------------------------------------------------------------------
    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    if( csmp_error.Verbose() )
    {
        cout <<"\n\treading "<< entries <<" finite-element type specifiers from 'pelement'..."<< endl;
        cout.flush();
    }
    int32_t*   pelmt = new int32_t[ entries ];
    fread( (void*) pelmt, ibytes, entries, fp );
    // checking the validity of the element types (valid range 2-23)
    for ( size_t i{0U}; i<entries; i++ )
      if ( pelmt[i] < 2 || pelmt[i] > 23 )
        throw csmp::Exception( ERROR,
                              "ANSYS_Interface::ReadPelementBinary",
                             "'pelement' value out of range ANSYS-TYPE range (2-23).");

    // adding element types to vset
    vector<int8_t> elmt_types;
    elmt_types.assign( pelmt, pelmt + entries );
    vset.ElementTypes( elmt_types );

    delete[] pelmt;

    return true;
}

template bool ANSYS_Interface::ReadPelementBinary( FILE*, VSet<1U>& );
template bool ANSYS_Interface::ReadPelementBinary( FILE*, VSet<2U>& );
template bool ANSYS_Interface::ReadPelementBinary( FILE*, VSet<3U>& );



// assumes that Pelmt data has already been read
template<uint32_t dim>
bool ANSYS_Interface::ReadPlistBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  uibytes = sizeof( uint32_t);
    size_t        entries(0);

    // setting up the storage for 'plist' in VSet
    assert( vset.Elements() > 0 );
    const size_t  nelements(vset.TotalNumberOfCells());

    // now the vset can be resized according to the new information
    deque<uint32_t>  ndele(nelements);
    size_t              n_plist_entries_expected{0};
    for ( size_t i{0U}; i<nelements; ++i ) {
         ndele[i] = csmp_elmt_specs::NodesPerElementOfType( vset.ElementType(i) );
         n_plist_entries_expected += ndele[i];
      }
    vset.ResizePlist( ndele );

    // reading nodes connected to elements 'plist' (unsigned int)
    // from the element types array build the plist
    // --------------------------------------------------------------
    // reading the size of the plist array
    fread( (void*) &entries, uibytes, 1U, fp );
    
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    assert( entries == n_plist_entries_expected );
    
    if ( csmp_error.Verbose() ) {
          cout <<"\n\treading "<< nelements <<" nodes-per-element records from 'plist' (size="<< entries <<")..."<< endl;
          cout.flush();
      }
  
    // reading the plist
    uint32_t*  plist = new uint32_t[ entries ];
    fread( (void*) plist, uibytes, entries, fp );

    deque<vector<int64_t> >::iterator  it(vset.PlistBegin());
    size_t  nentry(0U);

    // the elements of the plist (node ids) are assigned
    for ( size_t i{0U}; i<nelements; i++, it++ )
      for ( size_t j{0U}; j<ndele[i]; j++ )
          (*it)[j] = plist[nentry++];

    delete[] plist;

    return true;
}

template bool ANSYS_Interface::ReadPlistBinary( FILE*, VSet<1U>& );
template bool ANSYS_Interface::ReadPlistBinary( FILE*, VSet<2U>& );
template bool ANSYS_Interface::ReadPlistBinary( FILE*, VSet<3U>& );





template<uint32_t dim>
bool ANSYS_Interface::ReadPfvertsBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  uibytes = sizeof( uint32_t);
    int64_t       entries(0);

    // setting up the storage for 'pfverts' in VSet
    assert( vset.Elements() > 0 );
    assert( vset.HybridElementTypeMesh() );
    const size_t   nelements(vset.TotalNumberOfCells());

    // making an array of with the number of neighbors for each element
    deque<uint32_t>  nbors( nelements );
    size_t           n_pfverts_entries_expected{0};
    for ( size_t i{0U}; i<nelements; ++i ) {
         assert( vset.ElementType(i) >= -128 );
         assert( vset.ElementType(i) <=  128 );
         nbors[i] = csmp_elmt_specs::NeighborsPerElementOfType( vset.ElementType(i) );
         n_pfverts_entries_expected += nbors[i];
//         cout << nbors[i] <<":"<< parseAbbreviated_FE_Type( vset.ElementType(i) ) <<" ";
      }
    assert( nbors.size() == vset.Elements() );
    vset.ResizePfverts( nbors );

    // reading neighbors connected to elements 'pfverts' (int)
    // -----------------------------------------------------------
    // size of pfverts array
    fread( (void*) &entries, uibytes, 1U, fp );
    
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    
    if ( entries != n_pfverts_entries_expected ) {
         cerr <<"\nneighbor records "<< entries <<" vs expected: "<< n_pfverts_entries_expected;
         csmp_error.Note( WARNING, "ANSYS_Interface::ReadPfvertsBinary",
                           "neighbor element ('pfvert') record in binary file is corrupt and needs to be replaced");
      }
    if ( csmp_error.Verbose() ) {
         cout <<"\n\treading "<< nelements <<" neighbor-list records from 'pfverts' (size="<< entries <<")..."<< endl;
         cout.flush();
      }
    if ( entries >= 2147483647 )
      csmp_error.Note( ERROR, "ANSYS_Interface::ReadPfvertsBinary", "too many elements in file to be read by this reader");
      
    int32_t*  pfverts = new int32_t[ entries ];
    fread( (void*) pfverts, ibytes, entries, fp );


// DEBUGGING - what is actually been read
/*
cerr <<"\npfverts read from ANSYS file:\n";
for ( auto i{0}; i<n_pfverts_entries_expected; ++i )
  cerr << pfverts[i] <<" ";
cerr << endl;
*/

    // reading the C array into the resized pfverts deque inside VData
    // ---------------------------------------------------------------
    deque<vector<int64_t> >::iterator it(vset.PfvertsBegin());
    size_t  n_entry(0U);
    for ( size_t i{0U}; i<nelements; i++, ++it ) {
        // minimum number of neighbors per element
        assert( nbors[i] >= 2 );
        for ( size_t j{0U}; j<nbors[i]; j++ ) {
             // ignoring extra entries if ANSYS pfverts array is too short
             // later uses VData::EstablishNeighborConnectivity3D() to fix things up
             if ( n_entry >= entries ) break;
             (*it)[j] = pfverts[n_entry++];
          }
      }

    delete[] pfverts;

    return true;
}

template bool ANSYS_Interface::ReadPfvertsBinary( FILE*, VSet<1U>& );
template bool ANSYS_Interface::ReadPfvertsBinary( FILE*, VSet<2U>& );
template bool ANSYS_Interface::ReadPfvertsBinary( FILE*, VSet<3U>& );

// TESTING of previous function

// test: deque is OK, although XCode debugger shows it with a size of 0
//cerr <<"\nReadPfvertsBinary: neighbor info\n:";
//for ( auto it=nbors.begin(); it!=nbors.end(); ++it )
//  cerr << (*it) <<" ";

// TESTING
/*
cerr <<"\nReadPfvertsBinary: neighbor info\n:";
size_t i(0);
 for ( deque<vector<int64_t> >::const_iterator
       ft=vset.PfvertsBegin(); ft!=vset.PfvertsEnd(); ft++, i++ )
   {
      cout << i <<": \t";
      for ( size_t j{0U}; j<(*ft).size(); j++ ) cout << (*ft)[j] <<"\t ";
      cout << endl;
   }
*/





template<uint32_t dim>
bool ANSYS_Interface::ReadPmaterialBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  uibytes = sizeof( uint32_t);
    size_t        entries(0);

    // reading material information 'pmtrl'
    // ------------------------------------
    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    if( csmp_error.Verbose() )
    {
        cout <<"\n\t"<<"reading "<< entries <<" material-type specifiers for the elements from 'pmtrl'..."<< endl;
        cout.flush();
    }
    int32_t*  mtrls = new int32_t[ entries ];
    fread( (void*) mtrls, ibytes, entries, fp );
    vector<int32_t> elmt_mtrls;
    elmt_mtrls.assign( mtrls, mtrls + entries );
    
    vset.AddPmtrl( elmt_mtrls.begin(), elmt_mtrls.end() );
    
    delete[] mtrls;

    return true;
}

template bool ANSYS_Interface::ReadPmaterialBinary( FILE*, VSet<1U>& );
template bool ANSYS_Interface::ReadPmaterialBinary( FILE*, VSet<2U>& );
template bool ANSYS_Interface::ReadPmaterialBinary( FILE*, VSet<3U>& );







// ANSYS MODEL INTERFACE

ANSYS_ModelSettings::ANSYS_ModelSettings( const string& mesh_file_prefix )
    : mesh_file_prefix_     ( mesh_file_prefix ),
      irregular_mesh_       ( true ),
      binary_file_          ( true ),
      create_boundaries_    ( true )
{
  if ( doesDomainsFileExist( mesh_file_prefix.c_str() ) ){
        regions_.clear();
        readDesiredDomains( mesh_file_prefix.c_str(), regions_ );
    }
}

ANSYS_ModelSettings::ANSYS_ModelSettings( const ANSYS_ModelSettings& s )
    : mesh_file_prefix_     ( s.mesh_file_prefix_ ),
      regions_              ( s.regions_ ),
      irregular_mesh_       ( s.irregular_mesh_ ),
      binary_file_          ( s.binary_file_ ),
      create_boundaries_    ( s.create_boundaries_ )
{
}

ANSYS_ModelSettings& ANSYS_ModelSettings::operator=( const ANSYS_ModelSettings& s )
{
    if( &s != this )
    {
        mesh_file_prefix_     = s.mesh_file_prefix_;
        regions_              = s.regions_;
        irregular_mesh_       = s.irregular_mesh_;
        binary_file_          = s.binary_file_;
        create_boundaries_    = s.create_boundaries_;
    }
    return *this;
}

ANSYS_ModelSettings::~ANSYS_ModelSettings()
{
}

void ANSYS_ModelSettings
::MeshSetup( const string& regions_file_prefix,
             bool irregular_mesh,        /* true = non-box shaped model, false = box shaped model */
             bool binary_file,           /* true = binary, false = ascii */
             bool create_boundaries )    /* true = creates boundaries around model, false = does not create boundaries */
{
  if( doesDomainsFileExist( regions_file_prefix.c_str() ) )
    {
       regions_.clear();
       readDesiredDomains( regions_file_prefix.c_str(), regions_ );
    }
    irregular_mesh_     = irregular_mesh;
    binary_file_        = binary_file;
    create_boundaries_  = create_boundaries;
}

void ANSYS_ModelSettings
::MeshSetup( const set<string>& regions, // desired regions
             bool irregular_mesh,        /* true = non-box shaped model, false = box shaped model */
             bool binary_file,           /* true = binary, false = ascii */
             bool create_boundaries )    /* true = creates boundaries around model, false = does not create boundaries */
{
    regions_            = regions;
    irregular_mesh_     = irregular_mesh;
    binary_file_        = binary_file;
    create_boundaries_  = create_boundaries;
}


} // end namespace csmp
