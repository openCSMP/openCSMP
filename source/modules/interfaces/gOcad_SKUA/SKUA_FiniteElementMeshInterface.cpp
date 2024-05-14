#include "SKUA_FiniteElementMeshInterface.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "ModelTopology.h"
#include "VSet.h"
#include "Model.h"
#include "TextFileIO.h"

using namespace std;

namespace csmp {

/** Constructs the SKUA_FiniteElementMeshInterface interface object.

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
SKUA_FiniteElementMeshInterface::SKUA_FiniteElementMeshInterface( bool create_isoparametric_element_mesh )
   : file_header_("not initialized"),
     boundary_tag_( "BOUNDARY" ),
     isoparametric_(create_isoparametric_element_mesh)
 {
 }

 
SKUA_FiniteElementMeshInterface::~SKUA_FiniteElementMeshInterface()
 {
 }



void SKUA_FiniteElementMeshInterface::Clear()
{
    if ( !object_specs_.empty() )
      object_specs_.erase( object_specs_.begin(), object_specs_.end() );
    if ( !object_elements_.empty() )
      object_elements_.erase( object_elements_.begin(), object_elements_.end() );
}



/// Builds Model after it was constructed with the default constructor.
template<uint32_t dim>
void SKUA_FiniteElementMeshInterface::Read_SKUA_Mesh( const string& filename,
                                                      VSet<dim>&         vset,
                                                      ModelTopology&     mesh_topology,
                                                      bool               binary_input_file )
{
   /// read ansys mesh
   if ( binary_input_file )
       ReadMeshBinary( filename, vset, mesh_topology );
   else
       ReadMeshASCII( filename, vset, mesh_topology );
}


template void SKUA_FiniteElementMeshInterface::Read_SKUA_Mesh( const string&,VSet<1U>&,ModelTopology&,bool);
template void SKUA_FiniteElementMeshInterface::Read_SKUA_Mesh( const string&,VSet<2U>&,ModelTopology&,bool);
template void SKUA_FiniteElementMeshInterface::Read_SKUA_Mesh( const string&,VSet<3U>&,ModelTopology&,bool);





/**
 
Reads SKUA finite element meshes created by SKUA's finite element mesher. 
If box shaped model geometry is assumed
the method will attempt to set flags that will later be used
to assign boundary counditions in CSMP (see method documentation for
Model::AssignBoundaryValues3D() ).  

@section arguments Input Arguments 

The method needs the filename (without extension) of the '*.asc',
'*.dat' file pair that is output by the Input utility in the
SKUA mesher.  

@param vset The correctly read mesh and the family definitions are returned into
the user supplied vset and model topology objects, respectively.  

@param mesh_topology  the elements that make up the different model regions

@section implementation Implementation

Reads output which is expected to have the following number type formats.

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

NB: Element IDs in the SKUA dataset are expected to be numbered 0..n-1 by default. 

@section application Application

Creates an interface the SKUA suite of meshing tools. This interface
can also be used for mixed (poly-element type meshes).  

@section messages Messages 

To protect the
user against the erratic reading of transferred binary files a 
detailed reporting of the reading progress is output and 
errors are raised is encountered.  

*/
template<uint32_t dim>
void SKUA_FiniteElementMeshInterface::ReadMeshBinary( const string&  meshfile,
                                                      VSet<dim>&  vset,
                                                      ModelTopology&  mesh_topology )
 {
    ErrorHandler& csmp_error(ErrorHandler::Instance());

    // 1. Reading the '*.asc' file with the mesh topology description
    // ---------------------------------------------------------------
    string   asc_name  = meshfile; asc_name += ".asc";
    ifstream ifs_asc( asc_name.c_str() );
    if ( !ifs_asc.is_open() )
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                       "ASCII geometry input file with extension '.asc' could not be opened");
    
    if ( !ReadTitleASCII( ifs_asc, file_header_ ) )     // O.K.
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                    "File header not read correctly");

    if ( !ReadRegionsAndElementTypesASCII( ifs_asc ) ) { // O.K. - produces object_specs_
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                    "Region and element type information not read correctly");
      }
   
    ifs_asc.close(); // '*.asc' geometry file

    // 2.0 Reading the binary '*.dat' file with the mesh data
    // -------------------------------------------------------
    FILE*  ifs_dat(0);
    string dat_name = meshfile; dat_name += ".dat";

    if ( (ifs_dat=fopen( dat_name.c_str(), "rb" )) == nullptr )
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                         "DAT input file with extension '.dat' could not be opened");

    if ( csmp_error.Verbose() ) {
         cout <<"\n\nSKUA_FiniteElementMeshInterface::ReadMeshBinary: Reading connectivity file: '"<< dat_name;
         cout <<"'"<< endl << endl;
      }

    if ( !ReadNodeCoordinatesBinary( ifs_dat, vset ) )
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                 "Node coordinates not read correctly");

    if ( !ReadBoundaryFlagsAndConditionsBinary( ifs_dat, vset ) )
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                     "Boundary flags and conditions not read correctly");

    if ( !ReadPelementBinary( ifs_dat, vset ) )
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                     "Element types not read correctly");
                                       
    else { // if 'pelmt' was read correctly, CSMP element-type identifiers are created from SKUA integer identifiers
         assert( vset.HybridElementTypeMesh() );
         for ( size_t i{0U}; i<vset.Cells(); ++ i )
           vset.ElementType( i, convertSKUA_ElementType( vset.ElementType(i), isoparametric_ ) );
       }
    if ( !ReadPlistBinary( ifs_dat, vset ) )
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                   "Nodes per element information not read correctly");
                                       
    if ( !ReadPfvertsBinary( ifs_dat, vset ) )
      csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                   "Element neighbor information not read correctly");
                                       
    if ( !ReadPmaterialBinary( ifs_dat, vset ) )
       csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary",
                                   "Material property identifiers for elements not read correctly");
      
    // TODO: read property data  
      
    fclose( ifs_dat ); // '*.dat' pdata file

    if ( object_specs_.size() != object_elements_.size() )
    {
      csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshBinary():",
                                "Mismatch in object names and object specifiers");
    }
    else if( csmp_error.Verbose() )
    {
         cout <<"\n\nSKUA_FiniteElementMeshInterface::ReadMeshBinary(): Input (asc & dat) files '";
         cout << meshfile <<"' read successfully."<< endl;
    }

    mesh_topology.ModelName( meshfile.c_str() );

 } // ReadMeshBinary

template void SKUA_FiniteElementMeshInterface::ReadMeshBinary( const string&, VSet<1U>&, ModelTopology& );
template void SKUA_FiniteElementMeshInterface::ReadMeshBinary( const string&, VSet<2U>&, ModelTopology& );
template void SKUA_FiniteElementMeshInterface::ReadMeshBinary( const string&, VSet<3U>&, ModelTopology& );






















/** Reads SKUA '*.asc' and '*.dat' files as created by the Input module
in the MED interface. The dat files which store the mesh connectivity
must be in ASCII format.

This method reads models of an arbitrary shape and returns them into
the supplied VSet and mesh topology file. In order to assign boundary
conditions to irregular shaped models one can use the outer surfaces
if they are named appropriately. This must be done using the CAD
and SKUA tools.

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
void SKUA_FiniteElementMeshInterface::ReadMeshASCII( const string& meshfile,
                                                     VSet<dim>&  vset,
                                                     ModelTopology&  mesh_topology )
 {
    ErrorHandler& csmp_error(ErrorHandler::Instance());

    // 1. Reading the '*.asc' file with the mesh topology description
    // --------------------------------------------------------------
    string    asc_name = meshfile; asc_name += ".asc";
    ifstream  ifs_asc( asc_name.c_str() );
    if ( !ifs_asc.is_open() )
       csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                 "ASCII geometry input file with extension '.asc' could not be opened");

    if ( !ReadTitleASCII( ifs_asc, file_header_ ) ) {     // O.K.
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "File header not read correctly");
      }
    if ( !ReadRegionsAndElementTypesASCII( ifs_asc ) ) { // O.K.
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Region and element type information not read correctly");
      }
    ifs_asc.close(); // '*.asc' geometry file


    // 2. Reading the '*.dat' file with the mesh data
    // -------------------------------------------------------------
    string    dat_name = meshfile; dat_name += ".dat";
    ifstream  ifs_dat( dat_name.c_str() );

    if ( !ifs_dat.is_open() )
         csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                   "DAT pdata input file with extension '.dat' could not be opened");

    if( csmp_error.Verbose() )
    {
        cout <<"\n\nSKUA_FiniteElementMeshInterface::ReadMeshASCII: Reading connectivity file: '"<< dat_name;
        cout <<"'"<< endl << endl;
    }

    if ( !ReadNodeCoordinatesASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Node coordinates not read correctly");
      }
    if ( !ReadBoundaryFlagsAndConditionsASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Boundary flags and conditions not read correctly");
      }
    if ( !ReadPelementASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Element types not read correctly");
      }
    else { // if 'pelmt' was read correctly, CSMP element-type identifiers are created from SKUA integer identifiers
         assert( vset.HybridElementTypeMesh() );
         for ( size_t i{0U}; i<vset.Cells(); ++ i )
           vset.ElementType( i, convertSKUA_ElementType( vset.ElementType(i), isoparametric_ ) );
      }
    if ( !ReadPlistASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Nodes per element information not read correctly");
      }
    if ( !ReadPfvertsASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Element neighbor information not read correctly");
      }
    if ( !ReadPmaterialASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII",
                                    "Material property identifiers for elements not read correctly");
      } // OK
    // reading potential property data (SKM 2/12/20)
    if ( !ReadPropertyRecordsASCII( ifs_dat, vset ) ) {
         csmp_error.Note(  ERROR, "SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII",
                                    "Property data were not read correctly");
      }
    ifs_dat.close(); // '*.data' pdata file

    if ( object_specs_.size() != object_elements_.size() )
      {
        csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadMeshASCII():",
                                  "Mismatch in object names and object specifiers");
      }
    else if( csmp_error.Verbose() )
      {
           cout <<"\n\nSKUA_FiniteElementMeshInterface::ReadMeshASCII(): Input (asc & dat) files '";
           cout << meshfile <<"' read successfully."<< endl;
      }

    /// name topology
    mesh_topology.ModelName( meshfile.c_str() );

 } // ReadMeshASCII


template void SKUA_FiniteElementMeshInterface::ReadMeshASCII( const string&,VSet<1U>&,ModelTopology& );
template void SKUA_FiniteElementMeshInterface::ReadMeshASCII( const string&,VSet<2U>&,ModelTopology& );
template void SKUA_FiniteElementMeshInterface::ReadMeshASCII( const string&,VSet<3U>&,ModelTopology& );






// READING COMMENTS


/** If string commences with # sign true is returned, else
false is returned.
*/
bool SKUA_FiniteElementMeshInterface::IsCommentLine( char* str ) const
 {
    if ( str[0] == '#' )                  return true;
    if ( str[0] == ' ' && str[1] == '#' ) return true;

    return false;
 }



/** If a comment occurs in the current line everything after the comment
is skipped. The functions returns true if such a comment was found.
*/
bool SKUA_FiniteElementMeshInterface::SkipPotentialComment( ifstream& ifs ) const
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
        cout <<"\n\tSkipPotentialComment: Skipped comment: "<< comment << endl;

    return is_comment;

 } // end SkipPotentialComment




/** Advances the file stream to behind the comment line.
  */
void SKUA_FiniteElementMeshInterface::AdvancePastCommentLine( ifstream& ifs ) const
 {
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );
    if ( !ifs.is_open() or ifs.eof() ) {
         csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::AdvancePastCommentLine",
                            "end of file or file that is not open supplied as an argument." );
         return;
      }

    char  text_line[1024];
    do {
          ifs.getline( text_line, 1024 );
       }
//    while ( !IsCommentLine(text_line) && !(string(text_line).find('\n')==string::npos || string(text_line).find('\n')==string::npos) );
//    while ( !IsCommentLine(text_line) && text_line[0] != '\0' );
    while ( !IsCommentLine(text_line) and isBlankLine(text_line) );

    if ( csmp_error.Verbose() ) {
        //   linux                   osx                      windows
        if ( text_line[0] == '\n' or text_line[0] == '\r' or (text_line[0] == '\r' and text_line[1] == '\n') )
          cout <<"\n\tAdvancePastCommentLine: advanced past line break.\n";
        else
          cout <<"\n\tAdvancePastCommentLine: Skipped comment line: "<< text_line << endl;
      }

 } // end AdvancePastCommentLine













/// ASCII REGION'S FILE


/** Reads the title line of the '*.asc' file and echoes it to screen.

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The file title is returned.  
*/
bool SKUA_FiniteElementMeshInterface::ReadTitleASCII( ifstream& ifs, string& title )
 {
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    char  text_line[256];
    
    // 1. Reading file header 
    ifs.getline( text_line, 256 ); // title line
    if( csmp_error.Verbose() )
        cout <<"\nSKUA_FiniteElementMeshInterface::ReadTitleASCII: File header: "<< text_line << endl;
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
bool SKUA_FiniteElementMeshInterface::ReadRegionsAndElementTypesASCII( ifstream& ifs )
 {
    char            text_line[INFO_STRING];
    const char* const delims =" ,\t,:,\n,\r";
    string           object_name;
    string           elmt_specifier;
    int32_t          n_regions(0), region(0);
    size_t           n, n_elements;
    int32_t          material;
    vector<size_t>   empty_list;
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
                                       "SKUA_FiniteElementMeshInterface::ReadRegionsAndElementTypesASCII", 
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
                                "SKUA_FiniteElementMeshInterface::ReadRegionsAndElementTypesASCII", 
                                "Region does not contain any finite elements:", 
                                 object_name.c_str() );
                                   
               else {
                    // 2.5 inserting the region descrition into 'object_specs'
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
                    throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadRegionsAndElementTypesASCII",
                                "Region name is not contained in earlier tabulation:",
                                 object_name.c_str() ); 
                 
               // 3.5 reading number of elements                   
               n_elements  = atol( strtok(NULL,delims) );
               if ( n_elements <= 0 ) 
                 throw csmp::Exception( ERROR, 
                                "SKUA_FiniteElementMeshInterface::ReadRegionsAndElementTypesASCII", 
                                "Region does not contain any finite elements:", 
                                 object_name.c_str() );
               
               // 3.6 adding region and filling the elements into a list 
               //     If the region does already exist a new name is created using the element                   
               else { 
                    multimap<string,vector<size_t> >::iterator
                    rit = object_elements_.insert( make_pair(object_name,empty_list) );
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
    multimap<string,vector<size_t> >:: const_iterator it2(object_elements_.begin());
    if ( !object_specs_.empty() ) {
        if( csmp_error.Verbose() )
        {
            cout <<"\nSKUA_FiniteElementMeshInterface::ReadRegionsAndElementTypesASCII ";
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
bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    size_t  n_nodes(0U);

    // getting the record size    
    ifs >> n_nodes;
    SkipPotentialComment( ifs );
    
    // Px recornd
    deque<double> X(n_nodes);
    if ( n_nodes == 0U ) {
         throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII", 
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
        cout <<"\nSKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII: ";
        cout <<"Read x,y,z coordinates of: "<< n_nodes <<" nodes."<< endl;
    }

    vset.AddXYZ( X, Y, Z );  
      
    if( csmp_error.Verbose() )
    {
        cout <<"\nSKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII: ";
        cout <<"Coordinates of "<< n_nodes <<" nodes read successfully."<< endl;
    }

    return true;
    
 } // ReadNodeCoordinatesASCII

template bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII( ifstream&,VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII( ifstream&,VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesASCII( ifstream&,VSet<3U>& );


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
bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    int32_t         flag;
    double      bvalue;
    vector<bool>  bconds(vset.Vertices(),false);
    
    // Pbflags record
    if ( bconds.size() == 0 ) {
         throw csmp::Exception( WARNING, "SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsASCII", 
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
      }
      
    AdvancePastCommentLine( ifs );

    // Pbvals record
    for ( size_t i{0U}; i<bconds.size(); i++ ) {
         ifs >> bvalue;
//         if ( bconds[i] ) vset.AddBValue( i+1, bvalue ); ignoring these values
      }

    if( csmp_error.Verbose() )
    {
        cout <<"\nSKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsASCII: ";
        cout <<" Boundary flags and boundary values read successfully."<< endl;
    }

    return true;
    
 } // ReadBoundaryFlagsAndConditionsASCII


template bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsASCII( ifstream&,VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsASCII( ifstream&,VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsASCII( ifstream&,VSet<3U>& );


/**

Reads element-type integer-based identifiers for each element which are
later parsed to determine the number of nodes per element and neighbors
per element for the reading process.

@section arguments Input Arguments

A reference to the initialized ASCII input file stream.
 */
template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadPelementASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    unsigned long  records(0U);
    int8_t         etype;

    // Getting record size
    ifs >> records;

    SkipPotentialComment( ifs );

    if ( records == 0U ) {
         csmp_error.Note( WARNING, "SKUA_FiniteElementMeshInterface::ReadPelementASCII",
                                  "Record of element types per element appears empty." );
         return false;
      }

    vector<int8_t> elmt_types;
    set<int8_t>    range_of_types;
    elmt_types.reserve(records);
    for ( size_t i{0U}; i<records; i++ ) {
        ifs >> etype;
        elmt_types.push_back( etype );
        range_of_types.insert( etype );
      }

    if ( elmt_types.size() < vset.Elements() ) {
         throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadPelementASCII",
                                    "Element type information could not be obtained for all elements" );
         return false;
      }

    // if there is only a single type there is just one entry needed
    if ( range_of_types.size() == 1U ) {
         vset.SingleElementType( elmt_types[0] );
      }
    else vset.ElementTypes( elmt_types );

    return true;

 } // ReadPelementASCII

template bool SKUA_FiniteElementMeshInterface::ReadPelementASCII( ifstream&,VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPelementASCII( ifstream&,VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPelementASCII( ifstream&,VSet<3U>& );






/**
      reads vector<vector> from filestream where the elements of the vector are sequential
      the number of entries per vector is deduced from the element type stored in  'pelmt'
*/   
template<typename T>   
void readVectorOfVectors( ifstream& ifs, const deque<uint32_t>& vector_sizes, size_t total_items, deque<vector<T> >& file_records )
 {
    assert( ifs.is_open() );
    assert( !vector_sizes.empty() );
    assert( total_items > vector_sizes.size() ); 
    file_records.clear();
    size_t item=0U; 
    
    // case 1: all stored vectors have the same size
    // ---------------------------------------------
    // getting the number of nodes of the element to be read
    // reading record when all entries have the same size
    if ( vector_sizes.size() == 1 ) {
        const size_t entries_per_vector( vector_sizes[0] );
        while ( item < total_items )
          {
             // node ID's in file range 0...nodes-1
             vector<T> data;
             data.reserve( entries_per_vector );
             T id;
             for ( size_t i{0U}; i<entries_per_vector; ++i ) {
                  ifs >> id;
                  // assumption: there are not more nodes that elements * nodes_per_element
                  assert( id < static_cast<int64_t>(total_items) ); 
                  data.push_back( id );
                  item++;
               }
             file_records.emplace_back( data );
          }
      }

    // case 2: each stored vector has a different size
    // -----------------------------------------------
    // getting the number of nodes of the element to be read
    // reading record 
    else {
        size_t element_data_read(0U);
        while ( item < total_items )
          {
             const size_t entries_per_vector( vector_sizes[element_data_read] );
             // node ID's in file range 0...nodes-1
             vector<T> data;
             data.reserve( entries_per_vector );
             T id;
             for ( size_t i{0U}; i<entries_per_vector; ++i ) {
                  ifs >> id;
                  assert( id < static_cast<int64_t>(total_items) );
                  data.push_back( id );
                  item++;
               }
             file_records.emplace_back( data );
             element_data_read++;
          }
      }
      
    if ( item != total_items )
         throw csmp::Exception( ERROR, "readVectorOfVectors", 
                               "File record of vector<vector<typename>> was not correctly read" );

 } // end read inlined vector

template void readVectorOfVectors( ifstream&, const deque<uint32_t>& vector_sizes, size_t, deque<vector<uint32_t> >& );
template void readVectorOfVectors( ifstream&, const deque<uint32_t>& vector_sizes, size_t, deque<vector<int64_t> >& );



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
bool SKUA_FiniteElementMeshInterface::ReadPlistASCII( ifstream& ifs, VSet<dim>& vset, bool test_for_consecutive_node_numbering )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    map<size_t,vector<int64_t> >  plist;
    size_t                        total_items;
        
    // reading number of data identifiers in the record
    ifs >> total_items;
    SkipPotentialComment( ifs );

    if ( total_items == 0 ) {
         csmp_error.Note( WARNING, "SKUA_FiniteElementMeshInterface::ReadPlistASCII", 
                                     "File record of nodes per element (plist) appears empty" );
         return false;
      }
    
    // now the element types can be deduced and the 'plist' can be resized according to the new information
    deque<uint32_t>  ndele;
    if ( !vset.HybridElementTypeMesh() ) {
         ndele.push_back( csmp_elmt_specs::NodesPerElementOfType( vset.ElementType(0U) ) );
         vset.ResizePlist( total_items / ndele[0] ); // number of elements
      }
    else {
         assert( vset.HybridElementTypeMesh() );
         for ( size_t i{0U}; i<vset.Cells(); i++ )
           ndele.push_back( csmp_elmt_specs::NodesPerElementOfType( vset.ElementType(i) ) );
         vset.ResizePlist( ndele ); 
      }
      
    // reading plist
    deque<vector<int64_t> >  file_records;
    readVectorOfVectors( ifs, ndele, total_items, file_records );   

    // testing whether the nodes are numbered consecutively from 0..n-1
    if ( test_for_consecutive_node_numbering ) {
         set<size_t> node_ids;
         for ( deque<vector<int64_t> >::const_iterator it=file_records.begin(); it!=file_records.end(); ++it )
           for ( vector<int64_t>::const_iterator nit=(*it).begin(); nit!=(*it).end(); ++nit )
             node_ids.insert( (*nit) );
         // does the record start with 0 and ends with n-1?
         if ( (*node_ids.begin()) != 0U ) {
             cerr <<"\nID of first node: "<< (*node_ids.begin());
             csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadPlistASCII", 
                                       "Node numbering does not start with zero" );
           }
         else if ( (*node_ids.rbegin()) != vset.Vertices()-1U ) {
             cerr <<"\nID of last node vs. nodes in VSet: "<< (*node_ids.rbegin()) <<" vs. "<< vset.Vertices();
             csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadPlistASCII", 
                                       "Node numbering does not end with number of nodes in VSet-1 (=non-consecutive)" );
           }
      }

    vset.AddPlist( file_records.begin(), file_records.end() );
   
    if ( csmp_error.Verbose() ) {
          cout <<"\nSKUA_FiniteElementMeshInterface::ReadPlistASCII: ";
          cout <<"Member node IDs read for: "<< file_records.size() <<" elements."<< endl;
      }

    return true;
    
 } // ReadPlistASCII


template bool SKUA_FiniteElementMeshInterface::ReadPlistASCII( ifstream&,VSet<1U>&, bool );
template bool SKUA_FiniteElementMeshInterface::ReadPlistASCII( ifstream&,VSet<2U>&, bool );
template bool SKUA_FiniteElementMeshInterface::ReadPlistASCII( ifstream&,VSet<3U>&, bool );


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
bool SKUA_FiniteElementMeshInterface::ReadPfvertsASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    // reading how many neighbor-element data identifiers are in the file record
    size_t total_items;
    ifs >> total_items;
    assert( total_items >= 1 );

    SkipPotentialComment( ifs );

    if ( total_items == 0 ) {
         csmp_error.Note( WARNING, "SKUA_FiniteElementMeshInterface::ReadPfvertsASCII", 
                                     "File record of neighbors per element (pfverts) appears empty" );
         return false;
      }
    
    // making an array of numbers of neighbors of each element
    deque<uint32_t>  nbors; 
    if ( !vset.HybridElementTypeMesh() ) {
         nbors.push_back( csmp_elmt_specs::NeighborsPerElementOfType( vset.ElementType(0U) ) ); 
         vset.ResizePfverts( total_items / nbors[0] );
      }
    else {
         assert( vset.HybridElementTypeMesh() );
         for ( size_t i{0U}; i<vset.Cells(); ++i )
           nbors.push_back( csmp_elmt_specs::NeighborsPerElementOfType( vset.ElementType(i) ) );
         vset.ResizePfverts( nbors );
      }
      
    // reading the pfvert file record
    deque<vector<int64_t> >  file_records; 
    readVectorOfVectors( ifs, nbors, total_items, file_records );   

    // negative element numbers must now be identified as model boundaries
    // according to the families of the object neighbors which have these
    // element id's
    vset.AddPfverts( file_records.begin(), file_records.end() );
   
    if ( csmp_error.Verbose() )
      {
          cout <<"\nSKUA_FiniteElementMeshInterface::ReadPfvertsASCII: ";
          cout <<"Neighbor IDs read for: "<< file_records.size() <<" elements."<< endl;
      }

    return true;

 } // ReadPfvertsASCII

template bool SKUA_FiniteElementMeshInterface::ReadPfvertsASCII( ifstream&,VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPfvertsASCII( ifstream&,VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPfvertsASCII( ifstream&,VSet<3U>& );


// DEBUGGING
/*
cerr <<"\n'pfverts'"<< endl;
for ( size_t i{0U}; i<file_records.size(); ++i ) {
     cerr <<"\n"<< i <<": ";
     for ( size_t j{0U}; j<file_records[i].size(); j++ )
     cerr << file_records[i][j] <<" ";
  }
*/





/**
 
Method reads the material ID record from the file stream and assigns 
the data to the argument VSet.  

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return The material information is assigned to the supplied VSet.  
*/
template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadPmaterialASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // Getting record size that should be equivalent to the number of elements in the model
    size_t records;
    ifs >> records;
    assert( records == vset.Elements() );

    SkipPotentialComment( ifs );

    if ( records == 0 ) {
         csmp_error.Note( WARNING, "SKUA_FiniteElementMeshInterface::ReadPmaterialASCII:",
                                  "Record of material types per element appears empty." );
         return false;
      }
    
    // only for the volumetric elements material data are read from file
    vector<int32_t> elmt_mtrls;
    elmt_mtrls.reserve( records );
    int32_t         emtrl;
    for ( size_t i{0U}; i<records; ++i ) {
         ifs >> emtrl;
         elmt_mtrls.push_back(emtrl);
      }

    if ( elmt_mtrls.size() < vset.Elements() ) {
         cerr <<"\nFor the "<< vset.Elements() <<" elements, material parameter values were provided for only "<< records <<"\n";
         csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadPmaterialASCII:",
                           "Element material information could not be obtained for all elements." );
         return false;
      }
      
    // Adding the material ID record to the VSet
    vset.AddPmtrl( elmt_mtrls.begin(), elmt_mtrls.end() );  

    return true;

 } // ReadPmaterialASCII


template bool SKUA_FiniteElementMeshInterface::ReadPmaterialASCII( ifstream&,VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPmaterialASCII( ifstream&,VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPmaterialASCII( ifstream&,VSet<3U>& );




/**
 
Method reads a single property record into a PropertyData object and adds 
the data to the argument VSet. 

Each property record consists of :

name(string)    placement(enum PLACEMENT)   type (VARIABLE_TYPE)    (array length (int32))
number of data entries (int32)
flag values (int32)
property values (double)

for variable placement and type, upper case or lower case spellings are accepted.

@attention While space in variable names is expressed by an underscore which gets removed in CSMP.

@attention At this point, no range check is performed on the values read, but parsing errors related to variable placement, type and flags are reported.

@section arguments Input Arguments 

A reference to the initialized ASCII input file stream.  

@return true when the record was read correctly and completely. 

TODO: potentially, we could test variable names with the PropertyDatabase already here, to rectify issues 

*/
template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII( ifstream& ifs, VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    int32_t          integer;
   
    AdvancePastCommentLine( ifs );
   
    // 0. checking how many properties are contained in file if any
    if ( ifs.eof() ) 
      csmp_error.Note( WARNING, "SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII",
                        "end of file reached; input fle does not contain any property data" );  
    ifs >> integer;
    if ( integer <= 0 ) {
         csmp_error.Note( INFO, "SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII",
                            "text file does not contain any property data"); 
         return true;
      }
    const int32_t number_of_properties(integer);

    // 1. looping over the properties and storing them in the VSet
    for ( int32_t n=0; n<number_of_properties; ++n )
      {
        // 1.1 reading variable name, placement and type
        string property_name, str;
        // swallowing potential comments
        ifs >> property_name;
        if ( property_name[0] == '#' ) {
              getline( ifs, property_name );
              ifs >> property_name;
           }
        // replacing underscores with ' '
        for ( size_t i{0U}; i<property_name.size(); ++ i )
          if ( property_name[i] == '_' ) property_name[i] = ' ';
        ifs >> str;
        PLACEMENT place = parsePlacement( str.c_str() );
        ifs >> str;
        VARIABLE_TYPE type = parseType( str.c_str() );
        // 1.2 reading array length, if variable is an array or flagged array
        int32_t  array_length(0);
        if ( type == ARRAY || type == FLAGGEDARRAY ) {
             ifs >> integer;
             if ( integer <= 0 ) {
                  cerr <<"\n\tarray length = "<< integer;
                  csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII",
                                    "negative or zero array length read for variable", property_name.c_str() );
                  return false;
               }
             array_length = integer;
          }
        
        // reading number of value entries, cross-checking them using the variable type 
        AdvancePastCommentLine( ifs );
        size_t records;
        ifs >> records;

        if ( records == 0 ) {
             csmp_error.Note( WARNING, "SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII:",
                                      "Record of variable values appears to be empty." );
             return false;
          }
        if   ( place == ELEMENT ) assert( records == vset.Elements() );
        else if ( place == NODE ) assert( records == vset.Vertices() );
     
        // property data are read from file, flags first, then values
        PropertyData pdata( place, type, dim, array_length );
        pdata.Reserve( records );
        double value;
        
        // flags and values
        if ( type == SCALAR || type == VECTOR || type == FLAGGEDARRAY )
          {
             if      ( type == VECTOR ) records *= dim;
             else if ( type == FLAGGEDARRAY ) records *= array_length;
             // raw insertions 
             for ( size_t i{0U}; i<records; ++i ) {
                  ifs >> integer;
                  pdata.PushBack( intToVARIABLE_FLAG(integer) );
               }
             for ( auto i{0}; i<records; ++i ) {
                  ifs >> value;
                  pdata.PushBack( value );
               }
          }
        // more complex variable types
        else if ( type == TENSOR )
          {
             // dim flags per tensor
             records *= dim;
             // raw insertions 
             for ( size_t i{0U}; i<records; ++i ) {
                  ifs >> integer;
                  pdata.PushBack( intToVARIABLE_FLAG(integer) );
               }
             // dim * dim values per tensor  
             records *= dim;
             for ( size_t i{0U}; i<records; ++i ) {
                  ifs >> value;
                  pdata.PushBack( value );
               }
          }
        else if ( type == ARRAY )
          {
             // raw insertions - one flag value per array
             for ( size_t i{0U}; i<records; ++i ) {
                 ifs >> integer;
                 pdata.PushBack( intToVARIABLE_FLAG(integer) );
               }
             // array_length  values per tensor  
             records *= array_length;
             for ( size_t i{0U}; i<records; ++i ) {
                 ifs >> value;
                 pdata.PushBack( value );
               }
          }
        // adding the property data to the VSet
        vset.AddData( property_name.c_str(), pdata );   
          
      } // end for n properties
      
    return true;

 } // ReadPropertyRecordsASCII


template bool SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII( ifstream&,VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII( ifstream&,VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPropertyRecordsASCII( ifstream&,VSet<3U>& );

/*
cerr <<"\nremaining data in file:\n";    
while ( !ifs.eof() ) {
     char c;
     ifs.get( c );
     cout << c;
  }    
*/    
    





// ----------------------------------------------------------------------------------------------

// BINARY FILE READING

// ----------------------------------------------------------------------------------------------



/**
      Node coordinates
      
      written in Tcl using the formating (sizes in int32)
             
      puts -nonewline $TCLID [binary format d $TCLLine] ; x
      puts -nonewline $TCLID [binary format d $TCLLine] ; y
      puts -nonewline $TCLID [binary format d $TCLLine] ; z

*/
template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( fp == nullptr )
      csmp_error.Note( FATAL_ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary", "file pointer is zero.");

    const size_t  uibytes = sizeof(int32_t);
    const size_t  dbytes  = sizeof(double);
    size_t        entries(0);

    // reading node coordinates
    // -------------------------------

    // reading number of nodes
    assert( uibytes <=sizeof(size_t) );
    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX ); // max size of unsigned 32 bit integer
    if ( csmp_error.Verbose() )
      {
          cout <<"\n\treading "<< entries <<" node coordinates 'px,py,pz'..."<< endl;
          cout.flush();
      }
    // reading node coordinates 'px', 'py', 'pz' (double)
    double* px = new double[ entries ];
    double* py = new double[ entries ];
    double* pz = new double[ entries ];
    if ( fread( (void*) px, dbytes, entries, fp ) != entries )
     throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary",
                              "'px' array did not read correctly");
    if ( fread( (void*) py, dbytes, entries, fp ) != entries )
     throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary",
                              "'py' array did not read correctly");
    if ( fread( (void*) pz, dbytes, entries, fp ) != entries )
     throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary",
                              "'pz' array did not read correctly");

    // resizing the node array
    vset.ResizeNodes( entries );

    // assign node coordinates
cerr <<"\nSKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary: node coordinates:\n";
    for ( size_t i{0U}; i<entries; i++ ) {
cerr <<"\n"<< i<<": "<< px[i] <<","<< py[i] <<","<< pz[i];
         vset.Px( i, px[i] );
         vset.Py( i, py[i] ); 
         vset.Pz( i, pz[i] ); // SKUA models will always have a three coordinate's
      }
    delete[] px;
    delete[] py;
    delete[] pz;

    if ( extra_checks_on_binary_file_ ) {
        // coordinate range checking
        double* xmin=min_element( px, px + entries );
        double* xmax=max_element( px, px + entries );
        double* ymin=min_element( py, py + entries );
        double* ymax=max_element( py, py + entries );
        double* zmin=min_element( pz, pz + entries );
        double* zmax=max_element( pz, pz + entries );

        if ( (*xmax) - (*xmin) < numeric_limits<double>::epsilon() )
          throw csmp::Exception( FATAL_ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary:",
                                              "all nodes have same X-coordinate, but this is a 3D SKUA model.");

        if ( (*ymax) - (*ymin) < numeric_limits<double>::epsilon() )
          throw csmp::Exception( FATAL_ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary:",
                                              "all nodes have same Y-coordinate, but this is a 3D SKUA model.");

        if ( (*zmax) - (*zmin) < numeric_limits<double>::epsilon() )
          throw csmp::Exception( FATAL_ERROR, "SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary:",
                                              "all nodes have same Z-coordinate, but this is a 3D SKUA model.");
      }

    return true;
    
} // end 

template bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary( FILE*,VSet<1U>&);
template bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary( FILE*,VSet<2U>&);
template bool SKUA_FiniteElementMeshInterface::ReadNodeCoordinatesBinary( FILE*,VSet<3U>&);




template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  dbytes  = sizeof(double);

    // 1. reading boundary flags 'pbflags' (int)
    // ------------------------------------------

    if ( csmp_error.Verbose() )
      {
          cout <<"\n\treading boundary flags 'pbflags'..."<< endl;
          cout.flush();
      }
    set<uint32_t>  bnodes;
    int32_t        ival;
    const int32_t  min28(-30), zero(0);
    uint32_t       counter(0);
    size_t nodes(vset.Vertices());
    for ( auto i{0U}; i<nodes; i++ ){
         // TODO: to read models with 10s of billions of cells, ibytes must be int64_t 
         fread( (void*) &ival, ibytes, 1U, fp );
         if ( ival < min28 || ival > zero )
           throw csmp::Exception( ERROR, "SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsBinary","'pbflag' value out of range.");
         // if this is a boundary node
         if ( ival < zero ) {
              vset.BFlag( i, ival );
              bnodes.insert(i);
              counter++;
           }
      }

    // 2. reading boundary condition values 'pbounds' (double) and ignoring them
    // --------------------------------------------------------------------------

    if( csmp_error.Verbose() )
    {
        cout <<"\n\treading "<< counter <<" boundary values 'pbounds'..."<< endl;
        cout.flush();
    }
    double dval;
    for ( size_t i=1; i<=nodes; i++ ) {
         fread( (void*) &dval, dbytes, 1U, fp );
         // only if a boundary flag was stored, a boundary value is stored as well
      }
    bnodes.erase( bnodes.begin(), bnodes.end() );
    return true;
}

template bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsBinary( FILE*,VSet<1U>&);
template bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsBinary( FILE*,VSet<2U>&);
template bool SKUA_FiniteElementMeshInterface::ReadBoundaryFlagsAndConditionsBinary( FILE*,VSet<3U>&);







template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadPelementBinary( FILE* fp, VSet<dim>& vset )
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
                              "SKUA_FiniteElementMeshInterface::ReadPelementBinary",
                             "'pelement' value out of range SKUA-TYPE range (2-23).");

    // adding element types to vset
    vector<int8_t> elmt_types;
    elmt_types.assign( pelmt, pelmt + entries );
    vset.ElementTypes( elmt_types );

    delete[] pelmt;

    return true;
}

template bool SKUA_FiniteElementMeshInterface::ReadPelementBinary( FILE*, VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPelementBinary( FILE*, VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPelementBinary( FILE*, VSet<3U>& );




template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadPlistBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  uibytes = sizeof( uint32_t);
    size_t        entries(0);

    // setting up the storage for 'plist' in VSet
    assert( vset.HybridElementTypeMesh() );
    const size_t  nelements(vset.Cells());

    // now the vset can be resized according to the new information
    deque<uint32_t>  ndele(nelements);
    for ( size_t i{0U}; i<nelements; ++i )
      ndele[i] = csmp_elmt_specs::NodesPerElementOfType( vset.ElementType(i) );
    vset.ResizePlist( ndele );

    // reading nodes connected to elements 'plist' (unsigned int)
    // from the element types array build the plist
    // --------------------------------------------------------------
    // reading the size of the plist array
    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    if( csmp_error.Verbose() )
    {
        cout <<"\n\treading "<< nelements <<" nodes-per-element records from 'plist' (size="<< entries <<")..."<< endl;
        cout.flush();
    }
  
    // reading the plist
    uint32_t*  plist = new uint32_t[ entries ];
    fread( (void*) plist, uibytes, entries, fp );

    deque<vector<int64_t> >::iterator  it(vset.PlistBegin());
    size_t                            nentry(0U);

    // the elements of the plist (node ids) are assigned
    for ( size_t i{0U}; i<nelements; i++, it++ )
      for ( size_t j{0U}; j<ndele[i]; j++ )
          (*it)[j] = plist[nentry++];

    delete[] plist;

    return true;
}

template bool SKUA_FiniteElementMeshInterface::ReadPlistBinary( FILE*, VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPlistBinary( FILE*, VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPlistBinary( FILE*, VSet<3U>& );





template<uint32_t dim>
bool SKUA_FiniteElementMeshInterface::ReadPfvertsBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  uibytes = sizeof( uint32_t);
    size_t        entries(0);

    // setting up the storage for 'pfverts' in VSet
    assert( vset.HybridElementTypeMesh() );
    const size_t   nelements(vset.Cells());

    // making an array of numbers of neighbors of each element
    deque<uint32_t>  nbors( nelements );
    for ( size_t i{0U}; i<nelements; ++i )
      nbors[i] = csmp_elmt_specs::NeighborsPerElementOfType( vset.ElementType(i) );
    vset.ResizePfverts( nbors );

    // reading neighbors connected to elements 'pfverts' (int)
    // -----------------------------------------------------------
    // size of pfverts array
    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    if ( csmp_error.Verbose() ) {
         cout <<"\n\treading "<< nelements <<" neighbor-list records from 'pfverts' (size="<< entries <<")..."<< endl;
         cout.flush();
      }
    // TODO: eventually this must be an array of 'int64_t ' records
    int32_t*  pfverts = new int32_t[ entries ];
    fread( (void*) pfverts, ibytes, entries, fp );

    // reading the C array into the resized pfverts deque inside VData
    // ---------------------------------------------------------------
    deque<vector<int64_t> >::iterator it(vset.PfvertsBegin());
    size_t  nentry(0U);
    for ( size_t i=0U; i<nelements; i++, it++ )
      if (nbors[i] > 2) // SKM: not sure anymore why the restriction was imposed => JC: check it later due to some errors without this restriction especially for fault_boundary_test in BoundaryInterface_Test.
        for ( size_t j{0U}; j<nbors[i]; j++ )
          (*it)[j] = pfverts[nentry++];

    delete[] pfverts;

    return true;
}

template bool SKUA_FiniteElementMeshInterface::ReadPfvertsBinary( FILE*, VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPfvertsBinary( FILE*, VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPfvertsBinary( FILE*, VSet<3U>& );

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
bool SKUA_FiniteElementMeshInterface::ReadPmaterialBinary( FILE* fp, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const size_t  ibytes  = sizeof(int32_t);
    const size_t  uibytes = sizeof( uint32_t);
    size_t        entries(0);

    // reading material information 'pmtrl' (unsigned int)
    // -------------------------------------------------------

    fread( (void*) &entries, uibytes, 1U, fp );
    assert( entries > 0 );
    assert( entries < ULONG_MAX );
    if( csmp_error.Verbose() )
    {
        cout <<"\n\treading "<< entries <<" material-type specifiers for the elements from 'pmtrl'..."<< endl;
        cout.flush();
    }
    int32_t*  mtrls = new int32_t[ entries ];
    fread( (void*) mtrls, ibytes, entries, fp );
    vector<int32_t> elmt_mtrls;
    elmt_mtrls.assign( mtrls, mtrls + entries );
    delete[] mtrls;

    return true;
}

template bool SKUA_FiniteElementMeshInterface::ReadPmaterialBinary( FILE*, VSet<1U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPmaterialBinary( FILE*, VSet<2U>& );
template bool SKUA_FiniteElementMeshInterface::ReadPmaterialBinary( FILE*, VSet<3U>& );




// PLAIN FUNCTIONS


/** returns the corresponding CSMP element type, taking into account whether an isoparametric FEM formulation is used
*/
CSMP_FEM_TYPE  convertSKUA_ElementType( int32_t etype, bool isoparametric )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( isoparametric ) {
         if ( etype == 4 ) return ISOPARAMETRIC_LINEAR_TETRAHEDRON; 					        // TETRA_4          = 4,
         if ( etype == 8 ) return ISOPARAMETRIC_LINEAR_TRIANGLE;  					          // TRI_3            = 8,
         if ( etype == 2 ) return ISOPARAMETRIC_LINEAR_BAR;    						            // BAR_2            = 2,
/*
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_BAR, 						            // BAR_3            = 3,
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_BAR, 
         if ( etype == 2 ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE, 		      // TRI_3_X          = 9,
         if ( etype == 2 ) return SOPARAMETRIC_QUADRATIC_TRIANGLE, 					        // 2D & 3D TRI_6    = 10,
         if ( etype == 2 ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE, 		  // TRI_6_X          = 11,
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_TRIANGLE,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_TETRAHEDRON, 			          // TETRA_10         = 5,
         if ( etype == 2 ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON,    // TETRA_11
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_TETRAHEDRON,
         if ( etype == 2 ) return ISOPARAMETRIC_LINEAR_PYRAMID,						            // PYRA_5           = 18,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_PYRAMID13, 					        // PYRA_13          = 24,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_PYRAMID14,     				      // PYRA_14          = 22,
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_PYRAMID,
         if ( etype == 2 ) return ISOPARAMETRIC_LINEAR_PRISM,  						            // PENTA_6          = 12,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_PRISM15,         			      // PENTA_15         = 13,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_PRISM18,           			    // PENTA_18         = 21,
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_PRISM,
         if ( etype == 2 ) return ISOPARAMETRIC_LINEAR_QUADRILATERAL,                 // QUAD_4           = 14,
         if ( etype == 2 ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL,     // QUAD_4_X         = 15,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL,              // QUAD_8           = 16,
         if ( etype == 2 ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL,  // QUAD_8_X         = 17,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9,  			      // QUAD_9           = 19,
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_QUADRILATERAL,
         if ( etype == 2 ) return ISOPARAMETRIC_LINEAR_HEXAHEDRON,              		  // HEXA_8           = 6,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20,       		    // HEXA_20          = 7,
         if ( etype == 2 ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27, 				      // HEXA_27          = XX,
         if ( etype == 2 ) return ISOPARAMETRIC_CUBIC_HEXAHEDRON,
      */       
         cerr <<"\n\telement type "<< etype << endl;
         csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ConvertSKUA_ElementType",
                           "isoparametric finite element type not recognised" );
      }

   if ( etype == 4 ) return LINEAR_TETRAHEDRON; 					        // TETRA_4          = 4,
   if ( etype == 8 ) return LINEAR_TRIANGLE;  					          // TRI_3            = 8,
   if ( etype == 2 ) return LINEAR_BAR;    						            // BAR_2            = 2,

   cerr <<"\n\telement type "<< etype << endl;
   csmp_error.Note( ERROR, "SKUA_FiniteElementMeshInterface::ConvertSKUA_ElementType",
                     "finite element type not recognised" );

   // types not yet treated
   /*
      UNKNOWN,
      LINEAR_BAR,    										                  // BAR_2            = 2,
      QUADRATIC_BAR, 										                  // BAR_3            = 3,
      CUBIC_BAR, 										                      // BAR_4  
      LINEAR_TRIANGLE,   
      LINEAR_TRIANGLE3D, 									                // TRI_3            = 8,
      LINEAR_QUADRILATERAL,
      LINEAR_CUBOID,
      LINEAR_RECTANGLE,
      BARYCENTRIC_LINEAR_TRIANGLE, 			    	           	// TRI_3_X          = 9,
      QUADRATIC_TRIANGLE, 								                // 2D & 3D TRI_6    = 10,
      BARYCENTRIC_QUADRATIC_TRIANGLE, 					          // TRI_6_X          = 11,
      CUBIC_TRIANGLE,
      LINEAR_TETRAHEDRON, 								                // TETRA_4          = 4,
      QUADRATIC_TETRAHEDRON, 								              // TETRA_10         = 5,
      BARYCENTRIC_QUADRATIC_TETRAHEDRON,                  // PYRA_5           = 18,
      CUBIC_TETRAHEDRON,
      ZERO_DIMENSIONAL_FACE,                              // the face of a line element 
      POINT_ELEMENT,
      POLYGONAL_ELEMENT,
      POLYHEDRAL_ELEMENT,
      EXPERIMENTAL_ELEMENT,
   */
   
   return UNKNOWN;

} // end 



} // end namespace csmp
