// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FracMan_Interface.h"

using namespace std;

namespace csmp {

FRACMAN_Interface::FRACMAN_Interface()
 : n_points(0),
   n_fractures(0),
   n_properties(0),
   scale_factor(1),
   dataset("not initialized")
 {
 }
 
 
FRACMAN_Interface::~FRACMAN_Interface()
 {
 }
 

/**

Reads:

BEGIN FORMAT
    Format = Ascii
    Scale = 100
    No_Fractures = 6
    No_Nodes = 24
    No_Properties = 3
END FORMAT

*/
// returns number of points
long FRACMAN_Interface::ReadFORMAT( ifstream& ifs, long& scale, 
                                    long& fracs, long& props, bool ascii )
 {
    char  text_line[256];
    char* result(0);
    char* token(0);
    const char* const bdelims = ":=#%<>";
    const char* const delims = " ,=#%<>"; 
    string   text;
    long     points;
   
    // 1. searching the BEGIN FORMAT
    while ( ifs.getline( text_line, 256 ) ) {
         token = strtok( text_line, bdelims ); 
         if ( (result=strstr( text_line, "BEGIN FORMAT" )) != NULL ) break; 
      }
    if ( result == NULL ) {
         cout <<"\nFRACMAN_Interface::ReadFORMAT: Could not find format specs."<< endl;
         return 0;
      }

    // 2. reading FORMAT specifications
    // format
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims );
    text  = token;
    assert( text == "Format" );
    token = strtok( NULL, delims );
    text  = token;
    if ( text == "ascii" ) ascii = true;
    else                   ascii = false;
    // XAxis - direction
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims );
    text  = token;
    assert( text == "XAxis" );
    token = strtok( NULL, delims );
    string direc = token; // NOT USED
    // Scale
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims );
    text  = token;
    assert( text == "Scale" );
    token = strtok( NULL, delims );
    scale = atoi(token);
    // No_Fractures
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims );
    text  = token;
    assert( text == "No_Fractures" );
    token = strtok( NULL, delims );
    fracs = atoi(token);
    // No_Nodes
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims );
    text  = token;
    assert( text == "No_Nodes" );
    token = strtok( NULL, delims );
    points = atol(token);
    // No_Properties
    ifs.getline( text_line, 256 );
    token = strtok( text_line, delims );
    text  = token;
    assert( text == "No_Properties" );
    token = strtok( NULL, delims );
    props = atoi(token);

    // END FORMAT
    ifs.getline( text_line, 256 );
    token = strtok( text_line, bdelims );
    if ( (result=strstr( text_line, "END FORMAT" )) == NULL ) {
         cout <<"\nFRACMAN_Interface::ReadFORMAT: Could not find END FORMAT."<< endl;
         return 0;
      }
    
    return points;
    
 } // end 
    
    
/**

Reads:

BEGIN PROPERTIES
    Prop1    =    (Real*4)    "Permeability"
    Prop2    =    (Real*4)    "Compressibility"
    Prop3    =    (Real*4)    "Aperture"
END PROPERTIES
*/
bool  FRACMAN_Interface::ReadPROPERTIES( ifstream& ifs, long n_props, list<string>& props )
 {
    char  text_line[256];
    char* result(0);
    char* token(0);
    const char* const  bdelims = ":,=#%<>\"";
    const char* const  delims = " ,=#%<>\""; 
    string   text;
   
    assert( n_props > 0 );
   
    // 1. searching the BEGIN PROPERTIES
    while ( ifs.getline( text_line, 256 ) ) {
         token = strtok( text_line, bdelims ); 
         if ( (result=strstr( text_line, "BEGIN PROPERTIES" )) != NULL ) break; 
      }
    if ( result == NULL ) {
         cout <<"\nFRACMAN_Interface::ReadPROPERTIES: Could not find property specs."<< endl;
         return 0;
      }

    // 2. reading FORMAT specifications
    if ( !props.empty() ) props.erase( props.begin(), props.end() );
    for ( long i=0; i<n_props; i++ ) {
         ifs.getline( text_line, 256 );
         // propname is ignored
         token = strtok( text_line, delims );
         if ( token == NULL ) return false;
         // number format is ignored
         token = strtok( NULL, delims );
         if ( token == NULL ) return false;
         // name of property is read
         token = strtok( NULL, delims );
         if ( token == NULL ) return false;
         text  = token;
         props.push_back( text );
      }
    
    // 3. reading closing statement of property block
    // END PROPERTIES
    ifs.getline( text_line, 256 );
    strtok( text_line, bdelims );
    if ( (result=strstr( text_line, "END PROPERTIES" )) == NULL ) {
         cout <<"\nFRACMAN_Interface::ReadPROPERTIES: Could not find END PROPERTIES."<< endl;
         return false;
      }
    
    return true;
    
 } // end


void FRACMAN_Interface::InitializeFrom_FRACMAN_File( const char* ffb_file )
 {
    char  text_line[256], infile[200];
    char* result(0); 
    char* token(0);
    const char* const bdelims = ":,=#%<>";
    bool           ascii(true);
    string         frac_name;
    FRACMAN_Fracture  frac;
 
    dataset = ffb_file;
    strcpy( infile, ffb_file );
    strcat( infile, ".ffb" );

    ifstream ifs( infile );
    assert( ifs.is_open() );
   
    // 1. reading the file header 
    n_points = ReadFORMAT( ifs, scale_factor, n_fractures, n_properties, ascii );
    
    // 2. reading property specs
    if ( !ReadPROPERTIES( ifs, n_properties, properties ) ) {
         cout <<"\nFRACMAN_Interface::Read_FFR_File: Could not read PROPERTIES block."<< endl;
         return;
      }

    // 3. reading fracture data
    while ( ifs.getline( text_line, 256 ) ) {
         token = strtok( text_line, bdelims ); 
         if ( (result=strstr( text_line, "BEGIN FRACTURE" )) != NULL ) break; 
      }
    if ( result == NULL ) {
         cout <<"\nFRACMAN_Interface::Read_FFR_File: Could not find FRACTURE specs."<< endl;
         return;
      }
    
    // fracture data
    for ( long i=0; i<n_fractures; i++ ) {
         frac.InitializeFrom( n_properties, ifs );
         frac_name  = "fracture";
         frac_name += frac.TextSetID();
         frac_name += "_";
         frac_name += frac.TextID();
         fractures[ frac_name ] = frac;
         frac_name.erase( frac_name.begin(), frac_name.end() );
      } 

    // 4. reading closing statement of fracture block
    // END FRACTURE
    ifs.getline( text_line, 256 );
    strtok( text_line, bdelims );
    if ( (result=strstr( text_line, "END FRACTURE" )) == NULL )
      cout <<"\nFRACMAN_Interface::Read_FFR_File: Could not find END FRACTURE."<< endl;

    ifs.close();
      
    cout <<"\nFRACMAN_Interface::Read_FFR_File: file '"<< infile;
    cout <<"' read successfully..." << endl;
 }
 
 
 
/// returns volume of box
double FRACMAN_Interface::BoundingBox( Point<3>& pmin, Point<3>& pmax ) const
 {
    map<string,FRACMAN_Fracture,less<string> >::const_iterator  it = fractures.begin();
    list<Point<3>>::const_iterator  pit = (*it).second.Begin();
    
    if ( fractures.empty() ) {
         cout <<"\nFRACMAN_Interface::BoundingBox: Currently no data are stored."<< endl;
         return 0.0;
      }
      
    // 1. Finding extrema of points
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = xmax = (*pit)[0];
    ymin = ymax = (*pit)[1];
    zmin = zmax = (*pit)[2];
    
    for ( it=fractures.begin(); it!=fractures.end(); it++ )
      for ( pit=(*it).second.Begin(); pit!=(*it).second.End(); pit++ ) {
           // X
           if (   (*pit)[0] < xmin  ) xmin = (*pit)[0];
           if (   (*pit)[0] > xmax  ) xmax = (*pit)[0];
           // Y
           if (   (*pit)[1] < ymin  ) ymin = (*pit)[1];
           if (   (*pit)[1] > ymax  ) ymax = (*pit)[1];
           // Z
           if (   (*pit)[2] < zmin  ) zmin = (*pit)[2];
           if (   (*pit)[2] > zmax  ) zmax = (*pit)[2];
        }
    
    // 2. Getting volume of bounding box
    pmin.Set( xmin, ymin, zmin );
    pmax.Set( xmax, ymax, zmax );
    Point<3>  pt = pmax - pmin;
    
    return pt[0] * pt[1] * pt[2];
        
 } // end

 
 
/// places minimum xyz coordinate of all fractures into
/// the origin of the coordinate system
void FRACMAN_Interface::MoveGeometryToOrigin()
 {
    Point<3> pmin, pmax;
 
    BoundingBox( pmin, pmax );

    // displacing all fractures
    map<string,FRACMAN_Fracture,less<string> >::iterator  it;
    
    for ( it=fractures.begin(); it!=fractures.end(); it++ )
      (*it).second.Move( -pmin[0], -pmin[1], -pmin[2] );
 }
 


/// scales geometry by supplied factors
void FRACMAN_Interface::ScaleGeometry( double xfac, double yfac, double zfac )
 {
    map<string,FRACMAN_Fracture,less<string> >::iterator  it;
    
    for ( it=fractures.begin(); it!=fractures.end(); it++ )
      (*it).second.Scale( xfac, yfac, zfac );
 }
 
 
 
/// writes entities only DXF file compatible with DXF 12 or earlier
void FRACMAN_Interface::OutputToDXF( const char* dxf_file ) const
 {
    if ( fractures.empty() ) {
         cout <<"\nFRACMAN_Interface::OutputToDXF: Currently no data to output."<< endl;
         return;
      }

    map<string,FRACMAN_Fracture>::const_iterator  it;
    list<Point<3>>::const_iterator  pit;
    char                               name[100];
    strcpy ( name, dxf_file );
    strcat( name, ".dxf" );

    ofstream  ofs( name );
    assert( ofs );
    ofs.precision(6);
    
    // write first two lines
//    ofs << "999" << endl;      // 0 signifies that something starts
//    ofs << "DXF created by 'FRACMAN_Interface::OutputToDXF'" << endl;
    ofs << "  0" << endl;      // 0 signifies that something starts
    ofs << "SECTION" << endl;
    ofs << "  2" << endl;      // 2 -> block name will follow
    ofs << "ENTITIES" << endl; // only named entities will be written
    
    for ( it=fractures.begin(); it!=fractures.end(); it++ )
      {
         // write beginning of object
         ofs << "  0"        << endl;    
         ofs << "POLYLINE"   << endl;
         ofs << "  8"        << endl;    
         ofs << "DEFAULT"    << endl;
         ofs << "  5"        << endl;    
         ofs << "1E"         << endl;   // color number  
         ofs << " 66"        << endl;    
         ofs << "     1"     << endl;    
         ofs << " 10"        << endl;     
         ofs << "  0.0"      << endl;    
         ofs << " 20"        << endl;    
         ofs << "  0.0"      << endl;    
         ofs << " 30"        << endl;    
         ofs << "  0.0"      << endl;    
         ofs << " 70"        << endl;     
         ofs << "     8"     << endl;     
         
         for ( pit=(*it).second.Begin(); pit!=(*it).second.End(); pit++ )
           {
              ofs << "  0"     << endl;    
              ofs << "VERTEX"  << endl;
              ofs << "  8"     << endl;    
              ofs << "DEFAULT" << endl;
              ofs << "  5"     << endl;    
              ofs << "1F"      << endl;   // color number  
              // writing coordinates
              ofs << " 10"     << endl;  // x-coordinate 
              ofs << (*pit)[0] << endl;
              ofs << " 20"     << endl;  // y-coordinate 
              ofs << (*pit)[1] << endl;
              ofs << " 30"     << endl;  // z-coordinate 
              ofs << (*pit)[2] << endl;  
           }
        
         // repeat first node of closed polyline
         ofs << "  0"     << endl;    
         ofs << "VERTEX"  << endl;
         ofs << "  8"     << endl;    
         ofs << "DEFAULT" << endl;
         ofs << "  5"     << endl;    
         ofs << "1F"      << endl;   // color number  
         // writing coordinates of first node again
         ofs << " 10"     << endl;  // x-coordinate 
         ofs << (*(*it).second.Begin())[0] << endl;
         ofs << " 20"     << endl;  // y-coordinate 
         ofs << (*(*it).second.Begin())[1] << endl;
         ofs << " 30"     << endl;  // z-coordinate 
         ofs << (*(*it).second.Begin())[2] << endl;  
        
         // terminate current face  
         ofs << "  0"        << endl;    
         ofs << "SEQEND"     << endl;  // new thingy
         ofs << "  8"        << endl;  // new thingy
         ofs << "DEFAULT"    << endl;  // new thingy
         ofs << "  6"        << endl;  // new thingy
         ofs << "1E"         << endl;   // color number  
      }
      
    // terminating file
    ofs << "  0"    << endl;    
    ofs << "ENDSEC" << endl;    // end of file
    ofs << "  0"    << endl;    // 2 -> block name will follow
    ofs << "EOF"    << endl;    // end of file

    ofs.close();

    cout <<"\nFRACMAN_Interface::OutputToDXF: file '"<< name;
    cout <<"' written successfully." << endl;

} // end OutputToDXF                   




void FRACMAN_Interface::OutputSelectedRegions( const char* sregions ) const
 {
    if ( fractures.empty() ) {
         cout <<"\nFRACMAN_Interface::OutputSelectedRegions: ";
         cout <<"Currently no regions to output."<< endl;
         return;
      }
    if ( properties.empty() ) {
         cout <<"\nFRACMAN_Interface::OutputSelectedRegions: ";
         cout <<"Currently no property data available for output."<< endl;
         return;
      }

    map<string,FRACMAN_Fracture>::const_iterator  it;
    list<string>::const_iterator               pit;
    list<double>::const_iterator               dit;
    char                                       name[100];
    strcpy( name, sregions );
    strcat( name, "-regions.txt" );

    ofstream  ofs( name );
    assert( ofs );
    ofs.precision(6);
    
    // 1. writing header line
    ofs <<"FRED fracture model: '"<< dataset <<"' output file: '"<< sregions;
    ofs <<"' specifying the properties: "<< endl;
    
    // 2. writing properties
    for ( pit=properties.begin(); pit!=properties.end(); pit++ )
      ofs << (*pit) <<"  ";
    ofs << endl;

    // 3. writing the region=fracture names and the associated property values
    for ( it=fractures.begin(); it!=fractures.end(); it++ )
      {
         ofs << (*it).first <<"  ";
         for ( dit=(*it).second.PropertiesBegin(); 
               dit!=(*it).second.PropertiesEnd(); dit++ )
           ofs << (*dit) <<"  ";
         ofs << endl;
      }

    ofs.close();

    cout <<"\nFRACMAN_Interface::OutputSelectedRegions: file '"<< name;
    cout <<"' written successfully." << endl;

} // end OutputSelectedRegions



/// outputs fracture numbers and diameters to text
void FRACMAN_Interface::OutputSelectedFractureDiameters( const char* sregions ) const
 {
    if ( fractures.empty() ) {
         cout <<"\nFRACMAN_Interface::OutputSelectedFractureDiameters: ";
         cout <<"Currently no regions to output."<< endl;
         return;
      }
    if ( properties.empty() ) {
         cout <<"\nFRACMAN_Interface::OutputSelectedFractureDiameters: ";
         cout <<"Currently no property data available for output."<< endl;
         return;
      }

    map<string,FRACMAN_Fracture>::const_iterator  it;
    list<string>::const_iterator               pit;
    list<double>::const_iterator               dit;
    char                                       name[100];
    strcpy( name, sregions );
    strcat( name, "-diameters.txt" );

    ofstream  ofs( name );
    assert( ofs );
    ofs.precision(6);
    
    // 1. writing header line
    ofs <<"FRED fracture model: '"<< dataset <<"' output file: '"<< sregions;
    ofs <<"' N fracture, N perimeter points, perimeter and diameter assuming circular shape. "<< endl;
    
    // 2. writing the region=fracture names and the associated property values
    for ( it=fractures.begin(); it!=fractures.end(); it++ )
      {
         ofs << (*it).second.ID() <<"\t"<< (*it).second.PolygonPoints() <<"\t";
         ofs << (*it).second.Perimeter() <<"\t"<< (*it).second.Diameter() << endl;
      }

    ofs.close();

    cout <<"\nFRACMAN_Interface::OutputSelectedFractureDiameters: file '"<< name;
    cout <<"' written successfully." << endl;

} // end OutputSelectedFractureDiameters




void FRACMAN_Interface::Out() const
 {
    map<string,FRACMAN_Fracture,less<string> >::const_iterator  it;
    list<string>::const_iterator                             pit;
    
    // 1. writing header line
    cout <<"\nFRACMAN_Interface::Out: FRED fracture model: '"<< dataset;
    cout <<"' specifying the properties: "<< endl;
    
    // 2. writing properties
    for ( pit=properties.begin(); pit!=properties.end(); pit++ )
      cout << (*pit) <<"  ";
    cout << endl;

    // 3. writing the region=fracture names and the associated property values
    for ( it=fractures.begin(); it!=fractures.end(); it++ ) {
         cout <<"\n\nFracture: "<< (*it).first << endl;
         (*it).second.Out();
      }
    cout << endl;

 } // end Out

} // end namespace csmp





