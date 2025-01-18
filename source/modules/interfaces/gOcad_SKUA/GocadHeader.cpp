#include "GocadHeader.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;


namespace csmp {

GocadHeader::GocadHeader()  // default
  : gocad_type("TSolid"),
    name("CSP_object"),
    geological_type("rock"),
    cn(true),
    solid(false),
    sample_size(0),
    no_data_values(1.0e-30)
//    properties;
 {
 }


/**
    Checks header for specifications.
    
    'gocad_type' is the specific gocad type name, e.g. TSolid or so.
    'geological_type' is user-defined special name
*/
GocadHeader::GocadHeader( const char* gocad_name, const char* name, const char* geo_type, bool issolid )
 : gocad_type(gocad_name), name(name), geological_type(geo_type), solid(issolid)
 {
    if ( string(gocad_name) != "TSolid" &&  string(gocad_name) != "TSurf" && string(gocad_name) != "PLine" )
    cout <<"\nGocadHeader::(ct): WARNING gocad_name does not refer to a gocad object..."<< endl;
 }




GocadHeader::GocadHeader( const GocadHeader& h )  // default
 {
    *this = h;
 }




GocadHeader::~GocadHeader() 
 {
 }
 
 
 
 
GocadHeader&  GocadHeader::operator=( const GocadHeader& h )
 {
    if ( &h == this ) return *this;
    gocad_type      = h.gocad_type;
    name            = h.name;
    geological_type = h.geological_type;
    cn              = h.cn;
    solid           = h.solid;
    sample_size     = h.sample_size;
    no_data_values  = h.no_data_values;
    properties      = h.properties;
    
    return *this;
    
 } // end assignment    



/// assumes that "HEADER" has already been read
GocadHeader::GocadHeader( ifstream& ifn )
 {
    char  text[256], *token;
    const char* const delimiter1 = "\n,\r,:";
    
    // look at line remembering the line position
    streampos  pos = ifn.tellg();
    ifn.getline( text, 256 );

    if ( strstr( text, "HEADER {" ) == NULL )
      {
         cout <<"\nGocadHeader::InitializeFrom: no header information found..."<< endl;
         ifn.seekg( pos );
         return;
      }

    while ( ifn.getline( text, 256 ) )
      {
         token=strtok( text, delimiter1 );
         if ( !strcmp( token, "}" ) ) break;

         while ( token != NULL )
           {
              if ( !strcmp( token, "name" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   name = token;
                }
              if ( !strcmp( token, "*solid" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   if ( !strcmp( token, "true" ) ) solid = true;
                   else                            solid = false;
                }
             token=strtok( NULL, delimiter1 );
            }
      }
 } // ct






/// name string length comparison
bool  GocadHeader::operator<( const GocadHeader& h ) const
 {
     if ( name < h.name ) return true;
     return false;
 }
 
 
bool  GocadHeader::InitializeFrom( ifstream& ifn )
 {
    char  text[256], *token;
    const char* const delimiter1 = " ,\n,\r,:";
    
    // look at line remembering the line position
    streampos  pos = ifn.tellg();
    ifn.getline( text, 256 );

    if ( strstr( text, "HEADER {" ) == NULL )
      {
         cout <<"\nGocadHeader::InitializeFrom: no header information found..."<< endl;
         ifn.seekg( pos );
         return false;
      }

    while ( ifn.getline( text, 256 ) ) {
         token=strtok( text, delimiter1 );
         if ( !strcmp( token, "}" ) ) break;

         while ( token != NULL )
           {
              if ( !strcmp( token, "name" ) ) {
                   token=strtok( NULL, delimiter1 );
                   name = token;
                }
              if ( !strcmp( token, "*solid" ) ) {
                   token=strtok( NULL, delimiter1 );
                   if ( !strcmp( token, "true" ) ) solid = true;
                   else                            solid = false;
                }
             token=strtok( NULL, delimiter1 );
          }
      }
    return true;   
      
 } // end
 
 
 
void  GocadHeader::AddProperty( const char* pname, double lclip, double hclip, int32_t pclip, int32_t comp )
 {
   string temp(pname);
   replaceWhiteSpaceBy( temp, '_' );

   properties[ temp ] = GocadPropertyClassHeader( temp.c_str(), lclip, hclip, pclip, comp );

 } // end AddProperty
    


void GocadHeader::WriteToText( ofstream& ofs ) const
 {
    // 1. the first header block
    // -------------------------
    const double TOLERANCE = 0.0001,
                 NO_VALUE  = 1.0e-30;

    //              e.g. TSolid
    ofs <<"GOCAD" <<" "<< gocad_type <<" "<< TOLERANCE << endl;
    ofs <<"HEADER {"<< endl;
    ofs <<"name:"<< name << endl;
    if ( solid ) ofs <<"*solid:true"<<  endl;
    else         ofs <<"*solid:false"<< endl;

    // attributes that determine how the object will be drawn initially
    // ----------------------------------------------------------------
    ofs <<"*painted:true"          << endl;
    ofs <<"*painted*variable:"<< (*properties.begin()).first << endl;
    ofs <<"*mesh*color:gray30"     << endl;
    ofs <<"*cn:true"               << endl;
    ofs <<"*cn*size:0.2"           << endl;
    ofs <<"*cn*color:white"        << endl;
    ofs <<"*under_threshold:false" << endl;
    ofs <<"*border_only:false"     << endl;
    ofs <<"*color_coded:false"     << endl;
    ofs <<"*normals:true"          << endl;
    ofs <<"*solid*transparency:0"  << endl;
    ofs <<"*shrink_coef:1"         << endl;
    ofs <<"*solid*specular:gray70" << endl;
    ofs <<"}"<< endl;


    // 2. properties information
    // -------------------------
    map<string,GocadPropertyClassHeader>::const_iterator  it;

    ofs <<"PROPERTIES";
    for ( it=properties.begin(); it!=properties.end(); it++ )
      ofs <<" "<< (*it).first;
    ofs << endl;

    ofs <<"NO_DATA_VALUES";
    for ( size_t i{0U}; i<properties.size(); i++ ) ofs <<" "<< NO_VALUE;
    ofs << endl;

    ofs <<"PROPERTY_CLASSES";
    for ( it=properties.begin(); it!=properties.end(); it++ )
      ofs <<" "<< (*it).first;
    ofs << endl;

    // the dimensionality of the property variables
    ofs <<"ESIZES"; // so far only scalars
    for ( it=properties.begin(); it!=properties.end(); it++ ) 
      ofs <<" "<< (*it).second.Components();
    ofs << endl;

    // 3. Property class headers
    // -------------------------
    for ( it=properties.begin(); it!=properties.end(); it++ )
      (*it).second.WriteToText( ofs );

    // 4. Label of the type section
    // ----------------------------
    if      ( !strcmp("TSolid", gocad_type.c_str()) ) ofs <<"TVOLUME"<< endl;
    else if ( !strcmp("TSurf",  gocad_type.c_str()) ) ofs <<"TFACE"<<   endl;
    else if ( !strcmp("PLine",  gocad_type.c_str()) ) ofs <<"ILINE"<<   endl;
    else
    cout <<"\nGocadHeader::WriteToText: type label (name) not implemented yet."<< endl;

 } // end WriteToText



void GocadHeader::Out() const
 {
    cout <<"\nGocadHeader::Out: name: "<< name << endl << endl;
    cout <<"Geological type of object: "<< geological_type << endl;
    cout <<"With control nodes:        "<< cn << endl;   
    cout <<"Solid:                     "<< solid << endl;
    cout <<"Sample size:               "<< sample_size << endl;
    cout <<"No data value threshold:   "<< no_data_values << endl;
    
    if ( !properties.empty() ) cout <<"\nAssociated properties:" << endl;
    map<string,GocadPropertyClassHeader>::const_iterator  it;
    for ( it=properties.begin(); it!=properties.end(); it++ )
      {
         cout << (*it).first << endl;
         (*it).second.Out();
      }
    cout.flush();
      
 } // end 
 
 
} // csp
 
 
 
 
