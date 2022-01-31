#include "GocadPropertyClassHeader.h"
#include "Exception.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp {

GocadPropertyClassHeader::GocadPropertyClassHeader()
 : property_name("undefined"), property_class("undefined"),
   components(1), no_data_value(1.0e-30),
   low_clip(-1.0e-30), high_clip(1.0e+30), 
   p_clip(99),
   colormap("africa"),
   colormap_contrast(0.2)
 {
 }
 
GocadPropertyClassHeader::GocadPropertyClassHeader( const char* name, 
                                                    double lclip, double hclip, 
                                                    long pclip, size_t comp )
 : property_name(name), property_class(name),
   components(comp), no_data_value(1.0e-30),
   low_clip(lclip), high_clip(hclip),  p_clip(pclip),
   colormap("africa"),
   colormap_contrast(0.2)
 {
 } // end ct



GocadPropertyClassHeader&  GocadPropertyClassHeader::operator=( const GocadPropertyClassHeader& ch )
 {
    if ( &ch == this ) return *this;
    property_name     = ch.property_name;
    property_class    = ch.property_class;
    components        = ch.components;
    no_data_value     = ch.no_data_value;
    low_clip          = ch.low_clip;
    high_clip         = ch.high_clip;
    p_clip            = ch.p_clip;  
//    struct colormap:africa
//    struct colormap*contrast:0.2
    colormap          = ch.colormap;
    colormap_contrast = ch.colormap_contrast;
    
    return *this;
 }



 
GocadPropertyClassHeader::GocadPropertyClassHeader( const GocadPropertyClassHeader& ch )
 {
    *this = ch;
 }
 

GocadPropertyClassHeader::~GocadPropertyClassHeader()
 {
 }
 
 
bool GocadPropertyClassHeader::InitializeFrom( ifstream& ifn )
 {
    char  text[256], *token;
    const char* const  delimiter1 = " ,=,\n,\r,:,\t";
     
    // reading rest of line including {
    streampos  pos = ifn.tellg();
    ifn.getline( text, 256 );
    
    if ( strstr( text, "PROPERTY_CLASS_HEADER" ) == NULL )
      {
         ifn.seekg( pos );
         throw csmp::Exception( FATAL_ERROR, "GocadPropertyClassHeader::InitializeFrom", 
                         "No header information found..." );
         return false;
      }
    
    // reading property name before {
    token=strtok( text, delimiter1 );
    token=strtok( NULL, delimiter1 );
    property_name = token;

    // reading all property attributes
    while ( ifn.getline( text, 256 ) )
      {
         token=strtok( text, delimiter1 );
         if ( !strcmp( token, "}" ) ) break;

         while ( token != NULL )
           {
              if ( !strcmp( token, "*low_clip" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   low_clip = atof( token );
                }
              if ( !strcmp( token, "*high_clip" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   high_clip = atof( token );
                }
              if ( !strcmp( token, "*pclip" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   p_clip = atol( token );
                }
              if ( !strcmp( token, "*colormap" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   colormap = token;
                }
              if ( !strcmp( token, "*colormap*contrast" ) )
                {
                   token=strtok( NULL, delimiter1 );
                   colormap_contrast = atof( token );
                }
              token=strtok( NULL, delimiter1 );
           }
      }
    return true;   
      
 } // end


void GocadPropertyClassHeader::WriteToText( ofstream& ofs ) const
 {
    string name(property_name);
    replaceWhiteSpaceBy( name, '_');
    ofs <<"PROPERTY_CLASS_HEADER "<< name <<" {"<< endl;
    ofs <<"*low_clip:"  << low_clip  << endl;
    ofs <<"*high_clip:" << high_clip << endl;
    ofs <<"*pclip:"     << p_clip    << endl;

    // additional information regarding colormap
    // ------------------------------------------
    ofs <<"*colormap:rainbow"                     << endl;
    ofs <<"*colormap*transparency:true"           << endl;
    ofs <<"*colormap*low_clip_transparent:false"  << endl;
    ofs <<"*colormap*transparency_min:0"          << endl;         
    ofs <<"*colormap*transparency_max:1"          << endl;
    ofs <<"*colormap*transparency_power:1"        << endl;
    ofs <<"*colormap*high_clip_transparent:false" << endl;
    ofs <<"*colormap*nodata:false"                << endl;
    ofs <<"*cnp*symbol:point"                     << endl;

    ofs <<"}" << endl;    

 } // end WriteToText




void GocadPropertyClassHeader::Out() const
 {
    cout <<"\nGocadPropertyClassHeader::Out: property name: "<< property_name << endl;
    cout <<"Property class:      "<< property_class << endl;
    cout <<"Property components: "<< components << endl;
    cout <<"Low clip limit:      "<< low_clip   << endl;
    cout <<"High clip limit:     "<< high_clip  << endl;
    cout <<"P Clip limit:        "<< p_clip     << endl;  
    cout <<"Colormap:            "<< colormap   << endl;
    cout <<"Colormap contrast:   "<< colormap_contrast << endl;
    cout <<"ESIZE:               "<< components << endl;
    cout <<"NO_DATA_VALUES:      "<< no_data_value << endl;
    cout.flush();
    
 } // end Out


} // csmp












