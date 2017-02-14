#ifndef GOCAD_HEADER_H
#define GOCAD_HEADER_H

#include "CSMP_definitions.h"
#include "GocadPropertyClassHeader.h"

namespace csmp {

class GocadHeader {
    std::string  gocad_type;
    std::string  name;
    std::string  geological_type;
    bool         cn;    // control nodes y/n
    bool         solid; 
    long         sample_size;
    double64    no_data_values;
    std::map<std::string,GocadPropertyClassHeader>  properties;

  public:
    GocadHeader(); 
    GocadHeader( const char* gocad_type, const char* name, const char* geological_name, 
                 bool issolid=true ); 
                 
    GocadHeader( const GocadHeader& h );
    explicit GocadHeader( std::ifstream& ifn );
    GocadHeader&  operator=( const GocadHeader& h );
    bool          operator<( const GocadHeader& h ) const;
    ~GocadHeader();
    
    // methods
    std::string  Name() const { return name; };
    void        Name( const std::string& s ) { name = s; };
    void        Name( char* s ) { name = s; };
    
    bool        InitializeFrom( std::ifstream& ifs );
    void        WriteToText( std::ofstream& ofs ) const;
   
    void        AddProperty( const char* pname, double64 loclip, double64 hiclip, 
                             int32 pclip=99, int32 comp=1 );
    
    void Out() const { Out(std::cout); }
    void Out(std::ostream& os) const;
};      

} // csp

#endif

// potential extra attributes
/*

*painted
*vectors3d
*vectors3d*arrow:false
*painted*variable:Z
*shaded:true
*links:true
*solid:true
*mesh:false

*/
