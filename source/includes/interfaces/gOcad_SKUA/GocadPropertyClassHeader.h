#ifndef GOCAD_PROPERTY_CLASS_HEADER_H
#define GOCAD_PROPERTY_CLASS_HEADER_H

#include "CSMP_definitions.h"

namespace csmp {

class  GocadPropertyClassHeader {
    std::string  property_name, 
                 property_class;
    size_t    components;  // ESIZE
    double64    no_data_value;
    double64    low_clip;    // 1e-30
    double64    high_clip;   // 274.403
    long         p_clip;      // :99
    std::string  colormap;
    double64    colormap_contrast;
    
  public:
    GocadPropertyClassHeader();
    GocadPropertyClassHeader( const char* name, 
                              double64 lclip, double64 hclip, long pclip=99, 
                              size_t comp=1 );
                              
    GocadPropertyClassHeader( const GocadPropertyClassHeader& p );
    ~GocadPropertyClassHeader();
    GocadPropertyClassHeader&  operator=( const GocadPropertyClassHeader& p );
    
    
    std::string Property() const { return property_name; };
    void        Property( const std::string& s ) { property_name = s; };
    void        Property( const char* s ) { property_name = s; };
    
    std::string PropertyClass() const { return property_class; };
    void        PropertyClass( const std::string& s ) { property_class = s; };
    void        PropertyClass( const char* s ) { property_class = s; };
    
    double64   LowClip() const { return low_clip; };
    void        LowClip( double64 lc ) { low_clip = lc; };
    
    double64   HighClip() const { return high_clip; };
    void        HighClip( double64 hc ) { high_clip = hc; };
    
    long        PClip() const { return p_clip; };
    void        PClip( long pc ) { p_clip = pc; };
    
    double64   ColorContrast() const { return colormap_contrast; };
    void        ColorContrat( double64 ct ) { colormap_contrast = ct; };
    
    std::string ColorMap() const { return colormap; };
    void        ColorMap( const std::string& s ) { colormap = s; };
    void        ColorMap( char* s ) { colormap = s; };
    
    size_t   Components() const { return components; };
    
    size_t   ESize() const { return components; };
    void        ESize( size_t esize ) { components = esize; };
    
    double64   NoDataValue() const { return no_data_value; };
    void        NoDataValue( double64 val ) { no_data_value = val;  };
    
    bool        InitializeFrom( std::ifstream& ifn );
    void        WriteToText( std::ofstream& ofs ) const;

    void Out(std::ostream& os) const;
};

} // csp

#endif


// Let this class write info directly to file

/*
GEOLOGICAL_TYPE top
PROPERTIES thickness
NO_DATA_VALUES 1e-30 
PROPERTY_CLASSES thickness
ESIZES 1 
PROPERTY_CLASS_HEADER Z {
*low_clip:731.52
*high_clip:731.52
*pclip:99
}
PROPERTY_CLASS_HEADER thickness {
*low_clip:1e-30
*high_clip:274.403
*pclip:99
*colormap:africa
*colormap*contrast:0.2
}
*/
