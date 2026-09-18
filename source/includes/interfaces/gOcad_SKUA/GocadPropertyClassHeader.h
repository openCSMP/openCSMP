// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GOCAD_PROPERTY_CLASS_HEADER_H
#define GOCAD_PROPERTY_CLASS_HEADER_H

#include "CSMP_definitions.h"

namespace csmp {

class  GocadPropertyClassHeader {
    std::string  property_name, 
                 property_class;
    size_t    components;  // ESIZE
    double    no_data_value;
    double    low_clip;    // 1e-30
    double    high_clip;   // 274.403
    long         p_clip;      // :99
    std::string  colormap;
    double    colormap_contrast;
    
  public:
    GocadPropertyClassHeader();
    GocadPropertyClassHeader( const char* name, 
                              double lclip, double hclip, long pclip=99, 
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
    
    double   LowClip() const { return low_clip; };
    void        LowClip( double lc ) { low_clip = lc; };
    
    double   HighClip() const { return high_clip; };
    void        HighClip( double hc ) { high_clip = hc; };
    
    long        PClip() const { return p_clip; };
    void        PClip( long pc ) { p_clip = pc; };
    
    double   ColorContrast() const { return colormap_contrast; };
    void        ColorContrat( double ct ) { colormap_contrast = ct; };
    
    std::string ColorMap() const { return colormap; };
    void        ColorMap( const std::string& s ) { colormap = s; };
    void        ColorMap( char* s ) { colormap = s; };
    
    size_t   Components() const { return components; };
    
    size_t   ESize() const { return components; };
    void        ESize( size_t esize ) { components = esize; };
    
    double   NoDataValue() const { return no_data_value; };
    void        NoDataValue( double val ) { no_data_value = val;  };
    
    bool        InitializeFrom( std::ifstream& ifn );
    void        WriteToText( std::ofstream& ofs ) const;

    void Out() const;
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
