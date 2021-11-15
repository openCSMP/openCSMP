#ifndef CSMP_COLOR_PALETTE_H
#define CSMP_COLOR_PALETTE_H

#include "CSMP_definitions.h"

namespace csmp {

/**

@brief Implementation of an RGB color table that is used in the creation
of JPG output files. Used either to manage an array of RGBA values
or as a runtime calculater of RGB colors.

16777216 (256^3) colors.
This class should hold 256 RGB values of the template type // WHY 256?????
It can take and give these values from the program or a file
The palette can be composed with different saturations in
different ranges of the spectrum. There is nothing it cannot do.


@section goal Goal
To serve OpenGL command: glColor4fv( array ); // or glColor4dv(a);

@section vars Main variables

- 2D vector 256x4: 256 colors 3 RGB and one for alpha blending
   muss ev array sein statt vector damit es direct an OpenGL
   uebergeben werden kann. Or one may be able to give an iterator?
- file name for color palette file
- saturation, default 1
- spectrum range, default all
- selected color
- number of colors in palette


@section methodas Main methods

- Constructor creates standard rainbow with alpha 1
- Calculate palette: default rainbow spectrum with full saturation
- Read in color palette from file
- Read in color palette from program
- Give a pointer to the the palette array
- Give a pointer or iterator to the xxx. color
- Give a pointer to the the closest color to a value between 0..255
- Set or read artificial range instead of 0..255 (e.g. -200..+50000)
- Give a pointer to closest color in artifical range
- Set or read saturation
- Set or read spectrum range
- Set or read every x color y number of colors = black
- Reverse the palette
- Graphically select a color        and what belongs to it...
- Graphically select a color range
- Change number of colors of palette
- Transform color models: RgbToHsv(); HsvToRgb();
- Inverse an RGB color

*/
class ColorPalette {
  public:
    ColorPalette();
    /// for autoscaling
    ColorPalette( float low, float max );  
    ~ColorPalette();
    ColorPalette( const ColorPalette& a ); 
    ColorPalette& operator=( const ColorPalette& a ); 

    void  GiveRgb( float val, float* vessel ); 
    void  GiveHsv( float val, float* vessel ); 
    void  Saturation( float sat );  
    float Saturation() const; 
    void  Lightness( float light ); 
    float Lightness() const;
    void  Blend( float alpha ); 
    float Blend() const; 
    void  PaletteSize( size_t size ); 
    size_t  PaletteSize() const; 
    void  MakeRainbowPalette();
    void  MakeGreyPalette();
    bool  ReadColorPaletteFile( const char* file, const std::string& what );
    bool  ReadRgbColorPaletteFile( const char* file ); 
    bool  ReadHsvColorPaletteFile( const char* file );
    bool  WriteColorPaletteFile( const char* file, const std::string& what ); 
    bool  WriteRgbColorPaletteFile( const char* file ); 
    bool  WriteRgbColorPaletteFile(); 
    void  WriteRgbColorPaletteToStdout();
    void  WriteRgbColorPaletteToStdoutWithLineNumbers(); 
    bool  WriteHsvColorPaletteFile();
    bool  WriteHsvColorPaletteFile( const char* file );
    void  WriteHsvColorPaletteToStdout(); 
    void  WriteHsvColorPaletteToStdoutWithLineNumbers(); 
    void  RgbToHsv( const std::vector<float>& rgb, std::vector<float>& hsv );
    void  HsvToRgb( const std::vector<float>& hsv, std::vector<float>& rgb );
    void  RgbToHsv();
    void  HsvToRgb();
    void  InverseRgb( std::vector<float>& rgb );
    void  InversePalette( std::vector<std::vector<float> >& rgb );
    void  AdjustPaletteSize();
    void  ScaleColorRangeTo( float new_min, float new_max );
    void  ColorRangeRGB() const;

  private:
    void RGB_To_HSV( float r, float g, float b, float& h, float& s, float& v );
    void HSV_To_RGB( float& r, float& g, float& b, float h, float s, float v );

    std::vector<std::vector<float> >  rgbColors, hsvColors;     ///< x*{R,G,B},x*{H,S,V}
    float                   blend;                              ///< alpha value
    float                   saturation;
    float                   lightness;
    const char*               inputPaletteFile;
    const char*               outputPaletteFile;
    size_t                    paletteSize;                        ///< how many colors
    bool                      autoScaling;                        ///< different than 0-255
    float                   autoScalingLow;  
    float                   autoScalingHigh;
    float                   upperPaletteBound;
    float                   lowerPaletteBound;
};

} // csmp

#endif /* ColorPalette.h */
