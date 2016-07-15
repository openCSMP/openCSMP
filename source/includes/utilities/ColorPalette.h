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
    ColorPalette( float32 low, float32 max );  
    ~ColorPalette();
    ColorPalette( const ColorPalette& a ); 
    ColorPalette& operator=( const ColorPalette& a ); 

    void  GiveRgb( float32 val, float32* vessel ); 
    void  GiveHsv( float32 val, float32* vessel ); 
    void  Saturation( float32 sat );  
    float32 Saturation() const; 
    void  Lightness( float32 light ); 
    float32 Lightness() const;
    void  Blend( float32 alpha ); 
    float32 Blend() const; 
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
    void  RgbToHsv( const std::vector<float32>& rgb, std::vector<float32>& hsv );
    void  HsvToRgb( const std::vector<float32>& hsv, std::vector<float32>& rgb );
    void  RgbToHsv();
    void  HsvToRgb();
    void  InverseRgb( std::vector<float32>& rgb );
    void  InversePalette( std::vector<std::vector<float32> >& rgb );
    void  AdjustPaletteSize();
    void  ScaleColorRangeTo( float32 new_min, float32 new_max );
    void  ColorRangeRGB() const;

  private:
    void RGB_To_HSV( float32 r, float32 g, float32 b, float32& h, float32& s, float32& v );
    void HSV_To_RGB( float32& r, float32& g, float32& b, float32 h, float32 s, float32 v );

    std::vector<std::vector<float32> >  rgbColors, hsvColors;     ///< x*{R,G,B},x*{H,S,V}
    float32                   blend;                              ///< alpha value
    float32                   saturation;
    float32                   lightness;
    const char*               inputPaletteFile;
    const char*               outputPaletteFile;
    size_t                    paletteSize;                        ///< how many colors
    bool                      autoScaling;                        ///< different than 0-255
    float32                   autoScalingLow;  
    float32                   autoScalingHigh;
    float32                   upperPaletteBound;
    float32                   lowerPaletteBound;
};







// --------------------------------------------------------------t----
inline void  ColorPalette::Saturation( float32 sat )
{
  saturation = sat;
}



// ---------------------------------------------------------------t---
inline float32     ColorPalette::Saturation() const
{
  return saturation;
}



// ---------------------------------------------------------------t---
inline void  ColorPalette::Lightness( float32 light )
{
  lightness = light;
}



// ---------------------------------------------------------------t---
inline float32     ColorPalette::Lightness()  const
{
  return lightness;
}



// ---------------------------------------------------------------t---

inline void  ColorPalette::Blend( float32 alpha )
{
  blend = alpha;
}



// ---------------------------------------------------------------t---

inline float32     ColorPalette::Blend()  const
{
  return blend;
}



// ---------------------------------------------------------------t---

inline void ColorPalette::PaletteSize( size_t size )
{
  paletteSize = size;
  AdjustPaletteSize();
}



// ---------------------------------------------------------------t---

inline size_t ColorPalette::PaletteSize()  const
{
  return paletteSize;
}






// --------------------------------------------------------------t----

inline bool ColorPalette::ReadRgbColorPaletteFile( const char* file )
{
  std::string what("rgb");
  return ReadColorPaletteFile( file, what );
}



// ---------------------------------------------------------------t---

inline bool ColorPalette::ReadHsvColorPaletteFile( const char* file )
{
  std::string what("hsv");
  return ReadColorPaletteFile( file, what );
}






// WriteRgbColorPaletteFile()
// --------------------------------------------------------------t----

inline bool ColorPalette::WriteRgbColorPaletteFile()
{
  // Somehow select a filename
  //....
  char file[] = "TestRgbPalette";
  return WriteRgbColorPaletteFile( file ); // "file" should exist so long
}




// ---------------------------------------------------------------t---

inline bool ColorPalette::WriteRgbColorPaletteFile( const char* file )
{
  std::string what("rgb");
  return WriteColorPaletteFile( file, what );
}




// ---------------------------------------------------------------t---

inline bool ColorPalette::WriteHsvColorPaletteFile()
{
  // Somehow select a filename
  //....
  char file[] = "TestHsvPalette";
  return WriteHsvColorPaletteFile( file );
}




// ---------------------------------------------------------------t---

inline bool ColorPalette::WriteHsvColorPaletteFile( const char* file )
{
  std::string what("hsv");
  return WriteColorPaletteFile( file, what );
}






// ---------------------------------------------------------------t---
inline void ColorPalette::RgbToHsv( const std::vector<float32>& rgb, 
                                    std::vector<float32>&       hsv )
{
   RGB_To_HSV( rgb[0], rgb[1], rgb[2], 
               hsv[0], hsv[1], hsv[2] );
}




// --------------------------------------------------------------t----
inline void ColorPalette::HsvToRgb( const std::vector<float32>& hsv, 
                                    std::vector<float32>&       rgb )
{
   HSV_To_RGB( rgb[0], rgb[1], rgb[2], 
               hsv[0], hsv[1], hsv[2] );
}



// --------------------------------------------------------------t----

inline void ColorPalette::RgbToHsv()
{
   std::vector<std::vector<float32> >::const_iterator rgb; 
   std::vector<std::vector<float32> >::iterator       hsv; 

   for( rgb=rgbColors.begin(), hsv=hsvColors.begin(); 
        rgb!=rgbColors.end();  rgb++, hsv++ )
     RgbToHsv( (*rgb), (*hsv) );
}



// ---------------------------------------------------------------t---

inline void ColorPalette::HsvToRgb()
{
   std::vector<std::vector<float32> >::const_iterator hsv; 
   std::vector<std::vector<float32> >::iterator       rgb; 

   for( hsv=hsvColors.begin(), rgb=rgbColors.begin(); 
        hsv!=hsvColors.end(); rgb++, hsv++ )
     HsvToRgb( (*hsv), (*rgb) );
}



// --------------------------------------------------------------t----

inline void ColorPalette::RGB_To_HSV( float32 r, float32 g, float32 b, 
                                      float32& h, float32& s, float32& v )
{
  // Computer Graphics p. 592 ( rgb each in 0..1 )
  float32 max(1.0e30f), min(-1.0e30f);
  if(      r >= g && r >= b ) max = r;
  else if( g >= r && g >= b ) max = g;
  else if( b >= r && b >= g ) max = b;
  if(      r <= g && r <= b ) min = r;
  else if( g <= r && g <= b ) min = g;
  else if( b <= r && b <= g ) min = b;
  float32 delta = max - min;
  v = max;
  s = ( max != 0.0 ) ? ( delta / max ) : 0.0F;
  if( s == 0.0 ) h = 0.0F;
  else
    {
      if(      r == max ) h = ( g - b ) / delta;
      else if( g == max ) h = 2.0F + ( b - r ) / delta;
      else if( b == max ) h = 4.0F + ( r - g ) / delta;
      h *= 60.0F;
      while( h < 0.0 ) h += 360.0F;
    }
}



// HSV_To_RGB()
// ---------------------------------------------------------------t---

inline void ColorPalette::HSV_To_RGB( float32& r, float32& g, float32& b, 
                                      float32 h, float32 s, float32 v )
{
  // Computer Graphics p. 593 ( h in 0..360, s and v in 0..1 )
  float32 f,p,q,t;
  size_t i;
  if( s == 0.0 )
    {
      if( h == 0.0 )
        {
          r = v; g = v; b = v;
        }
      else
        {
          std::cout << "ColorPalette::HSV_To_RGB: if s is zero h must be NULL"<< std::endl;
          return;
        }
    }
  else
    {
      if( h == 360.0F ) h = 0.0F;
      h /= 60.0F;
      i = static_cast<size_t>(std::floor(h));
      f = h - i;
      p = v * ( 1.0F - s );
      q = v * ( 1.0F - ( s * f ) );
      t = v * ( 1.0F - ( s * ( 1.0F -f ) ) );
      switch( i )
        {
	case 0 : r=v; g=t; b=p; break;
        case 1 : r=q; g=v; b=p; break;
        case 2 : r=p; g=v; b=t; break;
        case 3 : r=p; g=q; b=v; break;
        case 4 : r=t; g=p; b=v; break;
        case 5 : r=v; g=p; b=q; break;
        }
    }
}







// GiveRgb( float32 val, float32* vessel )
// ------------------------------------------------------------------
inline void ColorPalette::GiveRgb( float32 val, float32* vessel )
{
  // For OpenGL one needs an array with four values with R,G,B,alpha
  // Provide the array as pointer in the second parameter
  // Use as: glColor4Tv( GiveRgb( someValue, array[for 4] ) );

  // Cast data value to size_t for indexing vector
  size_t which = static_cast<size_t>(val);

  // Safety
  if( which >= paletteSize )
    {
      std::cout << "ColorPalette::GiveRgb\n";
      std::cout << "    (size_t)val was: " << which << " instead of >=0 and < ";
      std::cout << paletteSize << std::endl;
      std::cout << "\nSetting palette values to zero !"<< std::endl;
      for( size_t i=0; i<3; i++ ) vessel[i] = 0.0F; // black
   }

  // Fill in the array
  for( size_t i=0; i<3; i++ ) vessel[i] = rgbColors[ which ][i];
  vessel[3] = blend;
}




// GiveHsv( float32 val, float32* vessel )
// ------------------------------------------------------------------
inline void ColorPalette::GiveHsv( float32 val, float32* vessel )
{
   size_t which = static_cast<size_t>(val);

  // Safety
  if( which >= paletteSize )
    {
      std::cout  << "ColorPalette::GiveRgb\n";
      std::cout  << "    (size_t)val was: " << which << " instead of >=0 and < ";
      std::cout  << paletteSize << std::endl;
      std::cout  << "\nSetting palette values to zero !"<< std::endl;
      for( size_t i=0; i<3; i++ ) vessel[i] = 0.0F; // black
   }
  
   for( size_t i=0; i<3; i++ ) vessel[i] = hsvColors[ which ][i];
   vessel[3] = blend;
}

} // csp

#endif /* ColorPalette.h */
