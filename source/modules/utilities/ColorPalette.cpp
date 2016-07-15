#include "ColorPalette.h"

using namespace std;

namespace csmp {

// ColorPalette()
// ------------------------------------------------------------------
ColorPalette::ColorPalette()
: inputPaletteFile(NULL),
  outputPaletteFile(NULL),
  autoScalingLow( 0.0F ),
  autoScalingHigh( 255.0F ),
  autoScaling( false ),
  paletteSize(256),
  upperPaletteBound( 240.0F ),
  lowerPaletteBound(   0.0F ),
  blend(1.0F),
  saturation(1.0F),
  lightness(255.0F),
  rgbColors( 256, vector<float32>(3) ), //  Don't use rgbColors( 256, 3 ), 3=float is not defined !
  hsvColors( 256, vector<float32>(3) )  
{
  // Make a default color palette
  cout <<"\nColorPalette::ColorPalette: building palettes"<< endl;
  cout.flush();
  paletteSize = 256;
  MakeRainbowPalette();
}


// ColorPalette( float32 low, float32 high )
// ---------------------------------------------------------------t---
ColorPalette::ColorPalette( float32 low, float32 high )
: inputPaletteFile(NULL),
  outputPaletteFile(NULL),
  autoScalingLow( low ),
  autoScalingHigh( high ),
  autoScaling( true ),
  paletteSize(256),
  upperPaletteBound( 240.0F ),
  lowerPaletteBound(   0.0F ),
  blend(1.0F),
  saturation(1.0F),
  lightness(255.0F),
  rgbColors( 256, vector<float32>(3) ), 
  hsvColors( 256, vector<float32>(3) ) 
{
  // Make a default color palette
  paletteSize = 256;
  cout <<"\nColorPalette::ColorPalette(float32 low, float32 high): building palettes"<< endl;
  cout.flush();
  MakeRainbowPalette();
  //MakeGreyPalette();
}


// ---------------------------------------------------------------t---
ColorPalette::ColorPalette( const ColorPalette& clr )
 : inputPaletteFile(NULL),
   outputPaletteFile(NULL)
{
   *this = clr;
}

// ---------------------------------------------------------------t---
ColorPalette::~ColorPalette()
{
}




// assignment operator
// --------------------------------------------------------------t----
ColorPalette& ColorPalette::operator=( const ColorPalette& a )
{
  if ( &a != this ) {
  // copy commands ...
  blend              = a.blend;                // GL alpha value
  saturation         = a.saturation;
  lightness          = a.lightness;
  inputPaletteFile   = NULL;
  outputPaletteFile  = NULL;
  paletteSize        = a.paletteSize;          // how many colors
  autoScaling        = a.autoScaling;          // different than 0-255
  autoScalingLow     = a.autoScalingLow;  
  autoScalingHigh    = a.autoScalingHigh;
  upperPaletteBound  = a.upperPaletteBound;
  lowerPaletteBound  = a.lowerPaletteBound;
  rgbColors          = a.rgbColors;
  hsvColors          = a.hsvColors;    
    }
  return *this;
}



// AdjustPaletteSize()
// ---------------------------------------------------------------t---
void ColorPalette::AdjustPaletteSize()
{
  // Only the size of the palette is adjusted. But later on one should
  //   also spread the old spectrum over the new range

  // Adapt the size of the palette
  size_t  newSize = paletteSize;

  if( newSize > hsvColors.size() )
    {
      // Enlarge vector
      hsvColors.reserve( newSize );
      rgbColors.reserve( newSize );
      while( newSize > hsvColors.size() ) hsvColors.push_back( vector<float32>(3) );
      while( newSize > rgbColors.size() ) rgbColors.push_back( vector<float32>(3) );
      assert( newSize == hsvColors.size() );
      assert( newSize == rgbColors.size() );
    }
  else if( newSize < hsvColors.size() )
    {
      // Shrink vector
      hsvColors.erase( hsvColors.begin()+static_cast<long>(newSize), hsvColors.end() );
      rgbColors.erase( rgbColors.begin()+static_cast<long>(newSize), rgbColors.end() );
      assert( newSize == hsvColors.size() );
      assert( newSize == rgbColors.size() );
    }
}



// MakeRainbowPalette()
// --------------------------------------------------------------t----
void ColorPalette::MakeRainbowPalette()
{
  vector<vector<float32> >::iterator  hsvIt;
  // Rainbow is easier to make in HSV
  size_t i;

  // Define range of palette

  // defaults: upperPaletteBound 240.0, lowerPaletteBound 0.0
  //upperPaletteBound = 359.8;
  //lowerPaletteBound = 0.0;

  float32 paletteRangeInDegree = upperPaletteBound - lowerPaletteBound;
  assert( paletteSize != 0 ); // to prevent division by zero
  float32   step = paletteRangeInDegree / paletteSize;

  // maybe a check for the range would be handy...

  for( hsvIt  = hsvColors.begin(), i=0;
       hsvIt != hsvColors.end(); hsvIt++, i++ )
    {
       (*hsvIt)[0] = upperPaletteBound - i * step; // hue between 0.0 and 360.0 degrees
       (*hsvIt)[1] = saturation; // saturation
       (*hsvIt)[2] = lightness;  // value
    }

  // Switch the palette to rgb
  HsvToRgb();

  // Scale Palette
  if( autoScaling == true )
      ScaleColorRangeTo( autoScalingLow, autoScalingHigh );
}



// MakeGreyPalette()
// --------------------------------------------------------------t----
void ColorPalette::MakeGreyPalette()
{
  vector<vector<float32> >::iterator  rgbIt;
  // Grey scaling is easier to make in RGB
  size_t i,j;
  
  assert( paletteSize != 0 );
  float32   step = 254.8F / paletteSize;
  for( rgbIt  = rgbColors.begin(), i=0;
       rgbIt != rgbColors.end(); rgbIt++, i++ )
    {
      j = static_cast<size_t>(i * step);
      (*rgbIt)[0] = j;         // red
      (*rgbIt)[1] = j;         // green
      (*rgbIt)[2] = j;         // blue    
    }
  RgbToHsv();

  // Scale Palette
  if( autoScaling == true )
      ScaleColorRangeTo( autoScalingLow, autoScalingHigh );
}




// ReadColorPaletteFile()
// ------------------------------------------------------------------
bool ColorPalette::ReadColorPaletteFile( const char* file, const string& what )
{
  // Read either an RGB or HSV file
  // Format: Comment lines are any with something additional to 0-9 and "."
  //         Empty lines
  //         Lines with 3 numbers of type size_t or float32
  //         The number of lines with values will be the size of the
  //             color palette
  // Parameters: char* file = "filename", String what = "rgb" OR "hsv"


  inputPaletteFile = file;
  const size_t    maxIntext = 1024;   // length of line accepted
  size_t          maxLines = 16384;   // = max palette size
  size_t          minLines = 256;     // = min palette size
  float32            minAllowed = 0.0F;  // for rgb
  float32            maxAllowed = 255.0F; // for hsv
  ifstream           ifs;                // palette file stream
  vector<float32> color(3); 
  size_t          i = 0;              // i = i :-)
  char*        p;
  char         intext[maxIntext];  // for input line
  size_t    item = 0;
  string       field[3];           // for the three values per line
  size_t    where = 0;          // 0 in space, 1 in word
  size_t    colorCounter = 0;   // number of colors in file
  size_t    lineCounter = 0;    // physical lines

  vector<vector<float32> >  allColors;   // tmp vector for all colors
  allColors.reserve( minLines );
  for( i=0; i<3; i++ ) color[i] = 0.0F;


  // Check for rgb or hsv in second parameter
  if( what != "rgb" && what != "hsv" )
    {
      cout << "ColorPalette::ReadColorPaletteFile:\n"
           << "    Second parameter must be rgb or hsv and not:"
           << what << endl;
      return false;
    } 

  // Opening the input file
  ifs.open ( inputPaletteFile );
  if ( ! ifs )
    { 
      cout << "ColorPalette::ReadColorPaletteFile:\n"
           << "     Can't find color input file: " << inputPaletteFile << endl; 
      return false;
    }    



  // Parse the input file linewise
  NEXTLINE:
  while( ifs.getline( intext, maxIntext ) )
    {
      lineCounter++; // first line is now line 1

      // Discard empty lines
      size_t len = strlen( intext );
      if( len == 0 ) goto NEXTLINE;

      // Understand if line too long
      if( len >= maxIntext - 1 )
        {
          cout << "ColorPalette::ReadColorPaletteFile:"
               << "    Line " << lineCounter 
               << " too long in palette file" << endl;
          return false;
        }

      // Pick up the 3 numbers out of the line
      item  = 0;                                // for value 1 of 3
      where = 0;                                // where? in space...
      for( i=0; i<3; i++ ) field[i] = "";       // set the Strings to ""

      for( p = intext; *p != '\0'; p++ )   // iterate over characters in line
        {

          // It's a space ////////////////////////////////////////////
          if( isspace(*p) )
            {
              switch( where )
                {
		case 0 :
                         continue;
                case 1 :
                         where = 0; 
                         item++;
                         break;
                default: cout << "ColorPalette::ReadColorPaletteFile:"
                              << "   \"where\" error. It is: " << where
                              << endl;
                        return false;
                }
            }

          // It's a number or a dot ///////////////////////////////////////
          else if( isdigit(*p) || *p == '.' )
            {
              where = 1;
              if( item > 2 )
                {
                  cout << "ColorPalette::ReadColorPaletteFile:\n"
                       << "    More than 3 values on line: "
                       << intext << endl;
                       return false;
                }
              field[item] += static_cast<char>(*p);
            }

          // It's neither a number nor a dot /////////////////////////
          else
            {
              goto NEXTLINE;
            }
        } // end of processing the line with for()


      // Check whether we really got 3 items
      for( i=0;i<3;i++ )
        {
          if( field[i] == "" )
            {
              cout << "ColorPalette::ReadColorPaletteFile:\n"
                   << "    Did not get three values. Missing value: "
                   << i << " on line: " << lineCounter << endl;
              return false;
            }
        }


      // Now transform String to T
      for( i=0; i<3; i++ )
        {
          // Make a float32 out of String
          color[i] = static_cast<float32>(atof( field[i].c_str() ));
          // Check for value in range
          if( what == "rgb" )
            {
              if( maxAllowed < color[i]  || minAllowed > color[i]  )
                {
                  cout << "ColorPalette::ReadColorPaletteFile:\n"
                       << "   Palette value out of range: " 
                       << minAllowed << " - " << maxAllowed << endl
                       << "    It was: " << color[i] 
                       << " on line: " << lineCounter
                       << endl << "    line was: " << intext << endl;
                  return false;
                }
            }
        }
      colorCounter++;

      // Realize if file too large
      if( colorCounter > maxLines )
        {
          cout << "ColorPalette::ReadColorPaletteFile:"
               << "    Palette file too large. Max: " << maxLines << endl;
          return false;
        }


      // Now put the three fields into the main tmp vector
      allColors.push_back( color );

    } // end of while( ifs.getline )

  ifs.close(); // File processed:-)

  // Put tmp vector to main vector
  if( what == "rgb" ) rgbColors = allColors;
  else                hsvColors = allColors;

  // Determine the palette size
  paletteSize = allColors.size();

  // Create the complementary palette
  if( what == "rgb" ) RgbToHsv();
  else                HsvToRgb();

  return true;

} // end of: ReadColorPaletteFile




// WriteColorPaletteFile()
// ---------------------------------------------------------------t---
bool ColorPalette::WriteColorPaletteFile( const char* file, const string& what )
{
  vector<vector<float32> >::iterator  rgbIt, hsvIt;
  // Check parameter 
  if( what != "rgb" && what != "hsv" )
    {
      cout << "ColorPalette::WriteColorPaletteFile:\n"
           << "    Variable 'what' can only be: rgb or hsv" << endl;
      return false;
    }

  // Open file for output
  ofstream ofs( file );
  if( ! ofs )
    {
      cout << "ColorPalette::WriteColorPaletteFile:\n"
           << "     Could not open file: " << file << endl;
      return false;
    }

  // Write header
  ofs << "Palette type: " << what        << endl;
  ofs << "Palette size: " << paletteSize << endl;

  // Write palette
  if( what == "rgb" )
    {
      for( rgbIt  = rgbColors.begin();
           rgbIt != rgbColors.end(); rgbIt++ )
        {
          for( size_t i=0; i<3; i++ )
            {
              ofs << '\t' << (*rgbIt)[i];
            }
          ofs << endl;
        }
    }
  else
    {
      for( hsvIt  = hsvColors.begin();
           hsvIt != hsvColors.end(); hsvIt++ )
        {
          for( size_t i=0; i<3; i++ ) ofs << '\t' << (*hsvIt)[i];
          ofs << endl;
        }
    }
  ofs.close();

  ifstream ifs(file);
  if(ifs) return true;
  else    return false;

} // end of WriteColorPaletteFile




// WriteRgbColorPaletteToStdout()
// ---------------------------------------------------------------t---
void  ColorPalette::WriteRgbColorPaletteToStdout()
{
  vector<vector<float32> >::iterator  rgbIt;
  
  for( rgbIt  = rgbColors.begin();
       rgbIt != rgbColors.end(); rgbIt++ )
    {
      for( size_t i=0; i<3; i++ ) cout << '\t' << ((*rgbIt)[i]);
      cout << endl;
    }
}



// WriteRgbColorPaletteToStdoutWithLineNumbers()
// ---------------------------------------------------------------t---
void  ColorPalette::WriteRgbColorPaletteToStdoutWithLineNumbers()
{
   vector<vector<float32> >::iterator  rgbIt;
  size_t lineCounter = 1;
  for( rgbIt  = rgbColors.begin();
       rgbIt != rgbColors.end(); rgbIt++ )
    {
      cout << lineCounter++;
      for( size_t i=0; i<3; i++ ) cout << '\t' << ((*rgbIt)[i]);
      cout << endl;
    }
}



// WriteHsvColorPaletteToStdout()
// ---------------------------------------------------------------t---
void  ColorPalette::WriteHsvColorPaletteToStdout()
{
  vector<vector<float32> >::iterator  hsvIt;
  
  for( hsvIt  = hsvColors.begin();
       hsvIt != hsvColors.end(); hsvIt++ )
    {
      for( size_t i=0; i<3; i++ ) cout << '\t' << (*hsvIt)[i];
      cout << endl;
    }
}




// WriteHsvColorPaletteToStdoutWithLineNumbers()
// ---------------------------------------------------------------t---
void ColorPalette::WriteHsvColorPaletteToStdoutWithLineNumbers()
{
  vector<vector<float32> >::iterator  hsvIt;
  size_t lineCounter = 1;
  
  for( hsvIt  = hsvColors.begin();
       hsvIt != hsvColors.end(); hsvIt++ )
    {
      cout << lineCounter++;
      for( size_t i=0; i<3; i++ ) cout << '\t' << (*hsvIt)[i];
      cout << endl;
    }
}



// ScaleColorRangeTo( const T& new_min, const T& new_max )
// ------------------------------------------------------------------
// tested: SKM 13/8/98
void ColorPalette::ScaleColorRangeTo( float32 new_min, float32 new_max )
 {
    // find min and max value
    float32                             old_min, old_max;
    vector<vector<float32> >::iterator  it;
    vector<float32>::iterator           jt;
    
    // the min, max values may be different among the r,g, b
    // columns. We are however interested in the absolute 
    // range
    // --------
    old_min = old_max = (*(*rgbColors.begin()).begin());

    // finding min /max
    // ----------------
    for ( it=rgbColors.begin(); it!=rgbColors.end(); it++ )
      for ( jt=(*it).begin(); jt!=(*it).end(); jt++ )
        {
           if ( *jt > old_max ) old_max = *jt;
           if ( *jt < old_min ) old_min = *jt;
        }
    assert( old_max >= old_min );
    assert( old_max > 0.0 );

    // test whether we are already O.K.
    // --------------------------------
    if ( old_max == new_max && old_min == new_min ) return;

    // scaling the color range to the new values
    // -----------------------------------------
    float32  old_range = old_max - old_min;
    float32  new_range = new_max - new_min;
    
    for ( it=rgbColors.begin(); it!=rgbColors.end(); it++ )
      for ( jt=(*it).begin(); jt!=(*it).end(); jt++ )
        {
          *jt = new_min + ((*jt - old_min) / old_range) * new_range;
        }

 } // end ScaleColorRangeTo




// ColorRangeRGB() const
// --------------------------------------------------------------------
// S.K. Matthai 13/8/98
void ColorPalette::ColorRangeRGB() const
 {
    float32  minR, maxR, minG, maxG, minB, maxB;
    vector<vector<float32> >::const_iterator it;

    minR = maxR = (*rgbColors.begin())[0];
    minG = maxG = (*rgbColors.begin())[1];
    minB = maxB = (*rgbColors.begin())[2];

    for ( it=rgbColors.begin(); it!=rgbColors.end(); it++ )
      {
         // minimum
         if ( minR > (*it)[0] ) minR = (*it)[0];
         if ( minG > (*it)[1] ) minG = (*it)[1];
         if ( minB > (*it)[2] ) minB = (*it)[2];
         // maximum
         if ( maxR < (*it)[0] ) maxR = (*it)[0];
         if ( maxG < (*it)[1] ) maxG = (*it)[1];
         if ( maxB < (*it)[2] ) maxB = (*it)[2];
      }
    cout <<"\nColorPalette::ColorRangeRGB:" << endl;
    cout <<"  Range for red:    "<< minR <<"\t"<< maxR << endl;
    cout <<"  Range for green:  "<< minR <<"\t"<< maxR << endl;
    cout <<"  Range for blue:   "<< minR <<"\t"<< maxR << endl;

 } // end ColorRangeRGB

} // end namespace csmp
