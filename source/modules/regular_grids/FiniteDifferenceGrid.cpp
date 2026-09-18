// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FiniteDifferenceGrid.h"
#include "binaryReadWrite.h"
#include "CSMP_highLevelUtilities.h"
#include "ColorPalette.h"

#ifdef CSMP_WITH_IMAGE_OUTPUT
  #if defined(_WIN32)
      // Windows
      #include "jpeglib.h"
  #elif defined(__APPLE__)
      // macOS / XCode
     #if defined(__x86_64__)
      // Intel 64-bit build
      #include "jpeglib.h"
      #elif defined(__aarch64__) || defined(__arm64__)
          // Apple Silicon build
          // #include </opt/homebrew/include/jpeglib.h>
          // TODO: temp fix for Rosetta project
          #include "jpeglib.h"
      #endif
  #elif defined(__linux__)
      // Linux
      #include "jpeglib.h"
  #endif
#endif // CSMP_WITH_IMAGE_OUTPUT

using namespace std;

namespace csmp {


FiniteDifferenceGrid::FiniteDifferenceGrid()
  : xresolution(1.0), yresolution(1.0),
    x_min(0.0), y_min(0.0), x_max(0.0), y_max(0.0), z_min(0.0), z_max(0.0),
    size_x(0), size_y(0), size_z(0), rowlength(0),
    xfr(0), yfr(0), zfr(0)
{
//   cout <<"\nFiniteDifferenceGrid: called default constructor..."<< endl;
}




FiniteDifferenceGrid::FiniteDifferenceGrid( double x_dim,
                                            double y_dim,
                                            double xres, double yres, int32_t frame_width )
  {
    Initialize( x_dim, y_dim, xres, yres, frame_width );
 }






FiniteDifferenceGrid::FiniteDifferenceGrid( double x_dmin, double x_dmax,
                                                double y_dmin, double y_dmax,
                                                double xres, double yres, int32_t frame_width )
 {
    Initialize( x_dmin, x_dmax, y_dmin, y_dmax, xres, yres, frame_width );
 }



double&  FiniteDifferenceGrid::operator()( int32_t row, int32_t col )
  {
#ifndef NDEBUG
      if ( row+yfr < 0 || row > size_y+yfr )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::operator(int32_t,int32): "<< std::endl;
           std::cout <<"Row access violation in row: "<< row;
           std::cout <<" (max row index = "<< size_y+yfr <<")" << std::endl;
           std::cout.flush();
           throw std::range_error("FiniteDifferenceGrid::operator(row access)");
        }
      if ( col+xfr < 0 || col > size_x+xfr )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::operator(int32_t,int32): "<< std::endl;
           std::cout <<"Column access violation in column: "<< col;
           std::cout <<" (max column index = "<< size_x+xfr <<")" << std::endl;
           std::cout.flush();
           throw std::range_error("FiniteDifferenceGrid::operator(column access)");
        }
#endif
      return grid[ (row+yfr)*rowlength + (col+xfr) ];
  }




double  FiniteDifferenceGrid::operator()( int32_t row, int32_t col ) const
  {
#ifndef NDEBUG
      if ( row+yfr < 0 || row > size_y+yfr )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::operator(int32_t,int32): "<< std::endl;
           std::cout <<"Row access violation in row: "<< row;
           std::cout <<" (max row index = "<< size_y+yfr <<")" << std::endl;
           std::cout.flush();
           throw std::range_error("const FiniteDifferenceGrid::operator(row access)");
        }
      if ( col+xfr < 0 || col > size_x+xfr )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::operator(int32_t,int32): "<< std::endl;
           std::cout <<"Column access violation in column: "<< col;
           std::cout <<" (max column index = "<< size_x+xfr <<")" << std::endl;
           std::cout.flush();
           throw std::range_error("const FiniteDifferenceGrid::operator(column access)");
        }
#endif
      return grid[ static_cast<size_t>((row+yfr)*rowlength + (col+xfr)) ];
  }




double  FiniteDifferenceGrid::Value( int32_t row, int32_t col ) const
  {
#ifndef NDEBUG
      if ( row+yfr < 0 || row > size_y+yfr )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::Value(int32_t,int32): "<< std::endl;
           std::cout <<"Row access violation in row: "<< row;
           std::cout <<" (max row index = "<< size_y+yfr <<")" << std::endl;
           std::cout.flush();
           throw std::range_error("FiniteDifferenceGrid::Value(row access)");
        }
      if ( col+xfr < 0 || col > size_x+xfr )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::Value(int32_t,int32): "<< std::endl;
           std::cout <<"Column access violation in column: "<< col;
           std::cout <<" (max column index = "<< size_x+xfr <<")" << std::endl;
           std::cout.flush();
           throw std::range_error("FiniteDifferenceGrid::Value(column access)");
        }
#endif
      return grid[ (row+yfr)*rowlength + (col+xfr) ];
  }



double& FiniteDifferenceGrid::N( int32_t row, int32_t col )
  {
      return (*this)(row-1,col);
  }



double& FiniteDifferenceGrid::NW( int32_t row, int32_t col )
  {
      return (*this)(row-1,col-1);
  }



double& FiniteDifferenceGrid::NE( int32_t row, int32_t col )
  {
      return (*this)(row-1,col+1);
  }



double& FiniteDifferenceGrid::S( int32_t row, int32_t col )
  {
      return (*this)(row+1,col);
  }



double& FiniteDifferenceGrid::SW( int32_t row, int32_t col )
  {
      return (*this)(row+1,col-1);
  }



double& FiniteDifferenceGrid::SE( int32_t row, int32_t col )
  {
      return (*this)(row+1,col+1);
  }



double& FiniteDifferenceGrid::W( int32_t row, int32_t col )
  {
      return (*this)(row,col-1);
  }




double& FiniteDifferenceGrid::E( int32_t row, int32_t col )
  {
      return (*this)(row,col+1);
  }



///returns x-coordinate for column c-style index 0...n
double  FiniteDifferenceGrid::X( int32_t i ) const
 {
#ifndef NDEBUG
    if ( i+xfr<0 || i>size_x+xfr ) {
         std::cout <<"\nFiniteDifferenceGrid<double>::X: Index out of range: "<< i;
         std::cout <<" versus ("<< -xfr <<"-"<< size_x+xfr <<")."<< std::endl;
         throw std::out_of_range("FiniteDifferenceGrid<double>::X");
      }
#endif
    return i * xresolution + x_min;
 }



///returns x-coordinate for row c-style index 0...n
double  FiniteDifferenceGrid::Y( int32_t j ) const
 {
#ifndef NDEBUG
    if ( j+yfr<0 || j>size_y+yfr ) {
         std::cout <<"\nFiniteDifferenceGrid<double>::Y: Index out of range: "<< j;
         std::cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< std::endl;
         throw std::out_of_range("FiniteDifferenceGrid<double>::Y");
      }
#endif
     return j * yresolution + y_min;
 }



/**

Computes the grid indices (i{0,n-1}, j{0..m-1}) of the grid point (x,y)
that is located the closest to the entered point

NB:  i,j are NOT numbered in the usual matrix sense
     - i refers to x-coordinate direction
     - j refers to y-coordinate direction

*/
void  FiniteDifferenceGrid::ClosestGridPointTo( double x, double y,
                                                       int32_t& row, int32_t& col )
 const
 {
#ifndef NDEBUG
      if ( x < x_min || x > x_max )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::ClosestGridPointTo: "<< std::endl;
           std::cout <<"x-coordinate lies outside of grid: "<< x;
           std::cout <<" (xmin="<< x_min <<", xmax="<< x_max <<")" << std::endl;
           std::cout.flush();
           throw std::domain_error("FiniteDifferenceGrid::ClosestGridPointTo(error x)");
        }
      if ( y < y_min || y > y_max )
        {
           std::cout <<"\nFiniteDifferenceGrid<double>::ClosestGridPointTo: "<< std::endl;
           std::cout <<"y-coordinate lies outside grid: "<< y;
           std::cout <<" (ymin="<< y_min <<", ymax="<< y_max <<")" << std::endl;
           std::cout.flush();
           throw std::domain_error("FiniteDifferenceGrid::ClosestGridPointTo(error y)");
        }
#endif
      col = rint((x-x_min) / xresolution);
      row = rint((y-y_min) / yresolution);
 }






bool FiniteDifferenceGrid::In( const char* fname )
 {
    ifstream ifs(fname);
    if ( !ifs )
      {
         cout <<"\nFiniteDifferenceGrid<double>::In: input file cannot be opened. ";
         cout <<" Nothing was done..."<< endl;
         return false;
      }
    char intext[255];
    double     fxmin, fxmax, fymin, fymax, fdx, fdy;
    int32_t  frows, fcols;

    // reading header
    ifs.getline( intext, 255 );
    cout <<"\nTitle of file: '"<< fname <<"': "<< intext << endl;

    // overall dimensions
    ifs >> fdx >> fdy;

     // x & y ranges
    ifs >> fxmin >> fxmax >> fymin >> fymax;

    // resolution x and y
    ifs >> xresolution >> yresolution;

    // rows x columns
    ifs >> frows >> fcols;

    // frame width
    ifs >> xfr >> yfr;
    if ( xfr != yfr ) {
         cout <<"\nFiniteDifferenceGrid<double>::In:  Cannot handle different X & Y frame sizes."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::In");
      }

    // rebuilding grid
    Initialize( fxmin, fxmax, fymin, fymax, xresolution, yresolution, xfr );

    // reading core data of grid from input file (not the frame)
    // extra row and column are used for the last values of the range
    for ( int32_t i=-yfr; i<frows+yfr; i++ )
      for ( int32_t j=-xfr; j<fcols+xfr; j++ ) ifs >> (*this)(i,j);
    ifs.close();

    cout <<"\nFiniteDifferenceGrid<double>::In: grid build successfully from text file." << endl;
    cout.flush();

    return true;
 }





void FiniteDifferenceGrid::Initialize( double x_dim, double y_dim,
				                       double xres, double yres, int32_t frame_width )
 {
     Initialize( 0.0, x_dim, 0.0, y_dim, xres, yres, frame_width );
 }




void FiniteDifferenceGrid::Initialize( double x_dmin, double x_dmax,
                                       double y_dmin, double y_dmax,
                                       double xres, double yres, int32_t frame_width )
 {
     // 0. removing old grid
     x_min       = x_dmin;
     x_max       = x_dmax;
     y_min       = y_dmin;
     y_max       = y_dmax;
     xresolution = xres;
     yresolution = yres;
     xfr = yfr = zfr = frame_width;

     // 0. checks
     z_min = z_max = 0.0;
     double swap;
     if ( x_min > x_max ) { swap=x_min; x_min=x_max; x_max=swap; }
     if ( y_min > y_max ) { swap=y_min; y_min=y_max; y_max=swap; }

     // 1. getting the resolution
     size_x = static_cast<int32_t>((x_max-x_min) / xresolution + 1.);
     size_y = static_cast<int32_t>((y_max-y_min) / yresolution + 1.);
     size_z = 0;

     // 2. extra data member for speed
     rowlength = size_x + 2 * xfr;

     // 3. setting up the vector-grid
     cout <<"\nFiniteDifferenceGrid: building new grid; allocating ";
     cout << (((size_x+2*xfr) * (size_y+2*yfr) * sizeof(double)) / 1.0e+6);
     cout <<" MByte of memory..." << endl;

     // 4. give initial value of zero
     grid.assign( (size_x+2*xfr) * (size_y+2*yfr), 0.0 );
 }



    double FiniteDifferenceGrid::ResolutionX() const { return xresolution; }


    double FiniteDifferenceGrid::ResolutionY() const { return yresolution; }


    int32_t   FiniteDifferenceGrid::Rows()        const { return size_y-1; }


    int32_t   FiniteDifferenceGrid::Columns()     const { return size_x-1; }


    bool   FiniteDifferenceGrid::HasFrame()    const { return (xfr!=0||yfr!=0||zfr!=0); }


    double FiniteDifferenceGrid::MinX()        const { return x_min; }


    double FiniteDifferenceGrid::MinY()        const { return y_min; }


    double FiniteDifferenceGrid::MinZ()        const { return z_min; }


    double FiniteDifferenceGrid::MaxX()        const { return x_max; }


    double FiniteDifferenceGrid::MaxY()        const { return y_max; }


    double FiniteDifferenceGrid::MaxZ()        const { return z_max; }


bool   FiniteDifferenceGrid::IsInitialized() const
  {
     return !(size_x==0 || size_y==0);
  }




FiniteDifferenceGrid& FiniteDifferenceGrid::operator=( const FiniteDifferenceGrid& g )
 {
    if ( &g == this ) return *this;
    xresolution = g.xresolution;
    yresolution = g.yresolution;
    size_x      = g.size_x;
    size_y      = g.size_y;
    size_z      = g.size_z;
    rowlength   = g.rowlength;
    x_min       = g.x_min;
    y_min       = g.y_min;
    z_min       = g.z_min;
    x_max       = g.x_max;
    y_max       = g.y_max;
    z_max       = g.z_max;
    xfr         = g.xfr;
    yfr         = g.yfr;
    zfr         = g.zfr;
    grid        = g.grid;

    return *this;
 }




FiniteDifferenceGrid& FiniteDifferenceGrid::operator=( double val )
 {
    std::fill( grid.begin(), grid.end(), val );
    return *this;
 }


void   FiniteDifferenceGrid::SetRowTo( int32_t row, double val )
 {
    if ( row+xfr<0 || row>size_x+xfr ) {
         cout <<"\nFiniteDifferenceGrid<double>::SetRowTo: Row index out of range: "<< row;
         cout <<" versus ("<< -xfr <<"-"<< size_x+xfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::SetRowTo");
      }
    for ( int32_t j=-xfr; j<(size_x+xfr); j++ ) (*this)(row,j) = val;
 }



void   FiniteDifferenceGrid::SetColumnTo( int32_t col, double val )
 {
    if ( col+yfr<0 || col>size_y+yfr ) {
         cout <<"\nFiniteDifferenceGrid<double>::SetColumnTo: column index out of range: "<< col;
         cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::SetColumnTo");
      }
    for ( int32_t i=-yfr; i<(size_y+yfr); i++ ) (*this)(i,col) = val;
 }




void   FiniteDifferenceGrid::SetFirstNRowsTo( int32_t n, double val )
 {
    if ( n<0 || n>size_y ) {
         cout <<"\nFiniteDifferenceGrid<double>::SetFirstNRowsTo: row index out of range: "<< n;
         cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::SetFirstNRowsTo");
      }
    for ( int32_t i=0; i<n; i++ ) SetRowTo( i, val );
 }




void   FiniteDifferenceGrid::SetLastNRowsTo( int32_t n, double val )
 {
    if ( n<0 || n>size_y ) {
         cout <<"\nFiniteDifferenceGrid<double>::SetLastNRowsTo: row index out of range: "<< n;
         cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::SetLastNRowsTo");
      }
    for ( int32_t i=size_y-1; i>(size_y-1-n); i-- ) SetRowTo( i, val );
 }



double FiniteDifferenceGrid::RowAverage( int32_t row ) const
 {
    double avg = 0.0;
    assert( row+yfr >= 0 && row < size_y+yfr );
    for ( int32_t j=-xfr; j<(size_x+xfr); j++ ) avg += Value(row,j);
    return avg /= size_x+2*xfr;
 }




double FiniteDifferenceGrid::ColumnAverage( int32_t col ) const
 {
    double avg = 0.0;
    assert( col+xfr >= 0 && col < size_x+xfr );
    for ( int32_t i=-yfr; i<(size_y+yfr); i++ ) avg += Value(i,col);
    return avg /= size_y+2*yfr;
 }



void   FiniteDifferenceGrid::SetRegionTo( double xmin, double xmax, double ymin, double ymax, double val )
 {
     int32_t i, j, imin, imax, jmin, jmax;

     assert( xmin < xmax );
     assert( ymin < ymax );
     assert( xmax <= x_max );
     assert( ymax <= y_max );

     ClosestGridPointTo( xmin, ymin, imin, jmin );
     ClosestGridPointTo( xmax, ymax, imax, jmax );

     for ( i=imin; i<=imax; i++ )
       for ( j=jmin; j<=jmax; j++ ) (*this)(i,j) = val;
 }




void   FiniteDifferenceGrid::Out() const
 {
    int32_t i, j;
    cout <<"\n\nFiniteDifferenceGrid<double>::Out(): printing grid of size: ";
    cout << x_max-x_min <<" by "<< y_max-y_min <<" m" << endl;
    cout <<"x-range:    "<< x_min <<"-"<< x_max <<" m" << endl;
    cout <<"y-range:    "<< y_min <<"-"<< y_max <<" m" << endl;
    cout <<"resolution (X): "<< xresolution <<" m"        << endl;
    cout <<"resolution (Y): "<< yresolution <<" m"        << endl;
    if ( HasFrame() ) cout <<"grid with frame of width: "<< xfr << endl;
    cout << size_y <<" rows, "<< size_x <<" columns" << endl;
    cout  << "\n\t";
    for ( j=-xfr; j<size_x+xfr; j++ ) cout << "col "<< j <<"\t";

    for ( i=-yfr; i<size_y+yfr; i++ )
      {
         cout <<"\nrow "<< i <<"\t";
         for ( j=-xfr; j<size_x+xfr; j++ )
           cout << (*this)(i,j) <<"\t";
      }
    cout <<"\n\nDone..." << endl;
 }




void   FiniteDifferenceGrid::Out( const char* fname, int32_t tstep, bool with_frame ) const
 {
    char name[200], num[30];
    strcpy( name, fname );
    snprintf( num, sizeof(num), "%d", tstep );
    strcat( name, num );

    ofstream ifs( name );
    assert( ifs.is_open() );

    ifs <<"FiniteDifferenceGrid<double>::Out: ASCII text file from double grid." << endl;
    // overall dimensions
    ifs << x_max-x_min <<"\t"<< y_max-y_min << endl;
     // x & y ranges
    ifs << x_min <<"\t"<< x_max <<"\t"<< y_min <<"\t"<< y_max << endl;
    // resolution x and y
    ifs << xresolution <<"\t"<< yresolution << endl;
    // rows x columns
    ifs << size_y <<"\t"<< size_x << endl;
    // frame width
    if ( with_frame ) ifs << xfr <<" "<< yfr << endl;
    else ifs << 0 <<" "<< 0 << endl;

    // central portion of grid to file
    if ( !with_frame )
      for ( int32_t i=0; i<size_y; i++ ) {
           for ( int32_t j=0; j<size_x; j++ ) ifs << (*this)(i,j) <<"\t";
           ifs << endl;
        }
    else
      for ( int32_t i=-yfr; i<size_y+yfr; i++ ) {
           for ( int32_t j=-xfr; j<size_x+xfr; j++ ) ifs << (*this)(i,j) <<"\t";
           ifs << endl;
        }

    ifs.close();
    cout <<"\n\n'" << name <<"' written successfully..." << endl;
 }



bool   FiniteDifferenceGrid::BinaryOut( const char* bin_name, int32_t tstep, bool with_frame ) const
 {
    char name[200], num[30], heading[200];
    strcpy( name, bin_name );
    snprintf( num, sizeof(num), "%d", tstep );
    strcat( name, num );
    strcpy( heading, "FiniteDifferenceGrid::BinaryOut: double grid as binary file");

    fstream fp (name, ios::out | ios::binary);
    if ( !fp.is_open() ) {
        cout <<"\nFiniteDifferenceGrid::BinaryOut: File: "<< bin_name << " could not be opened"<< endl;
        return false;
      }
    binaryFileWrite( fp, heading );

    // stores dimensions of grid
    std::vector<double>    dim_fT(8);
    // stores rows, columns, and frame
    std::vector<int32_t> dim_int(4);
    // grid data
    std::vector<double>    grid_data((static_cast<uint32_t>(size_x)+2*static_cast<uint32_t>(xfr)) * (static_cast<uint32_t>(size_y)+2*static_cast<uint32_t>(yfr)));
    size_t counter(0);

    // overall dimensions
    dim_fT[0] = x_max-x_min;
    dim_fT[1] = y_max-y_min;
     // x & y ranges
    dim_fT[2] = x_min;
    dim_fT[3] = x_max;
    dim_fT[4] = y_min;
    dim_fT[5] = y_max;
    // resolution x and y
    dim_fT[6] = xresolution;
    dim_fT[7] = yresolution;
    // rows x columns
    dim_int[0] = size_x;
    dim_int[1] = size_y;
    // frame width
    if ( with_frame ) {
        dim_int[2] = xfr;
        dim_int[3] = yfr;
      }
    else {
        dim_int[2] = 0;
        dim_int[3] = 0;
      }

    binaryFileWrite( fp, dim_fT );
    binaryFileWrite( fp, dim_int );

    // writing grid data to storage vector
    if ( !with_frame )
      for ( int32_t i=0; i<size_y; i++ ) {
           for ( int32_t j=0; j<size_x; j++ ) {
               grid_data[counter] = (*this)(i,j);
               counter++;
             }
        }
    else
      for ( int32_t i=-yfr; i<size_y+yfr; i++ ) {
           for ( int32_t j=-xfr; j<size_x+xfr; j++ ) {
               grid_data[counter] = (*this)(i,j);
               counter++;
             }
        }

   binaryFileWrite( fp, grid_data );
   fp.close();

   cout <<"\n\n'" << bin_name <<"' written successfully..." << endl;

   return true;
 }



bool FiniteDifferenceGrid::BinaryIn( const char* bin_name )
 {
    fstream fp(bin_name, ios::in | ios::binary);
    if ( !fp.is_open() ) {
        cout <<"\nFiniteDifferenceGrid::BinaryIn: File: "<< bin_name;
        cout <<" could not be opened"<< endl;
        return false;
     }
    char heading[200];
    binaryFileRead( fp, heading );
    cout <<"\nFiniteDifferenceGrid::BinaryIn: Reading: "<< heading << endl;

    // read dimensions of grid
    std::vector<double>  dim_fT;
    // read rows, columns, and frame
    std::vector<int32_t> dim_int;
    // read data
    std::vector<double>  grid_data;
    size_t counter(0);

    // read info about dimension
    binaryFileRead( fp, dim_fT );
    binaryFileRead( fp, dim_int );

    if ( dim_int[2] != dim_int[3] ) {
         cout <<"\nFiniteDifferenceGrid::BinaryIn:  Cannot handle different X & Y frame sizes."<< endl;
         throw length_error("FiniteDifferenceGrid::BinaryIn");
      }

    // rebuild grid
    Initialize( dim_fT[2], dim_fT[3], dim_fT[4], dim_fT[5], dim_fT[6], dim_fT[7], xfr );

    // read data and transfer to grid
    binaryFileRead( fp, grid_data );
    for ( int32_t i=-dim_int[3]; i<dim_int[1]+dim_int[3]; i++ ) {
        for ( int32_t j=-dim_int[2]; j<dim_int[0]+dim_int[2]; j++ ) {
            (*this)(i,j) = grid_data[counter];
            counter++;
          }
      }

 	  fp.close();

    cout <<"\nFiniteDifferenceGrid<double>::BinaryIn: grid build successfully from binary file." << endl;
    cout.flush();

    return true;
 }



void  FiniteDifferenceGrid::DataMinMaxWithoutNAN( double& dmin, double& dmax,
					                                        bool frame_included ) const
 {
    int32_t  i, j;
    double val;
    bool first_call(true);

    if ( frame_included )
      {
	// getting min and max of the data (including frame)
        for ( dmax=dmin=grid[0], i=1; i<((size_x+2*xfr) * (size_y+2*yfr)); i++ )
           {
             if ( !isnan( grid[i] ) ) {
                 if ( first_call ) {
                     dmin = dmax = grid[i];
                     first_call = false;
                   }
                 if ( grid[i] > dmax ) dmax = grid[i];
                 if ( grid[i] < dmin ) dmin = grid[i];
               }
           }
      }
    else
      {
	// getting min and max of the data (excluding frame)
        dmax=dmin=(*this)( static_cast<int32_t>(0), static_cast<int32_t>(0) );
        for ( i=0; i<size_y; i++ )
          for ( j=0; j<size_x; j++ )
           {
               if ( !isnan( (*this)(i,j) ) ) {
                   if ( first_call ) {
                       dmin = dmax = (*this)(i,j);
                       first_call = false;
                     }
                   val = (*this)(i,j);
                   if ( val > dmax ) dmax = val;
                   if ( val < dmin ) dmin = val;
                 }
           }
      }
 } // end



void  FiniteDifferenceGrid::DataMinMax( double& dmin, double& dmax,
					                                  bool frame_included ) const
 {
    int32_t  i, j;
    double val;

    if ( frame_included )
      {
	// getting min and max of the data (including frame)
        for ( dmax=dmin=grid[0], i=1; i<((size_x+2*xfr) * (size_y+2*yfr)); i++ )
           {
             if ( grid[i] > dmax ) dmax = grid[i];
             if ( grid[i] < dmin ) dmin = grid[i];
           }
      }
    else
      {
	// getting min and max of the data (excluding frame)
        dmax=dmin=(*this)( static_cast<int32_t>(0),static_cast<int32_t>(0) );
        for ( i=0; i<size_y; i++ )
          for ( j=0; j<size_x; j++ )
           {
             val = (*this)(i,j);
             if ( val > dmax ) dmax = val;
             if ( val < dmin ) dmin = val;
           }
      }
 } // end



///scales all grid data to range: dmin to dmax
void   FiniteDifferenceGrid::ScaleDataToRange( double dmin, double dmax )
 {
    double  old_min, old_max;
    DataMinMax( old_min, old_max, false );
    double  old_range = old_max - old_min;
    double  new_range = dmax - dmin;

    for ( int32_t i=0; i<((size_x+2*xfr) * (size_y+2*yfr)); i++ )
      grid[i] = dmin + ((grid[i] - old_min)/old_range) * new_range;

 } // end ScaleDataToRange











/// 4 point, third order polynomial extrapolation
void  FiniteDifferenceGrid::PolynomialInterpolation( double* xa, double* ya,
                                                     double  x,  double& y,
                                                     double& dy ) const
{
    // constraint point arrays are: xa[], double ya[],
    double  c[5], d[5];
    int32_t     i, m, ns=1;
    double  den, dif, dift, ho, hp, w;

	dif=fabs(x-xa[1]);
	for (i=1;i<=4;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<4;m++) {
		for (i=1;i<=4-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
			  {
			     cout <<"\nFiniteDifferenceGrid<double>::PolynomialInterpolation: ";
			     cout <<"Error. No interpolation was accomplished."<< endl;
			     cout <<"Method input: x=" << x <<", modified y="<< y <<", error="<< dy << endl;
			     cout <<"Constraint point arrays: "<< endl <<"position data: ";
			     for ( int32_t j=1; j<=4; j++ ) cout << xa[j] <<" ";
			     cout << endl <<"dependent variable data: ";
			     for ( int32_t j=1; j<=4; j++ ) cout << ya[j] <<" ";
			     cout << endl;
			     // setting output to zero
			     y  = 0.0;
			     dy = 0.0;
			     return;
			  }
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*ns < (4-m) ? c[ns+1] : d[ns--]));
	}

} // end PolynomialInterpolation


/**

given the point x,y outside of the grid use linear interpolation of the
values at the grid boundary to find the property value outside the grid.
The point which is extrapolated must be in the first square outside the grid.

A third order polynomial extrapolation is used which performs well if the grid
variable varies smoothly and monotonous.

*/
double  FiniteDifferenceGrid::ExtrapolateTo( double x, double y ) const
 {
    static double  xa[5], ya[5];
    int32_t        i;
    double         a = std::numeric_limits<double>::quiet_NaN();
   double          err{0.};

    // 1. outside normal boundaries
    // ----------------------------
    if ( x >= x_min && x <= x_max )
      {
         // BOTTOM OUT
         if ( y > y_max )
           {
              // getting 4 interpolation points for 3rd order polynomial extrapolation
              // outside of grid
              for ( i=1, a=y_max-3*yresolution; a<=y_max; a+=yresolution, i++ )
                {
                   xa[i] = a;
                   ya[i] = (*this)(x,a);
                }
           }
         else
         // TOP OUT
         if ( y < y_min )
           {
              for ( i=1, a=y_min; a<=y_min+3*yresolution; a+=yresolution, i++ )
                {
                   xa[i] = a;
                   ya[i] = (*this)(x,a);
                }
           }
         // doing a 3rd order interpolation, result returned into a
         PolynomialInterpolation( xa, ya, y, a, err );
         return a;
      }
    if ( y >= y_min && y <= y_max )
      {
         // LEFT OUT
         if ( x < x_min )
           {
              for ( i=1, a=x_min; a<=x_min+3*xresolution; a+=xresolution, i++ )
                {
                   xa[i] = a;
                   ya[i] = (*this)(a,y);
                }
           }
         else
         // RIGHT OUT
         if ( x > x_max )
           {
              for ( i=1, a=x_max-3*xresolution; a<=x_max; a+=xresolution, i++ )
                {
                   xa[i] = a;
                   ya[i] = (*this)(a,y);
                }
           }
         PolynomialInterpolation( xa, ya, x, a, err );

         return a;
      }

    // 2. corner cases
    // ---------------
    double m, len_fac;

    // BELOW ORIGIN
    if ( x < x_min && y < y_min )
      {
         // slope of line from interpolation point through origin to find on-grid
         // interpolation coordinates f(x) = mx + x0
         // NOTE: If the slope is very steep again a point may get interpolated outside
         //       of the grid in order to find the appropriate value.
         m       = (y-y_min) / (x-x_min);
         len_fac = sqrt( 1.0 + m*m );

         for ( i=1, a=x_min; a<=x_min+3*xresolution; a+=xresolution, i++ )
           {
              // scaling x to reflect pathlength aint32_t sloping line
              xa[i] = (i-1)*xresolution*len_fac;
              // interpolating on grid f(x)
              ya[i] = (*this)(a, (m * (a-x_min) + y_min) );
           }
         // transforming desired coordinate for extrapolation
         a = len_fac * (x-x_min);
      }
    // LOWER-LEFT CORNER
    if ( x < x_min && y > y_max )
      {
         m       = (y-y_max) / (x-x_min);
         len_fac = sqrt( 1.0 + m*m );

         for ( i=1, a=x_min; a<=x_min+3*xresolution; a+=xresolution, i++ )
           {
              xa[i] = (i-1)*xresolution*len_fac;
              ya[i] = (*this)(a, (m * (a-x_min) + y_max) );
           }
         a = len_fac * (x-x_min);
      }
    // LOWER RIGHT CORNER
    if ( x > x_max && y > y_max )
      {
         m       = (y-y_max) / (x-x_max);
         len_fac = sqrt( 1.0 + m*m );

         for ( i=1, a=x_max-3*xresolution; a<=x_max; a+=xresolution, i++ )
           {
              xa[i] = (4-i)*xresolution*len_fac;
              ya[i] = (*this)(a, (m * (a-x_max) + y_max) );
           }
         a = len_fac * (x_max-x);
      }
    // TOP-RIGHT CORNER
    if ( x > x_max && y < y_min )
      {
         m       = (y-y_min) / (x-x_max);
         len_fac = sqrt( 1.0 + m*m );

         for ( i=1, a=x_max-3*xresolution; a<=x_max; a+=xresolution, i++ )
           {
              xa[i] = (4-i)*xresolution*len_fac;
              ya[i] = (*this)(a, (m * (a-x_max) + y_min) );
           }
         a = len_fac * (x_max-x);
      }
    // extrapolating beyond grid corners
    PolynomialInterpolation( xa, ya, a, y, err );

//    cout <<"\nAccumulated interpolation function constraint values: "<< endl;
//    for ( i=1; i<=4; i++ )
//       cout <<"x: "<< xa[i] <<", f(x): "<< ya[i] << endl;
//    cout << endl;

    return y;

 } // end ExtrapolateTo


#ifdef CSMP_WITH_IMAGE_OUTPUT

// =============================================================================
// RAII helpers — local to this translation unit
// =============================================================================

namespace {

/// RAII wrapper for jpeg_compress_struct — ensures jpeg_destroy_compress
/// is always called on any exit path, including exceptions.
struct JpegCompressGuard
{
    jpeg_compress_struct* p;
    ~JpegCompressGuard() { jpeg_destroy_compress( p ); }
};

/// Lambda deleter for std::unique_ptr<FILE> — avoids the
/// -Wignored-attributes warning produced by decltype(&fclose).
auto make_file_deleter()
{
    return []( FILE* f ) { if ( f ) fclose( f ); };
}

using FilePtr = std::unique_ptr<FILE, decltype(make_file_deleter())>;

/// Opens a file for binary writing and returns a RAII-managed pointer.
/// Returns an empty unique_ptr if the file cannot be opened.
FilePtr open_jpeg_file( const std::string& name )
{
    return FilePtr( fopen( name.c_str(), "wb" ), make_file_deleter() );
}

/// Builds the output filename from base name and zero-padded timestep.
std::string make_jpeg_filename( const char* filename, int32_t timestep )
{
    char num[30];
    snprintf( num, sizeof(num), "%d", timestep );
    std::string padded( num );
    replaceWhiteSpaceBy( padded, '0' );
    return std::string( filename ) + padded + ".jpg";
}

/// Writes a uniform-value text file instead of a JPEG when min == max.
/// Returns true if the text file was written.
bool write_uniform_value_file( const std::string& jpeg_name,
                               double value, int32_t timestep )
{
    std::string txt_name = jpeg_name + ".txt";
    std::ofstream ofs( txt_name, std::ios::out | std::ios::trunc );
    if ( !ofs ) {
        std::cout << "\nOutput file: " << txt_name
                  << " could not be opened. Nothing was done..." << std::endl;
        return false;
    }
    ofs << txt_name << ": Uniform variable value: "
        << value << " at timestep: " << timestep << std::endl;
    std::cout << "\n\n'" << txt_name << "' written instead of JPG file..." << std::endl;
    return true;
}

/// Configures a jpeg_compress_struct for the given image dimensions,
/// component count and colour space.
void configure_jpeg( jpeg_compress_struct& cinfo,
                     int32_t width, int32_t height,
                     int32_t components, bool greyscale )
{
    cinfo.image_width      = static_cast<JDIMENSION>( width  );
    cinfo.image_height     = static_cast<JDIMENSION>( height );
    cinfo.input_components = components;
    cinfo.in_color_space   = greyscale ? JCS_GRAYSCALE : JCS_RGB;
    jpeg_set_defaults( &cinfo );
}

/// Writes the image buffer to the open JPEG compressor scanline by scanline.
void write_jpeg_scanlines( jpeg_compress_struct& cinfo,
                           const std::vector<unsigned char>& image_buffer,
                           int32_t row_stride )
{
    jpeg_start_compress( &cinfo, TRUE );
    while ( cinfo.next_scanline < cinfo.image_height )
    {
        JSAMPROW row = const_cast<JSAMPROW>(
            &image_buffer[ cinfo.next_scanline
                           * static_cast<JDIMENSION>( row_stride ) ] );
        jpeg_write_scanlines( &cinfo, &row, 1 );
    }
    jpeg_finish_compress( &cinfo );
}

} // anonymous namespace


// =============================================================================
// SaveToJPGWithoutNAN — auto data range
// =============================================================================

void FiniteDifferenceGrid::SaveToJPGWithoutNAN( const char* filename,
                                                int32_t     timestep,
                                                bool        greyscale,
                                                bool        sqrt_of ) const
{
    ColorPalette rgb_colorizer;
    float        colors[4];
    if ( greyscale ) rgb_colorizer.MakeGreyPalette();

    const std::string name = make_jpeg_filename( filename, timestep );

    double old_min, old_max;
    DataMinMaxWithoutNAN( old_min, old_max, false );
    std::cout.setf( std::ios::scientific );
    std::cout << "\nFiniteDifferenceGrid::SaveToJPGWithoutNAN: '" << filename
              << "' data range: " << old_min << " to " << old_max << std::endl;
    std::cout.unsetf( std::ios::scientific );

    if ( old_min == old_max ) {
        write_uniform_value_file( name, old_max, timestep );
        return;
    }

    if ( sqrt_of ) { old_min = std::sqrt( old_min ); old_max = std::sqrt( old_max ); }
    const double old_range = old_max - old_min;
    const double new_range = 255.0;

    const int32_t image_components = greyscale ? 1 : 3;
    const int32_t row_stride       = size_x * image_components;

    std::vector<unsigned char> image_buffer(
        static_cast<size_t>( size_x * size_y * image_components ) );

    size_t incr = 0;
    for ( int32_t i = size_y - 1; i >= 0; --i )
        for ( int32_t j = 0; j < size_x; ++j )
        {
            const double val = (*this)( i, j );
            double out_val;
            if ( std::isnan( val ) )
                out_val = 255.0;
            else if ( sqrt_of )
                out_val = ( ( std::sqrt( val ) - old_min ) / old_range ) * new_range;
            else
                out_val = ( ( val - old_min ) / old_range ) * new_range;

            if ( greyscale ) {
                image_buffer[incr++] = static_cast<unsigned char>( out_val );
            } else {
                if ( std::isnan( val ) )
                    colors[0] = colors[1] = colors[2] = 255.0f;
                else
                    rgb_colorizer.GiveRgb( static_cast<float>( out_val ), colors );
                image_buffer[incr++] = static_cast<unsigned char>( colors[0] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[1] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[2] );
            }
        }

    auto outfile = open_jpeg_file( name );
    if ( !outfile ) {
        std::cerr << "\nFiniteDifferenceGrid::SaveToJPGWithoutNAN: cannot open: "
                  << name << std::endl;
        return;
    }

    jpeg_compress_struct cinfo;
    jpeg_error_mgr       jerr;
    JpegCompressGuard    guard{ &cinfo };

    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    jpeg_stdio_dest( &cinfo, outfile.get() );
    configure_jpeg( cinfo, size_x, size_y, image_components, greyscale );
    write_jpeg_scanlines( cinfo, image_buffer, row_stride );

    std::cout << "\n'" << name << "' written successfully." << std::endl;
}


// =============================================================================
// SaveToJPGWithoutNAN — fixed data range
// =============================================================================

void FiniteDifferenceGrid::SaveToJPGWithoutNAN( const char* filename,
                                                int32_t     timestep,
                                                double      vmin,
                                                double      vmax,
                                                bool        greyscale,
                                                bool        sqrt_of ) const
{
    ColorPalette rgb_colorizer;
    float        colors[4];
    if ( greyscale ) rgb_colorizer.MakeGreyPalette();

    const std::string name = make_jpeg_filename( filename, timestep );

    double old_min, old_max;
    DataMinMaxWithoutNAN( old_min, old_max, false );
    std::cout.setf( std::ios::scientific );
    std::cout << "\nFiniteDifferenceGrid::SaveToJPGWithoutNAN: '" << filename
              << "' data range: " << old_min << " to " << old_max << std::endl;
    std::cout.unsetf( std::ios::scientific );

    // Clamp to prescribed range, warn if actual data exceeds it
    if ( vmin > old_min )
        std::cout << "\nData minimum smaller than prescribed minimum — using actual value.\n";
    else
        old_min = vmin;

    if ( vmax < old_max )
        std::cout << "\nData maximum larger than prescribed maximum — using actual value.\n";
    else
        old_max = vmax;

    if ( old_min == old_max ) {
        write_uniform_value_file( name, old_max, timestep );
        return;
    }

    if ( sqrt_of ) { old_min = std::sqrt( old_min ); old_max = std::sqrt( old_max ); }
    const double old_range = old_max - old_min;
    const double new_range = 255.0;

    const int32_t image_components = greyscale ? 1 : 3;
    const int32_t row_stride       = size_x * image_components;

    std::vector<unsigned char> image_buffer(
        static_cast<size_t>( size_x * size_y * image_components ) );

    size_t incr = 0;
    for ( int32_t i = size_y - 1; i >= 0; --i )
        for ( int32_t j = 0; j < size_x; ++j )
        {
            const double val = (*this)( i, j );
            double out_val;
            if ( std::isnan( val ) )
                out_val = 255.0;
            else if ( sqrt_of )
                out_val = ( ( std::sqrt( val ) - old_min ) / old_range ) * new_range;
            else
                out_val = ( ( val - old_min ) / old_range ) * new_range;

            if ( greyscale ) {
                image_buffer[incr++] = static_cast<unsigned char>( out_val );
            } else {
                if ( std::isnan( val ) )
                    colors[0] = colors[1] = colors[2] = 255.0f;
                else
                    rgb_colorizer.GiveRgb( static_cast<float>( out_val ), colors );
                image_buffer[incr++] = static_cast<unsigned char>( colors[0] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[1] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[2] );
            }
        }

    auto outfile = open_jpeg_file( name );
    if ( !outfile ) {
        std::cerr << "\nFiniteDifferenceGrid::SaveToJPGWithoutNAN: cannot open: "
                  << name << std::endl;
        return;
    }

    jpeg_compress_struct cinfo;
    jpeg_error_mgr       jerr;
    JpegCompressGuard    guard{ &cinfo };

    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    jpeg_stdio_dest( &cinfo, outfile.get() );
    configure_jpeg( cinfo, size_x, size_y, image_components, greyscale );
    write_jpeg_scanlines( cinfo, image_buffer, row_stride );

    std::cout << "\n'" << name << "' written successfully." << std::endl;
}


// =============================================================================
// SaveToJPG — auto data range
// =============================================================================

void FiniteDifferenceGrid::SaveToJPG( const char* filename,
                                      int32_t     timestep,
                                      bool        greyscale,
                                      bool        sqrt_of ) const
{
    ColorPalette rgb_colorizer;
    float        colors[4];
    if ( greyscale ) rgb_colorizer.MakeGreyPalette();

    const std::string name = make_jpeg_filename( filename, timestep );

    double old_min, old_max;
    DataMinMax( old_min, old_max, false );
    std::cout.setf( std::ios::scientific );
    std::cout << "\nFiniteDifferenceGrid::SaveToJPG: '" << filename
              << "' data range: " << old_min << " to " << old_max << std::endl;
    std::cout.unsetf( std::ios::scientific );

    if ( old_min == old_max ) {
        write_uniform_value_file( name, old_max, timestep );
        return;
    }

    if ( sqrt_of ) { old_min = std::sqrt( old_min ); old_max = std::sqrt( old_max ); }
    const double old_range = old_max - old_min;
    const double new_range = 255.0;

    const int32_t image_components = greyscale ? 1 : 3;
    const int32_t row_stride       = size_x * image_components;

    std::vector<unsigned char> image_buffer(
        static_cast<size_t>( size_x * size_y * image_components ) );

    size_t incr = 0;
    for ( int32_t i = size_y - 1; i >= 0; --i )
        for ( int32_t j = 0; j < size_x; ++j )
        {
            const double val     = (*this)( i, j );
            const double out_val = sqrt_of
                ? ( ( std::sqrt( val ) - old_min ) / old_range ) * new_range
                : ( ( val              - old_min ) / old_range ) * new_range;

            if ( greyscale ) {
                image_buffer[incr++] = static_cast<unsigned char>( out_val );
            } else {
                rgb_colorizer.GiveRgb( static_cast<float>( out_val ), colors );
                image_buffer[incr++] = static_cast<unsigned char>( colors[0] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[1] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[2] );
            }
        }

    auto outfile = open_jpeg_file( name );
    if ( !outfile ) {
        std::cerr << "\nFiniteDifferenceGrid::SaveToJPG: cannot open: "
                  << name << std::endl;
        return;
    }

    jpeg_compress_struct cinfo;
    jpeg_error_mgr       jerr;
    JpegCompressGuard    guard{ &cinfo };

    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    jpeg_stdio_dest( &cinfo, outfile.get() );
    configure_jpeg( cinfo, size_x, size_y, image_components, greyscale );
    write_jpeg_scanlines( cinfo, image_buffer, row_stride );

    std::cout << "\n'" << name << "' written successfully." << std::endl;
}


// =============================================================================
// SaveToJPG — fixed data range
// =============================================================================

void FiniteDifferenceGrid::SaveToJPG( const char* filename,
                                      int32_t     timestep,
                                      double      vmin,
                                      double      vmax,
                                      bool        greyscale,
                                      bool        sqrt_of ) const
{
    ColorPalette rgb_colorizer;
    float        colors[4];
    if ( greyscale ) rgb_colorizer.MakeGreyPalette();

    const std::string name = make_jpeg_filename( filename, timestep );

    double old_min, old_max;
    DataMinMax( old_min, old_max, false );
    std::cout.setf( std::ios::scientific );
    std::cout << "\nFiniteDifferenceGrid::SaveToJPG: '" << filename
              << "' data range: " << old_min << " to " << old_max << std::endl;
    std::cout.unsetf( std::ios::scientific );

    // Clamp to prescribed range, warn if actual data exceeds it
    if ( vmin > old_min )
        std::cout << "\nData minimum smaller than prescribed minimum — using actual value.\n";
    else
        old_min = vmin;

    if ( vmax < old_max )
        std::cout << "\nData maximum larger than prescribed maximum — using actual value.\n";
    else
        old_max = vmax;

    if ( old_min == old_max ) {
        write_uniform_value_file( name, old_max, timestep );
        return;
    }

    if ( sqrt_of ) { old_min = std::sqrt( old_min ); old_max = std::sqrt( old_max ); }
    const double old_range = old_max - old_min;
    const double new_range = 255.0;

    const int32_t image_components = greyscale ? 1 : 3;
    const int32_t row_stride       = size_x * image_components;

    std::vector<unsigned char> image_buffer(
        static_cast<size_t>( size_x * size_y * image_components ) );

    size_t incr = 0;
    for ( int32_t i = size_y - 1; i >= 0; --i )
        for ( int32_t j = 0; j < size_x; ++j )
        {
            const double val     = (*this)( i, j );
            const double out_val = sqrt_of
                ? ( ( std::sqrt( val ) - old_min ) / old_range ) * new_range
                : ( ( val              - old_min ) / old_range ) * new_range;

            if ( greyscale ) {
                image_buffer[incr++] = static_cast<unsigned char>( out_val );
            } else {
                rgb_colorizer.GiveRgb( static_cast<float>( out_val ), colors );
                image_buffer[incr++] = static_cast<unsigned char>( colors[0] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[1] );
                image_buffer[incr++] = static_cast<unsigned char>( colors[2] );
            }
        }

    auto outfile = open_jpeg_file( name );
    if ( !outfile ) {
        std::cerr << "\nFiniteDifferenceGrid::SaveToJPG: cannot open: "
                  << name << std::endl;
        return;
    }

    jpeg_compress_struct cinfo;
    jpeg_error_mgr       jerr;
    JpegCompressGuard    guard{ &cinfo };

    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    jpeg_stdio_dest( &cinfo, outfile.get() );
    configure_jpeg( cinfo, size_x, size_y, image_components, greyscale );
    write_jpeg_scanlines( cinfo, image_buffer, row_stride );

    std::cout << "\n'" << name << "' written successfully." << std::endl;
}

#endif // CSMP_WITH_IMAGE_OUTPUT




/**

interpolates its values on the target grid where the two grids overlap, elsewhere, the
target grid is left untouched. Only the area inside the frame is considered.
*/

void FiniteDifferenceGrid::LinearInterpolateOnTo( FiniteDifferenceGrid& grid2 ) const
 {
    if ( grid2.MinX() > MinX() ) {
         cout <<"\nFiniteDifferenceGrid<double>::LinearInterpolateOnTo: target grid min_x too large."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::LinearInterpolateOnTo");
      }
    if ( grid2.MinY() > MinY() ) {
         cout <<"\nFiniteDifferenceGrid<double>::LinearInterpolateOnTo: target grid min_y too large."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::LinearInterpolateOnTo");
      }
    if ( grid2.MaxX() < MaxX() ) {
         cout <<"\nFiniteDifferenceGrid<double>::LinearInterpolateOnTo: target grid max_x too small."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::LinearInterpolateOnTo");
      }
    if ( grid2.MaxY() < MaxY() ) {
         cout <<"\nFiniteDifferenceGrid<double>::LinearInterpolateOnTo: target grid max_y too small."<< endl;
         throw length_error("FiniteDifferenceGrid<double>::LinearInterpolateOnTo");
      }

    for ( int32_t i=0; i<grid2.Rows(); i++ )
      for ( int32_t j=0; j<grid2.Columns(); j++ ) {
           ClosestGridPointTo( grid2.X(j), grid2.Y(i), y1, x1 );
           grid2(i,j) = (*this)(x1,y1);
        }

 } // end LinearInterpolateOn



/**

outputs interpolated value of grid-point at the location x,y
in the domain of the regular grid 0.0...x_max, 0.0...y_max.
NB: if the queried x,y coordinates lie outside of the grid,
operator(x,y) returns the closest value on the grid. It does
not give a warning.
*/

double  FiniteDifferenceGrid::operator()( double x,  double y, bool rounded ) const
  {
      // 1. finding x and y indices of interpolation points
      if ( rounded ) {
           x1 = static_cast<int32_t>(rint( (x - x_min) / xresolution ));
           y1 = static_cast<int32_t>(rint( (y - y_min) / yresolution ));
        }
      else {
           x1 = static_cast<int32_t>((x-x_min) / xresolution);
           y1 = static_cast<int32_t>((y-y_min) / yresolution);
        }

      // 2. if one of the points lies outside of the central grid
      if ( x < x_min || x > x_max || y < y_min || y > y_max )
        return ExtrapolateTo( x, y );

      // 3. getting the values at the interpolation points
      p1 = (*this)( y1,   x1   );
      if ( xfr==0 && x1+1 >= size_x ) x1 = size_x - 2;
      p2 = (*this)( y1,   x1+1 );
      if ( yfr==0 && y1+1 >= size_y ) y1 = size_y - 2;
      p3 = (*this)( y1+1, x1+1 );
      p4 = (*this)( y1+1, x1   );

      // 4. if all points are the same, this shortcut is possible
      if ( p1 == p2 && p2 == p3 && p3 == p4 ) return p1;

      // 5. computing interpolation functions. Num. Recip. p. 105
      t = ((x-x_min) - x1*xresolution) / xresolution;
      u = ((y-y_min) - y1*yresolution) / yresolution;

      // 6. bi-linear interpolation
      return (1.0-t)*(1.0-u)*p1 + t*(1.0-u)*p2 + t*u*p3 + (1.0-t)*u*p4;
  }




double FiniteDifferenceGrid::InterpolateOutside( double x, double y ) const
{
    x1 = static_cast<int32_t>((x-x_min) / xresolution);
    y1 = static_cast<int32_t>((y-y_min) / yresolution);

  // 1. if smaller than the minimum of the y-range value is interpolated on the
  //    horizontal model boundary.
    if ( y <= y_min )
      {
         assert( x >= x_min && x <= x_max );
         delta = (y - y_min) / yresolution;
         if ( (*this)(1,x1)   == (*this)(0,x1)   &&
              (*this)(0,x1)   == (*this)(1,x1+1) &&
              (*this)(1,x1+1) == (*this)(0,x1+1) ) return (*this)(1,x1);

         //      interpolation slope                        y        b
         p1 = (((*this)(1,x1)-(*this)(0,x1))/yresolution)*delta + (*this)(0,x1);
         p2 = (((*this)(1,x1+1)-(*this)(0,x1+1))/yresolution)*delta + (*this)(0,x1+1);
         delta = ((x-x_min) - x1*xresolution) / xresolution;
         // interpolation parallel to x
         return p1 + delta * ((p2-p1)/xresolution);
      }

  // 2. if larger than the maximum of the y-range
    if ( y >= y_max )
      {
         assert( x >= x_min && x <= x_max );
         delta = (y - y_max) / yresolution;
         if ( (*this)(size_y-1,x1)   == (*this)(size_y-2,x1)   &&
              (*this)(size_y-2,x1)   == (*this)(size_y-1,x1+1) &&
              (*this)(size_y-1,x1+1) == (*this)(size_y-2,x1+1) ) return (*this)(size_y-1,x1);

         //      interpolation slope                        y        b
         p1 = (((*this)(size_y-1,x1)-(*this)(size_y-2,x1))/yresolution)*delta + (*this)(size_y-1,x1);
         p2 = (((*this)(size_y-1,x1+1)-(*this)(size_y-2,x1+1))/yresolution)*delta + (*this)(size_y-1,x1+1);
         delta = ((x-x_min) - x1*xresolution) / xresolution;
         // interpolation parallel to x
         return p1 + delta * ((p2-p1)/xresolution);
      }


  // 3. if smaller than the minimum of the x(horizontal) range
    if ( x <= x_min )
      {
         assert( y >= y_min && y <= y_max );
         delta = (x - x_min) / xresolution;
         if ( (*this)(y1,1)   == (*this)(y1,0)   &&
              (*this)(y1,0)   == (*this)(y1+1,1) &&
              (*this)(y1+1,1) == (*this)(y1+1,0) ) return (*this)(y1,1);

         //      interpolation slope                        x        b
         p1 = (((*this)(y1,1)-(*this)(y1,0))/xresolution)*delta + (*this)(y1,0);
         p2 = (((*this)(y1+1,1)-(*this)(y1+1,0))/xresolution)*delta + (*this)(y1+1,0);
         delta = ((y-y_min) - y1*yresolution) / yresolution;
         // interpolation parallel to y
         return p1 + delta * ((p2-p1)/yresolution);
      }

  // 4. if larger than the maximum of the x range
    if ( x >= x_max )
      {
         assert( y >= y_min && y <= y_max );
         delta = (x - x_max) / xresolution;
         if ( (*this)(y1,size_x-1)   == (*this)(y1,size_x-2)   &&
              (*this)(y1,size_x-2)   == (*this)(y1+1,size_x-1) &&
              (*this)(y1+1,size_x-1) == (*this)(y1+1,size_x-2) ) return (*this)(y1,size_x-1);

         //      interpolation slope                        x        b
         p1 = (((*this)(y1,size_x-1)-(*this)(y1,size_x-2))/xresolution)*delta + (*this)(y1,size_x-1);
         p2 = (((*this)(y1+1,size_x-1)-(*this)(y1+1,size_x-2))/xresolution)*delta + (*this)(y1+1,size_x-1);
         delta = ((y-y_min) - y1*yresolution) / yresolution;
         // interpolation parallel to y
         return p1 + delta * ((p2-p1)/yresolution);
      }
    std::cerr <<"\nFiniteDifferenceGrid:InterpolateOutside: could not handle point: "<< x <<","<< y << std::endl;
    return std::numeric_limits<double>::signaling_NaN();

} // end InterpolateOutside



/**

Interpolates strictly within the confines of the grid:
If y is smaller than the minimum of the y-range value is interpolated on the
horizontal model boundary x is set to the smallest x value. Then the desired value
is interpolated along the grid boundary between the marginal 2 points.
*/
double FiniteDifferenceGrid::InterpolateWithin( double& x, double& y ) const
{
    // 0. special corner cases
    if ( x < x_min && y < y_min ) return (*this)(0,0);
    if ( x < x_min && y > y_max ) return (*this)(size_y-1,0);
    if ( x > x_max && y < y_min ) return (*this)(0,size_x-1);
    if ( x > x_max && y > y_max ) return (*this)(size_y-1,size_x-1);

    x1 = static_cast<int32_t>((x-x_min) / xresolution);
    y1 = static_cast<int32_t>((y-y_min) / yresolution);

    if ( y <= y_min )
      {
         y     = y_min;
         delta = ((x-x_min) - x1 * xresolution) / xresolution;
         //  interpolation value y=0 -> y=1, x, x1 at boundary parallel to y-axis
         p1 = (*this)(0,x1);
         p2 = (*this)(0,x1+1);
         // interpolation parallel to x
         return p1 + delta * ((p2 - p1) / xresolution);
      }
    if ( y >= y_max )
      {
         y     = y_max;
         delta = ((x-x_min) - x1 * xresolution) / xresolution;
         p1 = (*this)(size_y-1,x1);
         p2 = (*this)(size_y-1,x1+1);
         return p1 + delta * ((p2 - p1) / xresolution);
      }

    if ( x <= x_min )
      {
         x     = x_min;
         delta = ((y-y_min) - y1 * yresolution) / yresolution;
         p1 = (*this)(y1,0);
         p2 = (*this)(y1+1,0);
         return p1 + delta * ((p2 - p1) / yresolution);
      }
    if ( x >= x_max )
      {
         x     = x_max;
         delta = ((y-y_min) - y1 * yresolution) / yresolution;
         p1 = (*this)(y1,size_x-1);
         p2 = (*this)(y1+1,size_x-1);
         return p1 + delta * ((p2 - p1) / yresolution);
      }

    // normal type interpolation
    return (*this)(x,y,false);

} // end InterpolateWithin


} // end namespace csmp
