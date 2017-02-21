#include "FiniteDifferenceGrid.h"
#include "binaryReadWrite.h"
#include "CSMP_highLevelUtilities.h"


#ifdef CSMP_WITH_IMAGE_OUTPUT
#include "ColorPalette.h"
#include "jpeg/jpeglib.h"
#endif

using namespace std;

namespace csmp {


FiniteDifferenceGrid::FiniteDifferenceGrid()
  : grid(0), 
    xresolution(1.0), yresolution(1.0),
    x_min(0.0), y_min(0.0), x_max(0.0), y_max(0.0), z_min(0.0), z_max(0.0),
    size_x(0), size_y(0), size_z(0), rowlength(0), 
    xfr(0), yfr(0), zfr(0)
{ 
//   cout <<"\nFiniteDifferenceGrid: called default constructor..."<< endl;
}




FiniteDifferenceGrid::FiniteDifferenceGrid( double64 x_dim,       
                                            double64 y_dim, 
                                            double64 xres, double64 yres, int32 frame_width )
  : grid(0)
 {
    Initialize( x_dim, y_dim, xres, yres, frame_width );
 }






FiniteDifferenceGrid::FiniteDifferenceGrid( double64 x_dmin, double64 x_dmax, 
                                                double64 y_dmin, double64 y_dmax, 
                                                double64 xres, double64 yres, int32 frame_width )
  : grid(0)
 {
    Initialize( x_dmin, x_dmax, y_dmin, y_dmax, xres, yres, frame_width );
 }






bool FiniteDifferenceGrid::In( const char* fname )
 {
    ifstream ifs(fname);
    if ( !ifs )
      {
         cout <<"\nFiniteDifferenceGrid<double64>::In: input file cannot be opened. ";
         cout <<" Nothing was done..."<< endl;
         return false;
      }
    char intext[255];
    double64     fxmin, fxmax, fymin, fymax, fdx, fdy;
    int32  frows, fcols;

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
         cout <<"\nFiniteDifferenceGrid<double64>::In:  Cannot handle different X & Y frame sizes."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::In");
      }
    
    // rebuilding grid
    Initialize( fxmin, fxmax, fymin, fymax, xresolution, yresolution, xfr );

    // reading core data of grid from input file (not the frame)
    // extra row and column are used for the last values of the range
    for ( int32 i=-yfr; i<frows+yfr; i++ )
      for ( int32 j=-xfr; j<fcols+xfr; j++ ) ifs >> (*this)(i,j);
    ifs.close();

    cout <<"\nFiniteDifferenceGrid<double64>::In: grid build successfully from text file." << endl;
    cout.flush();

    return true;
 }






FiniteDifferenceGrid::FiniteDifferenceGrid( const FiniteDifferenceGrid& g )
 : grid(0) // strangely this is an absolute must
 {
    *this = g;
 }





FiniteDifferenceGrid::~FiniteDifferenceGrid()
 {
    delete[] grid;
 }




void FiniteDifferenceGrid::Initialize( double64 x_dim, double64 y_dim,
				       double64 xres, double64 yres, int32 frame_width )
 {
     Initialize( 0.0, x_dim, 0.0, y_dim, xres, yres, frame_width );
 }




void FiniteDifferenceGrid::Initialize( double64 x_dmin, double64 x_dmax, 
                                           double64 y_dmin, double64 y_dmax, 
                                           double64 xres, double64 yres, int32 frame_width )
 {
     // 0. removing old grid
     delete[] grid; 
     x_min       = x_dmin;
     x_max       = x_dmax;
     y_min       = y_dmin;
     y_max       = y_dmax; 
     xresolution = xres;
     yresolution = yres;
     xfr = yfr = zfr = frame_width;

     // 0. checks
     z_min = z_max = 0.0;
     double64 swap;
     if ( x_min > x_max ) { swap=x_min; x_min=x_max; x_max=swap; }
     if ( y_min > y_max ) { swap=y_min; y_min=y_max; y_max=swap; }

     // 1. getting the resolution
     size_x = static_cast<int32>((x_max-x_min) / xresolution + 1.);
     size_y = static_cast<int32>((y_max-y_min) / yresolution + 1.);
     size_z = 0;

     // 2. extra data member for speed
     rowlength = size_x + 2 * xfr;
     
     // 3. setting up the vector-grid
     cout <<"\nFiniteDifferenceGrid: building new grid; allocating ";
     cout << (((size_x+2*xfr) * (size_y+2*yfr) * sizeof(double64)) / 1.0e+6);
     cout <<" MByte of memory..." << endl;
     
     grid  = new double64[ (size_x+2*xfr) * (size_y+2*yfr) ];
     
     // 4. give initial value of zero
     *this = 0.0;
 }



    double64 FiniteDifferenceGrid::ResolutionX() const { return xresolution; }
    

    double64 FiniteDifferenceGrid::ResolutionY() const { return yresolution; }
    

    int32   FiniteDifferenceGrid::Rows()        const { return size_y-1; }
    

    int32   FiniteDifferenceGrid::Columns()     const { return size_x-1; }
    

    bool   FiniteDifferenceGrid::HasFrame()    const { return (xfr!=0||yfr!=0||zfr!=0); }
    

    double64 FiniteDifferenceGrid::MinX()        const { return x_min; }
    

    double64 FiniteDifferenceGrid::MinY()        const { return y_min; }
    

    double64 FiniteDifferenceGrid::MinZ()        const { return z_min; }
    

    double64 FiniteDifferenceGrid::MaxX()        const { return x_max; }
    

    double64 FiniteDifferenceGrid::MaxY()        const { return y_max; }
    

    double64 FiniteDifferenceGrid::MaxZ()        const { return z_max; }


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
    delete[] grid;
    grid = new double64[ (size_x+2*xfr) * (size_y+2*yfr) ];
    for ( int32 i=0; i<((size_x+2*xfr) * (size_y+2*yfr)); i++ ) grid[i] = g.grid[i];

    return *this;    
 }
 
 


FiniteDifferenceGrid& FiniteDifferenceGrid::operator=( double64 val )
 {
    for ( int32 i=0; i<((size_x+2*xfr)*(size_y+2*yfr)); i++ ) grid[i] = val;
    return *this;    
 }


void   FiniteDifferenceGrid::SetRowTo( int32 row, double64 val )
 {
    if ( row+xfr<0 || row>size_x+xfr ) {
         cout <<"\nFiniteDifferenceGrid<double64>::SetRowTo: Row index out of range: "<< row;
         cout <<" versus ("<< -xfr <<"-"<< size_x+xfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::SetRowTo");
      }
    for ( int32 j=-xfr; j<(size_x+xfr); j++ ) (*this)(row,j) = val;
 }



void   FiniteDifferenceGrid::SetColumnTo( int32 col, double64 val )
 {
    if ( col+yfr<0 || col>size_y+yfr ) {
         cout <<"\nFiniteDifferenceGrid<double64>::SetColumnTo: column index out of range: "<< col;
         cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::SetColumnTo");
      }
    for ( int32 i=-yfr; i<(size_y+yfr); i++ ) (*this)(i,col) = val;
 }




void   FiniteDifferenceGrid::SetFirstNRowsTo( int32 n, double64 val )
 {
    if ( n<0 || n>size_y ) {
         cout <<"\nFiniteDifferenceGrid<double64>::SetFirstNRowsTo: row index out of range: "<< n;
         cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::SetFirstNRowsTo");
      }
    for ( int32 i=0; i<n; i++ ) SetRowTo( i, val );
 }




void   FiniteDifferenceGrid::SetLastNRowsTo( int32 n, double64 val )
 {
    if ( n<0 || n>size_y ) {
         cout <<"\nFiniteDifferenceGrid<double64>::SetLastNRowsTo: row index out of range: "<< n;
         cout <<" versus ("<< -yfr <<"-"<< size_y+yfr <<")."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::SetLastNRowsTo");
      }
    for ( int32 i=size_y-1; i>(size_y-1-n); i-- ) SetRowTo( i, val );
 }

 

double64 FiniteDifferenceGrid::RowAverage( int32 row ) const
 {
    double64 avg = 0.0;
    assert( row+yfr >= 0 && row < size_y+yfr );
    for ( int32 j=-xfr; j<(size_x+xfr); j++ ) avg += Value(row,j);
    return avg /= size_x+2*xfr;
 }




double64 FiniteDifferenceGrid::ColumnAverage( int32 col ) const
 {
    double64 avg = 0.0;
    assert( col+xfr >= 0 && col < size_x+xfr );
    for ( int32 i=-yfr; i<(size_y+yfr); i++ ) avg += Value(i,col);
    return avg /= size_y+2*yfr;
 } 



void   FiniteDifferenceGrid::SetRegionTo( double64 xmin, double64 xmax, double64 ymin, double64 ymax, double64 val )
 {
     int32 i, j, imin, imax, jmin, jmax;
     
     assert( xmin < xmax );
     assert( ymin < ymax );
     assert( xmax <= x_max );
     assert( ymax <= y_max );
     
     ClosestGridPointTo( xmin, ymin, imin, jmin ); 
     ClosestGridPointTo( xmax, ymax, imax, jmax ); 

     for ( i=imin; i<=imax; i++ )
       for ( j=jmin; j<=jmax; j++ ) (*this)(i,j) = val;
 }




void   FiniteDifferenceGrid::Out(std::ostream& os) const
 { 
    int32 i, j;
    os <<"\n\nFiniteDifferenceGrid<double64>::Out(): printing grid of size: ";
    os << x_max-x_min <<" by "<< y_max-y_min <<" m" << endl;
    os <<"x-range:    "<< x_min <<"-"<< x_max <<" m" << endl;
    os <<"y-range:    "<< y_min <<"-"<< y_max <<" m" << endl;
    os <<"resolution (X): "<< xresolution <<" m"        << endl;
    os <<"resolution (Y): "<< yresolution <<" m"        << endl;
    if ( HasFrame() ) os <<"grid with frame of width: "<< xfr << endl;
    os << size_y <<" rows, "<< size_x <<" columns" << endl;
    os  << "\n\t";
    for ( j=-xfr; j<size_x+xfr; j++ ) os << "col "<< j <<"\t";
    
    for ( i=-yfr; i<size_y+yfr; i++ )
      {
         os <<"\nrow "<< i <<"\t";
         for ( j=-xfr; j<size_x+xfr; j++ ) 
           os << (*this).Value(i,j) <<"\t";
      }
    os <<"\n\nDone..." << endl;
 }




void   FiniteDifferenceGrid::Out( const char* fname, int32 tstep, bool with_frame ) const
 {
    char name[200], num[30];
    strcpy( name, fname );
    sprintf( num, "%d", tstep );
    strcat( name, num );

    ofstream ifs( name );    
    assert( ifs.is_open() );
 
    ifs <<"FiniteDifferenceGrid<double64>::Out: ASCII text file from double64 grid." << endl;
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
      for ( int32 i=0; i<size_y; i++ ) {
           for ( int32 j=0; j<size_x; j++ ) ifs << (*this).Value(i,j) <<"\t";
           ifs << endl;
        }
    else
      for ( int32 i=-yfr; i<size_y+yfr; i++ ) {
           for ( int32 j=-xfr; j<size_x+xfr; j++ ) ifs << (*this).Value(i,j) <<"\t";
           ifs << endl;
        }

    ifs.close();
    cout <<"\n\n'" << name <<"' written successfully..." << endl;
 }



bool   FiniteDifferenceGrid::BinaryOut( const char* bin_name, int32 tstep, bool with_frame ) const
 {
    char name[200], num[30], heading[200];
    strcpy( name, bin_name );
    sprintf( num, "%d", tstep );
    strcat( name, num );
    strcpy( heading, "FiniteDifferenceGrid::BinaryOut: double64 grid as binary file");
 
    FILE*  fp;
    if ( (fp=fopen( name, "wb+")) == NULL ) {
        cout <<"\nFiniteDifferenceGrid::BinaryOut: File: "<< bin_name << " could not be opened"<< endl;
        return false;
      }
    skm_C_fwrite( fp, heading );  

    // stores dimensions of grid
    std::vector<double64>    dim_fT(8);
    // stores rows, columns, and frame
    std::vector<int32> dim_int(4);
    // grid data
    std::vector<double64>    grid_data((static_cast<size_t>(size_x)+2*static_cast<size_t>(xfr)) * (static_cast<size_t>(size_y)+2*static_cast<size_t>(yfr)));
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
    
    skm_C_fwrite( fp, dim_fT ); 
    skm_C_fwrite( fp, dim_int ); 
    
    // writing grid data to storage vector
    if ( !with_frame )
      for ( int32 i=0; i<size_y; i++ ) {
           for ( int32 j=0; j<size_x; j++ ) {
               grid_data[counter] = (*this).Value(i,j);
               counter++;
             }
        }
    else
      for ( int32 i=-yfr; i<size_y+yfr; i++ ) {
           for ( int32 j=-xfr; j<size_x+xfr; j++ ) {
               grid_data[counter] = (*this).Value(i,j);
               counter++;
             }
        }
    
   skm_C_fwrite( fp, grid_data ); 
   fclose( fp );

   cout <<"\n\n'" << bin_name <<"' written successfully..." << endl;
   
   return true;
 }



bool FiniteDifferenceGrid::BinaryIn( const char* bin_name )
 {
    FILE*  fp;  
    if ( (fp=fopen( bin_name, "rb")) == NULL ) {
        cout <<"\nFiniteDifferenceGrid::BinaryIn: File: "<< bin_name;
        cout <<" could not be opened"<< endl;
        return false;
     }
    char heading[200];
    skm_C_fread( fp, heading ); 
    cout <<"\nFiniteDifferenceGrid::BinaryIn: Reading: "<< heading << endl;

    // read dimensions of grid
    std::vector<double64>  dim_fT;
    // read rows, columns, and frame
    std::vector<int32> dim_int;
    // read data
    std::vector<double64>  grid_data;
    size_t counter(0);

    // read info about dimension
    skm_C_fread( fp, dim_fT ); 
    skm_C_fread( fp, dim_int ); 

    if ( dim_int[2] != dim_int[3] ) {
         cout <<"\nFiniteDifferenceGrid::BinaryIn:  Cannot handle different X & Y frame sizes."<< endl;
         throw length_error("FiniteDifferenceGrid::BinaryIn");
      }

    // rebuild grid
    Initialize( dim_fT[2], dim_fT[3], dim_fT[4], dim_fT[5], dim_fT[6], dim_fT[7], xfr );

    // read data and transfer to grid
    skm_C_fread( fp, grid_data ); 
    for ( int32 i=-dim_int[3]; i<dim_int[1]+dim_int[3]; i++ ) {
        for ( int32 j=-dim_int[2]; j<dim_int[0]+dim_int[2]; j++ ) {
            (*this)(i,j) = grid_data[counter];
            counter++;
          }
      }  

    fclose( fp );
    
    cout <<"\nFiniteDifferenceGrid<double64>::BinaryIn: grid build successfully from binary file." << endl;
    cout.flush();

    return true;
 }



void  FiniteDifferenceGrid::DataMinMaxWithoutNAN( double64& dmin, double64& dmax,
					                                  bool frame_included ) const
 {
    int32  i, j;
    double64 val;
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
        dmax=dmin=(*this).Value(0L,0L);
        for ( i=0; i<size_y; i++ )
          for ( j=0; j<size_x; j++ ) 
           {
               if ( !isnan( (*this).Value(i,j) ) ) {
                   if ( first_call ) { 
                       dmin = dmax = (*this).Value(i,j);
                       first_call = false;
                     }
                   val = (*this).Value(i,j);
                   if ( val > dmax ) dmax = val;
                   if ( val < dmin ) dmin = val;
                 }
           }
      }    
 } // end



void  FiniteDifferenceGrid::DataMinMax( double64& dmin, double64& dmax,
					                                  bool frame_included ) const
 {
    int32  i, j;
    double64 val;
         
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
        dmax=dmin=(*this).Value(0L,0L);
        for ( i=0; i<size_y; i++ )
          for ( j=0; j<size_x; j++ ) 
           {
             val = (*this).Value(i,j);
             if ( val > dmax ) dmax = val;
             if ( val < dmin ) dmin = val;
           }
      }    
 } // end



///scales all grid data to range: dmin to dmax
void   FiniteDifferenceGrid::ScaleDataToRange( double64 dmin, double64 dmax )
 {
    double64  old_min, old_max;
    DataMinMax( old_min, old_max, false );
    double64  old_range = old_max - old_min;
    double64  new_range = dmax - dmin;  

    for ( int32 i=0; i<((size_x+2*xfr) * (size_y+2*yfr)); i++ ) 
      grid[i] = dmin + ((grid[i] - old_min)/old_range) * new_range;

 } // end ScaleDataToRange











/// 4 point, third order polynomial extrapolation
void  FiniteDifferenceGrid::PolynomialInterpolation( double64* xa, double64* ya, 
                                                     double64  x,  double64& y, 
                                                     double64& dy ) const
{
    // constraint point arrays are: xa[], double64 ya[], 
    double64  c[5], d[5];
    int32     i, m, ns=1;
    double64  den, dif, dift, ho, hp, w;

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
			     cout <<"\nFiniteDifferenceGrid<double64>::PolynomialInterpolation: ";
			     cout <<"Error. No interpolation was accomplished."<< endl;
			     cout <<"Method input: x=" << x <<", modified y="<< y <<", error="<< dy << endl;
			     cout <<"Constraint point arrays: "<< endl <<"position data: ";
			     for ( int32 j=1; j<=4; j++ ) cout << xa[j] <<" ";
			     cout << endl <<"dependent variable data: ";
			     for ( int32 j=1; j<=4; j++ ) cout << ya[j] <<" ";
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
double64  FiniteDifferenceGrid::ExtrapolateTo( double64 x, double64 y ) const
 {
    static double64  xa[5], ya[5];
    int32        i;
    double64         a;
    double64         err;
    
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
    double64 m, len_fac;
    
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
              // scaling x to reflect pathlength aint32 sloping line
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

/* writes Pixmaps with the following header:
P3       // P3 signifies ascii format
# CREATOR: CSMP FiniteDifferenceGrid::SaveToPPM() Matthai & Roberts 1999
114 15   // rows x (columns * 3)
255      // colors
3  3  5    5  6  7  // R G B for each pixel
*/

void   FiniteDifferenceGrid::SaveToPPM( const char* filename, int32 timestep, bool greyscale ) const
 {
    // rgb color stuff (created only once)
    static ColorPalette  rgb_colorizer;  // default rainbow scale 1-255
    float                colors[4];      // color values RGB A
    if ( greyscale )     rgb_colorizer.MakeGreyPalette();
    
    // file name + extension
    char name[200], num[30];
    strcpy( name, filename );
    sprintf( num, "%d", timestep );
    strcat( name, num );
    strcat( name, ".ppm" );

    ofstream ifs( name );
    assert( ifs.is_open() );
    
    // obtain data range to scale the data to 0-255 for output
    double64  old_min, old_max;
    DataMinMax( old_min, old_max, false );
    double64  old_range = old_max - old_min;
    double64  new_range = 255.0; 
    double64  out_val; 
    int32    i, j;

    // data format tag 
    ifs <<"P3"<< endl;
    // binary output: ifs <<"P6"<< endl;

    // creator comment line
    ifs <<"# CREATOR: CSMP FiniteDifferenceGrid<double64>::SaveToPPM() Matthai & Roberts 1999." << endl;
    
    // writing the file name, timestep and data range into a comment line
    ifs <<"# Filename: "<< filename <<", timestep: "<< timestep;
    ifs <<", data value, min: "<< old_min <<", max: "<< old_max << endl;
    
    // output file size specifier
    ifs << size_x <<" "<< size_y << endl;
    
    // how many color values
    ifs << 255 << endl;    

    // creating RGB datablock for central grid portion without frame 
    for ( i=size_y-1; i>=0; i-- )
      {
        for ( j=0; j<size_x; j++ ) 
          {
             // scaling value to 0-256 scale
             out_val = (((*this).Value(i,j) - old_min)/old_range) * new_range;
             // getting red,green, blue color values
             rgb_colorizer.GiveRgb( static_cast<float>(out_val), colors );
             ifs.width(3);
             ifs << static_cast<short>(colors[0]) <<" ";
             ifs.width(3);
             ifs << static_cast<short>(colors[1]) <<" ";
             ifs.width(3);
             ifs << static_cast<short>(colors[2]) <<" ";
          }  
        ifs << endl;
      }

    ifs.close();
    cout <<"\n\n'" << name <<"' written successfully..." << endl;
          
 } // end SaveToPPM





/**

writes Grid to JPG image. Default is RGB in which case the method writes 
24 bit colors (3 1-byte floats for each pixel). If greyscale is chosen the 
writing becomes much faster and the image file will be only one-third
of the RGB image size. If a property value is NAN, the according color will be black.
*/

void FiniteDifferenceGrid::SaveToJPGWithoutNAN( const char* filename, int32 timestep, 
                                                    bool greyscale, bool sqrt_of ) const
 {
    // rgb color stuff (created only once)
    ColorPalette      rgb_colorizer;  // default rainbow scale 1-255
    float             colors[4];      // color values RGB A
    if ( greyscale )  rgb_colorizer.MakeGreyPalette();
    
    // file name + extension
    char num[30];
    sprintf( num, "%d", timestep );
    string  name(filename), padded_string( num );
    replaceWhiteSpaceBy( padded_string, '0' );
    name += padded_string;
    name +=".jpg";

    // obtain data range to scale the data to 0-255 for output
    double64  old_min, old_max;
    DataMinMaxWithoutNAN( old_min, old_max, false );
    cout.setf( ios::scientific );
    cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' Output data range: "<< old_min<<" to "<< old_max << endl;
    cout.unsetf( ios::scientific );
    
    // if old_min = old_max a textfile with the variable name and value is output instead
    // (there is no need to go through the whole JPEG procedure if only a single value would be output)
    ofstream ofs;

    if ( old_min == old_max )
      {
         name +=".txt";
         ofs.open ( name.c_str(), ios::out|ios::trunc );
         if ( !ofs ) 
           cout <<"\nOutput file: "<< name <<" could not be opened. Nothing was done..." << endl;
         else 
           {
              ofs << name <<": ";
              ofs <<" Uniform variable value: "<< old_max <<" at timestep: "<< timestep << endl;
              ofs.close();
              cout <<"\n\n'" << name <<"' written instead of JPG file..." << endl;
           }
         return;
      }
    
    // creating a scale factor for the data
    if ( sqrt_of )
      {
         old_min = sqrt( old_min );
         old_max = sqrt( old_max );
      }
    double64  old_range = old_max - old_min;
    double64  new_range = 255.0; 
    double64  out_val; 
    int32    i, j;

    // ----------------------------------------------------------------------------------
    // JPG Stuff
    // ----------------------------------------------------------------------------------
    // creating RGB scanline memory for central grid portion without frame
    JSAMPROW          row_pointer[1];
    int32             image_components = 3;
    if ( greyscale )  image_components = 1;
    int32             row_stride = size_x * image_components;
    unsigned char*    image_buffer;
    uint32            incr;
  
    // allocating image memory buffer and storing RGB values within it  
    image_buffer = new unsigned char[ (size_x * size_y * image_components) ];

      for ( incr=0, i=size_y-1; i>=0; i-- )
        for ( j=0; j<size_x; j++ ) 
          {
             // scaling value to 0-256 scale
             if ( sqrt_of ) {
                 if ( isnan( (*this).Value(i,j) ) ) out_val = 255.0;
                 else                               out_val = ((sqrt((*this).Value(i,j)) - old_min)/old_range) * new_range;
               }
             else {
                 if ( isnan( (*this).Value(i,j) ) ) out_val = 255.0;
                 else                               out_val = (((*this).Value(i,j) - old_min)/old_range) * new_range;
               }
             if ( greyscale ) image_buffer[incr++] = static_cast<unsigned char>(out_val);
             else
               {
                  // getting red,green, blue color values
                  if ( isnan( (*this).Value(i,j) ) ) colors[0] = colors[1] = colors[2] = colors[3] = 255;
                  else                               rgb_colorizer.GiveRgb( static_cast<float>(out_val), colors );
                  image_buffer[incr++] = static_cast<unsigned char>(colors[0]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[1]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[2]);
               }
          }
 
    // setting up the jpg storage structures
    jpeg_compress_struct  cinfo;
    jpeg_error_mgr        jerr;
    
    // initializing the error manager in the compression object
    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    
    FILE*  outfile;
    if ((outfile = fopen(name.c_str(), "wb")) == NULL ) 
      {
          cout <<"\nFiniteDifferenceGrid:SaveToJPG: cannot open outputfile: "<< name << endl;
          return;
      }
    jpeg_stdio_dest( &cinfo, outfile );
    
    // setting up image size and colorspace
    cinfo.image_width      = static_cast<uint32>(size_x);
    cinfo.image_height     = static_cast<uint32>(size_y);
    cinfo.input_components = 3;        // color values per pixel
    cinfo.in_color_space   = JCS_RGB;  // RGB or JCS_GRAY_SCALE (only 1 val per pixel)
    if ( greyscale ) 
      {
        cinfo.input_components = 1;        // color values per pixel
        cinfo.in_color_space   = JCS_GRAYSCALE; 
      }

    // assigning the values
    jpeg_set_defaults( &cinfo );
    
    // writing data to file: TRUE for complete jpg interchange datastream
    jpeg_start_compress( &cinfo, TRUE );
    while ( cinfo.next_scanline < cinfo.image_height )
      {
          row_pointer[0] = &image_buffer[ cinfo.next_scanline * row_stride ];
          jpeg_write_scanlines( &cinfo, row_pointer, 1 );
      }
    jpeg_finish_compress( &cinfo );
    fclose( outfile );

    // cleanup
    jpeg_destroy_compress( &cinfo );  
    delete[] image_buffer;

    cout <<"\n'" << name <<"' written successfully..." << endl;
          
 } // end SaveToJPG







/**

As method above but allows to specify a fixed value range, such that the
color range does not vary from timestep to timestep.  
*/

void FiniteDifferenceGrid::SaveToJPGWithoutNAN( const char* filename, int32 timestep, 
                                                    double64 vmin, double64 vmax, 
                                                    bool greyscale, bool sqrt_of ) const
 {
    // rgb color stuff (created only once)
    ColorPalette         rgb_colorizer;  // default rainbow scale 1-255
    float                colors[4];      // color values RGB A
    if ( greyscale )     rgb_colorizer.MakeGreyPalette();
    
    // file name + extension
    char num[30];
    sprintf( num, "%d", timestep );
    string  name(filename), padded_string( num );
    replaceWhiteSpaceBy( padded_string, '0' );
    name += padded_string;
    name +=".jpg";

    // obtain data range to scale the data to 0-255 for output
    double64  old_min, old_max;
    DataMinMaxWithoutNAN( old_min, old_max, false );
    cout.setf( ios::scientific );
    cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' Output data range: "<< old_min<<" to "<< old_max << endl;
    cout.unsetf( ios::scientific );
    
    if ( vmin > old_min ) {
         cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' ";
         cout <<"Data minimum is smaller than prescribed minimum; setting minimum to actual value."<< endl; 
      }
    // if the actual values are within the prescribed range, this range is used  
    else old_min = vmin;
    
    if ( vmax < old_max ) {
         cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' ";
         cout <<"Data maximum is larger than prescribed maximum; setting maximum to actual value."<< endl; 
      }
    else old_max = vmax;
    
    // if old_min = old_max a textfile with the variable name and value is output instead
    // (there is no need to go through the whole JPEG procedure if only a single value would be output)
    ofstream ofs;

    if ( old_min == old_max )
      {
         name +=".txt";
         ofs.open ( name.c_str(), ios::out|ios::trunc );
         if ( !ofs ) 
           cout <<"\nOutput file: "<< name <<" could not be opened. Nothing was done..." << endl;
         else 
           {
              ofs << name <<": ";
              ofs <<" Uniform variable value: "<< old_max <<" at timestep: "<< timestep << endl;
              ofs.close();
              cout <<"\n\n'" << name <<"' written instead of JPG file..." << endl;
           }
         return;
      }
    
    // creating a scale factor for the data
    if ( sqrt_of )
      {
         old_min = sqrt( old_min );
         old_max = sqrt( old_max );
      }
    double64  old_range = old_max - old_min;
    double64  new_range = 255.0; 
    double64  out_val; 
    int32    i, j;

    // ----------------------------------------------------------------------------------
    // JPG Stuff
    // ----------------------------------------------------------------------------------
    // creating RGB scanline memory for central grid portion without frame
    JSAMPROW          row_pointer[1];
    int               image_components = 3;
    if ( greyscale )  image_components = 1;
    int               row_stride = size_x * image_components;
    unsigned char*    image_buffer;
    uint32            incr;
  
    // allocating image memory buffer and storing RGB values within it  
    image_buffer = new unsigned char[ (size_x * size_y * image_components) ];

      for ( incr=0, i=size_y-1; i>=0; i-- )
        for ( j=0; j<size_x; j++ ) 
          {
             // scaling value to 0-256 scale
             if ( sqrt_of ) {
                 if ( isnan( (*this).Value(i,j) ) ) out_val = 255.0;
                 else                               out_val = ((sqrt((*this).Value(i,j)) - old_min)/old_range) * new_range;
               }
             else {
                 if ( isnan( (*this).Value(i,j) ) ) out_val = 255.0;
                 else                               out_val = (((*this).Value(i,j) - old_min)/old_range) * new_range;
               }
             if ( greyscale ) image_buffer[incr++] = static_cast<unsigned char>(out_val);
             else
               {
                  // getting red,green, blue color values
                  if ( isnan( (*this).Value(i,j) ) ) colors[0] = colors[1] = colors[2] = colors[3] = 255;
                  else                               rgb_colorizer.GiveRgb( static_cast<float>(out_val), colors );
                  image_buffer[incr++] = static_cast<unsigned char>(colors[0]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[1]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[2]);
               }

 
          }
 
    // setting up the jpg storage structures
    jpeg_compress_struct  cinfo;
    jpeg_error_mgr        jerr;
    
    // initializing the error manager in the compression object
    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    
    FILE*  outfile;
    if ((outfile = fopen(name.c_str(), "wb")) == NULL ) 
      {
          cout <<"\nFiniteDifferenceGrid:SaveToJPG: cannot open outputfile: "<< name << endl;
          return;
      }
    jpeg_stdio_dest( &cinfo, outfile );
    
    // setting up image size and colorspace
    cinfo.image_width      = static_cast<uint32>(size_x);
    cinfo.image_height     = static_cast<uint32>(size_y);
    cinfo.input_components = 3;        // color values per pixel
    cinfo.in_color_space   = JCS_RGB;  // RGB or JCS_GRAY_SCALE (only 1 val per pixel)
    if ( greyscale ) 
      {
        cinfo.input_components = 1;        // color values per pixel
        cinfo.in_color_space   = JCS_GRAYSCALE; 
      }

    // assigning the values
    jpeg_set_defaults( &cinfo );
    
    // writing data to file: TRUE for complete jpg interchange datastream
    jpeg_start_compress( &cinfo, TRUE );
    while ( cinfo.next_scanline < cinfo.image_height )
      {
          row_pointer[0] = &image_buffer[ cinfo.next_scanline * row_stride ];
          jpeg_write_scanlines( &cinfo, row_pointer, 1 );
      }
    jpeg_finish_compress( &cinfo );
    fclose( outfile );

    // cleanup
    jpeg_destroy_compress( &cinfo );  
    delete[] image_buffer;

    cout <<"\n'" << name <<"' written successfully..." << endl;
          
 } // end SaveToJPG (fixed range)


/**

writes Grid to JPG image. Default is RGB in which case the method writes 
24 bit colors (3 1-byte floats for each pixel). If greyscale is chosen the 
writing becomes much faster and the image file will be only one-third
of the RGB image size.

'sqrt' (default=false) allows you to store the square root of the value.
*/

void FiniteDifferenceGrid::SaveToJPG( const char* filename, int32 timestep, 
                                      bool greyscale, bool sqrt_of ) const
 {
    // rgb color stuff (created only once)
    ColorPalette  rgb_colorizer;  // default rainbow scale 1-255
    float                colors[4];      // color values RGB A
    if ( greyscale )     rgb_colorizer.MakeGreyPalette();
    
    // file name + extension
    char num[30];
    sprintf( num, "%d", timestep );
    string  name(filename), padded_string( num );
    replaceWhiteSpaceBy( padded_string, '0' );
    name += padded_string;
    name +=".jpg";

    // obtain data range to scale the data to 0-255 for output
    double64  old_min, old_max;
    DataMinMax( old_min, old_max, false );
    cout.setf( ios::scientific );
    cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' Output data range: "<< old_min<<" to "<< old_max << endl;
    cout.unsetf( ios::scientific );
    
    // if old_min = old_max a textfile with the variable name and value is output instead
    // (there is no need to go through the whole JPEG procedure if only a single value would be output)
    ofstream ofs;

    if ( old_min == old_max )
      {
         name +=".txt";
         ofs.open ( name.c_str(), ios::out|ios::trunc );
         if ( !ofs ) 
           cout <<"\nOutput file: "<< name <<" could not be opened. Nothing was done..." << endl;
         else 
           {
              ofs << name <<": ";
              ofs <<" Uniform variable value: "<< old_max <<" at timestep: "<< timestep << endl;
              ofs.close();
              cout <<"\n\n'" << name <<"' written instead of JPG file..." << endl;
           }
         return;
      }
    
    // creating a scale factor for the data
    if ( sqrt_of )
      {
         old_min = sqrt( old_min );
         old_max = sqrt( old_max );
      }
    double64  old_range = old_max - old_min;
    double64  new_range = 255.0; 
    double64  out_val; 
    int32    i, j;

    // ----------------------------------------------------------------------------------
    // JPG Stuff
    // ----------------------------------------------------------------------------------
    // creating RGB scanline memory for central grid portion without frame
    JSAMPROW          row_pointer[1];
    int32             image_components = 3;
    if ( greyscale )  image_components = 1;
    int32             row_stride = size_x * image_components;
    unsigned char*    image_buffer;
    uint32            incr;
  
    // allocating image memory buffer and storing RGB values within it  
    image_buffer = new unsigned char[ (size_x * size_y * image_components) ];

      for ( incr=0, i=size_y-1; i>=0; i-- )
        for ( j=0; j<size_x; j++ ) 
          {
             // scaling value to 0-256 scale
             if ( sqrt_of ) out_val = ((sqrt((*this).Value(i,j)) - old_min)/old_range) * new_range;
             else           out_val = (((*this).Value(i,j) - old_min)/old_range) * new_range;
             if ( greyscale ) image_buffer[incr++] = static_cast<unsigned char>(out_val);
             else
               {
                  // getting red,green, blue color values
                  rgb_colorizer.GiveRgb( static_cast<float>(out_val), colors );
                  image_buffer[incr++] = static_cast<unsigned char>(colors[0]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[1]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[2]);
               }
          }
 
    // setting up the jpg storage structures
    jpeg_compress_struct  cinfo;
    jpeg_error_mgr        jerr;
    
    // initializing the error manager in the compression object
    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    
    FILE*  outfile;
    if ((outfile = fopen( name.c_str(), "wb")) == NULL ) 
      {
          cout <<"\nFiniteDifferenceGrid:SaveToJPG: cannot open outputfile: "<< name << endl;
          return;
      }
    jpeg_stdio_dest( &cinfo, outfile );
    
    // setting up image size and colorspace
    cinfo.image_width      = static_cast<uint32>(size_x);
    cinfo.image_height     = static_cast<uint32>(size_y);
    cinfo.input_components = 3;        // color values per pixel
    cinfo.in_color_space   = JCS_RGB;  // RGB or JCS_GRAY_SCALE (only 1 val per pixel)
    if ( greyscale ) 
      {
        cinfo.input_components = 1;        // color values per pixel
        cinfo.in_color_space   = JCS_GRAYSCALE; 
      }

    // assigning the values
    jpeg_set_defaults( &cinfo );
    
    // writing data to file: TRUE for complete jpg interchange datastream
    jpeg_start_compress( &cinfo, TRUE );
    while ( cinfo.next_scanline < cinfo.image_height )
      {
          row_pointer[0] = &image_buffer[ cinfo.next_scanline * row_stride ];
          jpeg_write_scanlines( &cinfo, row_pointer, 1 );
      }
    jpeg_finish_compress( &cinfo );
    fclose( outfile );

    // cleanup
    jpeg_destroy_compress( &cinfo );  
    delete[] image_buffer;

    cout <<"\n'" << name <<"' written successfully..." << endl;
          
 } // end SaveToJPG







/**

As method above but allows to specify a fixed value range, such that the
color range does not vary from timestep to timestep.  
*/

void FiniteDifferenceGrid::SaveToJPG( const char* filename, int32 timestep, 
                                      double64 vmin, double64 vmax, 
                                      bool greyscale, bool sqrt_of ) const
 {
    // rgb color stuff (created only once)
    static ColorPalette  rgb_colorizer;  // default rainbow scale 1-255
    float                colors[4];      // color values RGB A
    if ( greyscale )     rgb_colorizer.MakeGreyPalette();
    
    // file name + extension
    char num[30];
    sprintf( num, "%d", timestep );
    string  name(filename), padded_string( num );
    replaceWhiteSpaceBy( padded_string, '0' );
    name += padded_string;
    name +=".jpg";

    // obtain data range to scale the data to 0-255 for output
    double64  old_min, old_max;
    DataMinMax( old_min, old_max, false );
    cout.setf( ios::scientific );
    cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' Output data range: "<< old_min<<" to "<< old_max << endl;
    cout.unsetf( ios::scientific );
    
    if ( vmin > old_min ) {
         cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' ";
         cout <<"Data minimum is smaller than prescribed minimum; setting minimum to actual value."<< endl; 
      }
    // if the actual values are within the prescribed range, this range is used  
    else old_min = vmin;
    
    if ( vmax < old_max ) {
         cout <<"\nFiniteDifferenceGrid:SaveToJPG: '"<< filename <<"' ";
         cout <<"Data maximum is larger than prescribed maximum; setting maximum to actual value."<< endl; 
      }
    else old_max = vmax;
    
    // if old_min = old_max a textfile with the variable name and value is output instead
    // (there is no need to go through the whole JPEG procedure if only a single value would be output)
    ofstream ofs;

    if ( old_min == old_max )
      {
         name +=".txt";
         ofs.open ( name.c_str(), ios::out|ios::trunc );
         if ( !ofs ) 
           cout <<"\nOutput file: "<< name <<" could not be opened. Nothing was done..." << endl;
         else 
           {
              ofs << name <<": ";
              ofs <<" Uniform variable value: "<< old_max <<" at timestep: "<< timestep << endl;
              ofs.close();
              cout <<"\n\n'" << name <<"' written instead of JPG file..." << endl;
           }
         return;
      }
    
    // creating a scale factor for the data
    if ( sqrt_of )
      {
         old_min = sqrt( old_min );
         old_max = sqrt( old_max );
      }
    double64  old_range = old_max - old_min;
    double64  new_range = 255.0; 
    double64  out_val; 
    int32    i, j;

    // ----------------------------------------------------------------------------------
    // JPG Stuff
    // ----------------------------------------------------------------------------------
    // creating RGB scanline memory for central grid portion without frame
    JSAMPROW          row_pointer[1];
    int               image_components = 3;
    if ( greyscale )  image_components = 1;
    int               row_stride = size_x * image_components;
    unsigned char*    image_buffer;
    uint32            incr;
  
    // allocating image memory buffer and storing RGB values within it  
    image_buffer = new unsigned char[ (size_x * size_y * image_components) ];

      for ( incr=0, i=size_y-1; i>=0; i-- )
        for ( j=0; j<size_x; j++ ) 
          {
             // scaling value to 0-256 scale
             if ( sqrt_of ) out_val = ((sqrt((*this).Value(i,j)) - old_min)/old_range) * new_range;
             else           out_val = (((*this).Value(i,j) - old_min)/old_range) * new_range;
             if ( greyscale ) image_buffer[incr++] = static_cast<unsigned char>(out_val);
             else
               {
                  // getting red,green, blue color values
                  rgb_colorizer.GiveRgb( static_cast<float>(out_val), colors );
                  image_buffer[incr++] = static_cast<unsigned char>(colors[0]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[1]);
                  image_buffer[incr++] = static_cast<unsigned char>(colors[2]);
               }
          }
 
    // setting up the jpg storage structures
    jpeg_compress_struct  cinfo;
    jpeg_error_mgr        jerr;
    
    // initializing the error manager in the compression object
    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress( &cinfo );
    
    FILE*  outfile;
    if ((outfile = fopen( name.c_str(), "wb")) == NULL ) 
      {
          cout <<"\nFiniteDifferenceGrid:SaveToJPG: cannot open outputfile: "<< name << endl;
          return;
      }
    jpeg_stdio_dest( &cinfo, outfile );
    
    // setting up image size and colorspace
    cinfo.image_width      = static_cast<uint32>(size_x);
    cinfo.image_height     = static_cast<uint32>(size_y);
    cinfo.input_components = 3;        // color values per pixel
    cinfo.in_color_space   = JCS_RGB;  // RGB or JCS_GRAY_SCALE (only 1 val per pixel)
    if ( greyscale ) 
      {
        cinfo.input_components = 1;        // color values per pixel
        cinfo.in_color_space   = JCS_GRAYSCALE; 
      }

    // assigning the values
    jpeg_set_defaults( &cinfo );
    
    // writing data to file: TRUE for complete jpg interchange datastream
    jpeg_start_compress( &cinfo, TRUE );
    while ( cinfo.next_scanline < cinfo.image_height )
      {
          row_pointer[0] = &image_buffer[ cinfo.next_scanline * row_stride ];
          jpeg_write_scanlines( &cinfo, row_pointer, 1 );
      }
    jpeg_finish_compress( &cinfo );
    fclose( outfile );

    // cleanup
    jpeg_destroy_compress( &cinfo );  
    delete[] image_buffer;

    cout <<"\n'" << name <<"' written successfully..." << endl;
          
 } // end SaveToJPG (fixed range)






#endif // end CSP_WITH_IMAGE_OUTPUT




/**

interpolates its values on the target grid where the two grids overlap, elsewhere, the
target grid is left untouched. Only the area inside the frame is considered.
*/

void FiniteDifferenceGrid::LinearInterpolateOnTo( FiniteDifferenceGrid& grid2 ) const
 {
    if ( grid2.MinX() > MinX() ) {
         cout <<"\nFiniteDifferenceGrid<double64>::LinearInterpolateOnTo: target grid min_x too large."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::LinearInterpolateOnTo");
      }
    if ( grid2.MinY() > MinY() ) {
         cout <<"\nFiniteDifferenceGrid<double64>::LinearInterpolateOnTo: target grid min_y too large."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::LinearInterpolateOnTo");
      }
    if ( grid2.MaxX() < MaxX() ) {
         cout <<"\nFiniteDifferenceGrid<double64>::LinearInterpolateOnTo: target grid max_x too small."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::LinearInterpolateOnTo");
      }
    if ( grid2.MaxY() < MaxY() ) {
         cout <<"\nFiniteDifferenceGrid<double64>::LinearInterpolateOnTo: target grid max_y too small."<< endl;
         throw length_error("FiniteDifferenceGrid<double64>::LinearInterpolateOnTo");
      }
 
    for ( int32 i=0; i<grid2.Rows(); i++ )
      for ( int32 j=0; j<grid2.Columns(); j++ ) {
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

double64  FiniteDifferenceGrid::operator()( double64 x,  double64 y, bool rounded ) const
  {
      // 1. finding x and y indices of interpolation points
      if ( rounded ) {
           x1 = rint((x - static_cast<double64>(x_min)) / xresolution);
           y1 = rint((y - static_cast<double64>(y_min)) / yresolution);
        }
      else {
           x1 = static_cast<int32>((x-static_cast<double64>(x_min)) / xresolution);
           y1 = static_cast<int32>((y-static_cast<double64>(y_min)) / yresolution);
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
      t = ((x-static_cast<double64>(x_min)) - x1*xresolution) / xresolution;
      u = ((y-static_cast<double64>(y_min)) - y1*yresolution) / yresolution;
         
      // 6. bi-linear interpolation
      return (1.0-t)*(1.0-u)*p1 + t*(1.0-u)*p2 + t*u*p3 + (1.0-t)*u*p4;   
  }




double64 FiniteDifferenceGrid::InterpolateOutside( double64 x, double64 y ) const
{
    x1 = static_cast<int32>((x-x_min) / xresolution);
    y1 = static_cast<int32>((y-y_min) / yresolution);

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
    return std::numeric_limits<double64>::quiet_NaN();
    
} // end InterpolateOutside



/**

Interpolates strictly within the confines of the grid:
If y is smaller than the minimum of the y-range value is interpolated on the
horizontal model boundary x is set to the smallest x value. Then the desired value
is interpolated along the grid boundary between the marginal 2 points.
*/
double64 FiniteDifferenceGrid::InterpolateWithin( double64& x, double64& y ) const
{
    // 0. special corner cases
    if ( x < x_min && y < y_min ) return (*this)(0,0);
    if ( x < x_min && y > y_max ) return (*this)(size_y-1,0);
    if ( x > x_max && y < y_min ) return (*this)(0,size_x-1);
    if ( x > x_max && y > y_max ) return (*this)(size_y-1,size_x-1);
  
    x1 = static_cast<int32>((x-x_min) / xresolution);
    y1 = static_cast<int32>((y-y_min) / yresolution);

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
