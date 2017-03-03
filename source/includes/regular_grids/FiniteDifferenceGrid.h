#ifndef CSMP_FINITE_DIFFERENCE_GRID_H
#define CSMP_FINITE_DIFFERENCE_GRID_H

#include "CSMP_definitions.h"
#include "CSMP_mathUtilities.h"

namespace csmp {

/**
 
@brief Two-dimensional regular grid with interpolation operations and an optional halo
across which interpolation still works; can be used for grid-point in finite-element cell transport computations.

@author S.K. Matthai
@date 28/7/97

FiniteDifferenceGrid class implements a regular grid i {0,i_max-1}
corresponding to the horizontal coordinate direction
and j {0,j_max-1} corresponding to the vertical coordinate direction.
Thus i=columns, j=rows.
an operator(x,y) is provided to interpolate values within the
grid domain.
X(), Y() give x and y values corresponding to column and row indices
of the regular grid, respectively.
Val(i,j) is used to input and output values at the exact grid locations
indicated by i(column=x) and j(row=y).

*/
class FiniteDifferenceGrid {
  public:
    FiniteDifferenceGrid();
    FiniteDifferenceGrid( const FiniteDifferenceGrid& g );
    
    FiniteDifferenceGrid( double64 x_dim, double64 y_dim, double64 xres, double64 yres, int32 frame_width=0 );
    
    FiniteDifferenceGrid( double64 x_min, double64 x_max, double64 y_min, double64 y_max, double64 xres, double64 yres, int32 frame_width=0 );
    
    ~FiniteDifferenceGrid();
    void Initialize( double64 x_dim, double64 y_dim, double64 xres, double64 yres, int32 frame_width=0 );
    
    void Initialize( double64 x_min, double64 x_max, double64 y_min, double64 y_max, double64 xres, double64 yres, int32 frame_width=0 );
    
    FiniteDifferenceGrid& operator=( const FiniteDifferenceGrid& g );
    FiniteDifferenceGrid& operator=( double64 val );
    void  LinearInterpolateOnTo( FiniteDifferenceGrid& grid ) const;

    double64    InterpolateOutside( double64 x, double64 y ) const;
    double64    InterpolateWithin( double64& x, double64& y ) const;
    double64    ExtrapolateTo( double64 x, double64 y ) const;
    double64    operator()( double64 x, double64 y, bool rounded=false ) const;
    double64&   operator()( int32 row, int32 col );
    double64    operator()( int32 row, int32 col ) const;
    double64    Value( int32 row, int32 col ) const; 
    double64&   N( int32 row, int32 col );
    double64&   S( int32 i, int32 j );
    double64&   E( int32 i, int32 j );
    double64&   W( int32 i, int32 j );
    double64&   NW( int32 i, int32 j );
    double64&   NE( int32 i, int32 j );
    double64&   SW( int32 i, int32 j );
    double64&   SE( int32 i, int32 j );
    double64    ResolutionX() const;
    double64    ResolutionY() const;
    int32 Rows()        const;
    int32 Columns()     const;
    bool  HasFrame()    const;
    double64    MinX()        const;
    double64    MinY()        const;
    double64    MinZ()        const;
    double64    MaxX()        const;
    double64    MaxY()        const;
    double64    MaxZ()        const;
    void  DataMinMax( double64& dmin, double64& dmax, bool incl_frame=false ) const;
    void  DataMinMaxWithoutNAN( double64& dmin, double64& dmax, bool incl_frame=false ) const;
    double64    X( int32 i ) const;
    double64    Y( int32 j ) const;
    void  ClosestGridPointTo( double64 x, double64 y, int32& index_row, int32& index_col ) const;
    double64    RowAverage( int32 row )    const;
    double64    ColumnAverage( int32 col ) const; 
    void  ScaleDataToRange( double64 dmin, double64 dmax );
    bool  IsInitialized() const;

  // setting grid values
    void  SetRowTo( int32 row, double64 val );
    void  SetColumnTo( int32 col, double64 val );
    void  SetFirstNRowsTo( int32 n, double64 val );
    void  SetLastNRowsTo( int32 n, double64 val );
    void  SetRegionTo( double64 xmin, double64 xmax, double64 ymin, double64 ymax, double64 val );
    void Out() const { Out(std::cout); }
    void  Out(std::ostream& os) const;
    void  Out( const char* fname, int32 tstep=0, bool with_frame=false ) const;
    bool  BinaryOut( const char* bin_name, int32 tstep=0, bool with_frame=false ) const;
    bool  In( const char* fname );
    bool  BinaryIn( const char* fname ); 

#ifdef CSMP_WITH_IMAGE_OUTPUT
    void  SaveToPPM( const char* filename, int32 timestep, bool greyscale=false ) const;
    void  SaveToJPG( const char* filename, int32 timestep, bool greyscale=false, bool sqrt_of=false ) const;
    void  SaveToJPG( const char* filename, int32 timestep, double64 vmin, double64 vmax, 
                     bool greyscale=false, 
                     bool sqrt_of=false ) const;
    void  SaveToJPGWithoutNAN( const char* filename, int32 timestep, bool greyscale=false, bool sqrt_of=false ) const;
    void  SaveToJPGWithoutNAN( const char* filename, int32 timestep, double64 vmin, double64 vmax, 
                               bool greyscale=false, 
                               bool sqrt_of=false ) const;
#endif
  
  private:
    double64*    grid;
    double64     xresolution, yresolution;
    double64     x_min, y_min, x_max, y_max, z_min, z_max;
    int32  size_x, size_y, size_z;
    int32  xfr, yfr, zfr, rowlength;

    // data only used by operator(double64,double64) and InterpolateOutside()
    mutable int32 x1, y1;
    mutable double64    delta, p1, p2, p3, p4, t, u;
    
    // (4 point, third order polynomial extrapolation)
    void   PolynomialInterpolation( double64* xa, double64* ya, double64 x, double64& y, double64& dy ) const;
 };

  

} // csmp


#endif




















