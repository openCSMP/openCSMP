// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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

    FiniteDifferenceGrid( double x_dim, double y_dim, double xres, double yres, int32_t frame_width=0 );
    
    FiniteDifferenceGrid( double x_min, double x_max, double y_min, double y_max, double xres, double yres, int32_t frame_width=0 );
    
    void Initialize( double x_dim, double y_dim, double xres, double yres, int32_t frame_width=0 );
    
    void Initialize( double x_min, double x_max, double y_min, double y_max, double xres, double yres, int32_t frame_width=0 );

    FiniteDifferenceGrid& operator=( const FiniteDifferenceGrid& g );

    FiniteDifferenceGrid& operator=( double val );

    void  LinearInterpolateOnTo( FiniteDifferenceGrid& grid ) const;

    double    InterpolateOutside( double x, double y ) const;
    double    InterpolateWithin( double& x, double& y ) const;
    double    ExtrapolateTo( double x, double y ) const;
    double    operator()( double x, double y, bool rounded=false ) const;
    double&   operator()( int32_t row, int32_t col );
    double    operator()( int32_t row, int32_t col ) const;
    double    Value( int32_t row, int32_t col ) const; 
    double&   N( int32_t row, int32_t col );
    double&   S( int32_t i, int32_t j );
    double&   E( int32_t i, int32_t j );
    double&   W( int32_t i, int32_t j );
    double&   NW( int32_t i, int32_t j );
    double&   NE( int32_t i, int32_t j );
    double&   SW( int32_t i, int32_t j );
    double&   SE( int32_t i, int32_t j );
    double    ResolutionX() const;
    double    ResolutionY() const;
    int32_t   Rows()        const;
    int32_t   Columns()     const;
    bool      HasFrame()    const;
    double    MinX()        const;
    double    MinY()        const;
    double    MinZ()        const;
    double    MaxX()        const;
    double    MaxY()        const;
    double    MaxZ()        const;
    void  DataMinMax( double& dmin, double& dmax, bool incl_frame=false ) const;
    void  DataMinMaxWithoutNAN( double& dmin, double& dmax, bool incl_frame=false ) const;
    double    X( int32_t i ) const;
    double    Y( int32_t j ) const;
    void  ClosestGridPointTo( double x, double y, int32_t& index_row, int32_t& index_col ) const;
    double    RowAverage( int32_t row )    const;
    double    ColumnAverage( int32_t col ) const; 
    void  ScaleDataToRange( double dmin, double dmax );
    bool  IsInitialized() const;

  // setting grid values
    void  SetRowTo( int32_t row, double val );
    void  SetColumnTo( int32_t col, double val );
    void  SetFirstNRowsTo( int32_t n, double val );
    void  SetLastNRowsTo( int32_t n, double val );
    void  SetRegionTo( double xmin, double xmax, double ymin, double ymax, double val );
    void  Out() const;
    void  Out( const char* fname, int32_t tstep=0, bool with_frame=false ) const;
    bool  BinaryOut( const char* bin_name, int32_t tstep=0, bool with_frame=false ) const;
    bool  In( const char* fname );
    bool  BinaryIn( const char* fname ); 

#ifdef CSMP_WITH_IMAGE_OUTPUT
    void  SaveToPPM( const char* filename, int32_t timestep, bool greyscale=false ) const;
    void  SaveToJPG( const char* filename, int32_t timestep, bool greyscale=false, bool sqrt_of=false ) const;
    void  SaveToJPG( const char* filename, int32_t timestep, double vmin, double vmax, 
                     bool greyscale=false, 
                     bool sqrt_of=false ) const;
    void  SaveToJPGWithoutNAN( const char* filename, int32_t timestep, bool greyscale=false, bool sqrt_of=false ) const;
    void  SaveToJPGWithoutNAN( const char* filename, int32_t timestep, double vmin, double vmax, 
                               bool greyscale=false, 
                               bool sqrt_of=false ) const;
#endif
  
  private:
    std::vector<double>  grid;
    double   xresolution, yresolution;
    double   x_min, y_min, x_max, y_max, z_min, z_max;
    int32_t  size_x, size_y, size_z;
    int32_t  xfr, yfr, zfr, rowlength;

    // data only used by operator(double,double) and InterpolateOutside()
    mutable int32_t x1, y1;
    mutable double    delta, p1, p2, p3, p4, t, u;
    
    // (4 point, third order polynomial extrapolation)
    void   PolynomialInterpolation( double* xa, double* ya, double x, double& y, double& dy ) const;
 };

  

} // csmp


#endif




















