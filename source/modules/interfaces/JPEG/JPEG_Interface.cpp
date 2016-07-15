#include "JPEG_Interface.h"
#include "Region.h"
#include "Model.h"
#include "FemToGridVisitor.h"
#include "FemFromGridVisitor.h"

using namespace std;

#ifdef CSMP_WITH_IMAGE_OUTPUT

namespace csmp {

// default constructor
JPEG_Interface::JPEG_Interface()
 {
    cout <<"\nJPEG_Interface(default constructor): Building JPEG file output facility."<< endl;
 }


// destructor
JPEG_Interface::~JPEG_Interface()
 {
    cout <<"\nJPEG_Interface(destructor): Removing JPEG file output facility."<< endl;
 }



/**
 
OutputDataToJPG() uses a public domain implementation of the JPEG 
library to output the values of a property to a JPEG, 256 color or 
greyscale image. As a default the colorscheme will be rainbow. As 
options, the file can be output in greyscale and / or the the square root
of the variable values can be displayed in order to have a better visual
resolution of variables whose values vary over many orders of magnitude. 

@section arguments Input Arguments 

The first argument specifies the name of the output file to which the 
extension '.jpg' will be appended automatically. The second argument 
specifies the name of the output variable, and the third argument gives
the output timestep which will also be appended to the output
filename before the extension. The last two arguments are optional and
allow the user to export as greyscale image (a third of the size of the
RGB image), and as an image which visualizes the sqare root of the 
property values, respectively. The default is a rainbow coloring of the
actual variable values. 

@return If OutputDataToJPG() successfully creates the JPEG file it will return
the boolean variable 'true', else it returns 'false'.

@section implementation Implementation

The method first maps the distributed variable onto a FiniteDifferenceGrid
in order to obtain regular-gridded data for the image. The 
minimum 'inner radius' of the finite elements in the current mesh is used
to define the resolution of the grid. 'inner radius' must therfore be defined
as a scalar physical element variable. The grid is 
declared as a static variable, such that it does not need to be 
re-allocated if the method is used multiple times. The method then calls 
the SaveToJPG() interface of the FiniteDifferenceGrid object which uses
a ColorPalette object to convert the variable values to a 256-RGB rainbow
palette color dataset. 

@note OutputDataToJPG() only works for two-dimensional models.

@section application Application

OutputDataToJPG() is an attractive option to directly obtain image data 
from a run. Also, the resulting JPEG images, are very small even for
spatially highly-resolved models and their output is accompanied by output 
to stdout of the original values of the output variable. Thus, the JPEG files
can be used to conviniently monitor a variable over the duration of a run. 
Also, the JPEG images can be strung together later into a movie, using 
one of the standard JPEG to MPEG converters. 

@section messages Messages 

If the model is not two-dimensional, the variable 'inner radius' is not
defined or the variable is a not a scalar, an error is reported and the 
variable is not output to file. 

 */
bool JPEG_Interface::OutputDataToJPG( Model<2U>& sg,
                                      const char* file_name, 
                                      const char* var_name, 
                                      long        timestep, 
                                      bool        gray,  
                                      bool        sqrt_of_value )
 {
    csmp::Index  a_key = sg.Database().StorageKey(var_name);

    // 1. testing variable A for suitability
    // -------------------------------------
    if ( !sg.Database().IsDefined("inner radius") )
      {
         cout << "Model<dim>::OutputDataToJPG: "<< endl;
         cout << "Element property 'inner radius' does not exist but is needed by this method. Nothing was done..."<< endl;
         return false;
      }
    if ( a_key.type != SCALAR )  
      {
         cout << "Model<dim>::OutputDataToJPG: "<< endl;
         cout <<"This method outputs only scalar properties. Use Interrelation to compute scalar."<< endl;
         return false;
      }

    // 2. building static FiniteDifferenceGrid for repeated data output
    // ----------------------------------------------------------------
    double64  rmin, rmax;
    
    if ( !regular_grid.IsInitialized() ) 
      {
         sg.AssignElementCharacteristicsTo("inner radius", "inner radius");
         sg.MinMaxOf("inner radius", rmin, rmax );       
         Point<2U>  xyz_min, xyz_max;
         sg.MinMaxCoordinates( xyz_min, xyz_max );
         //                       xmin        xmax        ymin        ymax       
         regular_grid.Initialize( xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], rmin, rmin );  
      }

    // 3. setting up static visitor to write repeatedly write data to grid
    // -------------------------------------------------------------------
    FemToGridVisitor<2U> writer( sg.Database(),  regular_grid, 
                                 var_name, sg.Region("Model").Elements() );
    writer.OverWrite( true );
    writer.OutputProperty( var_name );
     
    // 4. Writing data onto grid 
    // -------------------------
    sg.Accept( writer );
    
    // 5. Writing grid to JPG file
    // ---------------------------
    //                                      grayscale output 
    regular_grid.SaveToJPG( file_name, timestep, gray, sqrt_of_value );

    return true;

 } // end OutputDataToJPG


/**
 
As method above, but the user can specify the value range to which the
color spectrum will be mapped to. If the actual values are outside of the
user-specified range, the method will adopt to the range given by the 
values to prevent missing a part of the spectrum.  

@section arguments Input Arguments 

As above, except for 2 double64 arguments permitting the user to 
specify the output range.  
*/
bool JPEG_Interface::OutputDataToJPG( Model<2U>& sg,
                                          const char* file_name, 
                                          const char* var_name, 
                                          long   timestep, 
                                          double64    data_min, 
                                          double64    data_max,
                                          bool        gray,  
                                          bool        sqrt_of_value )
 {
    csmp::Index  a_key = sg.Database().StorageKey(var_name);

    // 1. testing variable A for suitability
    // -------------------------------------
    if ( !sg.Database().IsDefined("inner radius") )
      {
         cout << "Model<dim>::OutputDataToJPG: "<< endl;
         cout << "Element property 'inner radius' does not exist but is needed by this method. Nothing was done..."<< endl;
         return false;
      }
    if ( a_key.type != SCALAR )  
      {
         cout << "Model<dim>::OutputDataToJPG: "<< endl;
         cout <<"This method outputs only scalar properties. Use Interrelation to compute scalar."<< endl;
         return false;
      }

    // 2. building static FiniteDifferenceGrid for repeated data output
    // ----------------------------------------------------------------
    double64  rmin, rmax;
    
    if ( !regular_grid.IsInitialized() ) 
      {
         sg.AssignElementCharacteristicsTo("inner radius", "inner radius");
         sg.MinMaxOf("inner radius", rmin, rmax );       
         Point<2U>  xyz_min, xyz_max;
         sg.MinMaxCoordinates( xyz_min, xyz_max );
         //                       xmin        xmax        ymin        ymax       
         regular_grid.Initialize( xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], rmin, rmin );  
      }

    // 3. setting up static visitor to write repeatedly write data to grid
    // -------------------------------------------------------------------
    FemToGridVisitor<2U> writer( sg.Database(), regular_grid, 
                                 var_name, sg.Region("Model").Elements() );
    writer.OverWrite( true );
    writer.OutputProperty( var_name );
     
    // 4. Writing data onto grid 
    // -------------------------
    sg.Accept( writer );
    
    // 5. Writing grid to JPG file
    // ---------------------------
    regular_grid.SaveToJPG( file_name, timestep, data_min, data_max, gray, sqrt_of_value );

    return true;

 } // end OutputDataToJPG (version with fixed value range)

} // csmp

#endif //  CSMP_WITH_IMAGE_OUTPUT




