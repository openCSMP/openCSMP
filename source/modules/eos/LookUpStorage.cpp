#include "LookUpStorage.h"
#include "binaryReadWrite.h"

using namespace std;

namespace csmp {

template<typename fT>
LookUpStorage<fT>::LookUpStorage()
  : grid(0), 
    xresolution(1.0), yresolution(1.0),
    x_max(0.0), y_max(0.0),
    size_x(0), size_y(0), rowlength(0)
{ 
//   cout <<"\nLookUpStorage: called default constructor..."<< endl;
}



template<typename fT>
LookUpStorage<fT>::LookUpStorage( fT x_dim, fT y_dim, 
                                  fT xres, fT yres )
  : grid(0)
 {
    Initialize( x_dim, y_dim, xres, yres );
 }






template<typename fT>
LookUpStorage<fT>::LookUpStorage( const LookUpStorage<fT>& g )
 : grid(0)
 {
    *this = g;
 }




template<typename fT>
LookUpStorage<fT>::~LookUpStorage()
 {
 }



template<typename fT>
void LookUpStorage<fT>::Initialize( fT x_dmax, fT y_dmax,  fT xres, fT yres )
 {
     // 0. removing old grid
     x_max       = x_dmax;
     y_max       = y_dmax; 
     xresolution = xres;
     yresolution = yres;

     // 1. getting the resolution
     size_x = static_cast<int32>(x_max / xresolution + 1.);
     size_y = static_cast<int32>(y_max / yresolution + 1.);

     // 2. extra data member for speed
     rowlength = size_x;
     
     // 3. set vector size
     grid.resize(static_cast<size_t>(size_x * size_y));
    
     // 4. output memory requirements
     cout <<"\nLookUpStorage: building new grid; allocating ";
     cout << ((size_x * size_y * sizeof(fT)) / 1.0e+6);
     cout <<" MByte of memory..." << endl;
     
     
 }



template<typename fT>
LookUpStorage<fT>& LookUpStorage<fT>::operator=( const LookUpStorage<fT>& g )
 {
    if ( &g == this ) return *this;
    xresolution = g.xresolution;
    yresolution = g.yresolution;
    size_x      = g.size_x;
    size_y      = g.size_y;
    rowlength   = g.rowlength;
    x_max       = g.x_max;
    y_max       = g.y_max;
    grid        = g.grid;

    return *this;    
 }
 
 


template<typename fT>
bool   LookUpStorage<fT>::BinaryOut( const char* bin_name ) const
 {
    char name[200], heading[200];
    strcpy( name, bin_name );
    strcpy( heading, "LookUpStorage<fT>::BinaryOut: fT grid as binary file");
 
    fstream fp (name, ios::in | ios::binary);
    if ( !fp.is_open() ) {
        cout <<"\nLookUpStorage<fT>::BinaryOut: File: "<< bin_name << " could not be opened"<< endl;
        return false;
      }
    skm_C_fwrite( fp, heading );  

    // stores dimensions of grid
    std::vector<fT>    dim_fT(4);
    
    // overall dimensions
    dim_fT[0] = x_max;
    dim_fT[1] = y_max;
    // resolution x and y 
    dim_fT[2] = xresolution;
    dim_fT[3] = yresolution; 
    
    skm_C_fwrite( fp, dim_fT ); 
    
    skm_C_fwrite( fp, grid ); 
    
    fp.close();

   cout <<"\n\n'" << bin_name <<"' written successfully..." << endl;
   
   return true;
 }


template<typename fT>
bool LookUpStorage<fT>::BinaryIn( const char* bin_name )
 {
	 fstream fp(bin_name, ios::in | ios::binary);
    if ( !fp.is_open() ) {
        cout <<"\nLookUpStorage<fT>::BinaryIn: File: "<< bin_name;
        cout <<" could not be opened"<< endl;
        return false;
     }
    char heading[200];
    skm_C_fread( fp, heading ); 
    cout <<"\nLookUpStorage<fT>::BinaryIn: Reading: "<< heading << endl;

    // read dimensions of grid
    std::vector<fT>    dim_fT;
    // read data
    std::vector<fT>    input_data;

    // read info about dimension
    skm_C_fread( fp, dim_fT ); 

    // rebuild grid
    Initialize( dim_fT[0], dim_fT[1], dim_fT[2], dim_fT[3] );
    
    // read data and transfer to grid
    skm_C_fread( fp, input_data ); 
    grid = input_data;

    fp.close();
    
    cout <<"\nLookUpStorage<double64>::BinaryIn: grid build successfully from binary file." << endl;
    cout.flush();

    return true;
 }



template class LookUpStorage<double64>;


} // end namespace csmp
