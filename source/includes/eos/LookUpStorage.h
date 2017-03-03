#ifndef LOOK_UP_STORAGE_H
#define LOOK_UP_STORAGE_H

#include "CSMP_definitions.h"

namespace csmp {

// 2D grid
template<typename fT>
class LookUpStorage {
  public:
    LookUpStorage();
    LookUpStorage( const LookUpStorage<fT>& );
    
    LookUpStorage( fT x_dim, fT y_dim, fT xres, fT yres );
    
    ~LookUpStorage();
    void Initialize( fT x_dim, fT y_dim, fT xres, fT yres );
    
    LookUpStorage& operator=( const LookUpStorage& g );
    fT&   operator()( int32 row, int32 col );
    int32 Rows() const;
    int32 Columns() const;

    // setting grid values
    bool  BinaryOut( const char* bin_name ) const;
    bool  BinaryIn( const char* fname ); 

  private:
    std::vector<fT> grid;
    fT              xresolution, yresolution;
    fT              x_max, y_max;
    int32           size_x, size_y, rowlength;

 };

template<typename fT>
inline fT&  LookUpStorage<fT>::operator()( int32 row, int32 col ) 
  {
      return grid[ static_cast<size_t>(row*rowlength + col) ];
  }
template<typename fT>
inline int32   LookUpStorage<fT>::Rows() const { return size_y-1; }
    
template<typename fT>
inline int32   LookUpStorage<fT>::Columns() const { return size_x-1; }



} // csmp


#endif




















