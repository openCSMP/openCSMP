#ifndef CSMP_JPEG_INTERFACE_H
#define CSMP_JPEG_INTERFACE_H

#include "FiniteDifferenceGrid.h"

namespace csmp {

template<size_t> class Model;

class JPEG_Interface {
  public:
    JPEG_Interface();
    ~JPEG_Interface();
  
    bool OutputDataToJPG ( Model<2U>&, 
                           const char* file_name, const char* var_name, 
                           long timestep, bool gray=false, bool sqrt_of_value=false );
                           
    bool OutputDataToJPG ( Model<2U>&,
                           const char* file_name, const char* var_name, long timestep, 
                           double64 data_min, double64 data_max, bool gray=false, 
                           bool sqrt_of_value=false );  
  
  private:
    FiniteDifferenceGrid  regular_grid;
};

} // csmp

#endif
