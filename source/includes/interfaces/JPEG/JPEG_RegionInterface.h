#ifndef CSMP_JPEG_REGION_INTERFACE_H
#define CSMP_JPEG_REGION_INTERFACE_H

#include "FiniteDifferenceGrid.h"

namespace csmp {

template<size_t> class Model;

class JPEG_RegionInterface {
  public:
    JPEG_RegionInterface( Model<2U>&, const char* group_name );
    ~JPEG_RegionInterface();
  
    bool OutputRegionDataToJPG ( Model<2U>&, 
                                 const char* file_name, const char* var_name, 
                                 long timestep, bool gray=false, bool sqrt_of_value=false );
                           
    bool OutputRegionDataToJPG ( Model<2U>&,
                                 const char* file_name, const char* var_name, long timestep, 
                                 double64 data_min, double64 data_max, bool gray=false, 
                                 bool sqrt_of_value=false );  
  
  private:
    FiniteDifferenceGrid  regular_grid;
    std::string           name_;
};

}

#endif
