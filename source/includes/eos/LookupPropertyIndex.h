#ifndef LOOKUP_PROPERTY_INDEX_H
#define LOOKUP_PROPERTY_INDEX_H

namespace csmp
{
  // always insert new property indices before "max_index"
  enum lookup_property_index{temperature_index, 
       pressure_index, 
       composition_index, 
       density_index, 
       enthalpy_index,
       heatcapacity_index,
       compressibility_index,
       viscosity_index,
       max_index };
}
#endif
