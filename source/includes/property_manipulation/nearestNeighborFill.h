//
//  nearestNeighborFill.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 15/4/2022.
//

#ifndef CSMP_NEAREST_NEIGHBOR_FILL
#define CSMP_NEAREST_NEIGHBOR_FILL

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;

/// replaces no-data values of target variable with nearest-neighbor values until there are none left, by default NAN's are no-data values
template<uint32_t dim>
void nearestNeighborFill( Model<dim>&, const char* target_region, const char* variable, double no_data_value );

}

#endif /* CSMP_NEAREST_NEIGHBOR_FILL */
