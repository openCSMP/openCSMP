#ifndef PARALLEL_PLATE_PERMEABILITY_FROM_CHANNEL_WIDTH_H
#define PARALLEL_PLATE_PERMEABILITY_FROM_CHANNEL_WIDTH_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;

/**

Computes the parallel-plate permeability and the channel width in the
channel region. The aperture amounts to the maximum width of the channel 
at any particular segment. The results are stored in the variables 
'permeability' (m2) and 'channel width' (m).  

The calculated permeability is    k = a^2/12
where a is the channel width.

*/
template<uint32_t dim>
void parallelPlatePermeabilityFromChannelWidth( Model<dim>&, 
                                                const char* channel_region,
                                                const char* channel_width="channel width",
                                                double minimum_channel_width=1.0e-9 /* 10 Angstroem */ );
} // end namespace csmp

#endif
