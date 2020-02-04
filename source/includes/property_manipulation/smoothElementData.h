//
//  smoothData.hpp
//  CSMP_FECFVM_Simulator
//
//  Created by Stephan Matthai on 4/2/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_SMOOTH_ELEMENT_DATA_H
#define CSMP_SMOOTH_ELEMENT_DATA_H

#include "CSMP_definitions.h"

namespace csmp {

enum SMOOTHING_TYPE { MEAN, HARMONIC_MEAN, LOG10_MEAN, HORIZONTAL_MEAN };

/**
    smoothElementData - reads CSMP binary file  and performs various kinds of smoothing on the target element variable.
    
    @parameter smoothing_passes tells functions how often the same operation is to be performed
    
    @parameter log10smoothing  takes the decadic log value of the properties before the smoothing is performed.
    
    @parameter horizontal_direction_only restricts operation to the horizontal
*/
void smoothElementData( const std::string& model, 
                        const std::string& region_to_be_smoothed, 
                        const std::string& variable_name,
                        int smoothing_passes = 1, 
                        bool log10_smoothing = false, bool in_plane_smoothing = false );

} // end csmp

#endif /* CSMP_SMOOTH_ELEMENT_DATA_H */
