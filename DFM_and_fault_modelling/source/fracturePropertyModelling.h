//
//  fracturePropertyModelling.h
//
//  Created by Stephan Matthai on 1/11/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#ifndef FRACTURE_PROPERTY_MODELLING_H
#define FRACTURE_PROPERTY_MODELLING_H

#include <iostream>

namespace csmp {

/**  itf - ISF Tools: Flow-based upscaling via DFM models

Calculation of fracture aperture, permeability and porosity from in situ stress and mechanical properties of the rock


@todo  ! move stress tensor inside the bounding box of model so that one can see it.

*/
template<uint32_t dim>
void fracturePropertyModelling( const char* model_name );

}

#endif
