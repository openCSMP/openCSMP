//
//  fracturePropertyModelling.h
//
//  Created by Stephan Matthai on 1/11/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#ifndef FRACTURE_PROPERTY_MODELLING_H
#define FRACTURE_PROPERTY_MODELLING_H

#include "CSMP_definitions.h"

namespace csmp {

class InSituStress;

/**  itf - ISF Tools: Flow-based upscaling via DFM models

Calculation of fracture aperture, permeability and porosity from in situ stress and mechanical properties of the rock


@todo  ! move stress tensor inside the bounding box of model so that one can see it.

*/
template<uint32_t dim>
void fracturePropertyModelling( const char* model_name );

/**
  Barton et al's 1985 model that estimates shear-displacement from shear stress and effective normal stress, using the mobilised version of JRC
*/
void fractureApertureFromShearAndNormalStress( Model<2>&, const InSituStress& );

} // end csmp

#endif
