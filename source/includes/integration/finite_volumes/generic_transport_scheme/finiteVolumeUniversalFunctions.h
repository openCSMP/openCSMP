//
//  finiteVolumeUniversalFunctions.h
//
//  Created by Stephan Matthai on 9/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#ifndef FINITE_VOLUME_UNIVERSAL_FUNCTIONS_H
#define FINITE_VOLUME_UNIVERSAL_FUNCTIONS_H

#include <iostream>

namespace csmp {

template<size_t> class Model;
template<size_t> class Region;

/// sector volume, finite volume, FV pore volume
template<size_t dim> void initializeFiniteVolumeProperties( Model<dim>&, Region<dim>&  );

} // csmp

#endif /* defined(FINITE_VOLUME_UNIVERSAL_FUNCTIONS_H) */
