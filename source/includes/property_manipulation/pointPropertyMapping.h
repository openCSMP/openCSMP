// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only
//
//  pointPropertyMapping.h
//
//

#ifndef POINT_PROPERTY_MAPPING_H
#define POINT_PROPERTY_MAPPING_H

#include <iostream>

namespace csmp {

/** 
    Maps point properties from point clouds tagged with region names
    to current 3D CSMP model. 
    
    Properties are mapped either to volumes and / or surfaces,
    but the latter must be assigned last to avoid overwriting
    of pre-assigned values.
    
    After successful property assignment, the model is output to VTU for examination.
    If appropriate the results can be safed to current CSMP binary.
    
    The mapping involves both, volume/surface averaging, when multiple points
    fall into the same element, interpolation using the SAMG solver when the data are sparse,
    or nearest neighbor interpolation if where there are isolated element without property values.
    
    See  Mosser (2013, SPE 167629\00\00\00\00-STU) for a full description of the implemented algorithm.
    
    @todo TODO: implement averaging of the log10 of property value where these vary dramatically like permeability etc.
    
      @author Lukas Mosser 
      @author Stephan Matthai
    
*/

void pointPropertyMapping( const std::string& model_name );

} // end csmp

#endif /* POINT_PROPERTY_MAPPING_H */
