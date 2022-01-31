//
//  namedPropertyValuesToRegions.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 9/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_NAMED_PROPERTY_VALUES_TO_REGIONS_H
#define CSMP_NAMED_PROPERTY_VALUES_TO_REGIONS_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Model;

/// reads property-value - region-name mapping from the supplied text file (strings and IDs), checks and creates corresponding regions 
template<size_t dim>
void  namedPropertyValuesToRegions( Model<dim>&, const std::string& prop_name, const std::string& region_identifier_file );

/// to remove NO_DATA values converted to NAN (not a number), but only for scalar element property values.
template<size_t dim>
size_t replaceElement_NAN_ValuesWith( Model<dim>&, const std::string& element_var, double replacement_val );

} // end csmp

#endif /* CSMP_NAMED_PROPERTY_VALUES_TO_REGIONS_H */
