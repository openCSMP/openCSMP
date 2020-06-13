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

/// transforms rational (-1,0..5..n) property values into material IDs stored on the elements
template<size_t dim>
void material_IDs_FromPropertyValues( Model<dim>& model, const std::string& elmt_prop_name );

/// reads property-value - region-name mapping from the supplied text file (strings and IDs), checks and creates corresponding regions 
template<size_t dim>
void  namedPropertyValuesToRegions( Model<dim>& model, const std::string& prop_name, const std::string& region_identifier_file );

/// to remove NO_DATA values converted to NAN (not a number).
template<size_t dim>
void replaceElement_NAN_ValuesWith( Model<dim>& model, const char* element_var, double64 replacement_val );



/// reads material ID - region-name mapping from the supplied text file (strings and IDs), checks and creates corresponding regions 
// TODO: add this method to Material_ID_Manager or class with similar name
//template<size_t dim>
//void  createRegionsFromMaterialIDs( Model<dim>& model, const std::string& property_name_file );


} // end csmp

#endif /* CSMP_NAMED_PROPERTY_VALUES_TO_REGIONS_H */
