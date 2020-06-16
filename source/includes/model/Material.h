//
//  Material.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 13/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_MATERIAL_H
#define CSMP_MATERIAL_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Model;

/// transforms rational (-1,0..5..n) property values into material IDs stored on the elements
template<size_t dim>
void material_IDs_FromPropertyValues( Model<dim>&, const std::string& elmt_prop_name );

/// transforms rational (-1,0..5..n) property values into material IDs stored on the elements
template<size_t dim>
void propertyValuesPFromMaterial_IDs( Model<dim>&, const std::string& elmt_prop_name );

/// removes stair-steps in material boundaries where possible; stair-steps are identified by two or more element faces on the outside of a material domain
template<size_t dim>
void smoothMaterialInterfaces( Model<dim>& );

} // end csmp

#endif /* CSMP_MATERIAL_H */
