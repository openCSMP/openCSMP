#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

#include "CSMP_definitions.h"

namespace csmp
{

template<uint32_t> class Model;

template<uint32_t dim>
void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<dim>&,
                                                              const char* element_property,
                                                              const char* node_property,
                                                              const char* region,
                                                              double multiplication_factor = 1. );

template<uint32_t dim>
void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<dim>&,
                                                              const char* element_property,
                                                              const char* node_property,
                                                              const char* volume_modifier,
                                                              const char* region,
                                                              double multiplication_factor = 1. );

template<uint32_t dim>
void areaWeightedDistributionOfNodePropertyToNodeAndAdd( Model<dim>&,
                                                         const char* nodal_flux,
                                                         const char* nodal_value,
                                                         const char* region = "Model",
                                                         double multiplication_factor = 1.);

template<uint32_t dim>
void setPropertyToZero( Model<dim>&, const char* property );

} // end csmp

#endif
