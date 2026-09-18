// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PressureSaturationInitializer.h"
#include "Model.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "ErrorHandler.h"

namespace csmp {



template<uint32_t dim>
PressureSaturationInitializer<dim>::PressureSaturationInitializer( Model<dim>& model,
                                                                   const char* water_saturation,
                                                                   const char* oil_saturation,
                                                                   const char* water_pressure,
                                                                   const char* oil_pressure,
                                                                   double water_density,
                                                                   double oil_density,
                                                                   double water_oil_contact,
                                                                   double reference_depth,
                                                                   double reference_pressure,
                                                                   double lambda,
                                                                   double entry_pressure,
                                                                   double swc,
                                                                   double sor)
 : Visitor<dim>(MODEL, NODE),
   model_ref_(model),
   prop_ref_(model.Database()),
   water_pressure_key_(model.Database().StorageKey(water_pressure)),
   oil_pressure_key_(model.Database().StorageKey(oil_pressure)),
   water_saturation_key_(model.Database().StorageKey(water_saturation)),
   oil_saturation_key_(model.Database().StorageKey(oil_saturation)),
   water_density_(water_density),
   oil_density_(oil_density),
   water_oil_contact_(water_oil_contact),
   reference_depth_(reference_depth),
   reference_pressure_(reference_pressure),
   lambda_(lambda),
   entry_pressure_(entry_pressure),
   swc_(swc),
   sor_(sor),
   gravity_acceleration_(9.80665)
{

    if (water_pressure_key_.place != NODE || water_pressure_key_.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     water_pressure ," must be a scalar node property");

    if (oil_pressure_key_.place != NODE || oil_pressure_key_.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     oil_pressure ," must be a scalar node property");

    if (water_saturation_key_.place != NODE || water_saturation_key_.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     water_saturation ," must be a scalar node property");

    if (oil_saturation_key_.place != NODE || oil_saturation_key_.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     oil_saturation ," must be a scalar node property");

    if ( water_density <= oil_density )
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     "water density must be greater than oil density");

    if ( lambda < 0. )
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     "Brooks-Corey parameter must be greater than or equal to zero");

    if ( entry_pressure < 0. )
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     "entry pressure must be greater than or equal to zero");

    if ( swc + sor >= 1. )
        throw csmp::Exception(FATAL_ERROR, "PressureSaturationInitializer<dim>::(constructor) ",
                     "there is no mobile oil in the model");

} // end PressureSaturationInitializer




template<uint32_t dim>
PressureSaturationInitializer<dim>::~PressureSaturationInitializer()
{
} // end ~PressureSaturationInitializer




template<uint32_t dim>
double PressureSaturationInitializer<dim>::FreeWaterLevel()
{
    return ( water_oil_contact_ - entry_pressure_ /
           ( water_density_ - oil_density_ ) / gravity_acceleration_ );

} // end FreeWaterLevel




template<uint32_t dim>
double PressureSaturationInitializer<dim>::PressureAtFreeWaterLevel()
{

    double free_water_level = FreeWaterLevel();

    // if reference depth is higher than the water oil contact then the reference pressure is oil pressure
    // and we use oil gravity to calculate pressure at free water level
    if ( reference_depth_ > water_oil_contact_ )
        return ( reference_pressure_ + ( reference_depth_ - free_water_level ) *
                 gravity_acceleration_ * oil_density_ );
    // if reference depth is lower than the water oil contact then the reference pressure is water pressure
    // and we use water gravity to calculate pressure at free water level
    else
        return ( reference_pressure_ + ( reference_depth_ - free_water_level ) *
                 gravity_acceleration_ * water_density_ );

} // end PressureAtFreeWaterLevel




template<uint32_t dim>
double PressureSaturationInitializer<dim>::WaterSaturationFromPhasePressures(double oil_pressure,
                                                                               double water_pressure)
{
    if ( ( oil_pressure - water_pressure ) <= entry_pressure_ )
        return 1.;
    else
        return ( swc_ + ( 1. - swc_ + sor_ ) * std::pow( entry_pressure_ / ( oil_pressure - water_pressure ),
                                                         lambda_ ) );
} // end WaterSaturationFromPhasePressures



template<uint32_t dim>
void PressureSaturationInitializer<dim>::Visit(Node<dim>* node)
{

    double node_depth;
    if ( dim == 1U )
        node_depth = node->x();
    else
        node_depth = node->y();
    double free_water_level = FreeWaterLevel();
    double free_water_level_pressure = PressureAtFreeWaterLevel();
    ScalarVariable oil_pressure( PLAIN, free_water_level_pressure + oil_density_ * gravity_acceleration_ *
                                        ( free_water_level - node_depth ) );
    ScalarVariable water_pressure( PLAIN, free_water_level_pressure + water_density_ * gravity_acceleration_ *
                                          ( free_water_level - node_depth ) );
    ScalarVariable water_saturation( PLAIN, WaterSaturationFromPhasePressures(oil_pressure(), water_pressure()) );
    ScalarVariable oil_saturation( PLAIN, 1. - water_saturation() );

    node->Store(oil_pressure_key_, oil_pressure);
    node->Store(water_pressure_key_, water_pressure);
    node->Store(oil_saturation_key_, oil_saturation);
    node->Store(water_saturation_key_, water_saturation);

} // end Visit

template class PressureSaturationInitializer<1U>;
template class PressureSaturationInitializer<2U>;
template class PressureSaturationInitializer<3U>;

} // end namespace csmp
