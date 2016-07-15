#ifndef PRESSURE_SATURATION_INITIALIZER_VISITOR_H
#define PRESSURE_SATURATION_INITIALIZER_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Model;

/**
@author Sh. Bazrafkan
@date 2010
*/

/// Initialize the pressure and saturation for a model
template<size_t dim>
class PressureSaturationInitializer : public Visitor<dim> {
  
public:
      PressureSaturationInitializer(Model<dim>& model,
                                    const char* water_saturation,
                                    const char* oil_saturation,
                                    const char* water_pressure,
                                    const char* oil_pressure,
                                    double64 water_density,
                                    double64 oil_density,
                                    double64 water_oil_contact,
                                    double64 reference_depth,
                                    double64 reference_pressure,
                                    double64 lambda,
                                    double64 entry_pressure,
                                    double64 swc,
                                    double64 sor);
                 
      virtual ~PressureSaturationInitializer();
      
      double64 FreeWaterLevel();
      double64 PressureAtFreeWaterLevel();
      double64 WaterSaturationFromPhasePressures(double64 oil_pressure, 
                                                 double64 water_pressure);
    
      virtual void Visit(Node<dim>* node);
      virtual void Visit(Model<dim>* model);
    
private:
      const PropertyDatabase<dim>&  prop_ref_;
      Model< dim>&               model_ref_;
      csmp::Index                water_pressure_key_,
                                 oil_pressure_key_,
                                 water_saturation_key_,
                                 oil_saturation_key_;
      double64                   water_density_,
                                 oil_density_,
                                 water_oil_contact_,
                                 reference_depth_,
                                 reference_pressure_,
                                 lambda_,
                                 entry_pressure_,
                                 swc_,
                                 sor_;
      const double64             gravity_acceleration_;
                               
};


} // end namespace csmp

#endif
