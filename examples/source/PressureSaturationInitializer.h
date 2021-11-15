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
                                    double water_density,
                                    double oil_density,
                                    double water_oil_contact,
                                    double reference_depth,
                                    double reference_pressure,
                                    double lambda,
                                    double entry_pressure,
                                    double swc,
                                    double sor);
                 
      virtual ~PressureSaturationInitializer();
      
      double FreeWaterLevel();
      double PressureAtFreeWaterLevel();
      double WaterSaturationFromPhasePressures(double oil_pressure, 
                                                 double water_pressure);
    
      virtual void Visit( Model<dim>* m ) const 
        { std::cout <<"\nPressureSaturationInitializer:Visit(Model): "<< m->Name() <<"\n"; } 
            
      virtual void Visit( Node<dim>* );    
      // for all other targets the method stubs in the base class are used
      // this may will prompt some warnings  
    
private:
      const PropertyDatabase<dim>&  prop_ref_;
      Model< dim>&               model_ref_;
      csmp::Index                water_pressure_key_,
                                 oil_pressure_key_,
                                 water_saturation_key_,
                                 oil_saturation_key_;
      double                   water_density_,
                                 oil_density_,
                                 water_oil_contact_,
                                 reference_depth_,
                                 reference_pressure_,
                                 lambda_,
                                 entry_pressure_,
                                 swc_,
                                 sor_;
      const double             gravity_acceleration_;
                               
};


} // end namespace csmp

#endif
