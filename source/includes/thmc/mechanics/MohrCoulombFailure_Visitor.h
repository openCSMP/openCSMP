#ifndef MOHR_COULOMB_FAILURE_VISITOR_H
#define MOHR_COULOMB_FAILURE_VISITOR_H

#include "Visitor.h"
#include "TensorVariable.h"
#include "ScalarVariable.h"

namespace csmp {

class ScalarVariable;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

/** 
     Evaluates Mohr-Coulomb shear and tensile failure criteria that are output 
     to the element integration points.
     
     Input variables:
     - friction angle 
     - cohesion (scalar)
     - stress (tensor) = Cartesian stress
     - fluid pressure (if this variable is not specified its value is assumed to be zero)
     
     @todo SKM: this visitor should include the influence of pore pressure.
 
*/
template<uint32_t dim>
class MohrCoulombFailure_Visitor : public Visitor<dim> {
  public:
    explicit MohrCoulombFailure_Visitor( Model<dim>&,
                                         bool positive_compressive_stress_convention = true,
                                         bool verbose = false );

    ~MohrCoulombFailure_Visitor() = default;
    
    void Visit(Element<dim>* ) override final;

  private:
    PropertyDatabase<dim>&  stressref_;
    bool                    verbose_;
    csmp::Index             Stress_key_,      ///< stress
                            pf_key_,          ///< fluid pressure
                            Cohesion_key_,    ///< cohesion
                            Friction_key_,    ///< friction angle 0..90o
                            Failure_key_,     ///< failure
                            Failure01_key_;   ///< failure 0 1

    TensorVariable<dim>     Cartesian_stress_;

    const double degrees_to_radians_;
    const double sqrt3_;
    const double sqrt32_;
    double       fluid_pressure_;
    double       biot_coefficient_alpha_;
    double       sign_of_tensile_stress_;
};

} // end csmp

#endif
