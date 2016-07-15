#ifndef BRITTLE_FAILURE_VISITOR_H
#define BRITTLE_FAILURE_VISITOR_H

#include "Visitor.h"
#include "TensorVariable.h"
#include "ScalarVariable.h"
#include "MechanicalProperties.h"

namespace csmp {

class ScalarVariable;
template<size_t> class PropertyDatabase;
template<size_t> class Model;

/** 
     Evaluates (brittle) shear and tensile failure criteria that are output
     to the element integration points.
     
     Input variables:
     - friction angle 
     - cohesion (scalar)
     - stress (tensor) = Cartesian stress
     - fluid pressure (if this variable is not specified its value is assumed to be zero)
     - Biot coefficient
     
     @todo SKM: this visitor should include the influence of pore pressure.
 
*/
template<size_t dim>
class BrittleFailureVisitor : public Visitor<dim> {
  public:
    explicit BrittleFailureVisitor( Model<dim>&,
                                    bool verbose=false );

    virtual ~BrittleFailureVisitor();
    virtual void Visit( Element<dim>* );
  
  private:
    void InitializeInputProperties( Element<dim>* );

  private:
    csmp::Index           Youngs_key_,
                          Poissons_key_,
                          Stress_key_,     ///< stress
                          Biot_key_,       ///< Biot coefficient alpha
                          Cohesion_key_,   ///< cohesion
                          Friction_key_,   ///< friction angle 0..90o
                          Failure_key_,    ///< failure
                          Pressure_key_;   ///< fluid pressure

    TensorVariable<dim>   Cartesian_stress_;
    MechanicalProperties  mprops_;
    ScalarVariable        fluid_pressure_;
    bool                  verbose_;
};

} // end csmp

#endif /* BRITTLE_FAILURE_VISITOR_H */
