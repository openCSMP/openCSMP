#ifndef MOHR_COULOMB_FAILURE_CRITERION_H
#define MOHR_COULOMB_FAILURE_CRITERION_H

#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class MohrCoulombFailure : public Interrelation<dim> {
    Operand<dim>&  MS;        // mean stress
    Operand<dim>&  STRESS;
    Operand<dim>&  CRIT;
    Operand<dim>&  COH;       // cohesion of material
    double       phi;       // angle of friction
    double       Ts;        // tensile strength

    double    MeanStress( const TensorVariable<dim>& ts );
    double    DeviatoricStress( const TensorVariable<dim>& ts, double& t );
    double    Theta( const TensorVariable<dim>& ts, double t );
    double    G_OfTheta( double theta ); // Zienkiewitz II, p. 89
    
    TensorVariable<dim>  ts;
    ScalarVariable       ch;   

  public:
    MohrCoulombFailure( const PropertyDatabase<dim>& p,
                        double friction_angle );
                        
    ~MohrCoulombFailure();
    
    void Calculate();
};


} // csmp

#endif

