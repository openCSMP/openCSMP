#ifndef MOHR_COULOMB_FAILURE_CRITERION_H
#define MOHR_COULOMB_FAILURE_CRITERION_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class MohrCoulombFailure : public Interrelation<dim> {
    Operand<dim>&  MS;        // mean stress
    Operand<dim>&  STRESS;
    Operand<dim>&  CRIT;
    Operand<dim>&  COH;       // cohesion of material
    double64       phi;       // angle of friction
    double64       Ts;        // tensile strength

    double64    MeanStress( const TensorVariable<dim>& ts );
    double64    DeviatoricStress( const TensorVariable<dim>& ts, double64& t );
    double64    Theta( const TensorVariable<dim>& ts, double64 t );
    double64    G_OfTheta( double64 theta ); // Zienkiewitz II, p. 89
    
    TensorVariable<dim>  ts;
    ScalarVariable       ch;   

  public:
    MohrCoulombFailure( const PropertyDatabase<dim>& p,
                        double64 friction_angle );
                        
    ~MohrCoulombFailure();
    
    void Calculate();
};


} // csmp

#endif

