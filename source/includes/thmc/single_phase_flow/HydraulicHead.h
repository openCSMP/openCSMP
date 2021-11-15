#ifndef HYDRAULIC_HEAD_H
#define HYDRAULIC_HEAD_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class HydraulicHead : public Interrelation<dim> {
    Operand<dim>&   P;    // absolute fluid pressure
    Operand<dim>&   E;    // vertical elevation (zero at model base)
    Operand<dim>&   H;    // hydraulic head
    ScalarVariable height, pres;
    double                 rho0, g; // acceleration of gravity
    
  public:
    HydraulicHead( const PropertyDatabase<dim>& p, double ref_density=1000.0 );
    ~HydraulicHead() {};
    void Calculate();
};

} // csmp

#endif

