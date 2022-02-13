#ifndef QUARTZ_PRECIPITATION_H
#define QUARTZ_PRECIPITATION_H

#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class QuartzPrecipitation : public Interrelation<dim> {
    Operand<dim>&        F;   /// < fluid velocity
    Operand<dim>&        DS;  /// < quartz solubility gradient
    Operand<dim>&        FD;  /// < fluid density
    Operand<dim>&        Q;   /// < amount of quartz precipitated
    ScalarVariable      rho, reactant;
    VectorVariable<dim>  flux, dSdz, dSpr;
    double                      grad_S, dt, mass_flux;
    double                      sign, angle;
    
  public:
    QuartzPrecipitation( const PropertyDatabase<dim>& p, double time_increment );
    ~QuartzPrecipitation() {};
    void Calculate();
    void SetTimeIncrement( double time_increment ) { dt=time_increment; };
};

} // csmp


#endif



















