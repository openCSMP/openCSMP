#ifndef CONCENTRATION_FLUID_VOLUME_SOURCE_H
#define CONCENTRATION_FLUID_VOLUME_SOURCE_H

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class ConcentrationFluidVolumeSource : public Interrelation<dim> {
    Operand<dim>&  CONCP;    
    Operand<dim>&  CONCN;    
    Operand<dim>&  Q;    
    Operand<dim>&  PHI;    
    ScalarVariable concn, concp, phi;
    double dt, rho_increment, rho_zero;
    
  public:
    ConcentrationFluidVolumeSource( const PropertyDatabase<dim>&, double rho_max );
    virtual ~ConcentrationFluidVolumeSource() = default;
    void TimeIncrement( double time_increment );
    
    void Calculate() override final;
};

template<uint32_t dim>
inline void ConcentrationFluidVolumeSource<dim>::TimeIncrement( double time_increment ) 
{ dt = time_increment; }

}

#endif

