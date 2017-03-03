#ifndef TRANSPORT_STEP_SIZE_H
#define TRANSPORT_STEP_SIZE_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class TransportStepSize : public Interrelation<dim> {
  public:
    TransportStepSize( const PropertyDatabase<dim>& p, 
                       const char* n_velo, const char* e_velo, 
                       const char* i_radius );
                       
    ~TransportStepSize();
    
    double64 AdvectionTimeIncrement() const;
    
    void Calculate();
  
  private:
    Operand<dim>&        NV; /// < nodal velocity
    Operand<dim>&        EV; /// < element velocity
    Operand<dim>&        IR; /// < inner radius
    ScalarVariable       ir;
    VectorVariable<dim>  nv, ev;
    double64             advection_increment,
                         e_increment,
                         n_velo, e_velo, 
                         scalar_velocity;
    bool                 called;
};

} // csmp

#endif

