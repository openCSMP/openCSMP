#ifndef NEWTON_RAPHSON_INTEGRATOR_H
#define NEWTON_RAPHSON_INTEGRATOR_H

#include "IterativeIntegrator.h"

namespace csmp {

/** 
 
@brief Algorithm implementing a Newton iteration scheme for non-linear processes.

@author Adrian Burri
@date 2004

Based on the IterativeAlgorithm base class.

@section collaboration Collaboration

Models, Visitors and Interrelations. Calculations on single Region
objects is not supported.  

*/
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
class NewtonIntegrator : public IterativeIntegrator<dim,COMPUTATION_DOMAIN> {
public:
  NewtonIntegrator();
  
  virtual double Residual();
  
private:
  
};


} // end namespace csmp

#endif
