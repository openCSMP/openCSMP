#ifndef PICARD_INTEGRATOR_H
#define PICARD_INTEGRATOR_H

#include "IterativeIntegrator.h"

namespace csmp {

template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
class PicardIntegrator : public IterativeIntegrator<dim,COMPUTATION_DOMAIN> {
  public:
    PicardIntegrator();
    
    virtual double Residual();
    
  private:
    std::vector<double> resid_;
};




/** 
 
@class PicardIntegrator
@author Adrian Burri
@date 2004

Algorithm implementing a Picard iteration scheme for non-linear processes. 
Derived from the IterativeAlgorithm base class. 
 
*/

} // end namespace csmp
#endif
