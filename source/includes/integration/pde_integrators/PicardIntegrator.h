#ifndef PICARD_INTEGRATOR_H
#define PICARD_INTEGRATOR_H

#include "IterativeIntegrator.h"

namespace csmp {

template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
class PicardIntegrator : public IterativeIntegrator<dim,SIMPLICIAL_COMPLEX> {
  public:
    PicardIntegrator();
    
    virtual double64 Residual();
    
  private:
    std::vector<double64> resid_;
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
