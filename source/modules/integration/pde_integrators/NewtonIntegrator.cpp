#include "NewtonIntegrator.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

namespace csmp {

template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
NewtonIntegrator<dim,COMPUTATION_DOMAIN>::NewtonIntegrator() :
  IterativeIntegrator<dim,COMPUTATION_DOMAIN>() {}
  
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
double64 NewtonIntegrator<dim,COMPUTATION_DOMAIN>::Residual() {
  return std::sqrt(std::inner_product( this->rh_.begin(), 
                                       this->rh_.end(),
                                       this->rh_.begin(), 0.) );
}

template class NewtonIntegrator<1U,Region>;
template class NewtonIntegrator<2U,Region>;
template class NewtonIntegrator<3U,Region>;

template class NewtonIntegrator<1U,Boundary>;
template class NewtonIntegrator<2U,Boundary>;
template class NewtonIntegrator<3U,Boundary>;

template class NewtonIntegrator<1U,SplitBoundary>;
template class NewtonIntegrator<2U,SplitBoundary>;
template class NewtonIntegrator<3U,SplitBoundary>;

} // end namespace csmp
