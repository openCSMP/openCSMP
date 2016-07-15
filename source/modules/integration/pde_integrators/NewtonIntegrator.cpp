#include "NewtonIntegrator.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

namespace csmp {

template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
NewtonIntegrator<dim,SIMPLICIAL_COMPLEX>::NewtonIntegrator() :
  IterativeIntegrator<dim,SIMPLICIAL_COMPLEX>() {}
  
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
double64 NewtonIntegrator<dim,SIMPLICIAL_COMPLEX>::Residual() {
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
