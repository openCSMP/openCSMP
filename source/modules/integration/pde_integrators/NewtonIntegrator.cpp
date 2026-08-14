#include "NewtonIntegrator.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

namespace csmp {

template<uint32_t dim,template<uint32_t> class CELLTYPE>
NewtonIntegrator<dim,CELLTYPE>::NewtonIntegrator() :
  IterativeIntegrator<dim,CELLTYPE>() {}
  
template<uint32_t dim,template<uint32_t> class CELLTYPE>
double NewtonIntegrator<dim,CELLTYPE>::Residual() {
  return std::sqrt(std::inner_product( this->rh_.begin(), 
                                       this->rh_.end(),
                                       this->rh_.begin(), 0.) );
}


template class NewtonIntegrator<1U>;
template class NewtonIntegrator<2U>;
template class NewtonIntegrator<3U>;

template class NewtonIntegrator<1U,Face>;
template class NewtonIntegrator<2U,Face>;
template class NewtonIntegrator<3U,Face>;

} // end namespace csmp
