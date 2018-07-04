#include "PicardIntegrator.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

namespace csmp {

template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
PicardIntegrator<dim,COMPUTATION_DOMAIN>::PicardIntegrator() 
 : IterativeIntegrator<dim,COMPUTATION_DOMAIN>()
  {}
  
  
  
  
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
double64 PicardIntegrator<dim,COMPUTATION_DOMAIN>::Residual() 
{
  resid_.resize(this->G_.Rows());
  std::vector<double64>(resid_).swap(resid_);
  
  // Calculate residual vector
  this->G_.MultiplyWith(this->x_, resid_);
  std::transform(resid_.begin(),
                 resid_.end(), 
                 this->rh_.begin(),
                 resid_.begin(),
                 std::minus<double64>() );
  
  // Calculate euclidian vector norm
  return std::sqrt(std::inner_product(resid_.begin(), resid_.end(), resid_.begin(), 0.) );
}


template class PicardIntegrator<1U,Region>;
template class PicardIntegrator<2U,Region>;
template class PicardIntegrator<3U,Region>;

template class PicardIntegrator<1U,Boundary>;
template class PicardIntegrator<2U,Boundary>;
template class PicardIntegrator<3U,Boundary>;

template class PicardIntegrator<1U,SplitBoundary>;
template class PicardIntegrator<2U,SplitBoundary>;
template class PicardIntegrator<3U,SplitBoundary>;

} // end namespace csmp
