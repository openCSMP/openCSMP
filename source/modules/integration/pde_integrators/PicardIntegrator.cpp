#include "PicardIntegrator.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

namespace csmp {

template<uint32_t dim,template<uint32_t> class CELLTYPE>
PicardIntegrator<dim,CELLTYPE>::PicardIntegrator()
 : IterativeIntegrator<dim,CELLTYPE>()
  {}
  
  
  
  
template<uint32_t dim,template<uint32_t> class CELLTYPE>
double PicardIntegrator<dim,CELLTYPE>::Residual()
{
  resid_.resize(this->G_.Rows());
  std::vector<double>(resid_).swap(resid_);
  
  // Calculate residual vector
  this->G_.MultiplyWith(this->x_, resid_);
  std::transform(resid_.begin(),
                 resid_.end(), 
                 this->rh_.begin(),
                 resid_.begin(),
                 std::minus<double>() );
  
  // Calculate euclidian vector norm
  return std::sqrt(std::inner_product(resid_.begin(), resid_.end(), resid_.begin(), 0.) );
}


template class PicardIntegrator<1U>;
template class PicardIntegrator<2U>;
template class PicardIntegrator<3U>;

template class PicardIntegrator<1U,Face>;
template class PicardIntegrator<2U,Face>;
template class PicardIntegrator<3U,Face>;

} // end namespace csmp
