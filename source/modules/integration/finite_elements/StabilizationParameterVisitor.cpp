#include "StabilizationParameterVisitor.h"
#include "Model.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
StabilizationParameterVisitor<dim>::StabilizationParameterVisitor( Model<dim>& mdl,
                                                                   double mu, double coeff,
                                                                   const char* stab_param )
  : viscosity_(mu), coefficient_(coeff),
    stparam_key_(mdl.Database().StorageKey(stab_param))
{
    //extern ErrorHandler skm_err;
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(ELEMENT);
} // end constructor



template<uint32_t dim>
StabilizationParameterVisitor<dim>::StabilizationParameterVisitor( Model<dim>& mdl,
                                                                   double coeff,
                                                                   const char* stab_param)
  : viscosity_(1.6e-3), // water at room temperature
    coefficient_(coeff),
    visc_key_(mdl.Database().StorageKey("viscosity")),
    stparam_key_(mdl.Database().StorageKey(stab_param))
{
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(ELEMENT);

} // end constructor



template<uint32_t dim>
void StabilizationParameterVisitor<dim>::Visit( Element<dim>* e )
  {
     e->SegmentLengths( segments_ );
     const double min_segment = *min_element(segments_.begin(), segments_.end());
    
     // if the viscosity key has been defined
     if ( visc_key_.place != UNDEFINED ) viscosity_ = e->Read( visc_key_ );
    
     ScalarVariable sp( e->Status(stparam_key_), min_segment * min_segment / 12.0 / viscosity_ * coefficient_ );
     
     e->Store( stparam_key_, sp );
}
     
     
template class StabilizationParameterVisitor<1U>;
template class StabilizationParameterVisitor<2U>;
template class StabilizationParameterVisitor<3U>;
     
} // namespace csmp
