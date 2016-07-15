#include "StabilizationParameterVisitor.h"
#include "Model.h"

using namespace std;

namespace csmp {

template<size_t dim>
StabilizationParameterVisitor<dim>::StabilizationParameterVisitor( Model<dim>& mdl,
                                                                   double64 mu, double64 coeff,
                                                                   const char* stab_param)
  : pref(mdl.Database()),
    viscosity(mu), coefficient(coeff),
    stparam_key(mdl.Database().StorageKey(stab_param))
{
    //extern ErrorHandler skm_err;
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(ELEMENT);
} // end constructor


template<size_t dim>
StabilizationParameterVisitor<dim>::StabilizationParameterVisitor( Model<dim>& mdl,
                                                                   double64 coeff,
                                                                   const char* stab_param)
  : pref(mdl.Database()),
    viscosity(1.6e-3), // water at room temperature
    coefficient(coeff),
    visc_key_(mdl.Database().StorageKey("viscosity")),
    stparam_key(mdl.Database().StorageKey(stab_param))
{
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(ELEMENT);

} // end constructor



template<size_t dim>
StabilizationParameterVisitor<dim>::~StabilizationParameterVisitor() {}

template<size_t dim>
void StabilizationParameterVisitor<dim>::Visit( Model<dim>* m )
  {
  }

template<size_t dim>
void StabilizationParameterVisitor<dim>::Visit( Element<dim>* e )
  {
     e->SegmentLengths( segments );
     min_segment = *min_element(segments.begin(), segments.end());
    
     // if the viscosity key has been defined
     if ( visc_key_.place != UNDEFINED ) viscosity = e->Read( visc_key_ );
    
     sp() = min_segment*min_segment/12.0/viscosity*coefficient;
     
    /* for (stl_index i = 0; i <= dim; i++)
     {
       dir = 0.0;
       dir(i) = 1.0;
       min_segment = n->LengthOfElementInDirection(dir);
       sp(i) = min_segment*min_segment/12.0/viscosity;
     }
     */
     e->Store( stparam_key, sp );
}
     
     
template class StabilizationParameterVisitor<1U>;
template class StabilizationParameterVisitor<2U>;
template class StabilizationParameterVisitor<3U>;
     
} // namespace csmp
