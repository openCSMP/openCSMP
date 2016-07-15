#ifndef COMPUTETWOPHASEFLOWPROPERTIES_H
#define COMPUTETWOPHASEFLOWPROPERTIES_H
#include <cstddef>
namespace csmp
{
template<size_t, template<size_t> class>
class ModelSubDomain;

template<size_t>
class TwoPhaseModel;

template<size_t dim,template<size_t> class SIMPLEX>
void  compute2PhaseMobilityAtBaryCenter( ModelSubDomain<dim,SIMPLEX>& sg,
                                         TwoPhaseModel<dim>& relperm);
}
#endif // COMPUTETWOPHASEFLOWPROPERTIES_H
