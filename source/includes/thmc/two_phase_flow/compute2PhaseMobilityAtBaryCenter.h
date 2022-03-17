#ifndef COMPUTE_TWO_PHASE_FLOW_PROPERTIES_H
#define COMPUTE_TWO_PHASE_FLOW_PROPERTIES_H

#include "CSMP_definitions.h"

namespace csmp
{
template<uint32_t, template<uint32_t> class>
class ModelSubDomain;

template<uint32_t>
class TwoPhaseModel;

template<uint32_t dim,template<uint32_t> class CELL>
void  compute2PhaseMobilityAtBaryCenter( ModelSubDomain<dim,CELL>& sg,
                                         TwoPhaseModel<dim>& relperm);
}
#endif // COMPUTE_TWO_PHASE_FLOW_PROPERTIES_H
