#include "Event.h"


namespace csmp {

template<uint32_t dim>
Event<dim>::Event(Node<dim>* nd)
  : node_(nd),
    valid_(false),
    inPEPStack_(false),
    inQueue_(false),
    t_schedule_(0.)
{
}

template class Event<1U>;
template class Event<2U>;
template class Event<3U>;

} // end namespace csmp  
