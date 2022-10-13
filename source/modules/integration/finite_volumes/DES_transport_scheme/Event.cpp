#include "Event.h"


namespace csmp {

template<uint32_t dim>
Event<dim>::Event(Node<dim>* nd)
  : node_(nd),
    t_schedule_(0.),
    valid_(false),
    inPEPStack_(false),
    inQueue_(false)
{
}

//move constructor
  template<uint32_t dim>
  Event<dim>::Event(Event&& event)
  : node_(event.node_),
    t_schedule_(event.t_schedule_),
    valid_(event.valid_),
    inPEPStack_(event.inPEPStack_),
    inQueue_(event.inQueue_)
  {
  }

template class Event<1U>;
template class Event<2U>;
template class Event<3U>;

} // end namespace csmp  
