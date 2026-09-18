// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Event.h"


namespace csmp {

template<uint32_t dim>
Event<dim>::Event(Node<dim>* nd)
  : node_(nd),
    heap_node_(nullptr),
    t_schedule_(0.),
    valid_(false),
    inPEPStack_(false),
    inQueue_(false)
{
}

template class Event<1U>;
template class Event<2U>;
template class Event<3U>;

} // end namespace csmp  
