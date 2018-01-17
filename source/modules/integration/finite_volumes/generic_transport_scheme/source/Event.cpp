#include "Event.h"


namespace csmp {

//template<size_t dim>
Event::Event()
  : dC_cumulative_(0.),  //cumulative change of solution
    dC_tr_(0.),  //target increment of solution
    t_current_(0.),  //current time stamp
    t_next_(0.),  //next time stamp
    dt_CFL_(0.),  //CFL time increment
    dt_tr_(0.),  //target time increment
    valid_(false),
    inPEPStack_(false)
{
}

} // end namespace csmp  
