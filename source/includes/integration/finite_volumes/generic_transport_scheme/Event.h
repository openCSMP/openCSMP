#ifndef EVENT_H
#define EVENT_H

#include "CSMP_number_types.h"

namespace csmp {

struct Event {
    //constructor
    explicit Event();

    //function used for sorting event queue
    bool operator<(const Event &rhs) const { return t_next_ < rhs.t_next_; }; 

    //attributes relevant to event properties
    double64	dC_cumulative_; //cumulative change of solution
    double64	dC_tr_; //target increment of solution
    double64	t_current_; //current time stamp
    double64	t_next_; //next time stamp
    double64	dt_CFL_; //CFL time increment
    double64	dt_tr_; //target time increment

    bool valid_; //validity
    bool inPEPStack_; //whether in PEPStack
};

} // csmp

#endif









