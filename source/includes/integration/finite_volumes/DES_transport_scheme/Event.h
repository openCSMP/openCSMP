#ifndef CSMP_NEW_EVENT_H
#define CSMP_NEW_EVENT_H

#include "CSMP_number_types.h"
#include "Node.h"

namespace csmp {

template<size_t> class Node;

    
template<size_t dim>
class Event {
  public:
    //constructor
    Event(Node<dim>* nd);
    //getting node
    csmp::Node<dim>* getNode() const {return node_;}
    bool valid() const {return valid_;}
    void valid(bool validity) {valid_ = validity;}
    bool inPEPStack() const {return inPEPStack_;}
    void inPEPStack(bool boolean) {inPEPStack_ = boolean;}
    double64 t_schedule() const {return t_schedule_;}
    void t_schedule(double64 time) {t_schedule_ = time;}

  private:
    csmp::Node<dim>*	node_;
    bool valid_; //validity
    bool inPEPStack_; //whether in PEPStack       
    double64 t_schedule_; //scheduled time stamp used for sorting event queue
};

} // csmp

#endif









