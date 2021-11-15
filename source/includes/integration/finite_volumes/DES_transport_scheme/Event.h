#ifndef CSMP_NEW_EVENT_H
#define CSMP_NEW_EVENT_H

#include "Node.h"

namespace csmp {

template<size_t> class Node;

    
template<size_t dim>
class Event {
  public:
    //constructor
    explicit Event(Node<dim>* nd);
    Event() = delete;
    //getting node
    csmp::Node<dim>* getNode() const {return node_;}
    bool valid() const {return valid_;}
    void valid(bool validity) {valid_ = validity;}
    bool inPEPStack() const {return inPEPStack_;}
    void inPEPStack(bool boolean) {inPEPStack_ = boolean;}
    bool inQueue() const {return inQueue_;}
    void inQueue(bool boolean) {inQueue_ = boolean;}    
    double t_schedule() const {return t_schedule_;}
    void t_schedule(double time) {t_schedule_ = time;}
    std::vector<double> facetAreaCollection;
    std::vector<Point<dim>> facetNormalCollection;

  private:
    csmp::Node<dim>*	node_;
    bool valid_; //validity
    bool inPEPStack_; //whether in PEPStack  
    bool  inQueue_; //whether in event queue      
    double t_schedule_; //scheduled time stamp used for sorting event queue
};

} // csmp

#endif









