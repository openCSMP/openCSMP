#ifndef CSMP_NEW_EVENT_H
#define CSMP_NEW_EVENT_H

#include "Node.h"
#include "FibonacciHeap.h"

namespace csmp {

template<uint32_t> class Node;

    
template<uint32_t dim>
class Event {
  public:
    //constructor
    explicit Event(Node<dim>* nd);
    Event(Event&&); // TODO: complete rule of 5
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

    typedef ajb::detail::FibonacciHeap_Node<double,size_t> Heap_Node;
    void setHeapNode(Heap_Node* heap_node) {heap_node_ = heap_node;}
    Heap_Node* getHeapNode() {return heap_node_;}

  private:
    csmp::Node<dim>*	node_;
    Heap_Node* heap_node_;
    double t_schedule_; //scheduled time stamp used for sorting event queue
    //validity (false = event needs to be rescheduled after it is executed, becomes true after re-scheduling)
    bool valid_;
    bool inPEPStack_; //whether in PEPStack  
    bool  inQueue_; //whether in event queue
};

} // csmp

#endif









