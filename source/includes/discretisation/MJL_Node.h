#ifndef MJL_NODE_H
#define MJL_NODE_H

#include "MJL_geometry.h"

namespace mjl {

class Node {
  public:
    Node();
    virtual ~Node();
    Node* Next();
    Node* Prev();
    Node* Insert( Node* );
    Node* Remove();
    void  Splice( Node* );

  protected:
    Node* next_;
    Node* prev_;
};


inline Node::Node()
  : next_(this), prev_(this)
  {
  }
  

inline Node::~Node()
 {
 }


inline Node* Node::Next()
 {
    return next_;
 } 


inline Node* Node::Prev()
 {
    return prev_;
 } 


inline Node* Node::Remove()
 {
    prev_->next_  = next_;
    next_->prev_  = prev_;
    next_ = prev_ = this;
    
    return this; 
 } 

} // end namespace mjl

#endif 
