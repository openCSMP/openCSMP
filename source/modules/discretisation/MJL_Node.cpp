#include "MJL_Node.h"

namespace mjl {

// cannot be constant since b's connections are modified
Node* Node::Insert( Node* b )
 {
    Node* c = next_;
    b->next_    = c;
    b->prev_    = this;
    next_       = b;
    c->prev_    = b;
    
    return b;
 } 


// cannot be constant since b's connections are modified
void Node::Splice( Node* b )
 {
    Node* a = this;
    Node* an = a->next_;
    Node* bn = b->next_;
    a->next_     = bn;
    b->next_     = an;
    an->prev_    = b;
    bn->prev_    = a;
 } 

} // end namespace mjl
