#ifndef VISITOR_EXAMPLE_H
#define VISITOR_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// illustration of the 'Visitor' object-oriented design pattern that is used in CSMP
class  Visitor_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
};

} // csmp

#endif // VISITOR_EXAMPLE_H
