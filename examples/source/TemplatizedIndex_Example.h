#ifndef TEMPLATIZED_INDEX_EXAMPLE_H
#define TEMPLATIZED_INDEX_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// illustrates how Index has been extended into INDEX to support static polymorphism
class TemplatizedIndex_Example : public Example {
  public:
    virtual void Specifications();
    virtual void Run();

 };

 } // end csmp
 
#endif
