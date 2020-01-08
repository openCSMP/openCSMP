#ifndef QUADRILATERATOR_TO_CSMP_BINARY_EXAMPLE_H
#define QUADRILATERATOR_TO_CSMP_BINARY_EXAMPLE_H

#include "Example.h"

namespace csmp {

template<size_t> class Model;

bool is_NO_DATA_Value( double value );
void createInflowRegion( Model<2U>& );
  
class QuadrilateratorToCSMPbinary_Example : public Example {
   public:
     virtual void Run();
     virtual void Specifications();
};


} // csmp


#endif // QUADRILATERATOR_TO_CSMP_BINARY_EXAMPLE_H
