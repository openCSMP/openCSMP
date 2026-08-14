#ifndef CSMP_PROPERTY_DATA_TEST_H
#define CSMP_PROPERTY_DATA_TEST_H

#include "Test.h"
#include "PropertyData.h"

namespace csmp {

class PropertyData_Test : public Test {
  public:
    virtual void run();
    const static bool verbose_ = false;
};

} // namespace csmp

#endif // CSMP_PROPERTY_DATA_TEST_H
