#ifndef BOREHOLE_STABILITY2D_VVCase_H
#define BOREHOLE_STABILITY2D_VVCase_H

#include "Test.h"

namespace csmp {

template<size_t> class Model;

class  BoreHole_stability2D_VVCase : public Test{
public:
  BoreHole_stability2D_VVCase(const char* prefix);
  virtual void run();
};

} // csmp

#endif // BOREHOLE_STABILITY2D_VVCase_H
