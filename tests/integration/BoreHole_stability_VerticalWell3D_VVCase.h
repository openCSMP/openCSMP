#ifndef BOREHOLE_STABILITY_VERTICALWELL3D_VVCASE_H
#define BOREHOLE_STABILITY_VERTICALWELL3D_VVCASE_H

#include "Test.h"

namespace csmp {

class BoreHole_stability_VerticalWell3D_VVCase : public Test {
public:
   BoreHole_stability_VerticalWell3D_VVCase(const char* prefix);
   
   virtual void run();
};

} // csmp

#endif // BOREHOLE_STABILITY_VERTICALWELL3D_VVCASE_H
