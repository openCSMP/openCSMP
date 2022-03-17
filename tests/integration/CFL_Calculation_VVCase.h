#ifndef CFL_CALCULATION_VVCASE_H
#define CFL_CALCULATION_VVCASE_H

#include "Test.h"

namespace csmp{
template <uint32_t dim>
class CFL_Calculation_VVCase : public Test
{
public:
    CFL_Calculation_VVCase();
    CFL_Calculation_VVCase(const char* prefix);
    ~CFL_Calculation_VVCase();
    virtual void run();
};
}


#endif // CFL_CALCULATION_VVCASE_H
