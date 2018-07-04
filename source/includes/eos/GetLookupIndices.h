#ifndef GET_LOOKUP_INDICES_H
#define GET_LOOKUP_INDICES_H

#include "CSMP_definitions.h"

namespace csmp
{
    long     GetTemperatureIndex(const double64& t);
    long     GetPressureIndex(const double64& p);
    double64 GetTemperatureResolution(const double64& t);
    double64 GetPressureResolution(const double64& p);
}
#endif
