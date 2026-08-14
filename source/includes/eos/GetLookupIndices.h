#ifndef GET_LOOKUP_INDICES_H
#define GET_LOOKUP_INDICES_H

#include "CSMP_definitions.h"

namespace csmp
{
    long     GetTemperatureIndex(const double& t);
    long     GetPressureIndex(const double& p);
    double GetTemperatureResolution(const double& t);
    double GetPressureResolution(const double& p);
}
#endif
