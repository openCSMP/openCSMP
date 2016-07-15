#ifndef   ROCK_H
#define   ROCK_H

#include "CSMP_definitions.h"

namespace csmp
{

  class Rock
  {
  public:
    
    Rock(); // to be used with t-dependent heat capacity
    Rock( double64 heatcapacity ); // to be used with constant heat capacity
    ~Rock();

    double64 HeatCapacity( double64 t);
    double64 MinimumHeatCapacity();
    double64 Enthalpy( double64 t);
    
  private:

    //Rock();
    const double64  cp;
    const bool t_dependent;
  };
}

#endif
