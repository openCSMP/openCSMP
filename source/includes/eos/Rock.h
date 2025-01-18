#ifndef   ROCK_H
#define   ROCK_H

#include "CSMP_definitions.h"

namespace csmp
{

  class Rock
  {
  public:
    
    Rock(); // to be used with t-dependent heat capacity
    Rock( double heatcapacity ); // to be used with constant heat capacity
    ~Rock();

    double HeatCapacity( double t);
    double MinimumHeatCapacity();
    double Enthalpy( double t);
    
  private:

    //Rock();
    const double  cp;
    const bool t_dependent;
  };
}

#endif
