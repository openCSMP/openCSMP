#include <iostream>
#include <cmath>
#include "Rock.h"


using namespace std;

namespace csmp
{

  Rock::Rock()
    : cp(880.),
      t_dependent(true)
  {
    cout << "Rock::Rock( const double& temperature ) constructor:\n";
    cout << "      creating default Rock with temperature-dependent\n";
    cout << "      heat capacity\n\n";
  }

  Rock::Rock( double heatcapacity )
    : cp(heatcapacity), 
      t_dependent(false)
  {
    cout << "Rock::Rock( const double& temperature, double heatcapacity )\n";
    cout << "constructor:\n";
    cout << "      creating Rock with temperature-independent\n";
    cout << "      heat capacity of " << cp << "J/kg/K\n\n";
  }

  Rock::~Rock()
  {
  }

  double Rock::HeatCapacity( double t )
  {
    if(!t_dependent) return cp;
    else
      {
        if( t < 750.) return cp;
        else if( t < 800.) return cp+(t-750.)/50.*cp; // to avoid a jump function
        else return 2.0*cp; // doubled
      }
  }

  double Rock::MinimumHeatCapacity()
  {
    return cp;
  }


  double Rock::Enthalpy( double t )
  {
    if(!t_dependent) return cp*t; 
    else
      {
        if(t < 750.) return 880.*t;
        else if(t < 800.) return 
                            t*880. // base value 
                            + 0.5 * 880./50. * (t-750.)*(t-750.); // 0.5*slope*delta_t^2-0.5*slope*0^2;
        else return 
               t*880.
               + 0.5*880.*50.
               + (t-800.)*880.;
      }
  }
}
