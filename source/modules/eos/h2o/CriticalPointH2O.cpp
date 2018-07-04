#include "CriticalPointH2O.h"

#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

namespace csmp
{
  CriticalPointH2O::CriticalPointH2O(){}
  CriticalPointH2O::~CriticalPointH2O(){}


  double64 CriticalPointH2O::Temperature()      const { return 373.976e0; }
  double64 CriticalPointH2O::Pressure()         const { return 220.54915e5; }
  double64 CriticalPointH2O::MassFractionNaCl() const { return 0.0e0; }
  double64 CriticalPointH2O::Density()          const { return 321.89e0; }
  double64 CriticalPointH2O::MolarVolume()      const { return mh2o/321.89e0/1000.0; }
  double64 CriticalPointH2O::Enthalpy()         const { return 2086.0e3; }
}
