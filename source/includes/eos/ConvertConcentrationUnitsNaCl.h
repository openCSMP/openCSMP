#ifndef CONVERT_CONCENTRATION_UNITS_NACL_H
#define CONVERT_CONCENTRATION_UNITS_NACL_H

#include "CSMP_definitions.h"

namespace csmp
{
  const double64 mnacl = 58.443e0,mh2o=18.015e0,molal=55.5093e0;
  
  inline double64 MolarWeightSolution(double64 x)
  {
    return x*mnacl+(1.0e0-x)*mh2o;
  }

  inline double64 Weight2XNaCl(double64 weightpercent)
  {
    return weightpercent/mnacl/(weightpercent/mnacl+(100.0e0-weightpercent)/mh2o);
  }

  inline double64 XNaCl2Weight(double64 x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o)*100.0e0;
  }

  inline double64 XNaCl2Massfraction(double64 x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double64 Massfraction2XNaCl(double64 smf)
  {
    return smf/mnacl/(smf/mnacl+(1.0e0-smf)/mh2o);;
  }

  inline double64 XNaCl2WeightFraction(double64 x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double64 XNaCl2Molal(double64 x)
  {
    return molal*x/(1.0e0-x);
  }

  inline double64 Molal2XNaCl(double64 molality)
  {
    return molality/(molality+molal);
  }

  inline double64 DensityAndX2MolarVolume(double64 density, double64 x)
  {
    return 1.0e3/density*(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double64 DensityAndWt2MolarVolume(double64 density, double64 weightpercent)
  {
    return DensityAndX2MolarVolume(density,Weight2XNaCl(weightpercent));
  }

  inline double64 MolarVolumeAndX2Density(double64 molarvolume, double64 x)
  {
    // now in kg/m^3
    return (x*mnacl+(1.0-x)*mh2o)/molarvolume*1000.0;
  }

  inline double64 WtFractionNaClFromX(double64 x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double64 DXNaClDT2DWeightDT(double64 x, double64 dxdt)
  {
    return 100.0*dxdt*mnacl*mh2o/(x*mnacl+mh2o-mh2o*x)/(x*mnacl+mh2o-mh2o*x);
  }

  inline double64 DXNaClDP2DWeightDP(double64 x, double64 dxdp)
  {
    // caution, dxdp has unit of bar^-1
    return 100.0*dxdp*mnacl*mh2o/(x*mnacl+mh2o-mh2o*x)/(x*mnacl+mh2o-mh2o*x);
  }
}
#endif







