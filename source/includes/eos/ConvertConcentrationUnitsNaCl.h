// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONVERT_CONCENTRATION_UNITS_NACL_H
#define CONVERT_CONCENTRATION_UNITS_NACL_H

#include "CSMP_definitions.h"

namespace csmp
{
  const double mnacl = 58.443e0,mh2o=18.015e0,molal=55.5093e0;
  
  inline double MolarWeightSolution(double x)
  {
    return x*mnacl+(1.0e0-x)*mh2o;
  }

  inline double Weight2XNaCl(double weightpercent)
  {
    return weightpercent/mnacl/(weightpercent/mnacl+(100.0e0-weightpercent)/mh2o);
  }

  inline double XNaCl2Weight(double x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o)*100.0e0;
  }

  inline double XNaCl2Massfraction(double x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double Massfraction2XNaCl(double smf)
  {
    return smf/mnacl/(smf/mnacl+(1.0e0-smf)/mh2o);;
  }

  inline double XNaCl2WeightFraction(double x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double XNaCl2Molal(double x)
  {
    return molal*x/(1.0e0-x);
  }

  inline double Molal2XNaCl(double molality)
  {
    return molality/(molality+molal);
  }

  inline double DensityAndX2MolarVolume(double density, double x)
  {
    return 1.0e3/density*(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double DensityAndWt2MolarVolume(double density, double weightpercent)
  {
    return DensityAndX2MolarVolume(density,Weight2XNaCl(weightpercent));
  }

  inline double MolarVolumeAndX2Density(double molarvolume, double x)
  {
    // now in kg/m^3
    return (x*mnacl+(1.0-x)*mh2o)/molarvolume*1000.0;
  }

  inline double WtFractionNaClFromX(double x)
  {
    return x*mnacl/(x*mnacl+(1.0e0-x)*mh2o);
  }

  inline double DXNaClDT2DWeightDT(double x, double dxdt)
  {
    return 100.0*dxdt*mnacl*mh2o/(x*mnacl+mh2o-mh2o*x)/(x*mnacl+mh2o-mh2o*x);
  }

  inline double DXNaClDP2DWeightDP(double x, double dxdp)
  {
    // caution, dxdp has unit of bar^-1
    return 100.0*dxdp*mnacl*mh2o/(x*mnacl+mh2o-mh2o*x)/(x*mnacl+mh2o-mh2o*x);
  }
}
#endif







