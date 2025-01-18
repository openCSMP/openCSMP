#ifndef STRESSREGIME_H
#define STRESSREGIME_H

#include "CSMP_definitions.h"
#include "VectorVariable.h"

namespace csmp{

class StressRegime {
public:
  StressRegime( double Sv, double SH, double Sh, double trend );
  ~StressRegime(){};

  double            VerticalStressMagnitude() const;
  double            MaximumHorizontalStressMagnitude() const;
  double            MinimumHorizontalStressMagnitude() const;
  double            MaximumHorizontalStressTrend() const;
  VectorVariable<3U>  VerticalStressUnitVector() const;
  VectorVariable<3U>  MinimumHorizontalStressUnitVector() const;
  VectorVariable<3U>  MaximumHorizontalStressUnitVector() const;

private:
  void                CheckTrend();
  void                EstablishPrincipalStressUnitVectors();
  void                EstablishMinimumHorizontalStressVector();

  double            Sv_;
  double            SH_;
  double            Sh_;
  double            trend_;
  const double      PI_;
  VectorVariable<3U>  Svv_;
  VectorVariable<3U>  SHv_;
  VectorVariable<3U>  Shv_;

};

} //csmp

#endif // STRESSREGIME_H
