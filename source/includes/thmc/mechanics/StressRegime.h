#ifndef STRESSREGIME_H
#define STRESSREGIME_H

#include "CSMP_definitions.h"
#include "VectorVariable.h"

namespace csmp{

class StressRegime{
public:
  StressRegime( double64 Sv, double64 SH, double64 Sh, double64 trend );
  ~StressRegime(){};

  double64            VerticalStressMagnitude() const;
  double64            MaximumHorizontalStressMagnitude() const;
  double64            MinimumHorizontalStressMagnitude() const;
  double64            MaximumHorizontalStressTrend() const;
  VectorVariable<3U>  VerticalStressUnitVector() const;
  VectorVariable<3U>  MinimumHorizontalStressUnitVector() const;
  VectorVariable<3U>  MaximumHorizontalStressUnitVector() const;

private:
  void                CheckTrend();
  void                EstablishPrincipalStressUnitVectors();
  void                EstablishMinimumHorizontalStressVector();

  double64            Sv_;
  double64            SH_;
  double64            Sh_;
  double64            trend_;
  const double64      PI_;
  VectorVariable<3U>  Svv_;
  VectorVariable<3U>  SHv_;
  VectorVariable<3U>  Shv_;

};

} //csmp

#endif // STRESSREGIME_H
