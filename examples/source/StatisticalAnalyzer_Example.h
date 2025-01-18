#ifndef CSMP_STATISTICAL_ANALYZER_EXAMPLE_H
#define CSMP_STATISTICAL_ANALYZER_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  StatisticalAnalyzer_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();

private:
  double NormalDistributionGenerator( double mean, double sd,
                                      double minimum, double maximum );
};

} // csmp

#endif // CSMP_STATISTICAL_ANALYZER_EXAMPLE_H
