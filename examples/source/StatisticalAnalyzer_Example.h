#ifndef STATISTICALANALYZER_EXAMPLE_H
#define STATISTICALANALYZER_EXAMPLE_H

#include "CSMP_number_types.h"
#include "Example.h"

namespace csmp {

class  StatisticalAnalyzer_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();

private:
  double64 NormalDistributionGenerator( double64 mean, double64 sd,
                                        double64 minimum, double64 maximum );
};

} // csmp

#endif // STATISTICALANALYZER_EXAMPLE_H
