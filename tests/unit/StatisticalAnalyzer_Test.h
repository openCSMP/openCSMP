#ifndef CSMP_STATISTICAL_ANALYZER_TEST_H
#define CSMP_STATISTICAL_ANALYZER_TEST_H

#include "Test.h"
#include "StatisticalAnalyzer.h"
#include "Model1D.h"

namespace csmp 
{

class StatisticalAnalyzer_Test : public Test 
{
  public:
      StatisticalAnalyzer_Test();
      ~StatisticalAnalyzer_Test();
      void run(); // runs all the tests for the class (register other methods)
      void StatisticalAnalyzerDefineBins();
      void StatisticalAnalyzerRegionPropertyHistogramsElement();
      void StatisticalAnalyzerRegionPropertyHistograms();
      void StatisticalAnalyzer_TestRegionPropertyHistogramsElementProperty2BinningBasedOnProperty1();
  
  private:
    double64 fTolerance;
    Model1D<1U>* model_;
  
}; //end class

} //end csmp

#endif
