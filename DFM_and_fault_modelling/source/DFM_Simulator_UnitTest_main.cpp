//
//  DFM_Simulator_UnitTest_main.c
//  CSMP_DFM_Analyzer
//
//  Created by Stephan Matthai on 6/10/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include "DFM_PermeabilityAndFluxRatioAnalyzer_Test.h"
#include "TestSuite.h"

using namespace std;
using namespace csmp;

int main()
{
  TestSuite s("DFM equivalent permeability tests", &cout );

  s.addTest( new DFM_PermeabilityAndFluxRatioAnalyzer_Test() );

  s.run();
  s.report();
  
  return 0;
 
} // end testing

  
