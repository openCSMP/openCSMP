// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DFM_Simulator_UnitTest_main.c
//  CSMP_DFM_Analyzer
//
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

  
