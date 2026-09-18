// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DFM_PermeabilityAndFluxRatioAnalyzer_Test.h
//  CSMP_DFM_Analyzer
//
//

#ifndef DFM_PermeabilityAndFluxRatioAnalyzer_Test_H
#define DFM_PermeabilityAndFluxRatioAnalyzer_Test_H

#include "Test.h"
#include "DFM_PermeabilityAndFluxRatioAnalyzer.h"

namespace csmp {

class  DFM_PermeabilityAndFluxRatioAnalyzer_Test : public Test {
  public:
    DFM_PermeabilityAndFluxRatioAnalyzer_Test();
    virtual void run();

  private:
    DFM_PermeabilityAndFluxRatioAnalyzer<3U>  DFM_analyzer_;
  
};


} // end namespace


#endif /* defined(DFM_PermeabilityAndFluxRatioAnalyzer_Test_H) */
