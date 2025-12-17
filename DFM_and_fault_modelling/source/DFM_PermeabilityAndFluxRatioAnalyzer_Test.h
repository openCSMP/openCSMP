//
//  DFM_PermeabilityAndFluxRatioAnalyzer_Test.h
//  CSMP_DFM_Analyzer
//
//  Created by Stephan Matthai on 6/6/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
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
