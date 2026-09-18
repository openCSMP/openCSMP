//
//  DirichletPressureBoxModel_VVCase.hpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 1/02/2018.
//

#ifndef DIRICHLET_PRESSURE_BOX_MODEL_VV_CASE_H
#define DIRICHLET_PRESSURE_BOX_MODEL_VV_CASE_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

class DirichletPressureBoxModel_VVCase : public Test {
  public:
     DirichletPressureBoxModel_VVCase( const char* model="berea_sub100" );
     
     virtual void run();
  
  private:
    /// box model with standard isoparametric finite elements
    void TestModelFromANSYS();
  
    /// box model with analytically integrated elements
    void TestModelFromANSYS_AnalyticallyIntegrated();
  
    const std::string  model_name_;
    const bool         verbose_;
};


} // end csmp

#endif /* DIRICHLET_PRESSURE_BOX_MODEL_VV_CASE_H */
