//
//  SplitBoundaryPressureDiffusion_Test.h
//  CSMP integration tests
//
//  Created by Stephan Matthai on 14/10/2022.
//  Copyright (c) 2022 Stephan Matthai & Anne-Laure Tertois. All rights reserved.
//

#ifndef CSMP_SPLITBOUNDARY_PRESSURE_DIFFUSION_TEST_H
#define CSMP_SPLITBOUNDARY_PRESSURE_DIFFUSION_TEST_H

#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

/**
     Tests diffusive pressure transfer between a wavy fracture and the rock matrix
     employing Robin boundary conditions and a coupling coefficient that is based
     on the 1D analytic solution to the transient diffusion problem.
 */
class SplitBoundaryPressureDiffusion_Test : public Test {
  public:
    explicit SplitBoundaryPressureDiffusion_Test( const char* test_model="box_with_wavy_fracture",
                                                  const char* test_variables="SplitBoundaryPressureDiffusion_Test-variables.txt" );
                                     
    virtual ~SplitBoundaryPressureDiffusion_Test();
    
    virtual void run();
  
  private:
    
    /// permeability, porosity, total velocity
    void  ConfigureModel();
  
  private:
    const std::string test_model_, test_variable_file_;
    std::string fracture_name_;
    Model<3U>*  model_ptr_ = nullptr;
    double      model_length_, model_height_, model_width_;   ///< x-direction, y-direction, z-direction
    const bool  verbose_   = false;
};


// AUXILIARY STUFF

/// heat flow 1D (analytic)
double temperatureDiffusion1D( double T0, double kappa, double x, double t );

/// gradient of heat flow @ x0 multiplied with thermal conductivity K
double heatFlow1DAtX0( double T0, double kappa, double t, double K );

/// curvefitting the analytic solution to an actual temperature distribution



} // end csmp

#endif /* defined(CSMP_SPLITBOUNDARY_PRESSURE_DIFFUSION_TEST_H) */
