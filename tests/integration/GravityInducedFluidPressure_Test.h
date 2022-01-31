//
//  GravityInducedFluidPressure_Test.h
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 18/2/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#ifndef GRAVITY_INDUCED_FLUID_PRESSURE_TEST_H
#define GRAVITY_INDUCED_FLUID_PRESSURE_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<size_t> class Model;

/**
    Using a 1D model with 100 linear line elements,
    various ways of computing a fluid-static pressure are exercised using the
    compressible CO2 as the non-linear reference fluid.
 
    Vertical integration of fluid pressure to obtain the absolute pressure as a function
    of depth is tried in different ways, comparing the results with thos obtained by
    top down integration.
 
    @attention Test uses the model 'Greenshank_prism_mini' mini for the 3D testing.
 
    @author SKM
    @date 18/2/2019
 
    @discussion how to achieve both, maximum stability and high accuracy.
*/
class GravityInducedFluidPressure_Test : public Test {
  public:
    GravityInducedFluidPressure_Test();
    virtual ~GravityInducedFluidPressure_Test();
  
    virtual void run();
  
    void OutputResultsToText( const char* file_name ) const;
  
  private:
    void InitialiseModel1D( double model_height );
    void InitialiseModel3D( const char* model );
  
    void InitialiseTemperatureProfile1D( double T_top, double grad_T_K_per_m );
    void InitialiseTemperatureProfile3D( double T_top, double grad_T_K_per_m );
  
    void InitialisePressure3D( double pf_top, double fluid_density );

    void ReferencePressureByTopDownIntegration( double fluid_density );
    void ReferencePressureByTopDownIntegrationCO2( double pf_top );
    void ReferencePressureForFixedDensity( double ref_density );
  
    bool TestComputedWithReferencePressure();
  
    template<size_t dim>
    bool TestResultsByComparison( const Model<dim>* const, const char* test_variable, const char* reference_variable, double tolerance );
  
    /// like in vertical fluid pressure example
    void ComputeCO2Pressure_PDE_Integrator( double pf_top );
    // TODO: remove when errors are fixed
    void ComputeCO2Pressure_PDE_Integrator_CRM( double pf_top );

    /// total pressure like in simulator
    void ComputeCO2Pressure_PDE_Integrator2( double pf_top );

    /// reduced pressure based on deviations from a reference density
    void ComputeCO2PressureFromReducedPressure_PDE_Integrator2( Model<1U>* const, double pf_top );

    /// templatised version for 2D and 3D (small models, Gauss-Jordan solver)
    template<size_t dim>
    void ComputeCO2PressureFromReducedPressure_PDE_Integrator_GaussJordan( Model<dim>* const, double pf_top );

    /// templatised version for 2D and 3D (big models, SAMG)
    template<size_t dim>
    void ComputeCO2PressureFromReducedPressure_PDE_IntegratorCRM_SAMG( Model<dim>* const, double pf_top );

  private:
    Model<1U>* model1D_  = nullptr;
    Model<3U>* model3D_  = nullptr;
    const double patm_ = 100325.;  ///< Pa
    const double ptol_ = 5e3;      ///< tolerance for tests (5 kPa)
    double       top_;             ///< elevation of model top relevative to sea level, will be set in run
    const bool     verbose_;
};

} // end csmp

#endif /* GRAVITY_INDUCED_FLUID_PRESSURE_TEST_H */
