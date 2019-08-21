//
//  ExplicitTransport_Test.h
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 2/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_EXPLICIT_TRANSPORT_TEST_H
#define CSMP_EXPLICIT_TRANSPORT_TEST_H

#include "Test.h"

namespace csmp {

template<size_t> class Model;

/**
     Testing programme
 
     1. Conservation on the initialised model
        - prescribed velocity
        - computed velocity
        - tracer transport in steady-state pressure field, 1st order
        - tracer transport 2nd order
        - transport for transient pressure case
        - transport in the presence of sources and sinks
*/
class ExplicitTransport_Test : public Test {
  public:
    ExplicitTransport_Test( const char* test_model, const char* test_variables );
    virtual ~ExplicitTransport_Test();
    
    virtual void run();
  
  private:
    Model<3U>*  CreateModel( const char* ansys_input_data );
    Model<3U>*  CreateTetrahedralModel();
    Model<3U>*  CreateHexahedralModel();
  
    // UTILITIES
  
    /// permeability, porosity, total velocity
    void  AssignFlowProperties();
  
    /// computes fluid pressure and a divergence-free velocity field for LEFT-RIGHT fluid-pressure gradient
    void  DivergenceFreeTotalVelocityField( double64 delta_pf=1.0e7 );
  
    /// tests whether the normal of the facet is outward pointing relative to the FV sector; @return 1 or -1 (inward pointing)
    double64  FluxMultiplier( const Element<3U>* const eptr, size_t sector, size_t facet ) const;

    // TESTING
  
    /// testing the function that intialises the FV properties
    void  Test_initializeFiniteVolumeProperties( double64 tolerance_relaxation_factor = 1000. );
  
    /// is the flux balance equal to numeric zero in a divergence free velocity field ?
    void  TestInteriorFluxBalance( double64 tolerance_relaxation_factor = 100. );
  
    /// verifies that there is indeed no flux divergence in FVs that are located at a no-flow boundary
    void  TestNoFlowBoundaryFluxBalance( double64 tolerance_relaxation_factor = 100. );
  
    /// TVD? - first-order scheme during transport in a const velocity field and out of model
    void  TestFlowThroughModel( const char* model = "BOX40x3x10m", bool prescribed_velocity=false );
  
    // fluid sources and sinks
  
    // transient velocity field
  
    // conservation of a tracer plume over a certain transport distance
  
    // TVD? - second-order scheme during transport in a const velocity field and out of model

  
  private:
    const std::string test_model_, test_variable_file_;
    Model<3U>*  model_ptr_ = nullptr;
    double64    model_length_, model_height_, model_width_;   ///< x-direction, y-direction, z-direction
    const bool  verbose_   = false;
};


} // end csmp

#endif /* defined(CSMP_EXPLICIT_TRANSPORT_TEST_H) */
