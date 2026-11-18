// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_PDE_INTEGRATOR_COMPUTATION_TEST_H
#define CSMP_PDE_INTEGRATOR_COMPUTATION_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;

/** @brief Tests convergence, stability and performance of PDE_Integrator
 
    - compares performance of SparseMatrix and CompressedRowMatrix accumulation.
      - compares different solvers
      - compares accumulation performance for simplex and poly-element type meshes
      
    @section Input Models (all 3D)
    
    fracs3D VSet
    prism_test
    mechanical model of Shuiba formation
    
    @attention while convergence is tested, there is no comparison in here with analytical reference solutions
*/
class PDE_Integrator_Computation_Test : public Test {
  public:
    ~PDE_Integrator_Computation_Test();
    
    /// all the high-level of the testing here
    void run();
    
    static constexpr bool verbose_ = true;
   
  private:

   /// simplex model FracBox from 'vsetMakers'
   Model<3U>* BuildFracBoxModel3D();

   /// builds model and saves it to binary
   Model<3U>* BuildAnsysModel3D( const std::string& model_name );

    /// angle criteria tests for simplex and poly-element mesh returns potentially problematic elements for visualisation
    std::vector<Element<3>*> VerifyMeshQuality();
    
    /// Darcy's law computations
    void AssignProperties_PressureDiffusionAndFlow();
    
    void AssignEssentialConditions_PressureDiffusionAndFlow();
    
    /// Basic steady-state test including post-processing
    void TestSteadyState_PressureDiffusionAndFlow();

    /// Basic steady-state test including post-processing
    double TestTransient_PressureDiffusionAndFlow(); 

    Model<3>* model_ptr_ = nullptr;
                      
  };

} // end namespace csmp

#endif /* CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H */
