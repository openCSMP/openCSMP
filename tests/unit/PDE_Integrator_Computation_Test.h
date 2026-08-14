#ifndef CSMP_PDE_INTEGRATOR_COMPUTATION_TEST_H
#define CSMP_PDE_INTEGRATOR_COMPUTATION_TEST_H

//#include "Attorney.h"
#include "CSMP_definitions.h"
#include "Test.h"
//#include "CompressedRowMatrix.h"
//#include "PDE_Integrator.h"

namespace csmp {

/*
template<uint32_t dim,
         template<uint32_t> class CELLTYPE = Element, 
         class MATRIXTYPE = CompressedRowMatrix>
class PDE_Integrator_Attorney2 : public Attorney<PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>> {
    
    // Convenience alias to avoid repeating long template signatures
    using TargetIntegrator = PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>;

public:
    using Attorney<TargetIntegrator>::Attorney; // Inherit Attorney constructor

    // Expose member functions
    using TargetIntegrator::Accumulate;
    using TargetIntegrator::LateAccumulate;
    using TargetIntegrator::EstablishMatrixSetup;
    using TargetIntegrator::EliminateEssentialConditions;
    using TargetIntegrator::AssignEssentialConditions;
    using TargetIntegrator::AssignInitialConditions;

    // Expose member variables (operators and matrices)
    using TargetIntegrator::lhs_operators_;
    using TargetIntegrator::rhs_operators_;
    using TargetIntegrator::G_;
    using TargetIntegrator::rh_;
    using TargetIntegrator::DOF_indexes_;
    using TargetIntegrator::pivotVector_;
};
*/

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
