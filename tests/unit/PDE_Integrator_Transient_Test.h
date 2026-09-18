#ifndef CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H
#define CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H

#include "Attorney.h"
#include "CSMP_definitions.h"
#include "Test.h"
#include "CompressedRowMatrix.h"
#include "PDE_Integrator.h"

namespace csmp {

template<uint32_t dim,
         template<uint32_t> class CELLTYPE = Element, 
         class MATRIXTYPE = CompressedRowMatrix>
class PDE_Integrator_Attorney1 : public Attorney<PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>> {
    
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


class PDE_Integrator_Transient_Test : public Test {
  public:
    void run();
   
  private:

    void TestSteadyState( PDE_Integrator_Attorney1<2U>& attorney,
                          Region<2U>&              region,
                          Model<2U>&               model );
                      
    void TestLumpedMassMatrix( const PDE_Integrator_Attorney1<2U>& attorney,
                               const Region<2U>&                   region,
                               const Model<2U>&                    model,
                               double                              dt );
                            
    void TestSymmetry( const PDE_Integrator_Attorney1<2U>& attorney,
                       size_t                               n_free );
                     
    void TestInitialConditionRHS( const PDE_Integrator_Attorney1<2U>& attorney,
                                  const Region<2U>&                   region,
                                  const Model<2U>&                    model,
                                  double                              dt );
                                 
    void TestDirichletRHSModification( const PDE_Integrator_Attorney1<2U>& attorney,
                                       const Region<2U>&                   region,
                                       const Model<2U>&                    model );
    
     void OutVector(const std::vector<double>& vector, std::string file);

  };

} // end namespace csmp

#endif /* CSMP_PDE_INTEGRATOR_TRANSIENT_TEST_H */
