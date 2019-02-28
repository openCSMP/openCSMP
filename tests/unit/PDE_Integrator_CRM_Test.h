#ifndef PDE_INTEGRATOR_CRM_TEST_H
#define PDE_INTEGRATOR_CRM_TEST_H

#include "CSMP_definitions.h"

#include "VelocityAndVolumeFlux.h"

#include "Test.h"
#include "Model.h"
#include "VSet.h"

#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "PDE_Integrator_CRM.h"
#include "Integral_dNT_op_dN_dV.h"
#include "Integral_NT_op_N_dV.h"
#include "CSMP_highLevelUtilities.h"

#include "LinearSolver.h"

#include "PropertyHandle.h"
namespace csmp {

  class PDE_Integrator_CRM_Test : public Test {
  public:
    PDE_Integrator_CRM_Test();
    ~PDE_Integrator_CRM_Test();
    void run();

    void sameSolverTest();

    void constructorTest();
    void copyTest();

    void AddOperatorLHSTest(); // Add(MathOperatorLHS<dim>*);
    void AddOperatorRHSTest(); // Add(MathOperatorRHS<dim>*);
    void AddPostProcessTest();  //AddPostProcess(MathOperatorLHS<dim>*);
    void AddBoundaryIntegralsTest();	//  AddBoundaryIntegrals(MathOperatorRHS<dim>*);

    void TimeIncrementTest();           //TimeIncrement(double64 dt);
    void TransientTest();               //Transient() const;

    /// for a particular subregion of the model
    //void //IntegrateOver(INTEGRATION_DOMAIN<dim>&, bool debug = false);
    void  EstablishMatrixSetupTest();

    void  IntegrateOverTest();// IntegrateOver(Model<dim>& model, INTEGRATION_DOMAIN<dim>& domain, bool debug = false);

    bool  IdentifyShareBoundariesTest();// IdentifySharedBoundaries(const Model<dim>& model, const INTEGRATION_DOMAIN<dim>& subdomain, std::list<std::string>& shared_boundaries);

    void  AccumulateBoundaryIntegralsTest();// AccumulateBoundaryIntegrals(const INTEGRATION_DOMAIN<dim>& comp_domain, const Boundary<dim>& boundary);

    void  LateAccumulateBoundaryIntegralsTest(); //LateAccumulateBoundaryIntegrals(const INTEGRATION_DOMAIN<dim>& comp_domain, const Boundary<dim>& boundary);

    void  SetSolverTest();//        SetSolver(Solver* new_solver);
    void  GetSolverTest();





  private:

    Model<2U>* sg_;    
    PDE_Integrator_CRM<2U, Region>* alg_;
    Integral_dNT_op_dN_dV<2U, Element<2U> >* stiff_;
    Integral_NT_op_N_dV<2U, Element<2U> >* source_;

    //NumIntegral_dNT_op_dN_dV<2

    VelocityAndVolumeFlux<2U, Element<2U> >* fluid_velocity;
  };

} // end namespace csmp

#endif
