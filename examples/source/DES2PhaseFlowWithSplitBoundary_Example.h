#ifndef DES_2PHASE_FLOW_WITH_SPLIT_BOUNDARY_EXAMPLE_H
#define DES_2PHASE_FLOW_WITH_SPLIT_BOUNDARY_EXAMPLE_H

#include "Example.h"
#include "Model.h"
#include "FlowFunctionsModule.h"

namespace csmp {

class  DES2PhaseFlowWithSplitBoundary_Example : public Example
{

public:
    virtual void Run();
    virtual void Specifications();

private:

  template<uint32_t dim>
  void RunSimulation(Model<dim>& model);

  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k );

  template<uint32_t dim>
  void ComputeSteadyStatePressure(Model<dim>& mdl, bool with_gravity, bool with_tensor_k);

};

/*
  template<uint32_t dim>
  void RunSimulation(Model<dim>& model);

  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k );

  template<uint32_t dim>
  void ComputeSteadyStatePressure(Model<dim>& mdl, bool with_gravity, bool with_tensor_k);
*/

} // csmp

#endif // DES_2PHASE_FLOW_WITH_SPLIT_BOUNDARY_EXAMPLE_H
