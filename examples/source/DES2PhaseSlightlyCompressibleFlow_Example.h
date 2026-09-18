// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DES_2PHASE_SLIGHTLY_COMPRESSBILE_FLOW_EXAMPLE_H
#define DES_2PHASE_SLIGHTLY_COMPRESSBILE_FLOW_EXAMPLE_H

#include "Example.h"

#include "Model.h"
#include "TwoPhaseModel.h"
#include "FlowFunctionsModule.h"

namespace csmp {

class  DES2PhaseSlightlyCompressibleFlow_Example : public Example
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

} // csmp

#endif // DES_2PHASE_SLIGHTLY_COMPRESSBILE_FLOW_EXAMPLE_H
