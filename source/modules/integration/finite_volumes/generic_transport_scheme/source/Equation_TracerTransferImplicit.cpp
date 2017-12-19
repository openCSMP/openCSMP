//  Equation_TracerTransferImplicit.cpp
//
//  Created by Andrew Bromage on 8 Dec 2017.
//  Copyright (c) 2017 XXX. All rights reserved.
//

#include "Equation_TracerTransferImplicit.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "ImplicitTransport.h"
#include "finiteVolumeFunctions.h"
#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp {

template<size_t dim>
Equation_TracerTransferImplicit<dim>::Equation_TracerTransferImplicit( Model<dim>& model, const char* region )
    : accumulator(model,region),
      pv_op(model, "FV pore volume", "porosity"),
      flux_op(model, "facet flux"),
      // diffusion_op(model, "diffusivity"),
      pvrhs_op(model, "FV pore volume", "concentration", "porosity")
      // srclhs_op(model, "porosity", "total systems compressibility", "fluid pressure", "fluid pressure")
{
      pv_op.AccumulationMode(ADD_ACCUMULATE);
      pv_op.MultiplyWithTimeIncrement(true);
      accumulator.AddOperatorInterior(&pv_op);
      accumulator.AddOperatorPerimeter(&pv_op);

      flux_op.AccumulationMode(ADD_ACCUMULATE);
      accumulator.AddOperatorInterior(&flux_op);
      accumulator.AddOperatorPerimeter(&flux_op);

#if 0
      diffusion_op.AccumulationMode(ADD_ACCUMULATE);
      accumulator.AddOperatorInterior(&diffusion_op);
      accumulator.AddOperatorPerimeter(&diffusion_op);
#endif
      
      pvrhs_op.AccumulationMode(MULTIPLY_ACCUMULATE);
      pvrhs_op.MultiplyWithTimeIncrement(true);
      accumulator.AddOperatorInterior(&pvrhs_op);
      accumulator.AddOperatorPerimeter(&pvrhs_op);
      
#if 0
      srclhs_op.AccumulationMode(ADD_ACCUMULATE);
      accumulator.AddOperatorInterior(&srclhs_op);
      accumulator.AddOperatorPerimeter(&srclhs_op);
#endif

      accumulator.FinaliseOperators();
}


  
  template<size_t dim>
void
  Equation_TracerTransferImplicit<dim>::EquationTimeIncrement( double64 dt )
  {
    accumulator.TimeIncrement(1.0 / dt);
  }


  template struct Equation_TracerTransferImplicit<1U>;
  template struct Equation_TracerTransferImplicit<2U>;
  template struct Equation_TracerTransferImplicit<3U>;

} // end csmp
