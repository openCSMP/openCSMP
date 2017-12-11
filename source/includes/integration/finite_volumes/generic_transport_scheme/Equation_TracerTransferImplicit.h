#ifndef EQUATION_TRACER_TRANSFER_IMPLICIT_H
#define EQUATION_TRACER_TRANSFER_IMPLICIT_H

#include "LinearSystemAccumulator.h"

#include "CompressiblePressureSourceLHS.h"
#include "PoreVolume.h"
#include "DiffusionLHS.h"
#include "PoreVolumeRHS.h"
#include "SourceTermRHS.h"
#include "FluxLHS.h"


namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Model;
template<size_t> class Region;


template<size_t dim>
struct Equation_TracerTransferImplicit {
  Equation_TracerTransferImplicit( Model<dim>& model, const char* region );
  
  Equation_TracerTransferImplicit() = delete;
  Equation_TracerTransferImplicit(const Equation_TracerTransferImplicit&) = delete;
  Equation_TracerTransferImplicit(Equation_TracerTransferImplicit&&) = delete;

    void EquationTimeIncrement( double64 dt );

  LinearSystemAccumulator<dim> accumulator;
  PoreVolume<dim> pv_op;
  FluxLHS<dim> flux_op;
  // DiffusionLHS<dim> diffusion_op;
  PoreVolumeRHS<dim> pvrhs_op;
  // CompressiblePressureSourceLHS<dim> srclhs_op;
};



} // end csmp

#endif /* EQUATION_TRACER_TRANSFER_IMPLICIT_H */
