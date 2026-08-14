#ifndef PRESSURE_SOLVER_H
#define PRESSURE_SOLVER_H

#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "PointSource_rhsop.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

namespace csmp {

template<uint32_t> class Model;

template<uint32_t dim, template<uint32_t> class CELL>
class PressureSolver : public PDE_Integrator<dim,CELL>
  {
  public:
      PressureSolver( Model<dim>&, bool with_gravity, bool with_capillary );
      virtual ~PressureSolver();

#ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings& PressureSolverSettings();
#endif

  protected:
#ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings  settings_;
      SAMG_Solver    solver_;
#else
      CSMP_DEFAULT_LINEAR_SOLVER  solver_;
#endif
      NumIntegral_dNT_lhsop_dN_dV<dim,CELL>  total_mobility_;
      NumIntegral_dNT_op_dV<dim,CELL>        gravity_term_,
                                             capillary_term_;
      PointSource_rhsop<dim,CELL>            nodal_source_;
  };

} // end csmp

#endif
