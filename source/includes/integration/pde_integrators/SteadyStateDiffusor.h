// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STEADY_STATE_DIFFUSOR_H
#define STEADY_STATE_DIFFUSOR_H

#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "SAMG_Exception.h"
#else
#include "LinearSolver.h"
#endif

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;

template<uint32_t dim,template<uint32_t> class> class NumIntegral_NT_rhsop_N_dV;
template<uint32_t dim,template<uint32_t> class> class PointSource_rhsop;
template<uint32_t dim,template<uint32_t> class> class NumIntegral_dNT_op_dV;

/** Spatial integration of steady-state partial-differential diffusion eqns: k div^2=q+..., 

    @attention Gradient multiplier is currently not used and has no influence on results.
    @todo (3) Make generic so that it can take both numerically and analytically integrated FEs
*/
template<uint32_t dim,template<uint32_t> class CELLTYPE=Element>
class SteadyStateDiffusor : public PDE_Integrator<dim,CELLTYPE> {
  public:
	  SteadyStateDiffusor(  Model<dim>&,
                          const char* diffusivity,
                          const char* diffusing_variable,
                          const char* spatial_source_variable,
                          bool LumpedRHS );

    SteadyStateDiffusor( Model<dim>&, 
                          const char* diffusivity,
                          const char* diffusing_variable,
                          const char* spatial_source_variable );
                              
    SteadyStateDiffusor( Model<dim>&, 
                          const char* diffusivity,
                          const char* diffusing_variable,
                          const char* spatial_source_variable,  // on element or integration points
                          const char* point_source_variable );  // on the nodes

    /// integrating gradient terms applied to the righthand sides; no sources
    SteadyStateDiffusor( Model<dim>&, 
                          const char* diffusivity,
                          const char* diffusing_variable,
                          const char* gradient_variable, double gradient_multiplier );

    /// as above, but with FE source terms
    SteadyStateDiffusor( Model<dim>&, 
                          const char* diffusivity,
                          const char* diffusing_variable,
                          const char* spatial_source_variable,
                          const char* gradient_variable, double gradient_multiplier );

    SteadyStateDiffusor( Model<dim>&, 
                        const char* lhs_diffusivity,
                        const char* rhs_diffusivity,
                        const char* diffusing_variable,
                        const char* spatial_source_variable,
                        const char* gradient_variable, double gradient_multiplier );

    virtual ~SteadyStateDiffusor();

    /// K div^2 = q
    void ComputeSteadyState( ModelSubDomain<dim,CELLTYPE>& subDomain, bool verbose=false );

    virtual void AdjustSolverSettings();

#ifdef CSMP_WITH_SAMG_SOLVER
    /// Access to SAMG solver settings owned by SteadyStateDiffusor object
    SAMG_Settings& GetSolverSettings();
#else
    /// add extra functionality for alternative solver if needed
#endif
    
  protected:
    NumIntegral_dNT_lhsop_dN_dV<dim,CELLTYPE> conductance_;
    NumIntegral_NT_rhsop_N_dV<dim,CELLTYPE>*  source_ = nullptr;
    PointSource_rhsop<dim,CELLTYPE>*          nodal_source_;
    NumIntegral_dNT_op_dV<dim,CELLTYPE>*      gravity_ = nullptr;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings  settings_;
    SAMG_Solver    solver_;
#else
    CSMP_DEFAULT_LINEAR_SOLVER solver_;
#endif
    const double  grad_multiplier_;
    std::string   dep_var_name_;
    
  private:
    SteadyStateDiffusor();
    SteadyStateDiffusor( const SteadyStateDiffusor& );
#ifdef CSMP_WITH_SAMG_SOLVER
    void Adjust_SAMG_ForSubsequentSolves();
#else
    /// add extra functionality for alternative solver if needed
#endif

  private:
    bool   firstCall_;
};

} // end csmp

#endif 
