// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_SOLVER_CHOICE_H
#define CVFEM_SOLVER_CHOICE_H

/*  One place to pick the linear solver backend.

    All backends derive from csmp::Solver, which declares BOTH
    SolveMatrixEquation overloads (SparseMatrix and CompressedRowMatrix) as
    pure virtual — so every backend supports both matrix types and the choice
    is free at every call site, whatever the PDE_Integrator is instantiated
    with.

    Usage:
        SolverBundle bundle( SolverKind::PETSc );
        PDE_Integrator<dim, Element, SparseMatrix> P_FE( bundle.Get() );

    One bundle per system (pressure, temperature, well), NOT one shared —
    each system keeps its own solver object, as P_solver/T_solver did before.

    Note: SAMG's multiple-instance mechanism (SAMG_MULTIPLE_INSTANCES, the
    per-instance coarsening setup selected via SAMG_Settings) is legacy and is
    deliberately NOT used here. Every bundle leaves the instance at its
    default.

    Eigen is a DIRECT SparseLU solve. It suits the well Jacobian (small,
    dense-ish, nonsymmetric) and serves as a reference answer when an
    iterative solver is suspected of lying. It is NOT suitable for the
    reservoir: it refactorises every step and warns above 10k DOF.
    RequireReservoirCapable() enforces that.

    Benoit 14/09/2026
*/

#include <stdexcept>
#include <string>

#include "Solver.h"
#include "EigenSolver.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif

#if defined(CSMP_WITH_PETSC_SOLVER)
#include "PETSc_Solver.h"
#endif

namespace csmp {

enum class SolverKind { SAMG, PETSc, Eigen };

inline const char* SolverKindName( SolverKind k ) noexcept
{
    switch ( k ) {
        case SolverKind::SAMG:  return "SAMG";
        case SolverKind::PETSc: return "PETSc";
        case SolverKind::Eigen: return "Eigen";
    }
    return "unknown";
}

/// True if the backend was compiled into this build.
inline bool SolverKindAvailable( SolverKind k ) noexcept
{
    switch ( k ) {
        case SolverKind::SAMG:
#ifdef CSMP_WITH_SAMG_SOLVER
            return true;
#else
            return false;
#endif
        case SolverKind::PETSc:
#if defined(CSMP_WITH_PETSC_SOLVER)
            return true;
#else
            return false;
#endif
        case SolverKind::Eigen:
            return true;   // always built
    }
    return false;
}

/// Throws if the kind is unsuitable for a reservoir-sized system.
inline void RequireReservoirCapable( SolverKind k )
{
    if ( k == SolverKind::Eigen )
        throw std::runtime_error(
            "SolverBundle: Eigen is a direct SparseLU solve and is not intended "
            "for the reservoir system; use SAMG or PETSc." );
}

/** Owns one instance of each available backend and hands out the selected one
    as a Solver&. Settings members are declared BEFORE their solvers on
    purpose: members initialise in declaration order and each solver's
    constructor takes a pointer to its settings. */
class SolverBundle {
public:
    explicit SolverBundle( SolverKind kind )
        : kind_( kind )
#ifdef CSMP_WITH_SAMG_SOLVER
        , samg_( &samg_settings_ )
#endif
#if defined(CSMP_WITH_PETSC_SOLVER)
        , petsc_( &petsc_settings_ )
#endif
    {
        if ( !SolverKindAvailable( kind_ ) )
            throw std::runtime_error(
                std::string( "SolverBundle: backend '" ) + SolverKindName( kind_ ) +
                "' was not compiled into this build." );
    }

    SolverBundle( const SolverBundle& )            = delete;
    SolverBundle& operator=( const SolverBundle& ) = delete;

    /// The selected backend, for a PDE_Integrator or a direct Solve() call.
    Solver& Get()
    {
        switch ( kind_ ) {
#ifdef CSMP_WITH_SAMG_SOLVER
            case SolverKind::SAMG:  return samg_;
#endif
#if defined(CSMP_WITH_PETSC_SOLVER)
            case SolverKind::PETSc: return petsc_;
#endif
            case SolverKind::Eigen: return eigen_;
            default: break;
        }
        throw std::runtime_error( "SolverBundle::Get: backend not available" );
    }

    SolverKind  Kind() const noexcept { return kind_; }
    const char* Name() const noexcept { return SolverKindName( kind_ ); }

    /** Direct factorisation settings — sensible for a Newton Jacobian.
        No-op unless the selected backend is PETSc. */
    void UseDirectFactorisation()
    {
#if defined(CSMP_WITH_PETSC_SOLVER)
        if ( kind_ == SolverKind::PETSc ) {
            petsc_settings_.SetKSPType( "preonly" );
            petsc_settings_.SetPCType( "lu" );
        }
#endif
    }

    // ---- backend-specific tuning; only valid for the matching backend ----
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings&  SAMGSettings()  { return samg_settings_; }
#endif
#if defined(CSMP_WITH_PETSC_SOLVER)
    PETSc_Settings& PETScSettings() { return petsc_settings_; }
#endif

private:
    SolverKind      kind_;

#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings   samg_settings_;   // declared before samg_ — init order matters
    SAMG_Solver     samg_;
#endif

#if defined(CSMP_WITH_PETSC_SOLVER)
    PETSc_Settings  petsc_settings_;  // declared before petsc_
    PETSc_Solver    petsc_;
#endif

    EigenSolver     eigen_;           // no settings object
};

} // namespace csmp

#endif /* CVFEM_SOLVER_CHOICE_H */
