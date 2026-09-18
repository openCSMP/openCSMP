// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  PETSc_Settings.h
//  open-csmp-2024
//

#ifndef CSMP_PETSC_SETTINGS_H
#define CSMP_PETSC_SETTINGS_H

#include "CSMP_definitions.h"
#include "SolverSettings.h"

using namespace std;

namespace csmp {

/**
 @brief Defines PETSc linear solver settings, including KSP type, preconditioner type, convergence tolerances, and iteration limits. The optimal solver and preconditioner depend on matrix conditioning, mesh quality, problem size, etc.
 
 @author M. Dubnytska
 @date 2026
 
 @section overview Overview
 
 Recommended PETSc linear solver (KSP) and preconditioner (PC) combinations:
 
 | Problem / matrix type | Typical examples | Recommended KSP | Recommended PC | Comments |
 |------------------------------|-------------------------|-----------------------------|---------------------------|----------------|
 | Symmetric positive definite (SPD) | Pressure diffusion, Poisson equation, heat conduction, scalar diffusion | `"cg"` | `"gamg"` | Preferred default for large sparse SPD systems. CG requires both the matrix and preconditioner to be symmetric positive definite. |
 | SPD, very large parallel problem | Large diffusion or pressure problems on many MPI ranks | `"pipecg"` | `"gamg"` | Communication-reducing version of CG. Preferred only when communication/global reductions become a significant cost. Standard `"cg"` should normally be tried first. |
 | Symmetric indefinite | Some mixed formulations and saddle-point systems | `"minres"` | `"jacobi"` or suitable symmetric PC | MINRES is appropriate when the matrix is symmetric but not positive definite. The preconditioner must also be symmetric and positive definite. |
 | General nonsymmetric | Advection-diffusion, coupled transport, nonsymmetric discretizations | `"gmres"` | `"bjacobi"`, `"asm"` or problem-specific PC | Good general-purpose and robust choice when symmetry cannot be assumed. GMRES requires additional memory for Krylov vectors. |
 | General nonsymmetric with changing/iterative preconditioner | Nested solvers, nonlinear or variable preconditioning, FieldSplit with iterative subsolvers | `"fgmres"` | `"fieldsplit"` or other variable PC | Flexible GMRES permits the preconditioner to change between iterations. |
 | Nonsymmetric system where GMRES memory usage is too large | Large transport-type systems | `"bcgs"` | `"bjacobi"`, `"asm"` or problem-specific PC | BiCGStab uses less memory than GMRES but convergence can be less smooth and less robust. |
 | Small or moderate system / solver debugging | Any system for which a direct factorization is practical | `"preonly"` | `"lu"` | Applies an LU direct solve. Very useful for determining whether a problem is caused by matrix assembly or by iterative solver convergence. |
 | Small or moderate SPD system | Diffusion, Poisson, elasticity when the matrix is SPD | `"preonly"` | `"cholesky"` | Direct factorization specialized for symmetric positive-definite matrices. |
 | Least-squares / rectangular system | Overdetermined or underdetermined linear systems | `"lsqr"` | problem dependent | Intended for least-squares problems rather than the usual square FEM stiffness matrix. |
 | Multigrid smoother / special-purpose iteration | Multigrid levels, eigenvalue-bounded SPD operators | `"chebyshev"` | `"jacobi"` | Frequently useful as a smoother rather than as the main standalone linear solver. |
 |------------------------------|-------------------------|-----------------------------|---------------------------|----------------|
 
 @section recommendations Recommended solver settings
 
 @subsection scalar_diffusion Scalar diffusion problems
 
 For pressure diffusion, heat conduction, Poisson-type equations and similar scalar elliptic problems with positive coefficients, the assembled matrix is normally symmetric positive definite after appropriate Dirichlet boundary conditions have been imposed. A good default is:
 
 @code
 settings.SetKSPType("cg");
 settings.SetPCType("gamg");
 settings.SetRelativeTolerance(1.0e-8);
 @endcode
 
 @subsection nonsymmetric General nonsymmetric problems
 
 If the matrix is nonsymmetric, or its symmetry is not guaranteed, GMRES is a safer general-purpose starting point:
 
 @code
 settings.SetKSPType("gmres");
 settings.SetPCType("bjacobi");
 settings.SetRelativeTolerance(1.0e-8);
 @endcode
 
 For distributed-memory calculations, `"asm"` may also be considered instead of `"bjacobi"`, depending on the problem.
 
 @subsection debugging Debugging a linear system
 
 When an iterative method fails to converge, it might be useful to first test the same assembled system with a direct solver, if the problem size permits:
 
 @code
 settings.SetKSPType("preonly");
 settings.SetPCType("lu");
 @endcode
 
 If it converges with the direct solver, it is very likely that changing linear solver (KSP) and preconditioner (PC) combination will solve the problem.
 
 
 @section important_notes Some important notes
 
 `"cg"` should only be used when the matrix and preconditioner satisfy the symmetry and definiteness requirements of the conjugate-gradient method.
 
 `"minres"` permits a symmetric indefinite system matrix, but its preconditioner must remain symmetric and positive definite.
 
 `"gmres"` is usually the safest first choice for a genuinely nonsymmetric matrix.
 
 `"fgmres"` should be preferred to `"gmres"` when the preconditioner changes between iterations, for example when iterative nested solves are used as part of the preconditioner.
 
 `"gamg"` is particularly useful for large sparse elliptic problems. For vector problems such as elasticity, the matrix block structure and near-null-space information should be supplied correctly for optimal GAMG performance.
 
 `"preonly"` does not perform a Krylov iteration. It applies the selected preconditioner once and is therefore commonly combined with direct factorization preconditioners such as `"lu"` or `"cholesky"`.
 
 */
class PETSc_Settings : public SolverSettings {
public:
    PETSc_Settings() = default;
    ~PETSc_Settings() override = default;

    
    void SetKSPType( const string& value );
    void SetPCType(const string& value);
    void SetRelativeTolerance(double value);
    void SetAbsoluteTolerance(double value);
    void SetMaximumIterations(int value);
    void SetMonitorTrueResidual(bool value);
    void SetPrintConvergedReason(bool value);
    void SetUseInitialGuess(bool value);

    const string& GetKSPType() const;
    const string& GetPCType() const;
    double GetRelativeTolerance() const;
    double GetAbsoluteTolerance() const;
    int GetMaximumIterations() const;
    bool GetMonitorTrueResidual() const;
    bool GetPrintConvergedReason() const;
    bool GetUseInitialGuess() const;

private:
    
    std::string ksp_type_ = "gmres";
    std::string pc_type_ = "ilu";

    double relative_tolerance_ = 1.0e-10;
    double absolute_tolerance_ = 1.0e-50;

    int maximum_iterations_ = 1000;

    bool monitor_true_residual_ = false;
    bool print_converged_reason_ = false;
    
    bool use_initial_guess_ = false;
    
};

} // namespace csmp

#endif




