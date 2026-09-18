// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  PETSc_Solver.cpp
//  open-csmp-2024
//

#include "PETSc_Solver.h"

#include <chrono>
#include <stdexcept>

#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"
#include "Exception.h"

using namespace std;

namespace csmp {

PETSc_Solver::PETSc_Solver() :
      Solver(new PETSc_Settings()),
      settings_(static_cast<PETSc_Settings*>(solver_settings_)),
      owns_settings_(true)
{
}


PETSc_Solver::PETSc_Solver(PETSc_Settings* settings) :
      Solver(settings),
      settings_(settings),
      owns_settings_(false)
{
    if (settings_ == nullptr) {
        throw invalid_argument(
            "PETSc_Solver: settings pointer must not be null."
        );
    }
}



PETSc_Solver::~PETSc_Solver()
{
    if (owns_settings_) {
        delete settings_;
    }

    settings_ = nullptr;
    solver_settings_ = nullptr;
}



PETSc_Settings* PETSc_Solver::GetSolverSettings()
{
    return settings_;
}


void PETSc_Solver::InputSolverSettings(SolverSettings& settings)
{
    auto* petsc_settings =
        dynamic_cast<PETSc_Settings*>(&settings);

    if (petsc_settings == nullptr) {
        throw invalid_argument(
            "PETSc_Solver requires PETSc_Settings."
        );
    }

    if (owns_settings_) {
        delete settings_;
    }

    settings_ = petsc_settings;
    solver_settings_ = petsc_settings;
    owns_settings_ = false;
}




#if defined(CSMP_WITH_PETSC_SOLVER)

void PETSc_Solver::ConfigureKSP(KSP ksp)
{
    if (settings_ == nullptr) {
        throw runtime_error(
            "PETSc_Solver::ConfigureKSP: settings are null."
        );
    }

    PC pc;
    PetscCallAbort(PETSC_COMM_SELF, KSPGetPC(ksp, &pc));

    PetscCallAbort(
        PETSC_COMM_SELF,
        KSPSetType(ksp, settings_->GetKSPType().c_str())
    );

    PetscCallAbort(
        PETSC_COMM_SELF,
        PCSetType(pc, settings_->GetPCType().c_str())
    );

    PetscCallAbort(
        PETSC_COMM_SELF,
        KSPSetTolerances(
            ksp,
            settings_->GetRelativeTolerance(),
            settings_->GetAbsoluteTolerance(),
            PETSC_DEFAULT,
            settings_->GetMaximumIterations()
        )
    );

    if (settings_->GetMonitorTrueResidual()) {
        PetscCallAbort(
            PETSC_COMM_SELF,
            PetscOptionsSetValue(
                nullptr,
                "-ksp_monitor_true_residual",
                nullptr
            )
        );
    }

    if (settings_->GetPrintConvergedReason()) {
        PetscCallAbort(
            PETSC_COMM_SELF,
            PetscOptionsSetValue(
                nullptr,
                "-ksp_converged_reason",
                nullptr
            )
        );
    }
    
    PetscCallAbort(
        PETSC_COMM_SELF,
        KSPSetInitialGuessNonzero(
            ksp,
            settings_->GetUseInitialGuess() ? PETSC_TRUE : PETSC_FALSE
        )
    );

    /*
     * Options supplied externally override settings from the profile.
     */
    PetscCallAbort(PETSC_COMM_SELF, KSPSetFromOptions(ksp));
    
    if (settings_->GetMonitorTrueResidual()) {
        PetscCallAbort(
            PETSC_COMM_SELF,
            PetscOptionsClearValue(
                nullptr,
                "-ksp_monitor_true_residual"
            )
        );
    }

    if (settings_->GetPrintConvergedReason()) {
        PetscCallAbort(
            PETSC_COMM_SELF,
            PetscOptionsClearValue(
                nullptr,
                "-ksp_converged_reason"
            )
        );
    }
    
    
}



void PETSc_Solver::CheckPETSc_Initialized() const
{
    PetscBool initialized = PETSC_FALSE;

    PetscCallAbort(
        PETSC_COMM_SELF,
        PetscInitialized(&initialized)
    );

    if (!initialized) {
        throw runtime_error(
            "PETSc_Solver requires PetscInitialize() before Solve()."
        );
    }
}




void PETSc_Solver::SolveMatrixEquation( SparseMatrix& A,
                                        vector<double>& b,
                                        vector<double>& x,
                                        size_t /*no_unknowns */ )
 {
    // Petsc MPI initialisation check
    CheckPETSc_Initialized();
    
    
    PetscInt n = static_cast<PetscInt>(A.Rows());
     
    cout << "PETSc_Solver: Size of PetscInt: " << sizeof(PetscInt) << endl;
    
    auto start = chrono::high_resolution_clock::now();

    // 1. Create PETSc Vectors
    Vec petsc_b, petsc_x;
    VecCreateSeq(PETSC_COMM_SELF, n, &petsc_b);
    VecDuplicate(petsc_b, &petsc_x);
    
    
    if(settings_->GetUseInitialGuess()) {
         
         if(x.size() != static_cast<size_t>(n)){
             throw runtime_error(
                 "PETSc_Solver::SolveMatrixEquation: initial guess has incorrect size."
             );
         }
         
         for(PetscInt i = 0; i < n; i++) {
             PetscCallAbort(PETSC_COMM_SELF, VecSetValue(petsc_x, i, static_cast<PetscScalar>(x[static_cast<size_t>(i)]), INSERT_VALUES));
         }
         
         PetscCallAbort(PETSC_COMM_SELF, VecAssemblyBegin(petsc_x));
         PetscCallAbort(PETSC_COMM_SELF, VecAssemblyEnd(petsc_x));
     
    } else {
         
         PetscCallAbort(PETSC_COMM_SELF, VecSet(petsc_x, 0.0));
         
    }
    
    
    // Populate RHS
    for (PetscInt i = 0; i < n; ++i) {
        VecSetValue(petsc_b, i, b[static_cast<size_t>(i)], INSERT_VALUES);
    }
    VecAssemblyBegin(petsc_b);
    VecAssemblyEnd(petsc_b);

    // 2. Create PETSc Matrix
    Mat petsc_A;
    
    vector<PetscInt> nnz(static_cast<size_t>(n));
    vector<PetscBool> has_diag( static_cast<size_t>(n), PETSC_FALSE );

    for (PetscInt row = 0; row < n; ++row) {
        const auto rowEnd = A.RowEnd(static_cast<size_t>(row));

        for (auto it = A.RowBegin(static_cast<size_t>(row)); it != rowEnd; ++it) {
            ++nnz[static_cast<size_t>(row)];

            PetscInt col = static_cast<PetscInt>(it->first);
            if (col == row) {
                has_diag[static_cast<size_t>(row)] = PETSC_TRUE;
            }
        }

        // Reserve space for an explicit diagonal placeholder.
        if (!has_diag[static_cast<size_t>(row)]) {
            ++nnz[static_cast<size_t>(row)];

            cerr << "PETSc diagnostic: missing diagonal in CSMP matrix row "
                      << row;

            cerr << " ; if interleaved 2D displacement, node = "
                      << row / 2 << ", component = "
                      << (row % 2 == 0 ? "x" : "y");

            cerr << endl;
        }
    }

    MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 0, nnz.data(), &petsc_A);
    MatSetOption(petsc_A, MAT_FORCE_DIAGONAL_ENTRIES, PETSC_TRUE);
    

    // Populate the matrix
    for (size_t row = 0; row < A.Rows(); ++row) {
        const auto rowEnd = A.RowEnd(row);
        for (auto it = A.RowBegin(row); it != rowEnd; ++it) {
            PetscInt col = static_cast<PetscInt>(it->first);
            PetscScalar val = static_cast<PetscScalar>(it->second);
            MatSetValue(petsc_A, static_cast<PetscInt>(row), col, val, INSERT_VALUES);
        }
    }
    
    
    for (PetscInt row = 0; row < n; ++row) {
        if (!has_diag[static_cast<size_t>(row)]) {
            MatSetValue(petsc_A, row, row, 0.0, INSERT_VALUES);
        }
    }
    
    
    MatAssemblyBegin(petsc_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A, MAT_FINAL_ASSEMBLY);

    
    
    
    // 3. Setup Solver (KSP)
    KSP ksp;

    PetscCallAbort(PETSC_COMM_SELF, KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCallAbort(PETSC_COMM_SELF, KSPSetOperators(ksp, petsc_A, petsc_A));

    ConfigureKSP(ksp);
    
    

    // 4. Solve
    cout << "\nPETSc_Solver::SolveMatrixEquation: solving " << A.Rows() << " x " << A.Cols() << " system..." << endl;
     
    //cleanup
    auto cleanup = [&]() {
        KSPDestroy(&ksp);
        MatDestroy(&petsc_A);
        VecDestroy(&petsc_b);
        VecDestroy(&petsc_x);
    };
     
    
    // PETSc API error
    const PetscErrorCode error_code = KSPSolve(ksp, petsc_b, petsc_x);
    
    if (error_code != PETSC_SUCCESS) {
        const auto error = describe(static_cast<int>(error_code));

        cerr << "\tERROR: KSPSolve failed with PETSc error code "
            << static_cast<int>(error_code) << endl;
        printError(static_cast<int>(error_code));

        // If this error was caused by non-convergence (for example when -ksp_error_if_not_converged is enabled), print the KSP reason too.
        KSPConvergedReason error_reason = KSP_CONVERGED_ITERATING;
        if (KSPGetConvergedReason(ksp, &error_reason) == PETSC_SUCCESS && error_reason <= 0) {
            printKSPReason(static_cast<int>(error_reason));
        }

        cleanup();

        throw runtime_error(
            string("PETSc_Solver(SparseMatrix): KSPSolve failed: ") +
            string(error.symbol) + " (" + to_string(error.code) + ") - " +
            string(error.description)
        );
    }
    

    auto stop = chrono::high_resolution_clock::now();

    // --- Convergence & Statistics Reporting ---
    PetscInt its;
    PetscReal res_norm;
    PetscReal res_inf_norm;
    KSPConvergedReason reason;

    KSPGetIterationNumber(ksp, &its);
    KSPGetResidualNorm(ksp, &res_norm);
    KSPGetConvergedReason(ksp, &reason);
    
    // Positive reason = converged; zero/negative = no converged solution.
    if (reason <= 0) {
        const auto reason_info = describeKSPReason(static_cast<int>(reason));

        cerr << "\tERROR: PETSc_Solver(SparseMatrix) did not converge." << endl;
        printKSPReason(static_cast<int>(reason));

        cleanup();

        throw runtime_error(
            string("PETSc_Solver(SparseMatrix): solver did not converge: ") +
            string(reason_info.symbol) + " (" + to_string(reason_info.code) + ") - " +
            string(reason_info.description)
        );
    }
     
     
    // --- Calculate L-infinity Norm only for a converged solution ---
    Vec r;
    VecDuplicate(petsc_b, &r);
    MatMult(petsc_A, petsc_x, r);               // r = A * x
    VecAYPX(r, -1.0, petsc_b);                  // r = b - Ax
    VecNorm(r, NORM_INFINITY, &res_inf_norm);   // L-infinity norm of r
    VecDestroy(&r);

    cout << "\t#iterations:      " << its << endl;
    cout << "\tL-infinity norm: " << res_inf_norm << endl;

    if (settings_->GetPrintConvergedReason()) {
        printKSPReason(static_cast<int>(reason));
    }

    cout << "\tcompleted solution in " << chrono::duration_cast<chrono::seconds>(stop - start).count()
        << " seconds." << endl;

    // 5. Copy solution back to vector x
    PetscScalar *array_x;
    VecGetArray(petsc_x, &array_x);
    for ( size_t i{0}; i < static_cast<size_t>(n); ++i ) {
        x[i] = static_cast<double>( array_x[i] );
    }
    VecRestoreArray(petsc_x, &array_x);

    cleanup();
     
}




void PETSc_Solver::SolveMatrixEquation( CompressedRowMatrix& A,
                                        vector<double>& b,
                                        vector<double>& x,
                                        size_t /* no_unknowns */ )
{
    // Petsc MPI initialisation check
    CheckPETSc_Initialized();
    
    
    auto n = static_cast<PetscInt>(A.Rows());
    auto start = chrono::high_resolution_clock::now(); // Start Timer

    // 1. Create PETSc Vectors
    Vec petsc_b, petsc_x;
    VecCreateSeq(PETSC_COMM_SELF, n, &petsc_b);
    VecDuplicate(petsc_b, &petsc_x);
    
    if(settings_->GetUseInitialGuess()) {
        
        if(x.size() != static_cast<size_t>(n)){
            throw runtime_error(
                "PETSc_Solver::SolveMatrixEquation: initial guess has incorrect size."
            );
        }
        
        for(PetscInt i = 0; i < n; i++) {
            PetscCallAbort(PETSC_COMM_SELF, VecSetValue(petsc_x, i, static_cast<PetscScalar>(x[static_cast<size_t>(i)]), INSERT_VALUES));
        }
        
        PetscCallAbort(PETSC_COMM_SELF, VecAssemblyBegin(petsc_x));
        PetscCallAbort(PETSC_COMM_SELF, VecAssemblyEnd(petsc_x));
    
    } else {
        
        PetscCallAbort(PETSC_COMM_SELF, VecSet(petsc_x, 0.0));
        
    }
    
    
    // Populate RHS
    for (PetscInt i = 0; i < n; ++i) {
        VecSetValue(petsc_b, i, b[static_cast<size_t>(i)], INSERT_VALUES);
    }
    VecAssemblyBegin(petsc_b);
    VecAssemblyEnd(petsc_b);

    // 2. Create PETSc Matrix
    Mat petsc_A;

    const bool samg_format = A.IsFormattedForSAMG();

    auto rowBegin = [&A, samg_format](PetscInt row) -> PetscInt {
        return samg_format
            ? static_cast<PetscInt>(A.ia[static_cast<size_t>(row)] - 1)
            : static_cast<PetscInt>(A.ia[static_cast<size_t>(row)]);
    };

    auto rowEnd = [&A, samg_format](PetscInt row) -> PetscInt {
        return samg_format
            ? static_cast<PetscInt>(A.ia[static_cast<size_t>(row + 1)] - 1)
            : static_cast<PetscInt>(A.ia[static_cast<size_t>(row + 1)]);
    };

    auto colIndex = [samg_format](int32_t raw_col) -> PetscInt {
        return samg_format
            ? static_cast<PetscInt>(raw_col - 1)
            : static_cast<PetscInt>(raw_col);
    };

    // Preallocate, but add one slot for rows that are missing their diagonal.
    // This avoids PETSc ILU failing with "Matrix is missing diagonal entry ...".
    vector<PetscInt> nnz(static_cast<size_t>(n), 0);
    vector<PetscBool> has_diag(static_cast<size_t>(n), PETSC_FALSE);

    for (PetscInt row = 0; row < n; ++row) {
        const PetscInt begin = rowBegin(row);
        const PetscInt end   = rowEnd(row);

        if (begin < 0 || end < begin || end > static_cast<PetscInt>(A.ja.size())) {
            throw runtime_error(
                "PETSc_Solver::SolveMatrixEquation(CompressedRowMatrix): invalid CRS row pointer."
            );
        }

        for (PetscInt index = begin; index < end; ++index) {
            const PetscInt col = colIndex(A.ja[static_cast<size_t>(index)]);

            if (col < 0 || col >= n) {
                throw runtime_error(
                    "PETSc_Solver::SolveMatrixEquation(CompressedRowMatrix): invalid CRS column index."
                );
            }

            ++nnz[static_cast<size_t>(row)];

            if (col == row) {
                has_diag[static_cast<size_t>(row)] = PETSC_TRUE;
            }
        }

        if (!has_diag[static_cast<size_t>(row)]) {
            ++nnz[static_cast<size_t>(row)];

            cerr << "PETSc diagnostic: missing diagonal in CompressedRowMatrix row "
                      << row
                      << " ; if this is 2D interleaved displacement, node = "
                      << row / 2
                      << ", component = "
                      << ((row % 2 == 0) ? "x" : "y")
                      << endl;
        }
    }

    PetscCallAbort(PETSC_COMM_SELF, MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 0, nnz.data(), &petsc_A));

    // Ask PETSc to keep diagonal locations structurally present.
    PetscCallAbort(PETSC_COMM_SELF, MatSetOption(petsc_A, MAT_FORCE_DIAGONAL_ENTRIES, PETSC_TRUE));

    // Insert existing CRS entries.
    for (PetscInt row = 0; row < n; ++row) {
        const PetscInt begin = rowBegin(row);
        const PetscInt end   = rowEnd(row);

        for (PetscInt index = begin; index < end; ++index) {
            const PetscInt col = colIndex(A.ja[static_cast<size_t>(index)]);
            const PetscScalar val =
                static_cast<PetscScalar>(A.a[static_cast<size_t>(index)]);

            PetscCallAbort(PETSC_COMM_SELF, MatSetValue(petsc_A, row, col, val, INSERT_VALUES));
        }
    }

    // Insert explicit zero diagonals where the CRS structure lacks them.
    // This fixes the PETSc ILU structural error. If the real system is singular,
    // PETSc may still later report a zero pivot.
    for (PetscInt row = 0; row < n; ++row) {
        if (!has_diag[static_cast<size_t>(row)]) {
            PetscCallAbort(PETSC_COMM_SELF, MatSetValue(petsc_A, row, row, 0.0, INSERT_VALUES));
        }
    }

    PetscCallAbort(PETSC_COMM_SELF, MatAssemblyBegin(petsc_A, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(PETSC_COMM_SELF, MatAssemblyEnd(petsc_A, MAT_FINAL_ASSEMBLY));
    

    // 3. Setup Solver (KSP)
    KSP ksp;

    PetscCallAbort(PETSC_COMM_SELF, KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCallAbort(PETSC_COMM_SELF, KSPSetOperators(ksp, petsc_A, petsc_A));

    ConfigureKSP(ksp);

    // 4. Solve
    cout << "\nPETSc_Solver::SolveMatrixEquation: solving " << A.Rows() << " x " << A.Cols() << " system..." << endl;
    
    auto cleanup = [&]() {
        KSPDestroy(&ksp);
        MatDestroy(&petsc_A);
        VecDestroy(&petsc_b);
        VecDestroy(&petsc_x);
    };
    
    const PetscErrorCode error_code = KSPSolve(ksp, petsc_b, petsc_x);
    if (error_code != PETSC_SUCCESS) {
        const auto error = describe(static_cast<int>(error_code));

        cerr << "\tERROR: KSPSolve failed with PETSc error code " << static_cast<int>(error_code) << endl;
        printError(static_cast<int>(error_code));

        KSPConvergedReason error_reason = KSP_CONVERGED_ITERATING;
        if (KSPGetConvergedReason(ksp, &error_reason) == PETSC_SUCCESS && error_reason <= 0) {
            printKSPReason(static_cast<int>(error_reason));
        }

        cleanup();

        throw runtime_error(
            string("PETSc_Solver(CompressedRowMatrix): KSPSolve failed: ") +
            string(error.symbol) + " (" + to_string(error.code) + ") - " +
            string(error.description)
        );
    }

    auto stop = chrono::high_resolution_clock::now(); // Stop Timer

    // --- Convergence & Statistics Reporting ---
    PetscInt its;
    PetscReal res_inf_norm;
    KSPConvergedReason reason;

    KSPGetIterationNumber(ksp, &its);
    KSPGetConvergedReason(ksp, &reason);

    if (reason <= 0) {
        const auto reason_info = describeKSPReason(static_cast<int>(reason));

        cerr << "\tERROR: PETSc_Solver(CompressedRowMatrix) did not converge." << endl;
        printKSPReason(static_cast<int>(reason));

        cleanup();

        throw runtime_error(
            string("PETSc_Solver(CompressedRowMatrix): solver did not converge: ") +
            string(reason_info.symbol) + " (" + to_string(reason_info.code) + ") - " +
            string(reason_info.description)
        );
    }
    
    Vec r;
    VecDuplicate(petsc_b, &r);
    MatMult(petsc_A, petsc_x, r);        // r = A*x
    VecAYPX(r, -1.0, petsc_b);           // r = b - A*x
    VecNorm(r, NORM_INFINITY, &res_inf_norm);
    VecDestroy(&r);

    cout << "\t#iterations:      " << its << endl;
    cout << "\tL-infinity norm: " << res_inf_norm << endl;

    if (settings_->GetPrintConvergedReason()) {
        printKSPReason(static_cast<int>(reason));
    }

    cout << "\tcompleted solution in " << chrono::duration_cast<chrono::seconds>(stop - start).count()
        << " seconds." << endl;

    // 5. Copy solution back
    PetscScalar *array_x;
    VecGetArray(petsc_x, &array_x);
    for ( size_t i = 0; i < static_cast<size_t>(n); ++i ) x[i] = array_x[i];
    VecRestoreArray(petsc_x, &array_x);

    cleanup();
}

#endif

} // end namespace csmp
