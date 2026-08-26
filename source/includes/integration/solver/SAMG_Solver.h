#ifndef SAMG_SOLVER_H
#define SAMG_SOLVER_H

#include "Solver.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

namespace csmp {

class SAMG_Settings;

/**
 * @brief Algebraic multigrid solver wrapper for SAMG (SCAI / Fraunhofer, Germany).
 *
 * Solves the linear system @f$ A\,\mathbf{x} = \mathbf{b} @f$ using the SAMG
 * algebraic multigrid library obtainable from SCAPOS.
 *
 * ## Multiple solver instances
 *
 * When compiled with `SAMG_MULTIPLE_INSTANCES` and `LEGACY_SAMG`, up to six
 * named instances (`SAMG` … `SAMG5`) are available.  When compiled with
 * `SAMG_MULTIPLE_INSTANCES` but *without* `LEGACY_SAMG`, the modern
 * context-based API (`SAMG_CTX`, `SAMG_*_CTX`) is used instead, and the
 * instance index is passed as a string context — no explicit per-instance
 * headers are required.
 *
 * ## Ownership of SAMG_Settings
 *
 * | Constructor | Owns settings? |
 * |---|---|
 * | `SAMG_Solver()` (default) | **Yes** — allocates a default `SAMG_Settings` internally; destructor deletes it. |
 * | `SAMG_Solver(SAMG_Settings*)` | **No** — caller retains ownership. |
 * | `SAMG_Solver(const SAMG_Solver&)` | **Yes** — deep-copies the settings. |
 * | After `InputSolverSettings()` | **No** — caller retains ownership of the new object; any previously owned object is deleted first. |
 *
 * ## Stack vs. heap instantiation
 *
 * Prefer stack (automatic storage) instantiation when possible.  Dynamic
 * allocation combined with OpenMP first-touch initialisation has historically
 * caused invalid memory accesses on some Linux systems due to NUMA page
 * placement.
 *
 * @code
 * SAMG_Settings settings;
 * settings.Set_napproach(2);
 * SAMG_Solver solver(&settings);
 * PDE_Integrator<DIM, Region> pressure_diffusion(&solver);
 * @endcode
 *
 * @todo Verify that resetting `res_in_` to -1 before each solve is correct
 *       when a relative convergence criterion is active.
 * @todo Investigate runtime assertion when solving systems with a vector
 *       unknown on MSVC debug builds.
 */
class SAMG_Solver : public Solver
{
public:
    // ------------------------------------------------------------------
    // Construction / destruction
    // ------------------------------------------------------------------

    /**
     * @brief Default constructor.
     *
     * Allocates a default-constructed `SAMG_Settings` object internally.
     * The destructor will delete it.
     */
    SAMG_Solver();

    /**
     * @brief Settings-pointer constructor.
     *
     * @param settings  Non-owning pointer to an externally managed
     *                  `SAMG_Settings` instance.  The caller is responsible
     *                  for ensuring the object outlives this solver.
     */
    explicit SAMG_Solver(SAMG_Settings* settings);

    /// Copy constructor — performs a deep copy of the settings object.
    SAMG_Solver(const SAMG_Solver&);

    /// Copy-assignment operator — performs a deep copy of the settings object.
    SAMG_Solver& operator=(const SAMG_Solver&);

    /**
     * @brief Destructor.
     *
     * Deletes the `SAMG_Settings` object only if this instance owns it
     * (i.e. it was default- or copy-constructed, or the settings were
     * replaced via `InputSolverSettings`).
     */
    ~SAMG_Solver() override;

    // ------------------------------------------------------------------
    // Solver interface overrides
    // ------------------------------------------------------------------

    /**
     * @brief Replace the active solver settings.
     *
     * If this solver currently owns its settings object it is deleted before
     * the replacement takes place.  After this call the solver does **not**
     * own the new object — the caller retains responsibility for its lifetime.
     *
     * @param settings  Reference to the new settings (must be a `SAMG_Settings`).
     */
    void InputSolverSettings(SolverSettings& settings) override final;

    /// Returns a pointer to the active solver settings.
    SolverSettings* GetSolverSettings() override final;

    /// Returns the human-readable solver name.
    std::string Name() const override { return "SAMG_Solver"; }

    // ------------------------------------------------------------------
    // Diagnostics and file output
    // ------------------------------------------------------------------

    /**
     * @brief Write SAMG input data to a set of text files.
     *
     * Generates the following files (prefix = @p filename):
     * - `*.frm`  — format descriptor
     * - `*.amg`  — compressed row matrix (ia, ja, a)
     * - `*.rhs`  — right-hand side vector f
     * - `*.lhs`  — initial / final solution vector u
     * - `*.a`    — matrix entries with column indices (human-readable)
     * - `*.iu`   — unknown-type array (only for coupled systems)
     * - `*.ip`   — point-index array (only for point-based coupled systems)
     *
     * @param filename  Base name for the output files (no extension).
     * @return `true` if all files were written successfully.
     */
    bool Write_SAMG_TextInputFile(const char* filename) const;

    /**
     * @brief Enable or disable automatic text-file output after each solve.
     *
     * When enabled, `Write_SAMG_TextInputFile("csp_solution_data")` is called
     * automatically at the end of every `SolveMatrixEquation` call.
     *
     * @param write  `true` to enable, `false` (default) to disable.
     */
    void Write_SAMG_TextOutput(bool write);

    /// Returns the L2 residual of the last completed solve (@f$ r_\text{out} @f$).
    double LastSolverResidual() const;

    /// Writes the current u and f vectors to `single_output.txt`.
    void OutputVectors() const;

protected:
    // ------------------------------------------------------------------
    // SolveMatrixEquation overloads
    // ------------------------------------------------------------------

    /**
     * @brief Primary solve path — assembles from a `SparseMatrix`.
     *
     * The sparse matrix is converted to compressed-row storage internally
     * before being passed to SAMG.
     *
     * @param A            System matrix.
     * @param b            Right-hand side vector.
     * @param x            On entry: initial guess.  On exit: solution.
     * @param no_unknowns  Number of coupled unknowns per node (1 for scalar).
     */
    void SolveMatrixEquation(SparseMatrix& A,
                             std::vector<double>& b,
                             std::vector<double>& x,
                             size_t no_unknowns) override final;

    /**
     * @brief Legacy solve path — accepts a pre-built `CompressedRowMatrix`.
     *
     * Retained for backwards compatibility with code that assembles the
     * compressed-row matrix externally.  Prefer the `SparseMatrix` overload
     * for new code.
     *
     * @param A            Pre-built compressed-row matrix (will be converted
     *                     to SAMG format in-place via `ConvertToSAMGFormat()`).
     * @param b            Right-hand side vector.
     * @param x            On entry: initial guess.  On exit: solution.
     * @param no_unknowns  Number of coupled unknowns per node.
     */
    void SolveMatrixEquation(CompressedRowMatrix& A,
                             std::vector<double>& b,
                             std::vector<double>& x,
                             size_t no_unknowns) override final;

private:
    // ------------------------------------------------------------------
    // Private helpers — shared by both SolveMatrixEquation overloads
    // ------------------------------------------------------------------

    /**
     * @brief Rebuild iu, ip, iscale arrays and update nnu_/nna_/nsys_.
     *
     * Arrays are only reallocated when the problem size or number of
     * unknowns has changed since the last call, avoiding unnecessary work
     * in time-stepping loops.
     */
    void PrepareSystemArrays(int32_t new_nnu,
                             int32_t new_nna,
                             int32_t new_nsys);

    /**
     * @brief Copy x and b into the internal u_ and f_ vectors.
     *
     * Applies the point-based index permutation when
     * `settings_->UsePointBasedApproach()` is true.
     */
    void InitialiseSolutionVectors(const std::vector<double>& x,
                                   const std::vector<double>& b);

    /// Checks that crmat_, u_, and f_ are non-empty; throws on failure.
    void ValidateVectors() const;

    /**
     * @brief Retrieve settings, configure SAMG hidden parameters, and
     *        dispatch the solve.  Also performs convergence checking.
     */
    void CallSAMG();

    /**
     * @brief Copy the internal solution vector u_ back into x.
     *
     * Applies the inverse point-based permutation when required.
     */
    void CopySolutionBack(std::vector<double>& x) const;

    // ------------------------------------------------------------------
    // Private helpers — SAMG parameter configuration
    // ------------------------------------------------------------------

    /**
     * @brief Set SAMG hidden parameters for the active instance.
     *
     * Dispatches to the appropriate `SAMG_SET_*` / `SAMG_SET_*_CTX` macros
     * depending on the compile-time configuration:
     * - `!SAMG_MULTIPLE_INSTANCES`  — single instance, no suffix
     * - `SAMG_MULTIPLE_INSTANCES && LEGACY_SAMG` — explicit `SAMGn_*` calls
     * - `SAMG_MULTIPLE_INSTANCES && !LEGACY_SAMG` — context-based `*_CTX` calls
     */
    void SetSAMGHiddenParameters(
            int32_t ncg,
            int32_t levelx,
            int32_t clsolver_finest,
            int32_t ioform,
            int32_t ioform_length,
            int*    filnam_dump,
            int32_t filnam_dump_length,
            int32_t mode_mess,
            int32_t nrd,
            int32_t nru
#if defined(_OPENMP) && !defined(LEGACY_SAMG)
            ,
            int32_t icolor_omp,
            int32_t iordered_omp,
            int32_t irestriction_openmp,
            int32_t samg_omp_num_threads_external
#endif
    );

    /// Reset SAMG hidden parameters to library defaults for the active instance.
    void ResetSAMGHiddenParameters();

    /**
     * @brief Dispatch the actual `SAMG` / `SAMGn` / `SAMG_CTX` call.
     *
     * All parameters are passed by reference so that SAMG can modify
     * `iswtch` in-place (e.g. when the cycle criterion forces a new setup).
     */
    void DispatchSAMGSolve(int32_t& nsolve, int32_t& ifirst, double& eps,
                           int32_t& ncyc,   int32_t& iswtch,
                           double& a_cmplx, double& g_cmplx,
                           double& p_cmplx, double& w_avrge,
                           double& chktol,  int32_t& idump,
                           int32_t& iout,   int32_t& matrix);

    // ------------------------------------------------------------------
    // Private helpers — convergence and cycle criteria
    // ------------------------------------------------------------------

    /**
     * @brief Adjust `iswtch` to force a new coarsening setup if the number
     *        of iteration cycles has grown too large.
     *
     * Only active when `NO_PRIMARY_SOLVER_CONTROL` is defined.
     * Intended for IMPES time-stepping loops where reusing the coarsening
     * setup is beneficial but must be validated each step.
     */
    void CheckCycleCriterion(int32_t& iswtch) const;

    /**
     * @brief Update the stored best-cycle count after a successful solve.
     *
     * Only active when `NO_PRIMARY_SOLVER_CONTROL` is defined.
     */
    void UpdateCycleCriterion(int32_t iswtch);

    /// Calls `CheckCycleCriterion` only when `NO_PRIMARY_SOLVER_CONTROL` is defined.
    void CheckCycleCriterionIfNeeded(int32_t& iswtch) const;

    /// Calls `UpdateCycleCriterion` only when `NO_PRIMARY_SOLVER_CONTROL` is defined.
    void UpdateCycleCriterionIfNeeded(int32_t iswtch);

    /**
     * @brief Verify that the solver residual satisfies the convergence criterion.
     *
     * - Negative `eps`: absolute convergence — checks @f$ r_\text{out} \leq |eps| @f$.
     * - Positive `eps`: relative convergence — checks
     *   @f$ r_\text{out} / r_\text{in} \leq eps @f$.
     * - Zero `eps`: always returns `true` (no check performed).
     *
     * Issues a `WARNING` via `ErrorHandler` if the criterion is not met.
     *
     * @return `true` if converged (or `eps == 0`), `false` otherwise.
     */
    bool CheckConvergence(double eps) const;

    /**
     * @brief Check whether the matrix is too dense for single-level
     *        (transport-equation) mode.
     *
     * If `levelx == 1` and the ratio of non-zeros to unknowns exceeds 1.5,
     * `levelx` is reset to 25 (the SAMG default) and a warning is printed.
     * Only relevant when `RENOUNCE_COARSENING` is defined.
     */
    void CheckSparsityCriterion(int32_t& levelx) const;

    // ------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------

    SAMG_Settings* settings_ = nullptr; ///< Active solver settings (may or may not be owned).

    /// @name Problem dimensions
    /// Updated at the start of each `SolveMatrixEquation` call.
    ///@{
    int32_t nsys_ = 0;  ///< Number of coupled unknowns per node (1 = scalar).
    int32_t nnu_  = 0;  ///< Total number of degrees of freedom (matrix dimension).
    int32_t nna_  = 0;  ///< Total number of non-zero matrix entries.
    ///@}

    /// @name Solution and right-hand side vectors
    ///@{
    std::vector<double>  u_; ///< Initial guess on entry; solution on exit.
    std::vector<double>  f_; ///< Right-hand side vector.
    ///@}

    /// @name Coupled-system arrays
    /// Only populated when `nsys_ > 1`.
    ///@{
    std::vector<int32_t> iscale_; ///< Scaling flags per unknown (0 = no scaling).
    std::vector<int32_t> iu_;     ///< Unknown-type array: iu[i] in [1, nsys_].
    int32_t              ndiu_;   ///< Length of iu_ passed to SAMG (1 or nnu_).
    std::vector<int32_t> ip_;     ///< Node-index array (point-based approach only).
    int32_t              ndip_;   ///< Length of ip_ passed to SAMG (1 or nnu_).
    ///@}

    CompressedRowMatrix crmat_; ///< Internal compressed-row storage passed to SAMG.

    /// @name SAMG output parameters
    ///@{
    double  res_in_;      ///< Residual of the initial guess (@f$ r_\text{in} @f$).
    double  res_out_;     ///< Residual of the final approximation (@f$ r_\text{out} @f$).
    int32_t ncyc_done_;   ///< Number of multigrid cycles performed in the last solve.
    int32_t ncyc_best_;   ///< Lowest cycle count achieved across all solves (used by cycle criterion).
    int32_t ierr_;        ///< SAMG error code: 0 = success, >0 = fatal, <0 = warning.
    ///@}

    /// @name Ownership and output flags
    ///@{
    bool newed_SAMG_Settings_object;    ///< True if this object owns settings_ and must delete it.
    bool output_amg_data_to_text_files_; ///< If true, text files are written after each solve.
    ///@}
};

} // end namespace csmp

#endif // SAMG_SOLVER_H

