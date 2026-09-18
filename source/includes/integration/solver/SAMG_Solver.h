// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SAMG_SOLVER_H
#define SAMG_SOLVER_H

#include "Solver.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

namespace csmp {

class SAMG_Settings;

/**

@brief Subclass of Solver for linear algebra computations conducted with SAMG,
the algebraic multigrid solver for systems of equations from SCAPOS,
Fraunhofer Gesellschaft, Germany.

Solves the linear algebraic system A x = b using the algebraic multigrid
solver library from SCAI, Fraunhofer Gesellschaft (St. Augustin, Bonn),
which can be obtained from SCAPOS.

Up to 6 instances of the solver library can be used simultaneously when
compiled with SAMG_MULTIPLE_INSTANCES. Use separate instances if you have
multiple PDEs solved sequentially and want SAMG to remember its coarsening
setup between time steps.

@section segFault Segmentation Fault on Some Linux Systems

Invalid memory access has been observed when using dynamically allocated
SAMG_Solver objects. To avoid this, use static (stack) instantiation:

@code
SAMG_Settings settings;
settings.Set_napproach(2);
SAMG_Solver samgSolver(&settings);
PDE_Integrator<DIM,Region> pressure_diffusion(&samgSolver);
@endcode

@section licenceRetry Licence Server Retry

If the SAMG licence server is temporarily unavailable (e.g., during a
server restart), the solver will retry automatically. The number of
attempts and the delay between them are configurable via SAMG_Settings:

@code
settings.Set_license_retry_attempts(120);       // default: 120 attempts
settings.Set_license_retry_delay_seconds(60);   // default: 60 seconds
settings.Set_license_retry_log_path("samg_license_not_found.log");
@endcode

@todo (1) Check: is resetting res_in_ allowed when using relative convergence?
@todo (3) Convergence check: verify whether eps or fabs(eps) should be used
@todo (1) Investigate runtime debug error when solving systems with vector variables
@todo (3) Try mpl::is_volatile type trait to check whether pointers supplied
          in the constructor should be memory-managed by the solver

@author S.K. Matthai
@author various contributions
@date 1997

*/
class SAMG_Solver final : public Solver {
  public:

    // ========================================================================
    // CONSTRUCTORS AND DESTRUCTORS
    // ========================================================================

    SAMG_Solver();

    /// Creates a solver using the supplied settings object.
    /// @attention The SAMG_Solver does NOT take ownership of the settings pointer.
    explicit SAMG_Solver( SAMG_Settings* settings );

    SAMG_Solver( const SAMG_Solver& );
    SAMG_Solver& operator=( const SAMG_Solver& );
    virtual ~SAMG_Solver();

    // ========================================================================
    // SETTINGS INTERFACE
    // ========================================================================

    virtual void            InputSolverSettings( SolverSettings& settings ) override final;
    virtual SolverSettings* GetSolverSettings() override final;
    virtual std::string     Name() const override final { return "SAMG_Solver"; }

    // ========================================================================
    // DIAGNOSTIC OUTPUT
    // ========================================================================

    /// writes SAMG input data to text files for diagnostic purposes
    bool Write_SAMG_TextInputFile( const char* filename ) const;

    /// enables or disables writing of SAMG input/output data to text files
    void Write_SAMG_TextOutput( bool write );

    /// outputs u_ and f_ vectors to a text file
    void OutputVectors() const;

    // ========================================================================
    // CONVERGENCE AND CYCLE CRITERIA
    // ========================================================================

    /// returns the residual of the final approximation from the last solve
    double LastSolverResidual() const;

    /// resets levelx to 25 if the matrix is not nearly diagonal and levelx == 1
    void CheckSparsityCriterion( int32_t& levelx ) const;

    /// forces a new coarsening setup if the cycle count exceeds 1.5 * ncyc_best_
    void CheckCycleCriterion( int32_t& iswtch ) const;

    /// updates the best recorded cycle count after a successful solve
    void UpdateCycleCriterion( int32_t iswtch );

    /// checks whether the solver converged to within the requested tolerance
    bool CheckConvergence( double eps ) const;

  protected:

    // ========================================================================
    // SOLVE INTERFACE
    // ========================================================================

    /// solves A x = b where A is a CompressedRowMatrix (preferred path)
    virtual void SolveMatrixEquation( CompressedRowMatrix& A,
                              std::vector<double>& b,
                              std::vector<double>& x,
                                    size_t no_unknowns ) override;//Benoit add override

    /// solves A x = b where A is a CSMP SparseMatrix (converted internally)
    virtual void SolveMatrixEquation( SparseMatrix& A,
                              std::vector<double>& b,
                              std::vector<double>& x,
                                      size_t no_unknowns ) override;//Benoit add override

  private:

    // ========================================================================
    // INTERNAL HELPERS
    // ========================================================================

    /// initialises IU and IP index arrays for coupled systems
    void InitialiseIndexArrays( size_t new_nnu, size_t new_nsys );

    /// initialises f_ (RHS) and u_ (initial guess) from b and x
    void InitialiseVectors( const std::vector<double>& b,
                            const std::vector<double>& x );

    /// retrieves all SAMG parameters from the settings object
    void RetrieveSAMGParameters( int32_t& matrix,  int32_t& nsolve,
                                 int32_t& ifirst,  double&  eps,
                                 int32_t& ncyc,    int32_t& iswtch,
                                 double&  a_cmplx, double&  g_cmplx,
                                 double&  p_cmplx, double&  w_avrge,
                                 double&  chktol,  int32_t& idump,
                                 int32_t& iout,    int32_t& ncg,
                                 int32_t& levelx ) const;

    /// sets up SAMG hidden (secondary) parameters for the given solver instance
    void SetupHiddenParameters( int32_t idump, int32_t ncg,
                                int32_t levelx, int32_t instance );

    /// calls the appropriate SAMG solver instance
    void CallSAMG( int32_t instance,
                   int32_t nsolve, int32_t ifirst, double eps,
                   int32_t ncyc,   int32_t iswtch,
                   double a_cmplx, double g_cmplx,
                   double p_cmplx, double w_avrge,
                   double chktol,  int32_t idump,
                   int32_t iout,   int32_t matrix );

    /// calls SAMG with automatic retry on licence-not-found errors
    void CallSAMGWithLicenceRetry( int32_t instance,
                                   int32_t nsolve, int32_t ifirst, double eps,
                                   int32_t ncyc,   int32_t iswtch,
                                   double a_cmplx, double g_cmplx,
                                   double p_cmplx, double w_avrge,
                                   double chktol,  int32_t idump,
                                   int32_t iout,   int32_t matrix );

    /// appends a timestamped entry to the licence retry log file; thread-safe
    void LogLicenceRetry( int attempt, int maxAttempts );

    /// maps the SAMG solution vector u_ back into the user-supplied vector x
    void MapSolutionToOutput( std::vector<double>& x ) const;

    // ========================================================================
    // MEMBER VARIABLES
    // ========================================================================

    SAMG_Settings*  settings_;

    // --- System dimensions ---
    int32_t  nsys_;  ///< number of coupled solution variables (1 for scalar)
    int32_t  nnu_;   ///< number of unknowns (matrix size nnu x nnu)
    int32_t  nna_;   ///< number of non-zero matrix entries stored in a

    // --- Solution and RHS vectors ---
    std::vector<double>   u_;  ///< initial guess on entry; solution on exit
    std::vector<double>   f_;  ///< right-hand side vector

    // --- Coupled system index arrays ---
    std::vector<int32_t>  iscale_;  ///< scaling flags per unknown (0=no scaling, 1=scale)
    std::vector<int32_t>  iu_;      ///< unknown type array: iu[i] in [1, nsys]
    int32_t               ndiu_;    ///< size of iu_ (= nnu_ for coupled systems, 1 for scalar)
    std::vector<int32_t>  ip_;      ///< point number array for point-based approach
    int32_t               ndip_;    ///< size of ip_ (= nnu_ for point-based, 1 otherwise)

    // --- Internal matrix storage ---
    CompressedRowMatrix   crmat_;   ///< SAMG-format compressed row matrix

    // --- SAMG output parameters ---
    double   res_in_;      ///< residual of the initial guess
    double   res_out_;     ///< residual of the final approximation
    int32_t  ncyc_done_;   ///< total number of iteration cycles performed
    int32_t  ncyc_best_;   ///< lowest number of iterations achieved (for cycle criterion)
    int32_t  ierr_;        ///< error indicator: 0=ok, >0=fatal error, <0=warning

    // --- Ownership and output flags ---
    const bool  newed_SAMG_Settings_object;   ///< true if settings_ was allocated by this object
    bool        output_amg_data_to_text_files_;
};

} // end namespace csmp

// ============================================================================
// SAMG USER COORDINATE CALLBACK
// Required by the SAMG library; sets ndim=0 to indicate no coordinate data.
// ============================================================================

#ifndef CSP_PC_WINDOWS_NT_CODE_WARRIOR
void samg_user_coo ( int* i, int* ndim, double* x, double* y, double* z );
void samg_user_coo_( int* i, int* ndim, double* x, double* y, double* z );
#endif

#endif /* SAMG_SOLVER_H */
