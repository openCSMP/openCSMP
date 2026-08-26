#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "SAMG_Exception.h"
#include "CompressedRowMatrix.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"

// ---------------------------------------------------------------------------
// Header selection: legacy multi-instance vs. modern context-based API
// ---------------------------------------------------------------------------
#if defined(SAMG_MULTIPLE_INSTANCES) && defined(LEGACY_SAMG)
#  include "samg.h"
#  include "samg1.h"
#  include "samg2.h"
#  include "samg3.h"
#  include "samg4.h"
#  include "samg5.h"
#else
#  include "samg.h"
#endif

using namespace std;

namespace csmp {

// ===========================================================================
// Constructors / Destructor / Assignment
// ===========================================================================

SAMG_Solver::SAMG_Solver()
    : Solver(new SAMG_Settings()),
      settings_(dynamic_cast<SAMG_Settings*>(solver_settings_)),
      nsys_(0),
      ndiu_(1),
      ndip_(1),
      ip_(1),   // SKM FIX: MSVC needs minimum size of 1
      ierr_(0),
      res_in_(-1.),
      res_out_(-1.),
      ncyc_done_(0),
      ncyc_best_(0),
      output_amg_data_to_text_files_(false),
      newed_SAMG_Settings_object(true)
{
    assert(solver_settings_ != nullptr);
    assert(settings_        != nullptr);
}


SAMG_Solver::SAMG_Solver(SAMG_Settings* settings)
    : Solver(settings),
      settings_(settings),
      nsys_(0),
      ndiu_(1),
      ndip_(1),
      ip_(1),   // SKM FIX: MSVC needs minimum size of 1
      ierr_(0),
      res_in_(-1.),
      res_out_(-1.),
      ncyc_done_(0),
      ncyc_best_(0),
      output_amg_data_to_text_files_(false),
      newed_SAMG_Settings_object(false)
{
    assert(solver_settings_ != nullptr);
    assert(settings_        != nullptr);
}


SAMG_Solver::SAMG_Solver(const SAMG_Solver& solver)
    : Solver(new SAMG_Settings(*solver.settings_)),   // deep-copy settings
      settings_(dynamic_cast<SAMG_Settings*>(solver_settings_)),
      u_(solver.u_),
      f_(solver.f_),
      iscale_(solver.iscale_),
      nsys_(solver.nsys_),
      ndiu_(solver.ndiu_),
      ndip_(solver.ndip_),
      iu_(solver.iu_),
      ip_(solver.ip_),
      ierr_(solver.ierr_),
      res_in_(-1.),
      res_out_(-1.),
      ncyc_done_(0),
      ncyc_best_(0),
      output_amg_data_to_text_files_(solver.output_amg_data_to_text_files_),
      newed_SAMG_Settings_object(true)
{
    assert(solver_settings_ != nullptr);
    assert(settings_        != nullptr);
}


SAMG_Solver& SAMG_Solver::operator=(const SAMG_Solver& solver)
{
    if (this != &solver) {
        if (newed_SAMG_Settings_object) {
            delete settings_;
            settings_ = new SAMG_Settings(*solver.settings_);
            solver_settings_ = settings_;
        } else {
            settings_       = solver.settings_;
            solver_settings_ = solver.solver_settings_;
        }
        u_      = solver.u_;
        f_      = solver.f_;
        iscale_ = solver.iscale_;
        nsys_   = solver.nsys_;
        ndiu_   = solver.ndiu_;
        ndip_   = solver.ndip_;
        iu_     = solver.iu_;
        ip_     = solver.ip_;
        ierr_   = solver.ierr_;
        output_amg_data_to_text_files_ = solver.output_amg_data_to_text_files_;
    }
    assert(solver_settings_ != nullptr);
    assert(settings_        != nullptr);
    return *this;
}


SAMG_Solver::~SAMG_Solver()
{
    // Only delete if we own the settings object.
    // solver_settings_ and settings_ point to the same object, so delete once.
    if (newed_SAMG_Settings_object)
        delete settings_;
    // Do NOT delete solver_settings_ separately — it is the same pointer.
    settings_        = nullptr;
    solver_settings_ = nullptr;
}


// ===========================================================================
// Settings accessors
// ===========================================================================

/**
 * Replace the current solver settings.
 *
 * @warning Ownership semantics: if this solver owns its current settings
 *          (i.e. they were default-constructed or copy-constructed), the old
 *          object is deleted.  The caller retains ownership of @p settings and
 *          is responsible for its lifetime.  After this call
 *          newed_SAMG_Settings_object is set to false.
 */
void SAMG_Solver::InputSolverSettings(SolverSettings& settings)
{
    SAMG_Settings* cast = static_cast<SAMG_Settings*>(&settings);
    if (newed_SAMG_Settings_object && settings_ != nullptr)
        delete settings_;           // release previously owned object
    settings_                  = cast;
    solver_settings_           = cast;
    newed_SAMG_Settings_object = false; // caller owns the new object
}


SolverSettings* SAMG_Solver::GetSolverSettings()
{
    return settings_;
}


// ===========================================================================
// Internal helper: build a writable C-string from the solver-instance index
// (required by the modern context-based SAMG API).
// The caller is responsible for calling ReleaseContextString().
// ===========================================================================

#ifndef LEGACY_SAMG
static char* MakeContextString(int instance)
{
    std::string str = std::to_string(instance);
    char* buf = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), buf);
    buf[str.size()] = '\0';
    return buf;
}
#endif


// ===========================================================================
// Internal helper: set SAMG hidden parameters for a given instance.
// Centralises the per-instance preprocessor fan-out so that
// SolveMatrixEquation stays readable.
// ===========================================================================

void SAMG_Solver::SetSAMGHiddenParameters(
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
        )
{
    const int idmp = settings_->Get_idmp();

#ifdef RENOUNCE_COARSENING
    CheckSparsityCriterion(levelx);
#endif

#ifndef SAMG_MULTIPLE_INSTANCES
    // -----------------------------------------------------------------------
    // Single-instance path
    // -----------------------------------------------------------------------
    SAMG_SET_NCG(&ncg);
    SAMG_SET_LEVELX(&levelx);
    if (idmp > 1) {
        SAMG_ISET_IOFORM(&ioform, &ioform_length);
        SAMG_ISET_FILNAM_DUMP(filnam_dump, &filnam_dump_length);
    }
    SAMG_SET_MODE_MESS(&mode_mess);
    SAMG_SET_NRD(&nrd);
    SAMG_SET_NRU(&nru);

#else // SAMG_MULTIPLE_INSTANCES

#  ifdef LEGACY_SAMG
    // -----------------------------------------------------------------------
    // Legacy multi-instance path: explicit SAMGn_* calls
    // -----------------------------------------------------------------------
    const int inst = settings_->GetSolverInstance();

#   define SAMG_SET_HIDDEN_FOR_INSTANCE(PREFIX)                              \
        PREFIX##SET_NCG(&ncg);                                               \
        PREFIX##SET_LEVELX(&levelx);                                         \
        if (idmp > 1) {                                                      \
            PREFIX##ISET_IOFORM(&ioform, &ioform_length);                    \
            PREFIX##ISET_FILNAM_DUMP(filnam_dump, &filnam_dump_length);      \
        }                                                                    \
        PREFIX##SET_MODE_MESS(&mode_mess);                                   \
        PREFIX##SET_NRD(&nrd);                                               \
        PREFIX##SET_NRU(&nru);

    switch (inst) {
        case 0: { SAMG_SET_HIDDEN_FOR_INSTANCE(SAMG_)  break; }
        case 1: { SAMG_SET_HIDDEN_FOR_INSTANCE(SAMG1_) break; }
        case 2: { SAMG_SET_HIDDEN_FOR_INSTANCE(SAMG2_) break; }
        case 3: { SAMG_SET_HIDDEN_FOR_INSTANCE(SAMG3_) break; }
        case 4: { SAMG_SET_HIDDEN_FOR_INSTANCE(SAMG4_) break; }
        case 5: { SAMG_SET_HIDDEN_FOR_INSTANCE(SAMG5_) break; }
        default:
            throw csmp::Exception(ERROR,
                "SAMG_Solver::SetSAMGHiddenParameters",
                "Solver instance index out of range [0,5]");
    }
#   undef SAMG_SET_HIDDEN_FOR_INSTANCE

#  else // !LEGACY_SAMG — modern context-based API
    // -----------------------------------------------------------------------
    // Modern multi-instance path: context string passed to SAMG_*_CTX macros
    // -----------------------------------------------------------------------
    char* ctx = MakeContextString(settings_->GetSolverInstance());

    SAMG_SET_NCG_CTX(SAMG_CCTXT(ctx), &ncg);
    SAMG_SET_LEVELX_CTX(SAMG_CCTXT(ctx), &levelx);
    SAMG_SET_CLSOLVER_FINEST_CTX(SAMG_CCTXT(ctx), &clsolver_finest);

    if (idmp > 1) {
        SAMG_ISET_IOFORM_CTX(SAMG_CCTXT(ctx), &ioform, &ioform_length);
        SAMG_ISET_FILNAM_DUMP_CTX(SAMG_CCTXT(ctx), filnam_dump, &filnam_dump_length);
    }
    SAMG_SET_MODE_MESS_CTX(SAMG_CCTXT(ctx), &mode_mess);
    SAMG_SET_NRD_CTX(SAMG_CCTXT(ctx), &nrd);
    SAMG_SET_NRU_CTX(SAMG_CCTXT(ctx), &nru);

#    if defined(_OPENMP)
    SAMG_SET_ICOLOR_OMP_CTX(SAMG_CCTXT(ctx), &icolor_omp);
    SAMG_SET_IORDERED_OMP_CTX(SAMG_CCTXT(ctx), &iordered_omp);
    SAMG_SET_IRESTRICTION_OPENMP_CTX(SAMG_CCTXT(ctx), &irestriction_openmp);
    if (samg_omp_num_threads_external > 0)
        SAMG_SET_OMP_NUM_THREADS_EXTERNAL_CTX(SAMG_CCTXT(ctx), &samg_omp_num_threads_external);
#    endif

    delete[] ctx;
#  endif // LEGACY_SAMG
#endif   // SAMG_MULTIPLE_INSTANCES
}


// ===========================================================================
// Internal helper: reset SAMG hidden parameters for a given instance.
// ===========================================================================

void SAMG_Solver::ResetSAMGHiddenParameters()
{
#ifndef SAMG_OLD_INTERFACE

#  ifndef SAMG_MULTIPLE_INSTANCES
    SAMG_RESET_HIDDEN();

#  else // SAMG_MULTIPLE_INSTANCES

#    ifdef LEGACY_SAMG
    switch (settings_->GetSolverInstance()) {
        case 0: SAMG_RESET_HIDDEN();  break;
        case 1: SAMG1_RESET_HIDDEN(); break;
        case 2: SAMG2_RESET_HIDDEN(); break;
        case 3: SAMG3_RESET_HIDDEN(); break;
        case 4: SAMG4_RESET_HIDDEN(); break;
        case 5: SAMG5_RESET_HIDDEN(); break;
        default:
            throw csmp::Exception(ERROR,
                "SAMG_Solver::ResetSAMGHiddenParameters",
                "Solver instance index out of range [0,5]");
    }
#    else // modern context API
    char* ctx = MakeContextString(settings_->GetSolverInstance());
    SAMG_RESET_HIDDEN_CTX(SAMG_CCTXT(ctx));
    delete[] ctx;
#    endif

#  endif // SAMG_MULTIPLE_INSTANCES

#else // SAMG_OLD_INTERFACE
    SAMG_RESET_SECONDARY();
#endif
}


// ===========================================================================
// Internal helper: dispatch the actual SAMG solve call.
// ===========================================================================

void SAMG_Solver::DispatchSAMGSolve(
        int32_t& nsolve, int32_t& ifirst, double& eps,
        int32_t& ncyc,   int32_t& iswtch,
        double& a_cmplx, double& g_cmplx,
        double& p_cmplx, double& w_avrge,
        double& chktol,  int32_t& idump, int32_t& iout,
        int32_t& matrix)
{
#ifndef SAMG_MULTIPLE_INSTANCES
    // -----------------------------------------------------------------------
    // Single-instance path
    // -----------------------------------------------------------------------
    cout << "\n\n*** SAMG_Solver::SolveMatrixEquation: Calling SAMG"
         << "( nnu = " << nnu_ << ", nna = " << nna_ << ", ... ) ***\n\n";
    cout.flush();

#  ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
#    ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion(iswtch);
#    endif
#  endif

    SAMG(&nnu_, &nna_, &nsys_,
         &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
         &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
         &res_in_, &res_out_, &ncyc_done_, &ierr_,
         &nsolve, &ifirst, &eps, &ncyc, &iswtch,
         &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
         &chktol, &idump, &iout);

#  ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
#    ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion(iswtch);
#    endif
#  endif

#else // SAMG_MULTIPLE_INSTANCES

    if (Verbose())
        cout << "\n*** SAMG_Solver::SolveMatrixEquation: Calling SAMG instance: "
             << settings_->GetSolverInstance() << " ***\n\n";
    cout.flush();

#  ifdef LEGACY_SAMG
    // -----------------------------------------------------------------------
    // Legacy multi-instance path
    // -----------------------------------------------------------------------
#   define SAMG_SOLVE_INSTANCE(SAMGFN)                                       \
        do {                                                                  \
            CheckCycleCriterionIfNeeded(iswtch);                             \
            SAMGFN(&nnu_, &nna_, &nsys_,                                     \
                   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0],               \
                   &f_[0], &u_[0],                                           \
                   &iu_[0], &ndiu_, &ip_[0], &ndip_,                         \
                   &matrix, &iscale_[0],                                     \
                   &res_in_, &res_out_, &ncyc_done_, &ierr_,                 \
                   &nsolve, &ifirst, &eps, &ncyc, &iswtch,                   \
                   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,                  \
                   &chktol, &idump, &iout);                                  \
            UpdateCycleCriterionIfNeeded(iswtch);                            \
        } while(0)

    switch (settings_->GetSolverInstance()) {
        case 0: SAMG_SOLVE_INSTANCE(SAMG);  break;
        case 1: SAMG_SOLVE_INSTANCE(SAMG1); break;
        case 2: SAMG_SOLVE_INSTANCE(SAMG2); break;
        case 3: SAMG_SOLVE_INSTANCE(SAMG3); break;
        case 4: SAMG_SOLVE_INSTANCE(SAMG4); break;
        case 5: SAMG_SOLVE_INSTANCE(SAMG5); break;
        default:
            throw csmp::Exception(ERROR,
                "SAMG_Solver::DispatchSAMGSolve",
                "Desired solver instance is not available in current SAMG library");
    }
#   undef SAMG_SOLVE_INSTANCE

#  else // !LEGACY_SAMG — modern context-based API
    // -----------------------------------------------------------------------
    // Modern multi-instance path
    // -----------------------------------------------------------------------
    char* ctx = MakeContextString(settings_->GetSolverInstance());

    CheckCycleCriterionIfNeeded(iswtch);

    SAMG_CTX(SAMG_CCTXT(ctx),
             &nnu_, &nna_, &nsys_,
             &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
             &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
             &res_in_, &res_out_, &ncyc_done_, &ierr_,
             &nsolve, &ifirst, &eps, &ncyc, &iswtch,
             &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
             &chktol, &idump, &iout);

    UpdateCycleCriterionIfNeeded(iswtch);

    delete[] ctx;
#  endif // LEGACY_SAMG
#endif   // SAMG_MULTIPLE_INSTANCES
}


// ---------------------------------------------------------------------------
// Thin wrappers that honour the NO_PRIMARY_SOLVER_CONTROL guard
// ---------------------------------------------------------------------------

inline void SAMG_Solver::CheckCycleCriterionIfNeeded(int32_t& iswtch) const
{
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion(iswtch);
#else
    (void)iswtch;
#endif
}

inline void SAMG_Solver::UpdateCycleCriterionIfNeeded(int32_t iswtch)
{
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion(iswtch);
#else
    (void)iswtch;
#endif
}


// ===========================================================================
// SolveMatrixEquation — CompressedRowMatrix overload (retained from old version)
// ===========================================================================

void SAMG_Solver::SolveMatrixEquation(CompressedRowMatrix& A,
                                      vector<double>& b,
                                      vector<double>& x,
                                      size_t no_unknowns)
{
    const int32_t new_nnu  = static_cast<int32_t>(A.ia.size() - 1);
    const int32_t new_nna  = static_cast<int32_t>(A.ja.size());
    const int32_t new_nsys = static_cast<int32_t>(no_unknowns);

    u_.resize(static_cast<size_t>(new_nnu));
    f_.resize(static_cast<size_t>(new_nnu));

    crmat_ = A;
    crmat_.ConvertToSAMGFormat();

    if (Verbose())
        cout << "\nSAMG_Solver::SolveMatrixEquation (CRM): NNU=" << new_nnu
             << " NNA=" << new_nna << "\n";

    PrepareSystemArrays(new_nnu, new_nna, new_nsys);
    InitialiseSolutionVectors(x, b);
    ValidateVectors();
    CallSAMG();
    CopySolutionBack(x);

    if (output_amg_data_to_text_files_)
        Write_SAMG_TextInputFile("csp_solution_data");
}


// ===========================================================================
// SolveMatrixEquation — SparseMatrix overload (primary path)
// ===========================================================================

void SAMG_Solver::SolveMatrixEquation(SparseMatrix& A,
                                      vector<double>& b,
                                      vector<double>& x,
                                      size_t no_unknowns)
{
    const int32_t new_nnu  = static_cast<int32_t>(A.Rows());
    const int32_t new_nna  = static_cast<int32_t>(A.Entries());
    const int32_t new_nsys = static_cast<int32_t>(no_unknowns);

    u_.resize(static_cast<size_t>(new_nnu));
    f_.resize(static_cast<size_t>(new_nnu));

    if (Verbose())
        cout << "\n\n*** SAMG_Solver::SolveMatrixEquation: Initialising matrix ...\n\n";
    cout.flush();

    if (settings_->UsePointBasedApproach())
        crmat_.InitializePointBasedSAMG(A, static_cast<size_t>(new_nsys));
    else
        crmat_.Initialize(A);

    if (Verbose())
        cout << "\nSAMG_Solver::SolveMatrixEquation: NNU=" << new_nnu
             << " NNA=" << new_nna << "\n";

    PrepareSystemArrays(new_nnu, new_nna, new_nsys);
    InitialiseSolutionVectors(x, b);
    ValidateVectors();
    CallSAMG();
    CopySolutionBack(x);

    if (output_amg_data_to_text_files_)
        Write_SAMG_TextInputFile("csp_solution_data");
}


// ===========================================================================
// Private helpers shared by both SolveMatrixEquation overloads
// ===========================================================================

void SAMG_Solver::PrepareSystemArrays(int32_t new_nnu,
                                      int32_t new_nna,
                                      int32_t new_nsys)
{
    if (new_nsys > 1) {
        if (new_nnu != nnu_ || new_nsys != nsys_) {
            if (new_nnu != nnu_) {
                iu_.resize(static_cast<size_t>(new_nnu));
                ndiu_ = static_cast<int32_t>(iu_.size());
            }

            uint32_t k = 0U;
            if (settings_->UsePointBasedApproach()) {
                // Point-based: IU = [1,2,...,nsys, 1,2,...,nsys, ...]
                for (int32_t i = 0; i < new_nnu / new_nsys; ++i)
                    for (int32_t j = 0; j < new_nsys; ++j)
                        iu_[k++] = j + 1;
            } else {
                // Unknown-based: IU = [1,1,...,1, 2,2,...,2, ...]
                for (int32_t j = 0; j < new_nsys; ++j)
                    for (int32_t i = 0; i < new_nnu / new_nsys; ++i)
                        iu_[k++] = j + 1;
            }

            if (settings_->UsePointBasedApproach()) {
                if (new_nnu != nnu_) {
                    ip_.resize(static_cast<size_t>(new_nnu));
                    ndip_ = static_cast<int32_t>(ip_.size());
                }
                k = 0U;
                // IP = [1,1,...,1, 2,2,...,2, ...] (node index per unknown)
                for (int32_t i = 0; i < new_nnu / new_nsys; ++i)
                    for (int32_t j = 0; j < new_nsys; ++j)
                        ip_[k++] = i + 1;
            }
        }
    } else {
        ndiu_ = 1;
        ndip_ = 1;
        iu_.resize(1);
        ip_.resize(1);
    }

    if (new_nsys != nsys_) {
        iscale_.resize(static_cast<size_t>(new_nsys));
        fill(iscale_.begin(), iscale_.end(), 0);
    }

    nsys_ = new_nsys;
    nna_  = new_nna;
    nnu_  = new_nnu;
}


void SAMG_Solver::InitialiseSolutionVectors(const vector<double>& x,
                                            const vector<double>& b)
{
    if (settings_->UsePointBasedApproach()) {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int32_t i = 0; i < nnu_; ++i) {
            const size_t src = static_cast<size_t>(
                i % nsys_ * (nnu_ / nsys_) + i / nsys_);
            u_[static_cast<size_t>(i)] = x[src];
            f_[static_cast<size_t>(i)] = b[src];
        }
    } else {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int32_t i = 0; i < nnu_; ++i) {
            u_[static_cast<size_t>(i)] = x[static_cast<size_t>(i)];
            f_[static_cast<size_t>(i)] = b[static_cast<size_t>(i)];
        }
    }
}


void SAMG_Solver::ValidateVectors() const
{
    if (crmat_.a.empty())  cout << "\n a is empty.\n";
    if (crmat_.ja.empty()) cout << "\n ja is empty.\n";
    if (crmat_.ia.empty()) cout << "\n ia is empty.\n";
    if (u_.empty())        cout << "\n u is empty.\n";
    if (f_.empty())        cout << "\n f is empty.\n";

    if (crmat_.a.empty() || crmat_.ja.empty() || crmat_.ia.empty()
        || u_.empty() || f_.empty())
        throw csmp::Exception(FATAL_ERROR,
            "SAMG_Solver::SolveMatrixEquation",
            "Unable to allocate required memory for transfer arrays");
}


void SAMG_Solver::CallSAMG()
{
    if (Verbose()) {
        cout << "\nSAMG_Solver::CallSAMG: Calling SAMG ...\n";
        cout.flush();
    }

    res_in_  = -1.;
    res_out_ = -1.;

    // Retrieve all settings parameters
    int32_t matrix  = settings_->Get_matrix();
    int32_t nsolve  = settings_->Get_nsolve();
    int32_t ifirst  = settings_->Get_ifirst();
    double  eps     = settings_->Get_eps();
    int32_t ncyc    = settings_->Get_ncyc();
    int32_t iswtch  = settings_->Get_iswtch();
    double  a_cmplx = settings_->Get_a_cmplx();
    double  g_cmplx = settings_->Get_g_cmplx();
    double  p_cmplx = settings_->Get_p_cmplx();
    double  w_avrge = settings_->Get_w_avrge();
    double  chktol  = settings_->Get_chktol();
    int32_t idump   = settings_->Get_idump();
    int32_t iout    = settings_->Get_iout();
    int32_t mode_mess = settings_->Get_mode_mess();
    int32_t nrd     = settings_->Get_nrd();
    int32_t nru     = settings_->Get_nru();
    int32_t ncg     = settings_->Get_ncg();

    int32_t ioform             = settings_->Get_ioform();
    int32_t ioform_length      = settings_->Get_ioform_length();
    int*    filnam_dump        = settings_->Get_filnam_dump();
    int32_t filnam_dump_length = settings_->Get_filnam_dump_length();

#if defined(_OPENMP)
    int32_t icolor_omp                  = settings_->Get_icolor_omp();
    int32_t iordered_omp                = settings_->Get_iordered_omp();
    int32_t irestriction_openmp         = settings_->Get_irestriction_openmp();
    int32_t samg_omp_num_threads_external = settings_->Get_samg_omp_num_threads_external();
#endif

    // ------------------------------------------------------------------
    // Configure SAMG hidden / secondary parameters
    // ------------------------------------------------------------------
    if (settings_->ExplicitSecondary()) {

        int32_t levelx         = settings_->Get_levelx();
        int32_t clsolver_finest = settings_->Get_clsolver_finest();

        // Adaptive direct-solver threshold (new feature from new version)
#ifndef LEGACY_SAMG
        {
            int32_t nptmax = 0;
            char* ctx = MakeContextString(settings_->GetSolverInstance());
            SAMG_GET_NPTMAX_CTX(SAMG_CCTXT(ctx), &nptmax);
            nptmax = std::min(nptmax, settings_->Get_nptmax());
            SAMG_SET_NPTMAX_CTX(SAMG_CCTXT(ctx), &nptmax);
            delete[] ctx;

            if (nnu_ < std::abs(nptmax)) {
                levelx          = 1;
                clsolver_finest = 107;
            }
        }
#else
        {
            int32_t nptmax = 0;
            SAMG_GET_NPTMAX(&nptmax);
            nptmax = std::min(nptmax, settings_->Get_nptmax());
            SAMG_SET_NPTMAX(&nptmax);

            if (nnu_ < std::abs(nptmax)) {
                levelx          = 1;
                clsolver_finest = 107;
            }
        }
#endif

        SetSAMGHiddenParameters(
            ncg, levelx, clsolver_finest,
            ioform, ioform_length, filnam_dump, filnam_dump_length,
            mode_mess, nrd, nru
#if defined(_OPENMP) && !defined(LEGACY_SAMG)
            , icolor_omp, iordered_omp,
              irestriction_openmp, samg_omp_num_threads_external
#endif
        );

    } else {
        ResetSAMGHiddenParameters();
    }

    // ------------------------------------------------------------------
    // Dispatch the solve
    // ------------------------------------------------------------------
    ierr_ = 0;
    DispatchSAMGSolve(nsolve, ifirst, eps, ncyc, iswtch,
                      a_cmplx, g_cmplx, p_cmplx, w_avrge,
                      chktol, idump, iout, matrix);

    // ------------------------------------------------------------------
    // Error / warning reporting
    // ------------------------------------------------------------------
    ErrorHandler& csmp_error = ErrorHandler::Instance();

    if (ierr_ > 0) {
        csmp_error.Note(ERROR,
            "SAMG_Solver::SolveMatrixEquation",
            "SAMG solver returned with an error; error code: ",
            to_string(ierr_).c_str());
    } else if (ierr_ < 0 && ierr_ != -841) { // -841 = bicgstab restart
        if (Verbose())
            csmp_error.Note(WARNING,
                "SAMG_Solver::SolveMatrixEquation",
                "SAMG solver returned with a warning; code: ",
                to_string(ierr_).c_str());
    }

    CheckConvergence(eps);
}


void SAMG_Solver::CopySolutionBack(vector<double>& x) const
{
    if (settings_->UsePointBasedApproach()) {
        for (int32_t i = 0; i < nnu_; ++i)
            x[static_cast<size_t>(i % nsys_ * (nnu_ / nsys_) + i / nsys_)]
                = u_[static_cast<size_t>(i)];
    } else {
        x = u_;
    }
}


// ===========================================================================
// Miscellaneous public methods
// ===========================================================================

void SAMG_Solver::Write_SAMG_TextOutput(bool write)
{
    output_amg_data_to_text_files_ = write;
}


void SAMG_Solver::OutputVectors() const
{
    ofstream ofs("single_output.txt");
    ofs.setf(ios::scientific);
    ofs.precision(15);

    ofs << "\n\nu:\n";
    for (int32_t i = 0; i < nnu_; ++i)
        ofs << (i + 1) << ", " << u_[static_cast<size_t>(i)] << "\n";

    ofs << "\n\nf:\n";
    for (int32_t i = 0; i < nnu_; ++i)
        ofs << (i + 1) << ", " << f_[static_cast<size_t>(i)] << "\n";
}


/**
 * Outputs the input matrix and vectors into the file-based interface format
 * expected by SAMG.  Generates the following files:
 *
 * | Extension | Contents |
 * |-----------|----------|
 * | `.frm`    | Format descriptor (dimensions, matrix type, scaling flags) |
 * | `.amg`    | Compressed-row arrays ia, ja, a |
 * | `.rhs`    | Right-hand side vector f |
 * | `.lhs`    | Solution / initial-guess vector u |
 * | `.a`      | Human-readable matrix with column indices |
 * | `.iu`     | Unknown-type array (coupled systems only) |
 * | `.ip`     | Node-index array (point-based coupled systems only) |
 *
 * See the SAMG users guide for the precise meaning of each file.
 */
bool SAMG_Solver::Write_SAMG_TextInputFile(const char* file) const
{
    string   out_file(file);
    bool     complete_output(true);

    // ------------------------------------------------------------------
    // .frm — format descriptor
    // ------------------------------------------------------------------
    out_file += ".frm";
    ofstream ofs(out_file.c_str());

    ofs << out_file << "\nCSMP output file created for SAMG test.\n";
    ofs << "f         4\n";
    ofs << "# NNA  NNU  MATRIX  NSYS  NPNT\n";
    ofs << nna_ << "  " << nnu_ << "  "
        << settings_->Get_matrix() << "  " << nsys_ << "  ";
    ofs << (ndip_ > 1 ? 1 : 0) << "\n";
    ofs << "# ISCALE(1...NSYS)\n";
    for (int32_t i = 0; i < nsys_; ++i) ofs << 0 << "\n";
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
         << "' written successfully.\n";

    // ------------------------------------------------------------------
    // .amg — ia, ja, a arrays
    // ------------------------------------------------------------------
    assert(!crmat_.ia.empty());
    assert(!crmat_.ja.empty());
    assert(!crmat_.a.empty());

    out_file  = file;
    out_file += ".amg";
    ofs.open(out_file.c_str());

    for (size_t i = 0U; i < static_cast<size_t>(nnu_ + 1); ++i)
        ofs << crmat_.ia[i] << "\n";
    for (size_t i = 0U; i < static_cast<size_t>(nna_); ++i)
        ofs << crmat_.ja[i] << "\n";

    ofs.setf(ios::scientific);
    long prec = ofs.precision(15);
    for (size_t i = 0U; i < static_cast<size_t>(nna_); ++i)
        ofs << crmat_.a[i] << "\n";
    ofs.unsetf(ios::scientific);
    ofs.precision(prec);
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
         << "' written successfully.\n";

    // ------------------------------------------------------------------
    // .rhs — right-hand side vector f
    // ------------------------------------------------------------------
    assert(!f_.empty());

    out_file  = file;
    out_file += ".rhs";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    for (size_t i = 0U; i < static_cast<size_t>(nnu_); ++i)
        ofs << f_[i] << "\n";
    ofs.unsetf(ios::scientific);
    ofs.precision(prec);
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
         << "' written successfully.\n";

    // ------------------------------------------------------------------
    // .lhs — solution / initial-guess vector u
    // ------------------------------------------------------------------
    assert(!u_.empty());

    out_file  = file;
    out_file += ".lhs";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    for (size_t i = 0U; i < u_.size(); ++i)
        ofs << u_[i] << "\n";
    ofs.unsetf(ios::scientific);
    ofs.precision(prec);
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
         << "' written successfully.\n";

    // ------------------------------------------------------------------
    // .a — human-readable matrix with column indices
    // ------------------------------------------------------------------
    assert(!crmat_.a.empty());

    out_file  = file;
    out_file += ".a";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    {
        size_t entry = 0U;
        for (size_t row = 0U; row < static_cast<size_t>(nnu_); ++row) {
            ofs << "\n" << (row + 1) << "\t";
            for (int32_t j = crmat_.ia[row]; j != crmat_.ia[row + 1]; ++j, ++entry)
                ofs << crmat_.a[entry] << "\t[" << crmat_.ja[entry] << "]\t";
        }
    }
    ofs.unsetf(ios::scientific);
    ofs.precision(prec);
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
         << "' written successfully.\n";

    // ------------------------------------------------------------------
    // .iu — unknown-type array (coupled systems only)
    // ------------------------------------------------------------------
    if (nsys_ > 1 && ndiu_ > 1) {
        assert(!iu_.empty());
        out_file  = file;
        out_file += ".iu";
        ofs.open(out_file.c_str());
        for (size_t i = 0U; i < static_cast<size_t>(ndiu_); ++i)
            ofs << iu_[i] << "\n";
        ofs.close();
        cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
             << "' written successfully.\n";
    }

    // ------------------------------------------------------------------
    // .ip — node-index array (point-based coupled systems only)
    // ------------------------------------------------------------------
    if (nsys_ > 1 && ndip_ > 1) {
        assert(!ip_.empty());
        out_file  = file;
        out_file += ".ip";
        ofs.open(out_file.c_str());
        for (size_t i = 0U; i < static_cast<size_t>(ndip_); ++i)
            ofs << ip_[i] << "\n";
        ofs.close();
        cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: file '" << out_file
             << "' written successfully.\n";
    }

    cout.flush();
    return complete_output;
}


double SAMG_Solver::LastSolverResidual() const
{
    return res_out_;
}


void SAMG_Solver::CheckSparsityCriterion(int32_t& levelx) const
{
    if (levelx == 1 && nna_ / nnu_ > 1.5) {
        cout << "\n\n*** SAMG_Solver::CheckSparsityCriterion: "
                "Matrix sparsity criterion not fulfilled, levelx reset to default. ***\n"
             << "\t(Single-level setup should only be used for essentially diagonal "
                "matrices, e.g. transport equations.)\n"
             << "\tlevelx set to: ";
        levelx = 25;
        cout << levelx << "\n\n";
    }
}


void SAMG_Solver::CheckCycleCriterion(int32_t& iswtch) const
{
    if (ncyc_done_ > 1.5 * ncyc_best_) {
        iswtch += 1000000; // forces a new coarsening setup on the next call
        cout << "\n\n*** SAMG_Solver::CheckCycleCriterion: "
                "Iteration cycle count criterion failed — new coarsening setup forced! ***\n\n";
        cout.flush();
    }
}


void SAMG_Solver::UpdateCycleCriterion(int32_t iswtch)
{
    if (iswtch < 4000000)
        ncyc_best_ = std::min(ncyc_best_, ncyc_done_);
    else
        ncyc_best_ = ncyc_done_;
}



/**
 * Checks absolute or relative convergence against the user-supplied tolerance.
 *
 * - `eps == 0`: no check performed, returns `true` immediately.
 * - `eps < 0`:  absolute criterion — passes if @f$ r_\text{out} \leq |eps| @f$.
 * - `eps > 0`:  relative criterion — passes if
 *               @f$ r_\text{out} / r_\text{in} \leq eps @f$.
 *
 * A `WARNING` is issued via `ErrorHandler` if the criterion is not met, but
 * no exception is thrown — the caller decides how to handle non-convergence.
 */
bool SAMG_Solver::CheckConvergence(double eps) const
{
    if (fabs(eps) < 1.0e-30) return true;

    ErrorHandler& csmp_error = ErrorHandler::Instance();

    cout << "\n*** SAMG_Solver::CheckConvergence: eps = " << eps
         << ", res_out = " << res_out_
         << ", res_in = "  << res_in_ << " ***\n\n";
    cout.flush();

    // Absolute convergence
    if (eps < 0.) {
        if (res_out_ > fabs(eps)) {
            csmp_error.Note(WARNING,
                "SAMG_Solver::CheckConvergence",
                "Absolute error L2-norm solution criterion was not fulfilled");
            return false;
        }
    }

    // Relative convergence
    if (eps > 0.) {
        const double relative_residual =
            (res_in_ == 0.) ? 0. : res_out_ / res_in_;

        cout << "\n\trelative residual = " << relative_residual << "\n\n";

        if (relative_residual > eps) {
            csmp_error.Note(WARNING,
                "SAMG_Solver::CheckConvergence",
                "Relative solution criterion was not fulfilled");
            return false;
        }
        cout.flush();
    }

    return true;
}


} // end namespace csmp


// ---------------------------------------------------------------------------
// SAMG user coordinate callback — required by some SAMG builds.
// Setting *ndim = 0 tells SAMG that no geometric coordinates are available.
// ---------------------------------------------------------------------------
#ifndef CSP_PC_WINDOWS_NT_CODE_WARRIOR

void samg_user_coo (int*, int* ndim, double*, double*, double*);
void samg_user_coo_(int*, int* ndim, double*, double*, double*);

void samg_user_coo(int*, int* ndim, double*, double*, double*)
{
    cout << "\nSAMG_Solver::samg_user_coo: Setting *ndim to 0\n";
    *ndim = 0;
}

void samg_user_coo_(int*, int* ndim, double*, double*, double*)
{
    cout << "\nSAMG_Solver::samg_user_coo_: Setting *ndim to 0\n";
    *ndim = 0;
}

#endif

