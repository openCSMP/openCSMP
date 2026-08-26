#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "SAMG_Exception.h"
#include "CompressedRowMatrix.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include <fstream>
#include <thread>
#include <chrono>
#include <mutex>

// ============================================================================
// SAMG LIBRARY HEADERS
// ============================================================================

// Enforce dependency: SAMG_LEGACY requires SAMG_MULTIPLE_INSTANCES
#ifdef SAMG_LEGACY
#  ifndef SAMG_MULTIPLE_INSTANCES
#    define SAMG_MULTIPLE_INSTANCES
#  endif
#endif

// Include the appropriate SAMG headers based on build configuration
#ifndef SAMG_MULTIPLE_INSTANCES
    // Single instance: only the base SAMG header needed
#   include "samg.h"
#else
    // Multiple instances (SAMG_LEGACY): include all numbered instance headers
#   include "samg.h"
#   include "samg1.h"
#   include "samg2.h"
#   include "samg3.h"
#   include "samg4.h"
#   include "samg5.h"
#endif

// Compile-time verification that the headers expose what we need
#ifdef SAMG_MULTIPLE_INSTANCES
#  ifndef SAMG1_SET_NCG
#    error "samg1.h was included but SAMG1_SET_NCG is not defined. \
Check that the correct legacy SAMG library headers are being used."
#  endif
#endif

/// SAMG error code indicating that the licence server could not be reached
static constexpr int32_t SAMG_LICENSE_NOT_FOUND = 5;

using namespace std;

namespace csmp {

// ============================================================================
// CONSTRUCTORS AND DESTRUCTORS
// ============================================================================

SAMG_Solver::SAMG_Solver()
    : Solver(new SAMG_Settings()),
      settings_(dynamic_cast<SAMG_Settings*>(solver_settings_)),
      nsys_(0),
      ndiu_(1),
      ndip_(1),
      ip_(1),  ///< SKM FIX: MS Studio needs minimum size of 1
      ierr_(0),
      res_in_(-1.),
      res_out_(-1.),
      ncyc_done_(0),
      ncyc_best_(0),
      output_amg_data_to_text_files_(false),
      newed_SAMG_Settings_object(true)
{
    assert( solver_settings_ != NULL );
    assert( settings_ != NULL );
}


SAMG_Solver::SAMG_Solver( SAMG_Settings* settings )
    : Solver(settings),
      settings_(settings),
      nsys_(0),
      ndiu_(1),
      ndip_(1),
      ip_(1),  ///< SKM FIX: MS Studio needs minimum size of 1
      ierr_(0),
      res_in_(-1.),
      res_out_(-1.),
      ncyc_done_(0),
      ncyc_best_(0),
      output_amg_data_to_text_files_(false),
      newed_SAMG_Settings_object(false)
{
    assert( solver_settings_ != NULL );
    assert( settings_ != NULL );
}


SAMG_Solver::SAMG_Solver( const SAMG_Solver& solver )
    : Solver(new SAMG_Settings()),
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
    assert( solver_settings_ != NULL );
    assert( settings_ != NULL );
}


SAMG_Solver& SAMG_Solver::operator=( const SAMG_Solver& solver )
{
    if ( this != &solver ) {
        if ( newed_SAMG_Settings_object ) {
            delete settings_;
            settings_ = new SAMG_Settings(*solver.settings_);
            this->solver_settings_ = settings_;
        }
        else {
            settings_        = solver.settings_;
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
    assert( solver_settings_ != NULL );
    assert( settings_ != NULL );
    return *this;
}


SAMG_Solver::~SAMG_Solver()
{
    if ( newed_SAMG_Settings_object ) delete settings_;
}

// ============================================================================
// SETTINGS INTERFACE
// ============================================================================

void SAMG_Solver::InputSolverSettings( SolverSettings& settings )
{
    settings_ = static_cast<SAMG_Settings*>(&settings);
}

SolverSettings* SAMG_Solver::GetSolverSettings()
{
    return settings_;
}

// ============================================================================
// INTERNAL HELPERS
// ============================================================================

/**
    Initialises the IU and IP index arrays for coupled systems of equations.

    For a scalar system (nsys == 1), both arrays are set to a single dummy entry.
    For coupled systems, the arrays are built according to the point-based or
    unknown-based approach as configured in SAMG_Settings.

    @param new_nnu  number of unknowns in the system
    @param new_nsys number of coupled solution variables
*/
void SAMG_Solver::InitialiseIndexArrays( size_t new_nnu, size_t new_nsys )
{
    if ( new_nsys > 1U ) {
        // Only rebuild if sizes have changed
        if ( new_nnu != static_cast<size_t>(nnu_) ||
             new_nsys != static_cast<size_t>(nsys_) ) {

            iu_.resize( new_nnu );
            ndiu_ = static_cast<int32_t>( iu_.size() );

            uint32_t k{ 0U };
            if ( settings_->UsePointBasedApproach() ) {
                // POINT BASED (SAMG manual p.16):
                // solution vector [x1,y1,x2,y2,...,xn,yn]
                // IU pattern: [1,2,1,2,...,1,2]
                for ( int32_t i = 0; i < static_cast<int32_t>(iu_.size() / new_nsys); ++i )
                    for ( int32_t j = 0; j < static_cast<int32_t>(new_nsys); ++j )
                        iu_[k++] = j + 1;

                // IP pattern: [1,1,2,2,...,n,n]
                ip_.resize( new_nnu );
                ndip_ = static_cast<int32_t>( ip_.size() );
                k = 0U;
                for ( int32_t i = 0; i < static_cast<int32_t>(ip_.size() / new_nsys); ++i )
                    for ( int32_t j = 0; j < static_cast<int32_t>(new_nsys); ++j )
                        ip_[k++] = i + 1;
            }
            else {
                // UNKNOWN BASED (SAMG manual p.16):
                // solution vector [x1,x2,...,xn,y1,y2,...,yn]
                // IU pattern: [1,1,...,1,2,2,...,2]
                for ( int32_t j = 0; j < static_cast<int32_t>(new_nsys); ++j )
                    for ( int32_t i = 0; i < static_cast<int32_t>(iu_.size() / new_nsys); ++i )
                        iu_[k++] = j + 1;
            }
        }
    }
    else {
        // Scalar system: dummy arrays of size 1 with value 1
        // SKM FIX: MS Studio needs minimum size of 1
        ndiu_ = 1;
        ndip_ = 1;
        iu_.assign( 1, 1 );
        ip_.assign( 1, 1 );
    }
}


/**
    Initialises the RHS vector f_ and the initial guess vector u_ from
    the supplied vectors b and x, applying the point-based index
    transformation if required.

    @param b    right-hand side vector
    @param x    initial guess / solution vector
*/
void SAMG_Solver::InitialiseVectors( const vector<double>& b,
                                     const vector<double>& x )
{
    if ( settings_->UsePointBasedApproach() ) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for ( int32_t i = 0; i < nnu_; ++i ) {
            const int32_t mapped = i % nsys_ * (nnu_ / nsys_) + i / nsys_;
            u_[i] = x[mapped];
            f_[i] = b[mapped];
        }
    }
    else {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for ( int32_t i = 0; i < nnu_; ++i ) {
            u_[i] = x[i];
            f_[i] = b[i];
        }
    }
}


/**
    Retrieves all SAMG parameters from the settings object into local variables.
    These are passed by pointer to the Fortran SAMG interface.
*/
void SAMG_Solver::RetrieveSAMGParameters( int32_t& matrix,  int32_t& nsolve,
                                          int32_t& ifirst,  double&  eps,
                                          int32_t& ncyc,    int32_t& iswtch,
                                          double&  a_cmplx, double&  g_cmplx,
                                          double&  p_cmplx, double&  w_avrge,
                                          double&  chktol,  int32_t& idump,
                                          int32_t& iout,    int32_t& ncg,
                                          int32_t& levelx ) const
{
    matrix  = settings_->Get_matrix();
    nsolve  = settings_->Get_nsolve();
    ifirst  = settings_->Get_ifirst();
    eps     = settings_->Get_eps();
    ncyc    = settings_->Get_ncyc();
    iswtch  = settings_->Get_iswtch();
    a_cmplx = settings_->Get_a_cmplx();
    g_cmplx = settings_->Get_g_cmplx();
    p_cmplx = settings_->Get_p_cmplx();
    w_avrge = settings_->Get_w_avrge();
    chktol  = settings_->Get_chktol();
    idump   = settings_->Get_idump();
    iout    = settings_->Get_iout();
    ncg     = settings_->Get_ncg();
    levelx  = settings_->Get_levelx();
}


/**
    Sets up SAMG hidden parameters (secondary settings) for the given
    solver instance. If ExplicitSecondary() is false, the hidden parameters
    are reset to their defaults.

    @param idump    dump level (controls file output)
    @param ncg      number of coarse-grid levels
    @param levelx   maximum number of AMG levels
    @param instance solver instance index (0..5); only used with SAMG_MULTIPLE_INSTANCES
*/
void SAMG_Solver::SetupHiddenParameters( int32_t idump, int32_t ncg,
                                         int32_t levelx, int32_t instance )
{
    if ( settings_->ExplicitSecondary() ) {

#ifdef RENOUNCE_COARSENING
        CheckSparsityCriterion( levelx );
#endif

        // Retrieve secondary parameters
        int32_t  mode_mess  = settings_->Get_mode_mess();
        int32_t  nrd        = settings_->Get_nrd();
        int32_t  nru        = settings_->Get_nru();

        // Optional file output parameters
        int32_t  ioform     = settings_->Get_ioform();
        int32_t  ioform_len = settings_->Get_ioform_length();
        int*     filnam     = settings_->Get_filnam_dump();
        int32_t  filnam_len = settings_->Get_filnam_dump_length();

#ifndef SAMG_MULTIPLE_INSTANCES
        // ----------------------------------------------------------------
        // Single instance: only base SAMG_* macros available
        // ----------------------------------------------------------------
        SAMG_SET_NCG( &ncg );
        SAMG_SET_LEVELX( &levelx );
        if ( idump > 1 ) {
            SAMG_ISET_IOFORM( &ioform, &ioform_len );
            SAMG_ISET_FILNAM_DUMP( filnam, &filnam_len );
        }
        SAMG_SET_MODE_MESS( &mode_mess );
        SAMG_SET_NRD( &nrd );
        SAMG_SET_NRU( &nru );

#else
        // ----------------------------------------------------------------
        // Multiple instances: dispatch to the correct SAMG instance
        // SAMG1_* ... SAMG5_* macros only available with SAMG_MULTIPLE_INSTANCES
        // ----------------------------------------------------------------
        switch ( instance ) {
            case 0:
                SAMG_SET_NCG( &ncg );
                SAMG_SET_LEVELX( &levelx );
                if ( idump > 1 ) {
                    SAMG_ISET_IOFORM( &ioform, &ioform_len );
                    SAMG_ISET_FILNAM_DUMP( filnam, &filnam_len );
                }
                SAMG_SET_MODE_MESS( &mode_mess );
                SAMG_SET_NRD( &nrd );
                SAMG_SET_NRU( &nru );
                break;
#ifdef SAMG1_SET_NCG
            case 1:
                SAMG1_SET_NCG( &ncg );
                SAMG1_SET_LEVELX( &levelx );
                if ( idump > 1 ) {
                    SAMG1_ISET_IOFORM( &ioform, &ioform_len );
                    SAMG1_ISET_FILNAM_DUMP( filnam, &filnam_len );
                }
                SAMG1_SET_MODE_MESS( &mode_mess );
                SAMG1_SET_NRD( &nrd );
                SAMG1_SET_NRU( &nru );
                break;
#endif  // SAMG1_SET_NCG
#ifdef SAMG2_SET_NCG
            case 2:
                SAMG2_SET_NCG( &ncg );
                SAMG2_SET_LEVELX( &levelx );
                if ( idump > 1 ) {
                    SAMG2_ISET_IOFORM( &ioform, &ioform_len );
                    SAMG2_ISET_FILNAM_DUMP( filnam, &filnam_len );
                }
                SAMG2_SET_MODE_MESS( &mode_mess );
                SAMG2_SET_NRD( &nrd );
                SAMG2_SET_NRU( &nru );
                break;
#endif  // SAMG2_SET_NCG
#ifdef SAMG3_SET_NCG
            case 3:
                SAMG3_SET_NCG( &ncg );
                SAMG3_SET_LEVELX( &levelx );
                if ( idump > 1 ) {
                    SAMG3_ISET_IOFORM( &ioform, &ioform_len );
                    SAMG3_ISET_FILNAM_DUMP( filnam, &filnam_len );
                }
                SAMG3_SET_MODE_MESS( &mode_mess );
                SAMG3_SET_NRD( &nrd );
                SAMG3_SET_NRU( &nru );
                break;
#endif  // SAMG3_SET_NCG
#ifdef SAMG4_SET_NCG
            case 4:
                SAMG4_SET_NCG( &ncg );
                SAMG4_SET_LEVELX( &levelx );
                if ( idump > 1 ) {
                    SAMG4_ISET_IOFORM( &ioform, &ioform_len );
                    SAMG4_ISET_FILNAM_DUMP( filnam, &filnam_len );
                }
                SAMG4_SET_MODE_MESS( &mode_mess );
                SAMG4_SET_NRD( &nrd );
                SAMG4_SET_NRU( &nru );
                break;
#endif  // SAMG4_SET_NCG
#ifdef SAMG5_SET_NCG
            case 5:
                SAMG5_SET_NCG( &ncg );
                SAMG5_SET_LEVELX( &levelx );
                if ( idump > 1 ) {
                    SAMG5_ISET_IOFORM( &ioform, &ioform_len );
                    SAMG5_ISET_FILNAM_DUMP( filnam, &filnam_len );
                }
                SAMG5_SET_MODE_MESS( &mode_mess );
                SAMG5_SET_NRD( &nrd );
                SAMG5_SET_NRU( &nru );
                break;
#endif  // SAMG5_SET_NCG
            default:
                throw csmp::Exception( ERROR, "SAMG_Solver::SetupHiddenParameters",
                                       "Solver instance index out of range [0,5]" );
        }
#endif  // SAMG_MULTIPLE_INSTANCES

    }
    else {

        // Reset all hidden parameters to SAMG defaults
#ifndef SAMG_MULTIPLE_INSTANCES

#ifndef SAMG_OLD_INTERFACE
        SAMG_RESET_HIDDEN();
#else
        SAMG_RESET_SECONDARY();
#endif  // SAMG_OLD_INTERFACE

#else   // SAMG_MULTIPLE_INSTANCES

        switch ( instance ) {
            case 0:
                SAMG_RESET_HIDDEN();
                break;
#ifdef SAMG1_RESET_HIDDEN
            case 1: SAMG1_RESET_HIDDEN(); break;
#endif
#ifdef SAMG2_RESET_HIDDEN
            case 2: SAMG2_RESET_HIDDEN(); break;
#endif
#ifdef SAMG3_RESET_HIDDEN
            case 3: SAMG3_RESET_HIDDEN(); break;
#endif
#ifdef SAMG4_RESET_HIDDEN
            case 4: SAMG4_RESET_HIDDEN(); break;
#endif
#ifdef SAMG5_RESET_HIDDEN
            case 5: SAMG5_RESET_HIDDEN(); break;
#endif
            default:
                throw csmp::Exception( ERROR, "SAMG_Solver::SetupHiddenParameters",
                                       "Solver instance index out of range [0,5]" );
        }
#endif  // SAMG_MULTIPLE_INSTANCES

    }

} // end SetupHiddenParameters


/**
    Calls the appropriate SAMG solver instance with the supplied parameters.

    @param instance  solver instance index (0..5); only used with SAMG_MULTIPLE_INSTANCES
*/
void SAMG_Solver::CallSAMG( int32_t instance,
                             int32_t nsolve, int32_t ifirst, double eps,
                             int32_t ncyc,   int32_t iswtch,
                             double a_cmplx, double g_cmplx,
                             double p_cmplx, double w_avrge,
                             double chktol,  int32_t idump,
                             int32_t iout,   int32_t matrix )
{
    // Base call lambda: always available
    auto call_base = [&]() {
        SAMG( &nnu_, &nna_, &nsys_,
              &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
              &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
              &res_in_, &res_out_, &ncyc_done_, &ierr_,
              &nsolve, &ifirst, &eps, &ncyc, &iswtch,
              &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
              &chktol, &idump, &iout );
    };

#ifndef SAMG_MULTIPLE_INSTANCES
    // ----------------------------------------------------------------
    // Single instance: only base SAMG() call available
    // ----------------------------------------------------------------
    call_base();

#else
    // ----------------------------------------------------------------
    // Multiple instances: dispatch to the correct SAMG instance
    // SAMG1() ... SAMG5() only available with SAMG_MULTIPLE_INSTANCES
    // ----------------------------------------------------------------
    switch ( instance ) {
        case 0:
            call_base();
            break;
#ifdef SAMG1
        case 1:
            SAMG1( &nnu_, &nna_, &nsys_,
                   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
                   &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
                   &res_in_, &res_out_, &ncyc_done_, &ierr_,
                   &nsolve, &ifirst, &eps, &ncyc, &iswtch,
                   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
                   &chktol, &idump, &iout );
            break;
#endif  // SAMG1
#ifdef SAMG2
        case 2:
            SAMG2( &nnu_, &nna_, &nsys_,
                   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
                   &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
                   &res_in_, &res_out_, &ncyc_done_, &ierr_,
                   &nsolve, &ifirst, &eps, &ncyc, &iswtch,
                   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
                   &chktol, &idump, &iout );
            break;
#endif  // SAMG2
#ifdef SAMG3
        case 3:
            SAMG3( &nnu_, &nna_, &nsys_,
                   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
                   &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
                   &res_in_, &res_out_, &ncyc_done_, &ierr_,
                   &nsolve, &ifirst, &eps, &ncyc, &iswtch,
                   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
                   &chktol, &idump, &iout );
            break;
#endif  // SAMG3
#ifdef SAMG4
        case 4:
            SAMG4( &nnu_, &nna_, &nsys_,
                   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
                   &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
                   &res_in_, &res_out_, &ncyc_done_, &ierr_,
                   &nsolve, &ifirst, &eps, &ncyc, &iswtch,
                   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
                   &chktol, &idump, &iout );
            break;
#endif  // SAMG4
#ifdef SAMG5
        case 5:
            SAMG5( &nnu_, &nna_, &nsys_,
                   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
                   &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
                   &res_in_, &res_out_, &ncyc_done_, &ierr_,
                   &nsolve, &ifirst, &eps, &ncyc, &iswtch,
                   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
                   &chktol, &idump, &iout );
            break;
#endif  // SAMG5
        default:
            throw csmp::Exception( ERROR, "SAMG_Solver::CallSAMG",
                                   "Desired solver instance is not available "
                                   "in current SAMG library" );
    }
#endif  // SAMG_MULTIPLE_INSTANCES

} // end CallSAMG


/**
    Calls SAMG with automatic retry on licence-not-found errors.

    The number of retry attempts and the delay between them are
    configurable via SAMG_Settings. On each failed attempt, a
    timestamped entry is appended to the licence retry log file.
*/
void SAMG_Solver::CallSAMGWithLicenceRetry( int32_t instance,
                                            int32_t nsolve, int32_t ifirst, double eps,
                                            int32_t ncyc,   int32_t iswtch,
                                            double a_cmplx, double g_cmplx,
                                            double p_cmplx, double w_avrge,
                                            double chktol,  int32_t idump,
                                            int32_t iout,   int32_t matrix )
{
    const int maxAttempts  = settings_->Get_license_retry_attempts();
    const int delaySeconds = settings_->Get_license_retry_delay_seconds();

    for ( int attempt = 1; attempt <= maxAttempts; ++attempt )
    {
        try {
            CallSAMG( instance, nsolve, ifirst, eps, ncyc, iswtch,
                      a_cmplx, g_cmplx, p_cmplx, w_avrge,
                      chktol, idump, iout, matrix );

            if ( ierr_ == SAMG_LICENSE_NOT_FOUND )
            {
                cerr << "[Attempt " << attempt << "/" << maxAttempts
                     << "] SAMG licence not available (ierr=" << SAMG_LICENSE_NOT_FOUND
                     << "). Retrying in " << delaySeconds << "s...\n";

                LogLicenceRetry( attempt, maxAttempts );

                // Reset output variables before next attempt
                ierr_    = 0;
                res_in_  = -1.;
                res_out_ = -1.;

                if ( attempt < maxAttempts )
                    this_thread::sleep_for( chrono::seconds(delaySeconds) );
            }
            else
            {
                // Success or a different error — don't retry
                break;
            }
        }
        catch ( const exception& e ) {
            cerr << "[Attempt " << attempt << "] Exception during SAMG call: "
                 << e.what() << "\n";
            break;  // Don't retry on unexpected exceptions
        }
    }

} // end CallSAMGWithLicenceRetry


/**
    Appends a timestamped entry to the SAMG licence retry log file.
    Thread-safe via an internal mutex.

    @param attempt      current attempt number
    @param maxAttempts  maximum number of attempts configured
*/
void SAMG_Solver::LogLicenceRetry( int attempt, int maxAttempts )
{
    // Mutex ensures thread-safe file access when multiple solver
    // instances run concurrently (SAMG_MULTIPLE_INSTANCES)
    static mutex log_mutex;
    lock_guard<mutex> lock( log_mutex );

    const string logPath = settings_->Get_license_retry_log_path();

    ofstream log( logPath, ios::app );
    if ( !log.is_open() ) {
        cerr << "[SAMG_Solver] Warning: could not open licence retry log: "
             << logPath << "\n";
        return;
    }

    // Thread-safe timestamp
    const auto   now = chrono::system_clock::now();
    const time_t t   = chrono::system_clock::to_time_t( now );
    char timebuf[32];
#ifdef _WIN32
    struct tm tm_info;
    localtime_s( &tm_info, &t );
    strftime( timebuf, sizeof(timebuf), "%Y-%m-%d %H:%M:%S", &tm_info );
#else
    struct tm tm_info;
    localtime_r( &t, &tm_info );
    strftime( timebuf, sizeof(timebuf), "%Y-%m-%d %H:%M:%S", &tm_info );
#endif

    log << "[" << timebuf << "]"
        << " SAMG licence not found"
        << " | Attempt " << attempt << "/" << maxAttempts
        << "\n";

} // end LogLicenceRetry


/**
    Maps the SAMG solution vector u_ back into the user-supplied vector x,
    applying the point-based index transformation if required.

    @param x  solution vector to be updated
*/
void SAMG_Solver::MapSolutionToOutput( vector<double>& x ) const
{
    if ( settings_->UsePointBasedApproach() )
        for ( int32_t i = 0; i < nnu_; ++i )
            x[i % nsys_ * (nnu_ / nsys_) + i / nsys_] = u_[i];
    else
        x = u_;
}

// ============================================================================
// SOLVE: CompressedRowMatrix
// ============================================================================

/**
    Solves the linear system A x = b using SAMG, where A is supplied
    as a CompressedRowMatrix (already in SAMG-compatible format).

    @param A            solution matrix in compressed row storage
    @param b            right-hand side vector
    @param x            on entry: initial guess; on exit: solution
    @param no_unknowns  number of coupled solution variables (1 for scalar)
*/
void SAMG_Solver::SolveMatrixEquation( CompressedRowMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
{
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    const size_t new_nnu  = A.ia.size() - 1;
    const size_t new_nna  = A.ja.size();
    const size_t new_nsys = no_unknowns;

    u_.resize( new_nnu );
    f_.resize( new_nnu );

    crmat_ = A;
    crmat_.ConvertToSAMGFormat();

    if ( Verbose() )
        cout << "\nSAMG_Solver::SolveMatrixEquation (CompressedRowMatrix): "
             << "NNU=" << new_nnu << " NNA=" << new_nna << "\n";

    // Safety check — retained from original
    if ( crmat_.a.empty() || crmat_.ja.empty() || crmat_.ia.empty() ||
         u_.empty() || f_.empty() )
        throw csmp::Exception( FATAL_ERROR, "SAMG_Solver::SolveMatrixEquation",
                               "Unable to allocate required memory for transfer arrays" );

    InitialiseIndexArrays( new_nnu, new_nsys );

    if ( new_nsys != static_cast<size_t>(nsys_) )
        iscale_.assign( new_nsys, 0 );

    nsys_ = static_cast<int32_t>( new_nsys );
    nna_  = static_cast<int32_t>( new_nna );
    nnu_  = static_cast<int32_t>( new_nnu );

    InitialiseVectors( b, x );

    // Retrieve SAMG parameters
    int32_t matrix, nsolve, ifirst, ncyc, iswtch, idump, iout, ncg, levelx;
    double  eps, a_cmplx, g_cmplx, p_cmplx, w_avrge, chktol;
    RetrieveSAMGParameters( matrix, nsolve, ifirst, eps, ncyc, iswtch,
                            a_cmplx, g_cmplx, p_cmplx, w_avrge,
                            chktol, idump, iout, ncg, levelx );

    const int32_t instance = settings_->GetSolverInstance();
    SetupHiddenParameters( idump, ncg, levelx, instance );

    ierr_    = 0;
    res_in_  = -1.;
    res_out_ = -1.;

    if ( Verbose() )
        cout << "\n*** SAMG_Solver: Calling SAMG instance " << instance
             << " (nnu=" << nnu_ << ") ***\n" << endl;

#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    CallSAMGWithLicenceRetry( instance, nsolve, ifirst, eps, ncyc, iswtch,
                              a_cmplx, g_cmplx, p_cmplx, w_avrge,
                              chktol, idump, iout, matrix );

#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif

    if ( ierr_ > 0 )
        csmp_error.Note( ERROR, "SAMG_Solver::SolveMatrixEquation",
                         "SAMG solver returned with an error; error code: ",
                         to_string(ierr_).c_str() );
    else if ( ierr_ < 0 && ierr_ != -841 )  // -841 = BiCGStab restart, not an error
        csmp_error.Note( WARNING, "SAMG_Solver::SolveMatrixEquation",
                         "SAMG solver returned with a warning; code: ",
                         to_string(ierr_).c_str() );

    CheckConvergence( eps );
    MapSolutionToOutput( x );

    if ( output_amg_data_to_text_files_ )
        Write_SAMG_TextInputFile( "csp_solution_data" );

} // end SolveMatrixEquation (CompressedRowMatrix)

// ============================================================================
// SOLVE: SparseMatrix
// ============================================================================

/**
    Solves the linear system A x = b using SAMG, where A is supplied
    as a CSMP SparseMatrix. The matrix is converted internally to
    CompressedRowMatrix before being passed to SAMG.

    @param A            solution matrix in CSMP sparse format
    @param b            right-hand side vector
    @param x            on entry: initial guess; on exit: solution
    @param no_unknowns  number of coupled solution variables (1 for scalar)
*/
void SAMG_Solver::SolveMatrixEquation( SparseMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
{
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    const size_t new_nnu  = A.Rows();
    const size_t new_nna  = A.Entries();
    const size_t new_nsys = no_unknowns;

    u_.resize( new_nnu );
    f_.resize( new_nnu );

    // Convert CSMP SparseMatrix to SAMG-compatible CompressedRowMatrix
    if ( settings_->UsePointBasedApproach() )
        crmat_.InitializePointBasedSAMG( A, new_nsys );
    else
        crmat_.Initialize( A );

    if ( Verbose() )
        cout << "\nSAMG_Solver::SolveMatrixEquation (SparseMatrix): "
             << "NNU=" << new_nnu << " NNA=" << new_nna << "\n";

    // Safety check — retained from original
    if ( crmat_.a.empty() || crmat_.ja.empty() || crmat_.ia.empty() ||
         u_.empty() || f_.empty() )
        throw csmp::Exception( FATAL_ERROR, "SAMG_Solver::SolveMatrixEquation",
                               "Unable to allocate required memory for transfer arrays" );

    InitialiseIndexArrays( new_nnu, new_nsys );

    if ( new_nsys != static_cast<size_t>(nsys_) )
        iscale_.assign( new_nsys, 0 );

    nsys_ = static_cast<int32_t>( new_nsys );
    nna_  = static_cast<int32_t>( new_nna );
    nnu_  = static_cast<int32_t>( new_nnu );

    InitialiseVectors( b, x );

    // Retrieve SAMG parameters
    int32_t matrix, nsolve, ifirst, ncyc, iswtch, idump, iout, ncg, levelx;
    double  eps, a_cmplx, g_cmplx, p_cmplx, w_avrge, chktol;
    RetrieveSAMGParameters( matrix, nsolve, ifirst, eps, ncyc, iswtch,
                            a_cmplx, g_cmplx, p_cmplx, w_avrge,
                            chktol, idump, iout, ncg, levelx );

    const int32_t instance = settings_->GetSolverInstance();
    SetupHiddenParameters( idump, ncg, levelx, instance );

    ierr_    = 0;
    res_in_  = -1.;
    res_out_ = -1.;

    if ( Verbose() )
        cout << "\n*** SAMG_Solver: Calling SAMG instance " << instance
             << " (nnu=" << nnu_ << ") ***\n" << endl;

#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    CallSAMGWithLicenceRetry( instance, nsolve, ifirst, eps, ncyc, iswtch,
                              a_cmplx, g_cmplx, p_cmplx, w_avrge,
                              chktol, idump, iout, matrix );

#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif

    if ( ierr_ > 0 )
        csmp_error.Note( ERROR, "SAMG_Solver::SolveMatrixEquation",
                         "SAMG solver returned with an error; error code: ",
                         to_string(ierr_).c_str() );
    else if ( ierr_ < 0 && ierr_ != -841 )  // -841 = BiCGStab restart, not an error
        csmp_error.Note( WARNING, "SAMG_Solver::SolveMatrixEquation",
                         "SAMG solver returned with a warning; code: ",
                         to_string(ierr_).c_str() );

    CheckConvergence( eps );
    MapSolutionToOutput( x );

    if ( output_amg_data_to_text_files_ )
        Write_SAMG_TextInputFile( "csp_solution_data" );

} // end SolveMatrixEquation (SparseMatrix)

// ============================================================================
// DIAGNOSTIC OUTPUT
// ============================================================================

void SAMG_Solver::Write_SAMG_TextOutput( bool write )
{
    output_amg_data_to_text_files_ = write;
}


void SAMG_Solver::OutputVectors() const
{
    const string out( "single_output.txt" );
    ofstream ofs( out.c_str() );
    ofs.setf( ios::scientific );
    const long prec = ofs.precision(15);
    ofs.precision( prec );

    ofs << "\n\nu:\n";
    for ( int32_t i = 0; i < nnu_; ++i )
        ofs << i+1 << ", " << u_[i] << "\n";

    ofs << "\n\nf:\n";
    for ( int32_t i = 0; i < nnu_; ++i )
        ofs << i+1 << ", " << f_[i] << "\n";
}


/**
    Outputs the SAMG input data to a set of text files for diagnostic purposes.

    Files written:
    - .frm  format descriptor and matrix dimensions
    - .amg  ia, ja, a arrays (compressed row storage)
    - .rhs  right-hand side vector f
    - .lhs  initial guess / solution vector u
    - .a    matrix entries with column indices, row by row
    - .iu   unknown type array (only for coupled systems)
    - .ip   point number array (only for point-based coupled systems)

    @param file  base filename (without extension)
    @return true if all files were written successfully
*/
bool SAMG_Solver::Write_SAMG_TextInputFile( const char* file ) const
{
    string   out_file( file );
    bool     complete_output( true );

    // -------------------------------------------------------------------------
    // .frm  format descriptor and matrix dimensions
    // -------------------------------------------------------------------------
    out_file += ".frm";
    ofstream ofs( out_file.c_str() );
    ofs << out_file << "\nCSMP output file created for SAMG test.\n";
    ofs << "f         4\n";
    ofs << "# NNA  NNU  MATRIX  NSYS  NPNT\n";
    ofs << nna_ << "  " << nnu_ << "  " << settings_->Get_matrix()
        << "  " << nsys_ << "  ";
    ofs << ( ndip_ > 1 ? 1 : 0 ) << "\n";
    ofs << "# ISCALE(1...NSYS)\n";
    for ( int32_t i = 0; i < nsys_; ++i ) ofs << 0 << "\n";
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
         << "' written successfully.\n";

    // -------------------------------------------------------------------------
    // .amg  ia(1...nnu+1), ja(1...nna), a(1...nna)
    // -------------------------------------------------------------------------
    assert( !crmat_.ia.empty() );
    assert( !crmat_.ja.empty() );
    assert( !crmat_.a.empty() );
    out_file  = file;
    out_file += ".amg";
    ofs.open( out_file.c_str() );
    for ( size_t i{0U}; i < static_cast<size_t>(nnu_) + 1U; ++i )
        ofs << crmat_.ia[i] << "\n";
    for ( size_t i{0U}; i < static_cast<size_t>(nna_); ++i )
        ofs << crmat_.ja[i] << "\n";
    ofs.setf( ios::scientific );
    const long prec_amg = ofs.precision(15);
    for ( size_t i{0U}; i < static_cast<size_t>(nna_); ++i )
        ofs << crmat_.a[i] << "\n";
    ofs.unsetf( ios::scientific );
    ofs.precision( prec_amg );
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
         << "' written successfully.\n";

    // -------------------------------------------------------------------------
    // .rhs  right-hand side vector f(1...nnu)
    // -------------------------------------------------------------------------
    assert( !f_.empty() );
    out_file  = file;
    out_file += ".rhs";
    ofs.open( out_file.c_str() );
    ofs.setf( ios::scientific );
    const long prec_rhs = ofs.precision(15);
    for ( size_t i{0U}; i < static_cast<size_t>(nnu_); ++i )
        ofs << f_[i] << "\n";
    ofs.unsetf( ios::scientific );
    ofs.precision( prec_rhs );
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
         << "' written successfully.\n";

    // -------------------------------------------------------------------------
    // .lhs  initial guess / solution vector u(1...nnu)
    // -------------------------------------------------------------------------
    assert( !u_.empty() );
    out_file  = file;
    out_file += ".lhs";
    ofs.open( out_file.c_str() );
    ofs.setf( ios::scientific );
    const long prec_lhs = ofs.precision(15);
    for ( size_t i{0U}; i < u_.size(); ++i )
        ofs << u_[i] << "\n";
    ofs.unsetf( ios::scientific );
    ofs.precision( prec_lhs );
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
         << "' written successfully.\n";

    // -------------------------------------------------------------------------
    // .a  matrix entries with column indices, row by row
    // -------------------------------------------------------------------------
    assert( !crmat_.a.empty() );
    out_file  = file;
    out_file += ".a";
    ofs.open( out_file.c_str() );
    ofs.setf( ios::scientific );
    const long prec_a = ofs.precision(15);
    int32_t idx{ 0 };
    for ( size_t row{0U}; row < static_cast<size_t>(nnu_); ++row )
    {
        ofs << "\n" << (row+1) << "\t";
        for ( int32_t j = crmat_.ia[row]; j != crmat_.ia[row+1]; ++j, ++idx )
            ofs << crmat_.a[idx] << "\t[" << crmat_.ja[idx] << "]\t";
    }
    ofs.unsetf( ios::scientific );
    ofs.precision( prec_a );
    ofs.close();
    cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
         << "' written successfully.\n";

    // -------------------------------------------------------------------------
    // .iu  unknown type array (only for coupled systems)
    // -------------------------------------------------------------------------
    if ( nsys_ > 1 && ndiu_ > 1 ) {
        assert( !iu_.empty() );
        out_file  = file;
        out_file += ".iu";
        ofs.open( out_file.c_str() );
        for ( int32_t i = 0; i < ndiu_; ++i )
            ofs << iu_[i] << "\n";
        ofs.close();
        cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
             << "' written successfully.\n";
    }

    // -------------------------------------------------------------------------
    // .ip  point number array (only for point-based coupled systems)
    // -------------------------------------------------------------------------
    if ( nsys_ > 1 && ndip_ > 1 ) {
        assert( !ip_.empty() );
        out_file  = file;
        out_file += ".ip";
        ofs.open( out_file.c_str() );
        for ( int32_t i = 0; i < ndip_; ++i )
            ofs << ip_[i] << "\n";
        ofs.close();
        cout << "\nSAMG_Solver::Write_SAMG_TextInputFile: '" << out_file
             << "' written successfully.\n";
    }

    cout.flush();
    return complete_output;

} // end Write_SAMG_TextInputFile

// ============================================================================
// CONVERGENCE AND CYCLE CRITERIA
// ============================================================================

double SAMG_Solver::LastSolverResidual() const
{
    return res_out_;
}


/**
    Checks whether the matrix sparsity is compatible with the requested
    number of AMG levels. If levelx == 1 (single-level, i.e. Gauss-Seidel
    preconditioned BiCGStab without coarsening) but the matrix is not
    nearly diagonal, levelx is reset to its default value of 25.

    This check is only meaningful for advective transport problems.
*/
void SAMG_Solver::CheckSparsityCriterion( int32_t& levelx ) const
{
    if ( levelx == 1 && nna_ / nnu_ > 1.5 ) {
        cout << "\n\n*** SAMG_Solver::CheckSparsityCriterion: "
             << "Matrix sparsity criterion not fulfilled; levelx reset to default. ***\n"
             << "\t(Single-level setup should only be used for essentially diagonal matrices.)\n"
             << "\tlevelx set to: ";
        levelx = 25;
        cout << levelx << "\n";
    }
}


/**
    Checks whether the number of iteration cycles from the previous solve
    exceeds 1.5 times the best recorded count. If so, forces a new
    coarsening setup by incrementing iswtch.

    Only called when NO_PRIMARY_SOLVER_CONTROL is defined.
*/
void SAMG_Solver::CheckCycleCriterion( int32_t& iswtch ) const
{
    if ( ncyc_done_ > 1.5 * ncyc_best_ ) {
        iswtch += 1000000;
        cout << "\n\n*** SAMG_Solver::CheckCycleCriterion: "
             << "Iteration cycle criterion failed; new coarsening setup forced. ***\n";
        cout.flush();
    }
}


/**
    Updates the best recorded cycle count after a successful solve.

    Only called when NO_PRIMARY_SOLVER_CONTROL is defined.
*/
void SAMG_Solver::UpdateCycleCriterion( int32_t iswtch )
{
    if ( iswtch < 4000000 )
        ncyc_best_ = min( ncyc_best_, ncyc_done_ );
    else
        ncyc_best_ = ncyc_done_;
}


/**
    Checks whether the solver converged to within the requested tolerance.

    - If eps == 0: convergence is assumed (absolute criterion disabled).
    - If eps < 0:  absolute convergence criterion: res_out <= |eps|.
    - If eps > 0:  relative convergence criterion: res_out/res_in <= eps.

    @param eps  convergence tolerance as passed to SAMG
    @return true if converged, false otherwise
*/
bool SAMG_Solver::CheckConvergence( double eps ) const
{
    // Absolute criterion disabled
    if ( fabs(eps) < 1.0e-30 ) return true;

    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    cout << "\n*** SAMG_Solver::CheckConvergence: "
         << "eps=" << eps
         << ", res_out=" << res_out_
         << ", res_in=" << res_in_ << " ***\n";
    cout.flush();

    // Absolute convergence: negative eps
    if ( eps < 0. ) {
        if ( res_out_ > fabs(eps) ) {
            csmp_error.Note( WARNING,
                             "SAMG_Solver::CheckConvergence",
                             "Absolute error L2-norm solution criterion was not fulfilled" );
            return false;
        }
    }

    // Relative convergence: positive eps
    if ( eps > 0. ) {
        const double relative_residual = ( res_in_ == 0. ) ? 0. : res_out_ / res_in_;
        cout << "\n\trelative residual = " << relative_residual << "\n\n";
        if ( relative_residual > eps ) {
            csmp_error.Note( WARNING,
                             "SAMG_Solver::CheckConvergence",
                             "Relative solution criterion was not fulfilled" );
            return false;
        }
        cout.flush();
    }

    return true;

} // end CheckConvergence

} // end namespace csmp

// ============================================================================
// SAMG USER COORDINATE CALLBACK (required by SAMG library)
// ============================================================================

#ifndef CSP_PC_WINDOWS_NT_CODE_WARRIOR
void samg_user_coo ( int*, int* ndim, double*, double*, double* );
void samg_user_coo_( int*, int* ndim, double*, double*, double* );

void samg_user_coo( int*, int* ndim, double*, double*, double* ) {
    cout << "\nSAMG_Solver::samg_user_coo: Setting *ndim to 0\n";
    *ndim = 0;
}
void samg_user_coo_( int*, int* ndim, double*, double*, double* ) {
    cout << "\nSAMG_Solver::samg_user_coo_: Setting *ndim to 0\n";
    *ndim = 0;
}
#endif

