#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "SAMG_Exception.h"
#include "CompressedRowMatrix.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include <fstream>

#ifndef SAMG_MULTIPLE_INSTANCES
#include "samg.h"
#else
#include "samg.h"
#include "samg1.h"
#include "samg2.h"
#include "samg3.h"
#include "samg4.h"
#include "samg5.h"
#endif


using namespace std;

namespace csmp {


SAMG_Solver::SAMG_Solver() 
    : Solver(new SAMG_Settings()),
      settings_(dynamic_cast<SAMG_Settings*>(solver_settings_)),
      nsys_(0),
      ndiu_(1),
      ndip_(1),
      ip_(1), // SKM FIX: MS Studio needs minimum size of 1 else it throws assertion in vector-based approach as it calls &ip_[0]
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



SAMG_Solver::SAMG_Solver( SAMG_Settings* settings ) :
    Solver(settings),
    settings_(settings),
    nsys_(0),
    ndiu_(1),
    ndip_(1),
	ip_(1), // SKM FIX: MS Studio needs minimum size of 1 else it throws assertion in vector-based approach as it calls &ip_[0]
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
            this->solver_settings_=settings_;
          }
        else {
            settings_=solver.settings_;
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

/** Set new solver settings.

@param settings Pointer to a SAMG_Settings object.

@section implementation Implementation

The new setting is stored as a pointer variable.

@warning The SAMG_Solver object DOES NOT delete ANY pointer to SAMG_Settings
objects given to it. The user is therefore solely responsible for guaranteeing
the destruction of SAMG_Settings objects created with new
*/
void SAMG_Solver::InputSolverSettings( SolverSettings& settings )
 {
    settings_ = static_cast<SAMG_Settings*>(&settings);
 }

SolverSettings* SAMG_Solver::GetSolverSettings() {
    return settings_;
}

void SAMG_Solver::SolveMatrixEquation( CompressedRowMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  /// 1. Calculate number of non-zero elements in A
  size_t new_nnu = A.ia.size() - 1;
  size_t new_nna = A.ja.size();
  size_t new_nsys = no_unknowns;
  u_.resize( new_nnu );
  f_.resize( new_nnu );

  crmat_ = A;

  if ( Verbose() )
    cout << "\nSAMGp_Solver::SolveMatrixEquation: Sizes NNU(rows=cols): "
    << new_nnu << " NNA(total unknowns): " << new_nna << "   ";

  // set arrays for systems
  if ( new_nsys > 1U ) {
    // only rebuild arrays if settings have changed
    if ( new_nnu != nnu_ || new_nsys != nsys_ ) {
      if ( new_nnu != nnu_ ) {
        iu_.resize( new_nnu );
        ndiu_ = static_cast<int32_t>(iu_.size()); // update the size indicator for iu vector
      }

      uint32_t k( 0U );
      // POINT BASED APPROACH (SAMG manual page 16)
      // coupled systems with solution vector [x1, y1, x2, y2, x3, y3, ...., xn, yn]
      // IU must gave the form [1,2,1,2,1,2,...,1,2]
      if ( settings_->UsePointBasedApproach() ) {
        for ( int32_t i = 0; i<(iu_.size() / new_nsys); i++ )
          for ( int32_t j = 0; j<new_nsys; j++ )
            iu_[k++] = j + 1;
      }
      // UNKOWN BASED APPROACH (SAMG manual page 16)
      // coupled systems with solution vector [x1, x2, x3,..., xn; y1, y2, y3, ...., yn]
      // IU must have the form [1,1,1,...,n,2,2,2,...,n]
      else {
        for ( int32_t j = 0; j<new_nsys; j++ )
          for ( int32_t i = 0; i<(iu_.size() / new_nsys); i++ )
            iu_[k++] = j + 1;
      }
    }
    // POINT BASED APPROACH (SAMG manual page 16)
    // coupled systems with solution vector [x1, y1, x2, y2, x3, y3, ...., xn, yn]
    if ( settings_->UsePointBasedApproach() ) {
      // only rebuild arrays if settings have changed
      if ( new_nnu != nnu_ || new_nsys != nsys_ ) {
        if ( new_nnu != nnu_ ) {
          ip_.resize( new_nnu );
          ndip_ = static_cast<int32_t>(ip_.size()); // resizing the ip vector size indicator
        }

        uint32_t  k( 0 );
        // IP must be the corresponding node number format [1,1,2,2,3,3,...,n,n]
        for ( int32_t i = 0; i<(ip_.size() / new_nsys); i++ )
          for ( int32_t j = 0; j<new_nsys; j++ )
            ip_[k++] = i + 1;
      }
    }
  }

  else { // new_nsys < 1 or new_nsys=1
         //      if ( ndiu_ != 1 ) iu_.resize(1);
         //      if ( ndip_ != 1 ) ip_.resize(1);
         // These four lines below are left in case a solver might be used to solve more than one variable
         // and later it is used again to solve a single one (i.e. new_nsys=1).
    ndiu_ = 1;
    ndip_ = 1;
    iu_.resize( 1 );
    ip_.resize( 1 );
  }

  // To be done for scalar systems as well

  if ( new_nsys != nsys_ ) {
    iscale_.resize( new_nsys );
    // putting zero values into this array switches the scaling off
    fill( iscale_.begin(), iscale_.end(), 0 );
  }

  // assign new values to actual ones
  nsys_ = static_cast<int32_t>(new_nsys);
  nna_ = static_cast<int32_t>(new_nna);
  nnu_ = static_cast<int32_t>(new_nnu);

  // now the righthand and solution vectors are initialized
  if ( settings_->UsePointBasedApproach() ) {
#if defined(_OPENMP )
#pragma omp parallel for // algorithm has been done this way to complete idea of "first touch".
#endif
    for ( int32_t i = 0U; i < nnu_; i++ ) {
      u_[i] = x[i%nsys_*(nnu_ / nsys_) + i / nsys_]; // initial guess for the solution vector
      f_[i] = b[i%nsys_*(nnu_ / nsys_) + i / nsys_]; // right-hand side
    }
  }
  else {
#if defined(_OPENMP )
#pragma omp parallel for // algorithm has been done this way to complete idea of "first touch".
#endif
    for ( int32_t i = 0U; i<nnu_; i++ ) {
      u_[i] = x[i]; // initial guess for the solution vector
      f_[i] = b[i]; // right-hand side
    }
  }

  if ( crmat_.a.empty() ) cout << "\n a is empty." << endl;
  if ( crmat_.ja.empty() ) cout << "\n ja is empty." << endl;
  if ( crmat_.ia.empty() ) cout << "\n ia is empty." << endl;
  if ( x.empty() ) cout << "\n x is empty." << endl;
  if ( b.empty() ) cout << "\n b is empty." << endl;

  // test for empty vectors
  if ( crmat_.a.empty() || crmat_.ja.empty() || crmat_.ia.empty() || u_.empty() || f_.empty() )
    throw csmp::Exception( FATAL_ERROR,
                           "SAMG_Solver::SolveMatrixEquation",
                           "Unable to allocate required memory for transfer arrays" );

  if ( Verbose() ) {
    cout << "\nSAMG_Solver::SolveSAMG: Calling SAMG ..." << endl;
    cout.flush();
  }

  /// Initialize settings
  /// Output parameter
  res_in_ = -1.; // by default (ntake_res_in = 0 ) res_in is an output parameter
  res_out_ = -1.;

  /// Paramater from SAMG_Settings object
  int32_t matrix = settings_->Get_matrix();
  int32_t nsolve = settings_->Get_nsolve();
  int32_t ifirst = settings_->Get_ifirst();
  double eps = settings_->Get_eps();
  int32_t ncyc = settings_->Get_ncyc();
  int32_t iswtch = settings_->Get_iswtch();
  double a_cmplx = settings_->Get_a_cmplx();
  double g_cmplx = settings_->Get_g_cmplx();
  double p_cmplx = settings_->Get_p_cmplx();
  double w_avrge = settings_->Get_w_avrge();
  double chktol = settings_->Get_chktol();
  int32_t idump = settings_->Get_idump();
  int32_t iout = settings_->Get_iout();
  int32_t mode_mess = settings_->Get_mode_mess();
  int32_t nrd = settings_->Get_nrd();
  int32_t nru = settings_->Get_nru();
  int32_t ncg = settings_->Get_ncg();

  /// SAMG matrix output to file
  int32_t ioform = settings_->Get_ioform();                          // matrix output format parameter (ASCII characters)
  int32_t ioform_length = settings_->Get_ioform_length();            // matrix output format parameter lenght (number of ASCII characters)
  int* filnam_dump = settings_->Get_filnam_dump();                 // matrix output format filename (ASCII characters)
  int32_t filnam_dump_length = settings_->Get_filnam_dump_length();  // matrix output format filename lenght (number of ASCII characters)


                                                                   /// Define SAMG hidden parameters
  if ( settings_->ExplicitSecondary() ) {

    /// Define number of SAMG levels
    int32_t levelx = settings_->Get_levelx();

#ifndef SAMG_MULTIPLE_INSTANCES
    SAMG_SET_NCG( &ncg );

    /// Check applicability of renounced coarsening
#ifdef RENOUNCE_COARSENING
    CheckSparsityCriterion( levelx );
#endif
    SAMG_SET_LEVELX( &levelx );

    /// SAMG output to file
    if ( settings_->Get_idmp() > 1 ) {
      SAMG_ISET_IOFORM( &ioform, &ioform_length );
      SAMG_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
    }
    SAMG_SET_MODE_MESS( &mode_mess );
    SAMG_SET_NRD( &nrd );
    SAMG_SET_NRU( &nru );
#endif

#ifdef SAMG_MULTIPLE_INSTANCES
    if ( settings_->GetSolverInstance() == 0 ) {
      SAMG_SET_NCG( &ncg );

      /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
      CheckSparsityCriterion( levelx );
#endif
      SAMG_SET_LEVELX( &levelx );

      /// SAMG output to file
      if ( settings_->Get_idmp() > 1 ) {
        SAMG_ISET_IOFORM( &ioform, &ioform_length );
        SAMG_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
      }
      SAMG_SET_MODE_MESS( &mode_mess );
      SAMG_SET_NRD( &nrd );
      SAMG_SET_NRU( &nru );
    }
    else if ( settings_->GetSolverInstance() == 1 ) {
      SAMG1_SET_NCG( &ncg );

      /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
      CheckSparsityCriterion( levelx );
#endif
      SAMG1_SET_LEVELX( &levelx );

      /// SAMG output to file
      if ( settings_->Get_idmp() > 1 ) {
        SAMG1_ISET_IOFORM( &ioform, &ioform_length );
        SAMG1_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
      }
      SAMG1_SET_MODE_MESS( &mode_mess );
      SAMG1_SET_NRD( &nrd );
      SAMG1_SET_NRU( &nru );
    }
    else if ( settings_->GetSolverInstance() == 2 ) {
      SAMG2_SET_NCG( &ncg );

      /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
      CheckSparsityCriterion( levelx );
#endif
      SAMG2_SET_LEVELX( &levelx );

      /// SAMG output to file
      if ( settings_->Get_idmp() > 1 ) {
        SAMG2_ISET_IOFORM( &ioform, &ioform_length );
        SAMG2_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
      }
      SAMG2_SET_MODE_MESS( &mode_mess );
      SAMG2_SET_NRD( &nrd );
      SAMG2_SET_NRU( &nru );
    }
    else if ( settings_->GetSolverInstance() == 3 ) {
      SAMG3_SET_NCG( &ncg );

      /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
      CheckSparsityCriterion( levelx );
#endif
      SAMG3_SET_LEVELX( &levelx );

      /// SAMG output to file
      if ( settings_->Get_idmp() > 1 ) {
        SAMG3_ISET_IOFORM( &ioform, &ioform_length );
        SAMG3_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
      }
      SAMG3_SET_MODE_MESS( &mode_mess );
      SAMG3_SET_NRD( &nrd );
      SAMG3_SET_NRU( &nru );
    }
    else if ( settings_->GetSolverInstance() == 4 ) {
      SAMG4_SET_NCG( &ncg );

      /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
      CheckSparsityCriterion( levelx );
#endif
      SAMG4_SET_LEVELX( &levelx );

      /// SAMG output to file
      if ( settings_->Get_idmp() > 1 ) {
        SAMG4_ISET_IOFORM( &ioform, &ioform_length );
        SAMG4_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
      }
      SAMG4_SET_MODE_MESS( &mode_mess );
      SAMG4_SET_NRD( &nrd );
      SAMG4_SET_NRU( &nru );
    }
    else if ( settings_->GetSolverInstance() == 5 ) {
      SAMG5_SET_NCG( &ncg );

      /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
      CheckSparsityCriterion( levelx );
#endif
      SAMG5_SET_LEVELX( &levelx );

      /// SAMG output to file
      if ( settings_->Get_idmp() > 1 ) {
        SAMG5_ISET_IOFORM( &ioform, &ioform_length );
        SAMG5_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
      }
      SAMG5_SET_MODE_MESS( &mode_mess );
      SAMG5_SET_NRD( &nrd );
      SAMG5_SET_NRU( &nru );
    }
#endif

  }
  else {

#ifndef SAMG_OLD_INTERFACE
#ifndef SAMG_MULTIPLE_INSTANCES
    SAMG_RESET_HIDDEN();
#else
    if ( settings_->GetSolverInstance() == 0 ) {
      SAMG_RESET_HIDDEN();
    }
    else if ( settings_->GetSolverInstance() == 1 ) {
      SAMG1_RESET_HIDDEN();
    }
    else if ( settings_->GetSolverInstance() == 2 ) {
      SAMG2_RESET_HIDDEN();
    }
    else if ( settings_->GetSolverInstance() == 3 ) {
      SAMG3_RESET_HIDDEN();
    }
    else if ( settings_->GetSolverInstance() == 4 ) {
      SAMG4_RESET_HIDDEN();
    }
    else if ( settings_->GetSolverInstance() == 5 ) {
      SAMG5_RESET_HIDDEN();
    }

#endif
#else
    SAMG_RESET_SECONDARY();
#endif
  }

  ierr_ = 0;

#ifndef SAMG_MULTIPLE_INSTANCES
  cout << "\n\n*** SAMG_Solver::SolveMatrixEquation: 'Calling SAMG( nnu = " << nnu_ << ", nna = " << nna_ << ", ... ) ***\n\n";
  cout.flush();

  /// IMPES without SAMG Multiple Instances
#ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
  /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
  CheckCycleCriterion( iswtch );
#endif
#endif

  SAMG( &nnu_, &nna_, &nsys_,
        &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
        &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
        &res_in_, &res_out_, &ncyc_done_, &ierr_,
        &nsolve, &ifirst, &eps, &ncyc, &iswtch,
        &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
        &chktol, &idump, &iout );

  /// IMPES without SAMG Multiple Instances
#ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
  /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
  UpdateCycleCriterion( iswtch );
#endif
#endif
#endif

#ifdef SAMG_MULTIPLE_INSTANCES
  if ( Verbose() ) cout << "\n*** SAMG_Solver::SolveMatrixEquation: Calling SAMG instance: " << settings_->GetSolverInstance() << " ***\n\n";
  cout.flush();

  if ( settings_->GetSolverInstance() == 0 ) {

    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    /*
    if (ip_.empty()) {
    int* ip(0);
    SAMG(&nnu_, &nna_, &nsys_,
    &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
    &iu_[0], &ndiu_,
    -> !    ip,
    &ndip_, &matrix, &iscale_[0],
    &res_in_, &res_out_, &ncyc_done_, &ierr_,
    &nsolve, &ifirst, &eps, &ncyc, &iswtch,
    &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
    &chktol, &idump, &iout);
    }
    else
    */
    //crmat_.Out("test-compressed-row-matrix.txt");

    SAMG( &nnu_, &nna_, &nsys_,
          &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
          &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
          &res_in_, &res_out_, &ncyc_done_, &ierr_,
          &nsolve, &ifirst, &eps, &ncyc, &iswtch,
          &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
          &chktol, &idump, &iout );

    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif

  }
  else if ( settings_->GetSolverInstance() == 1 ) {

    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    SAMG1( &nnu_, &nna_, &nsys_,
           &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
           &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
           &res_in_, &res_out_, &ncyc_done_, &ierr_,
           &nsolve, &ifirst, &eps, &ncyc, &iswtch,
           &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
           &chktol, &idump, &iout );

    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif

  }
  else if ( settings_->GetSolverInstance() == 2 ) {

    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    SAMG2( &nnu_, &nna_, &nsys_,
           &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
           &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
           &res_in_, &res_out_, &ncyc_done_, &ierr_,
           &nsolve, &ifirst, &eps, &ncyc, &iswtch,
           &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
           &chktol, &idump, &iout );

    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif

  }
  else if ( settings_->GetSolverInstance() == 3 ) {
    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    SAMG3( &nnu_, &nna_, &nsys_,
           &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
           &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
           &res_in_, &res_out_, &ncyc_done_, &ierr_,
           &nsolve, &ifirst, &eps, &ncyc, &iswtch,
           &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
           &chktol, &idump, &iout );

    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif
  }
  else if ( settings_->GetSolverInstance() == 4 ) {
    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    SAMG4( &nnu_, &nna_, &nsys_,
           &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
           &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
           &res_in_, &res_out_, &ncyc_done_, &ierr_,
           &nsolve, &ifirst, &eps, &ncyc, &iswtch,
           &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
           &chktol, &idump, &iout );

    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif

  }
  else if ( settings_->GetSolverInstance() == 5 ) {
    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif

    SAMG5( &nnu_, &nna_, &nsys_,
           &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
           &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
           &res_in_, &res_out_, &ncyc_done_, &ierr_,
           &nsolve, &ifirst, &eps, &ncyc, &iswtch,
           &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
           &chktol, &idump, &iout );

    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif
  }
  else
    throw csmp::Exception( ERROR, "SAMG_Solver::SolveMatrixEquation:",
                           "Desired solver instance is not available in current SAMG library" );
#endif

  if ( ierr_ > 0 ) {
    csmp_error.notice( ERROR, "SAMG_Solver::SolveMatrixEquation: ",
                       "SAMG solver returned with an error; error code: ", (to_string( ierr_ )).c_str() );

  }
  else if ( ierr_ < 0 and ierr_ != -841 ) { // bicgstab restart
    csmp_error.notice( WARNING, "SAMG_Solver::SolveMatrixEquation:",
                       "SAMG solver returned with a warning; code: ", (to_string( ierr_ )).c_str() );
  }

  /// SAMG convergence check
  CheckConvergence( eps );

  // solver returned ok so lets place contents back into x
  if ( settings_->UsePointBasedApproach() )
    for ( auto i = 0U; i < nnu_; i++ )
      x[i%nsys_*(nnu_ / nsys_) + i / nsys_] = u_[i];
  else x = u_;

  // 2. writing SAMG solver input/output data to file
  // -----------------------------------------
  if ( output_amg_data_to_text_files_ ) {
    cout << "\nSAMG_Solver): SAMG TEXT FILE OUTPUT HAS BEEN ENABLED ! "
      << "Watch for '.frm', '.amg', '.rhs' and perhaps '.iu' and 'ip' files that will be written." << endl;
    Write_SAMG_TextInputFile( "csp_solution_data" );
  }
}


/**
    CSMP's sparse matrix is converted into a temporary
    CompressedRowStorage which is handed over to SAMG that solves it.
    The solution is returned into the vector x.
    
    Only 2 DLL instances (0->not numbered and (1) are supported.
*/
void SAMG_Solver::SolveMatrixEquation( SparseMatrix& A, 
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    /// 1. Calculate number of non-zero elements in A
    size_t new_nnu = A.Rows();
    size_t new_nna = A.Entries();
    size_t new_nsys = no_unknowns;
    u_.resize(new_nnu);
    f_.resize(new_nnu);

    if (Verbose()) cout << "\n\n*** SAMG_Solver::SolveMatrixEquation: 'Created CompressedRowMatrix ...' ***\n\n";
    cout.flush();

    if ( settings_->UsePointBasedApproach() )
        crmat_.InitializePointBased( A, new_nsys );
    else
        crmat_.Initialize( A );

    if ( Verbose() )
        cout << "\nSAMGp_Solver::SolveMatrixEquation: Sizes NNU(rows=cols): "
             << new_nnu << " NNA(total unknowns): " << new_nna << "   ";

    // set arrays for systems
    if ( new_nsys > 1U ) {
        // only rebuild arrays if settings have changed
        if ( new_nnu != nnu_ || new_nsys != nsys_ ) {
            if ( new_nnu != nnu_ ) {
                iu_.resize(new_nnu);
                ndiu_ = static_cast<int32_t>(iu_.size()); // update the size indicator for iu vector
            }

            uint32_t k(0U);
            // POINT BASED APPROACH (SAMG manual page 16)
            // coupled systems with solution vector [x1, y1, x2, y2, x3, y3, ...., xn, yn]
            // IU must gave the form [1,2,1,2,1,2,...,1,2]
            if ( settings_->UsePointBasedApproach() ) {
                for ( int32_t i=0; i<(iu_.size()/new_nsys); i++ )
                    for ( int32_t j=0; j<new_nsys; j++ )
                        iu_[k++] = j+1;
            }
            // UNKOWN BASED APPROACH (SAMG manual page 16)
            // coupled systems with solution vector [x1, x2, x3,..., xn; y1, y2, y3, ...., yn]
            // IU must have the form [1,1,1,...,n,2,2,2,...,n]
            else {
                for ( int32_t j=0; j<new_nsys; j++ )
                    for ( int32_t i=0; i<(iu_.size()/new_nsys); i++ )
                        iu_[k++] = j+1;
            }
        }
        // POINT BASED APPROACH (SAMG manual page 16)
        // coupled systems with solution vector [x1, y1, x2, y2, x3, y3, ...., xn, yn]
        if ( settings_->UsePointBasedApproach() ) {
            // only rebuild arrays if settings have changed
            if ( new_nnu != nnu_ || new_nsys != nsys_ ) {
                if ( new_nnu != nnu_ ) {
                    ip_.resize(new_nnu);
                    ndip_ = static_cast<int32_t>(ip_.size()); // resizing the ip vector size indicator
                }

                uint32_t  k(0);
                // IP must be the corresponding node number format [1,1,2,2,3,3,...,n,n]
                for ( int32_t i=0; i<(ip_.size()/new_nsys); i++ )
                    for ( int32_t j=0; j<new_nsys; j++ )
                        ip_[k++] = i+1;
            }
        }
    }

    else { // new_nsys < 1 or new_nsys=1
        //      if ( ndiu_ != 1 ) iu_.resize(1);
        //      if ( ndip_ != 1 ) ip_.resize(1);
        // These four lines below are left in case a solver might be used to solve more than one variable
        // and later it is used again to solve a single one (i.e. new_nsys=1).
        ndiu_=1;
        ndip_=1;
        iu_.resize(1);
        ip_.resize(1);
    }
    
    // To be done for scalar systems as well

    if (new_nsys != nsys_) {
        iscale_.resize(new_nsys);
        // putting zero values into this array switches the scaling off
        fill( iscale_.begin(), iscale_.end(), 0 );
    }

    // assign new values to actual ones
    nsys_ = static_cast<int32_t>(new_nsys);
    nna_  = static_cast<int32_t>(new_nna);
    nnu_  = static_cast<int32_t>(new_nnu);

    // now the righthand and solution vectors are initialized
    if ( settings_->UsePointBasedApproach() ) {
#if defined(_OPENMP )
#pragma omp parallel for // algorithm has been done this way to complete idea of "first touch".
#endif
        for ( int32_t i = 0U; i < nnu_; i++ ) {
            u_[i] = x[i%nsys_*(nnu_/nsys_)+i/nsys_]; // initial guess for the solution vector
            f_[i] = b[i%nsys_*(nnu_/nsys_)+i/nsys_]; // right-hand side
        }
    }
    else {
#if defined(_OPENMP )
#pragma omp parallel for // algorithm has been done this way to complete idea of "first touch".
#endif
        for ( int32_t i=0U; i<nnu_; i++ ) {
            u_[i] = x[i]; // initial guess for the solution vector
            f_[i] = b[i]; // right-hand side
        }
    }

    if ( crmat_.a.empty()) cout <<"\n a is empty." <<endl;
    if ( crmat_.ja.empty()) cout <<"\n ja is empty." <<endl;
    if ( crmat_.ia.empty()) cout <<"\n ia is empty." <<endl;
    if ( x.empty()) cout <<"\n x is empty." <<endl;
    if ( b.empty()) cout <<"\n b is empty." <<endl;

    // test for empty vectors
    if( crmat_.a.empty() || crmat_.ja.empty() || crmat_.ia.empty() || u_.empty() || f_.empty() )
        throw csmp::Exception( FATAL_ERROR,
                               "SAMG_Solver::SolveMatrixEquation",
                               "Unable to allocate required memory for transfer arrays");

    if ( Verbose() ) {
        cout << "\nSAMG_Solver::SolveSAMG: Calling SAMG ..." << endl;
        cout.flush();
    }

    /// Initialize settings
    /// Output parameter
    res_in_       = -1.; // by default (ntake_res_in = 0 ) res_in is an output parameter
    res_out_      = -1.;

    /// Paramater from SAMG_Settings object
    int32_t matrix = settings_->Get_matrix();
    int32_t nsolve = settings_->Get_nsolve();
    int32_t ifirst = settings_->Get_ifirst();
    double eps = settings_->Get_eps();
    int32_t ncyc = settings_->Get_ncyc();
    int32_t iswtch = settings_->Get_iswtch();
    double a_cmplx = settings_->Get_a_cmplx();
    double g_cmplx = settings_->Get_g_cmplx();
    double p_cmplx = settings_->Get_p_cmplx();
    double w_avrge = settings_->Get_w_avrge();
    double chktol = settings_->Get_chktol();
    int32_t idump = settings_->Get_idump();
    int32_t iout = settings_->Get_iout();
    int32_t mode_mess = settings_->Get_mode_mess();
    int32_t nrd = settings_->Get_nrd();
    int32_t nru = settings_->Get_nru();
    int32_t ncg = settings_->Get_ncg();

    /// SAMG matrix output to file
    int32_t ioform = settings_->Get_ioform();                          // matrix output format parameter (ASCII characters)
    int32_t ioform_length = settings_->Get_ioform_length();            // matrix output format parameter lenght (number of ASCII characters)
    int* filnam_dump = settings_->Get_filnam_dump();                 // matrix output format filename (ASCII characters)
    int32_t filnam_dump_length = settings_->Get_filnam_dump_length();  // matrix output format filename lenght (number of ASCII characters)


    /// Define SAMG hidden parameters
    if (settings_->ExplicitSecondary()) {

        /// Define number of SAMG levels
        int32_t levelx = settings_->Get_levelx();

#ifndef SAMG_MULTIPLE_INSTANCES
        SAMG_SET_NCG( &ncg );

        /// Check applicability of renounced coarsening
#ifdef RENOUNCE_COARSENING
        CheckSparsityCriterion( levelx );
#endif
        SAMG_SET_LEVELX( &levelx );

        /// SAMG output to file
        if ( settings_->Get_idmp() > 1 ){
            SAMG_ISET_IOFORM( &ioform, &ioform_length );
            SAMG_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
        }
        SAMG_SET_MODE_MESS(&mode_mess);
        SAMG_SET_NRD(&nrd);
        SAMG_SET_NRU(&nru);
#endif

#ifdef SAMG_MULTIPLE_INSTANCES
        if ( settings_->GetSolverInstance() == 0 ) {
            SAMG_SET_NCG( &ncg );

            /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
            CheckSparsityCriterion( levelx );
#endif
            SAMG_SET_LEVELX( &levelx );

            /// SAMG output to file
            if ( settings_->Get_idmp() > 1 ){
                SAMG_ISET_IOFORM(&ioform, &ioform_length);
                SAMG_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
            }
            SAMG_SET_MODE_MESS(&mode_mess);
            SAMG_SET_NRD(&nrd);
            SAMG_SET_NRU(&nru);
        }
        else if ( settings_->GetSolverInstance() == 1 ) {
            SAMG1_SET_NCG(&ncg);

            /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
            CheckSparsityCriterion( levelx );
#endif
            SAMG1_SET_LEVELX( &levelx );

            /// SAMG output to file
            if ( settings_->Get_idmp() > 1 ){
                SAMG1_ISET_IOFORM( &ioform, &ioform_length );
                SAMG1_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
            }
            SAMG1_SET_MODE_MESS(&mode_mess);
            SAMG1_SET_NRD(&nrd);
            SAMG1_SET_NRU(&nru);
        }
        else if ( settings_->GetSolverInstance() == 2 ) {
            SAMG2_SET_NCG(&ncg);

            /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
            CheckSparsityCriterion( levelx );
#endif
            SAMG2_SET_LEVELX( &levelx );

            /// SAMG output to file
            if ( settings_->Get_idmp() > 1 ){
                SAMG2_ISET_IOFORM( &ioform, &ioform_length );
                SAMG2_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
            }
            SAMG2_SET_MODE_MESS(&mode_mess);
            SAMG2_SET_NRD(&nrd);
            SAMG2_SET_NRU(&nru);
        }
        else if ( settings_->GetSolverInstance() == 3 ) {
            SAMG3_SET_NCG(&ncg);

            /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
            CheckSparsityCriterion( levelx );
#endif
            SAMG3_SET_LEVELX( &levelx );

            /// SAMG output to file
            if ( settings_->Get_idmp() > 1 ){
                SAMG3_ISET_IOFORM( &ioform, &ioform_length );
                SAMG3_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
            }
            SAMG3_SET_MODE_MESS(&mode_mess);
            SAMG3_SET_NRD(&nrd);
            SAMG3_SET_NRU(&nru);
        }
         else if ( settings_->GetSolverInstance() == 4 ) {
            SAMG4_SET_NCG(&ncg);

            /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
            CheckSparsityCriterion( levelx );
#endif
            SAMG4_SET_LEVELX( &levelx );

            /// SAMG output to file
            if ( settings_->Get_idmp() > 1 ){
                SAMG4_ISET_IOFORM( &ioform, &ioform_length );
                SAMG4_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
            }
            SAMG4_SET_MODE_MESS(&mode_mess);
            SAMG4_SET_NRD(&nrd);
            SAMG4_SET_NRU(&nru);
        }
        else if ( settings_->GetSolverInstance() == 5 ) {
            SAMG5_SET_NCG(&ncg);

            /// Check applicability of reused coarsening setup
#ifdef RENOUNCE_COARSENING
            CheckSparsityCriterion( levelx );
#endif
            SAMG5_SET_LEVELX( &levelx );

            /// SAMG output to file
            if ( settings_->Get_idmp() > 1 ){
                SAMG5_ISET_IOFORM( &ioform, &ioform_length );
                SAMG5_ISET_FILNAM_DUMP( filnam_dump, &filnam_dump_length );
            }
            SAMG5_SET_MODE_MESS(&mode_mess);
            SAMG5_SET_NRD(&nrd);
            SAMG5_SET_NRU(&nru);
        }
#endif

    } else {

#ifndef SAMG_OLD_INTERFACE
#ifndef SAMG_MULTIPLE_INSTANCES
        SAMG_RESET_HIDDEN();
#else
        if ( settings_->GetSolverInstance() == 0 ) {
            SAMG_RESET_HIDDEN();
        } else if ( settings_->GetSolverInstance() == 1 ) {
            SAMG1_RESET_HIDDEN();
        }
        else if ( settings_->GetSolverInstance() == 2 ) {
            SAMG2_RESET_HIDDEN();
        }
        else if ( settings_->GetSolverInstance() == 3 ) {
            SAMG3_RESET_HIDDEN();
        }
        else if ( settings_->GetSolverInstance() == 4 ) {
            SAMG4_RESET_HIDDEN();
        }
        else if ( settings_->GetSolverInstance() == 5 ) {
            SAMG5_RESET_HIDDEN();
        }

#endif
#else
        SAMG_RESET_SECONDARY();
#endif
    }

    ierr_ = 0;

#ifndef SAMG_MULTIPLE_INSTANCES
    cout <<"\n\n*** SAMG_Solver::SolveMatrixEquation: 'Calling SAMG( nnu = " << nnu_ << ", nna = " << nna_ << ", ... ) ***\n\n";
    cout.flush();

    /// IMPES without SAMG Multiple Instances
#ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
    /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    CheckCycleCriterion( iswtch );
#endif
#endif

    SAMG( &nnu_, &nna_, &nsys_,
          &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
          &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
          &res_in_, &res_out_, &ncyc_done_, &ierr_,
          &nsolve, &ifirst, &eps, &ncyc, &iswtch,
          &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
          &chktol, &idump, &iout );

    /// IMPES without SAMG Multiple Instances
#ifdef IMPES_WITHOUT_SAMG_MULTIPLE_INSTANCES
    /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
    UpdateCycleCriterion( iswtch );
#endif
#endif
#endif

#ifdef SAMG_MULTIPLE_INSTANCES
    if (Verbose()) cout <<"\n*** SAMG_Solver::SolveMatrixEquation: Calling SAMG instance: "<< settings_->GetSolverInstance() <<" ***\n\n";
    cout.flush();

    if ( settings_->GetSolverInstance() == 0 ) {

        /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        CheckCycleCriterion( iswtch );
#endif

/*
		if (ip_.empty()) {
			int* ip(0);
			SAMG(&nnu_, &nna_, &nsys_,
				&crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
				&iu_[0], &ndiu_, 
-> !    ip,
        &ndip_, &matrix, &iscale_[0],
				&res_in_, &res_out_, &ncyc_done_, &ierr_,
				&nsolve, &ifirst, &eps, &ncyc, &iswtch,
				&a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
				&chktol, &idump, &iout);
		  }
		else
*/
//crmat_.Out("test-compressed-row-matrix.txt");

			SAMG( &nnu_, &nna_, &nsys_,
              &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
              &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
              &res_in_, &res_out_, &ncyc_done_, &ierr_,
              &nsolve, &ifirst, &eps, &ncyc, &iswtch,
              &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
              &chktol, &idump, &iout );

        /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        UpdateCycleCriterion( iswtch );
#endif

    } else if ( settings_->GetSolverInstance() == 1 ) {

        /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        CheckCycleCriterion( iswtch );
#endif

        SAMG1( &nnu_, &nna_, &nsys_,
               &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
               &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
               &res_in_, &res_out_, &ncyc_done_, &ierr_,
               &nsolve, &ifirst, &eps, &ncyc, &iswtch,
               &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
               &chktol, &idump, &iout );

        /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        UpdateCycleCriterion( iswtch );
#endif

    } else if ( settings_->GetSolverInstance() == 2 ) {

        /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        CheckCycleCriterion( iswtch );
#endif

        SAMG2( &nnu_, &nna_, &nsys_,
               &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
               &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
               &res_in_, &res_out_, &ncyc_done_, &ierr_,
               &nsolve, &ifirst, &eps, &ncyc, &iswtch,
               &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
               &chktol, &idump, &iout );

        /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        UpdateCycleCriterion( iswtch );
#endif

    } else if ( settings_->GetSolverInstance() == 3 ) {
        /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        CheckCycleCriterion( iswtch );
#endif

        SAMG3( &nnu_, &nna_, &nsys_,
               &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
               &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
               &res_in_, &res_out_, &ncyc_done_, &ierr_,
               &nsolve, &ifirst, &eps, &ncyc, &iswtch,
               &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
               &chktol, &idump, &iout );

        /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        UpdateCycleCriterion( iswtch );
#endif
    }  else if ( settings_->GetSolverInstance() == 4 ) {
        /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        CheckCycleCriterion( iswtch );
#endif

        SAMG4( &nnu_, &nna_, &nsys_,
               &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
               &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
               &res_in_, &res_out_, &ncyc_done_, &ierr_,
               &nsolve, &ifirst, &eps, &ncyc, &iswtch,
               &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
               &chktol, &idump, &iout );

        /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        UpdateCycleCriterion( iswtch );
#endif

    } else if ( settings_->GetSolverInstance() == 5 ) {
        /// Check applicability of reusing previous coarsening setup unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        CheckCycleCriterion( iswtch );
#endif

        SAMG5( &nnu_, &nna_, &nsys_,
               &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
               &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
               &res_in_, &res_out_, &ncyc_done_, &ierr_,
               &nsolve, &ifirst, &eps, &ncyc, &iswtch,
               &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
               &chktol, &idump, &iout );

        /// Collect number of iteration cycles from SAMG output unless primary solver control is enabled by iswit(>5)
#ifdef NO_PRIMARY_SOLVER_CONTROL
        UpdateCycleCriterion( iswtch );
#endif
    }
    else
        throw csmp::Exception( ERROR, "SAMG_Solver::SolveMatrixEquation:",
                               "Desired solver instance is not available in current SAMG library" );
#endif

    if ( ierr_ > 0 ) {
        csmp_error.notice( ERROR, "SAMG_Solver::SolveMatrixEquation: ",
                          "SAMG solver returned with an error; error code: ", (to_string(ierr_)).c_str() );

    } else if ( ierr_ < 0 and ierr_ != -841 ) { // bicgstab restart
        csmp_error.notice( WARNING, "SAMG_Solver::SolveMatrixEquation:",
                           "SAMG solver returned with a warning; code: ", (to_string(ierr_)).c_str() );
    }

    /// SAMG convergence check
    CheckConvergence( eps );

    // solver returned ok so lets place contents back into x
    if ( settings_->UsePointBasedApproach() )
        for ( auto i = 0U; i < nnu_; i++ )
            x[i%nsys_*(nnu_/nsys_)+i/nsys_] = u_[i];
    else x = u_;

    // 2. writing SAMG solver input/output data to file
    // -----------------------------------------
    if ( output_amg_data_to_text_files_ ) {
        cout << "\nSAMG_Solver): SAMG TEXT FILE OUTPUT HAS BEEN ENABLED ! "
             << "Watch for '.frm', '.amg', '.rhs' and perhaps '.iu' and 'ip' files that will be written."<< endl;
        Write_SAMG_TextInputFile( "csp_solution_data" );
    }
}


void  SAMG_Solver::Write_SAMG_TextOutput( bool write ) {
    output_amg_data_to_text_files_ = write;
}


void SAMG_Solver::OutputVectors() const
{
    string out("single_output");
    out += ".txt";
    ofstream ofs(out.c_str());
    ofs.setf(ios::scientific);
    long prec = ofs.precision(15);
    ofs.precision(prec);

    // output u and f to screen
    ofs<<"\n\nu:"<< endl;
    for (unsigned int i=0;i!=nnu_;i++)
        ofs<<i+1<<", "<<u_[i]<< endl;

    ofs<<"\n\nf:"<< endl;
    for (unsigned int i=0;i!=nnu_;i++)
        ofs<<i+1<<", "<<f_[i]<< endl;
}

/**

Outputs the input SparseMatrix into the desired format for the 
AMG solver's file-based interface. This entails the generation of
a format file (*.frm) and the files "*.amg", "*.rhs", "*.lhs", "*.a",
as well as the optional files "*.iu" and "*.ip",
see SAMG users guide for their meaning.  

*/
bool  SAMG_Solver::Write_SAMG_TextInputFile( const char* file ) const
{
    string    out_file(file);
    bool      complete_output(true);
    
    out_file += ".frm";
    
    ofstream  ofs(out_file.c_str());
    
    // -----------------------------------------------------
    // .frm  type of data format and sizes of pointer arrays
    // -----------------------------------------------------
    // Hierbei sind nna, nnu, matrix, nsys und der Vektor iscale
    // vom samg(...) interface her bekannt.
    // NPNT ist 1 oder 0 je nachdem, ob Punktinformationen vorhanden
    // sind oder nicht (d.h. ob der Vektor pi_ip die Laenge m_ndip=nnu hat und mit
    // sinnvollen Punktnummern gefuellt ist, oder ob pi_iu ein dummy vector
    // der Laenge m_ndip=1 ist).
    ofs << out_file <<"\nCSMP output file created for SAMG test."<< endl;
    ofs <<"f         4"<< endl;
    // in the following matrix is a composite parameter which indicates whether
    // the solution matrix is symmetric (1vs2) and whether the rowsums are 0 (1vs2)
    ofs <<"# NNA  NNU  MATRIX  NSYS  NPNT"<< endl;
    ofs << nna_ <<"  "<< nnu_  <<"  "<< settings_->Get_matrix() << "  " << nsys_ <<"  ";
    // if point information is present
    if ( ndip_ > 1 ) ofs << 1 << endl;
    else             ofs << 0 << endl;
    ofs <<"# ISCALE(1...NSYS)"<< endl;
    // no-scaling of any of the unknowns shall occur
    for ( int32_t i=0; i<nsys_; i++ ) ofs << 0 << endl;
    ofs.close();
    cout <<"\nSAMG_Solver::Write_SAMG_TextInputFile: file '"<< out_file;
    cout <<" written successfully."<< endl;
    
    // -----------------------------------------------------
    // .amg  ia(1...nnu+1), ja(1...nna), a(1...nna) condensed
    //       matrix storage indices and matrix
    // -----------------------------------------------------
    // enthaelt hintereinanderweg die drei Vektoren ia,ja,a
    // (bei Dir pi_ia, pi_ja, pd_a). Immer nur 1 Eintrag pro Zeile,
    // dh. nnu+1 + nna +nna Zeilen.
    assert( !crmat_.ia.empty() );
    assert( !crmat_.ja.empty() );
    assert( !crmat_.a.empty() );
    out_file  = file;
    out_file += ".amg";
    ofs.open(out_file.c_str());
    // this order is O.K. because the indices are not printed
    for ( auto i{0}; i<nnu_+1U; i++ ) ofs << crmat_.ia[i] << endl;
    for ( auto i{0}; i<nna_;    i++ ) ofs << crmat_.ja[i] << endl;
    ofs.setf(ios::scientific);
    long prec = ofs.precision(15);
    for ( auto i{0}; i<nna_;   i++ ) ofs << crmat_.a[i] << endl;
    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Solver::Write_SAMG_TextInputFile: file '"<< out_file;
    cout <<" written successfully."<< endl;

    // -----------------------------------------------------
    // .rhs		pd_f(1...nnu)
    // -----------------------------------------------------
    // wenn vorhanden enthaelt die right hand side (d.h. nnu Zeilen).
    assert( !f_.empty() );
    out_file  = file;
    out_file += ".rhs";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    for ( size_t i=0; i<nnu_; i++ ) ofs << f_[i] << endl;
    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Solver::Write_SAMG_TextInputFile: file '"<< out_file;
    cout <<" written successfully."<< endl;

    // -----------------------------------------------------
    // .lhs		pd_f(1...nnu)
    // -----------------------------------------------------
    // wenn vorhanden enthaelt die right hand side (d.h. nnu Zeilen).
    assert( !u_.empty() );
    out_file  = file;
    out_file += ".lhs";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    for ( size_t i=0; i<u_.size(); i++ ) ofs << u_[i] << endl;
    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Solver::Write_SAMG_TextInputFile: file '"<< out_file;
    cout <<" written successfully."<< endl;

    // -----------------------------------------------------
    // .a		pd_f(1...nnu)
    // -----------------------------------------------------
    // wenn vorhanden enthaelt die right hand side (d.h. nnu Zeilen).
    assert( !crmat_.a.empty() );
    out_file  = file;
    out_file += ".a";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    auto i{0};
    for ( size_t row=0U; row!=nnu_; row++ )
    {
        ofs <<"\n"<<(row+1)<<"\t"; // node number
        for ( int32_t j=crmat_.ia[row]; j!=crmat_.ia[row+1]; j++, i++ )
            ofs <<crmat_.a[i]<<"\t["<<crmat_.ja[i]<<"]\t";
    }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Solver::Write_SAMG_TextInputFile: file '"<< out_file;
    cout <<" written successfully."<< endl;


    // -----------------------------------------------------
    // .iu   iu(1...nnu)	for the ith variable, iu(i)
    //                      is the identifier of the physical
    //                   	unknown (1...nsys)
    // -----------------------------------------------------
    // notwendig im Falle eines PDE-Systems
    // enthaelt den Vektor iu (bei Dir pi_iu),
    // also (m_)ndiu Zeilen
    // Falls die Matrix nur von einer skalaren PDE
    // herkommt, fallen die Dateien <name>.iu und <name.ip> weg.
    if ( nsys_ > 1 && ndiu_ > 1 ) {
        assert( !iu_.empty() );
        out_file  = file;
        out_file += ".iu";
        ofs.open(out_file.c_str());
        for ( i=0; i<ndiu_; i++ ) ofs << iu_[i] << endl;
        ofs.close();
        cout <<"\nSAMG_Solver::Write_SAMG_TextInputFile: file '"<< out_file;
        cout <<" written successfully."<< endl;
    }

    // -----------------------------------------------------
    // .ip	 ip(1...nnu)	for the ith variable, ip(i) is the identifier of the node where
    //						variable i is located (this can only be used if the variables are
    //						ordered point-wise and, for every point, unknown wise and ascending.
    // -----------------------------------------------------
    // enthaelt den Vektor ip (bei Dir pi_ip),
    // also (m_)ndip Zeilen
    if ( nsys_ > 1 && ndip_ > 1 ) {
        assert( !ip_.empty() );
        out_file  = file;
        out_file += ".ip";
        ofs.open(out_file.c_str());
        for ( i=0; i<ndip_; i++ ) ofs << ip_[i] << endl;
        ofs.close();
        cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file;
        cout <<" written successfully."<< endl;
    }

    cout.flush();

    return complete_output;
    
} // end Write_SAMG_TextInputFile



double SAMG_Solver::LastSolverResidual() const {
    return res_out_;
}



/// Check matrix sparsity
/**
@section arguments Input Arguments

Gauss-Seidel preconditioned BiCGstab without coarsening (levelx = 1) is permitted
only for nealy diagonal matrices. SAMG hidden parameter levelx defines number of
coarsening levels. If SAMG levelx is set to 1, and the sparsity criterion is not fulfilled
( ratio of non-diagonal matrix entries to unknowns > 2), levelx is reset to its
default parameter 25.

Within CSMP++ this optimization feature can be only applied when solving for advective transport

*/
void SAMG_Solver::CheckSparsityCriterion( int32_t& levelx ) const {
    if ( levelx == 1 && nna_ / nnu_ > 1.5 ) {
        cout <<"\n\n*** SAMG_Solver::CheckSparsityCriterion: Matrix sparsity criterion not fulfilled, levelx parameter set to default. ***\n\n";
        cout <<"\t(current single level setup should only be used for matrices that are essentially diagonal (transport eqn. etc.)).";
        cout <<"\n\n\tlevelx is set to: ";
        levelx = 25;
        cout << levelx <<"\n";
    }
}





/// SAMG iteration cycles benchmark criterion
/**
@section arguments Input Arguments

Below functions are called only if NO_PRIMARY_SOLVER_CONTROL is defined by the user.

When SAMG solver settings allow reuse of coarsening setups for subsequent time steps it must be
ensured that it is still applicable to the problem. The below function has access to the SAMG output
parameter ncyc_done_ and stores the lowest number of performed iteration cycles required to reach
convergence (ncyc_best_). If the criterion can not be fulfilled a new SAMG solver setting is enforce by
adjusting the SAMG paramters iswtch before it is passed on to the solver call. Subsequent solver calls
will query the unaffected iswit setting of the SAMG_Settings instance and hence resume reusing the new
solver setup.

Within CSMP++ this optimization feature can be only applied when solving for diffusion in an IMPES
numeric scheme or IMP-IMS in combination with SAMG Multiple Instances.
*/
void SAMG_Solver::CheckCycleCriterion( int32_t& iswtch ) const {
    if ( ncyc_done_ > 1.5 * ncyc_best_ ) {
        iswtch = iswtch + 1000000; // equals settings_->Set_iswit(4);
        cout <<"\n\n*** SAMG_Solver::CheckCycleCriterion: Number of iteration cycles criterion failed, new coarsening setup forced! ***\n\n";
        cout.flush();
    }
}





void SAMG_Solver::UpdateCycleCriterion( int32_t iswtch ) {
    if ( iswtch < 4000000 ) {
        ncyc_best_ = std::min( ncyc_best_, ncyc_done_ );
    }
    else {
        ncyc_best_ = ncyc_done_;
    }
}




/**
     Verbose checks for absolute or relative convergence as supplied by user.
     
     Returns true if the solution residual is smaller than the user
     supplied eps (this usually indicates convergence)
*/
bool SAMG_Solver::CheckConvergence( double eps ) const
{
    // if the (absolute) convergence criterion eps = 0.
    if ( fabs(eps) < 1.0e-30 ) return true;

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    cout << "\n*** SAMG_Solver::CheckConvergence: ";
    cout <<"eps = " << eps << ", res_out = " << res_out_ << ", versus res_in = " << res_in_ << "' ***\n\n";
    cout.flush();

    // Absolute convergence: negative values of eps indicate SAMG to use eps as an absolute convergence criterion
    if ( eps < 0. )
    {
        if ( res_out_ > fabs( eps ) ) {
            csmp_error.notice( WARNING,
                               "SAMG_Solver::CheckConvergence",
                               "Absolute error L2-norm solution criterion was not fulfilled");
            return false;
        }
    }

    // Relative convergence: positive values of eps prompt SAMG to use eps as a relative convergence criterion
    if ( eps > 0. )
    {
        const double  relative_residual = ( res_in_ == 0. ) ? 0. : res_out_/res_in_;
        
        cout <<"\n\trelative residual = " << relative_residual << "\n\n";
        // make sure that the target residual is strongly violated (best single prec. solution)
        if ( relative_residual > eps ) {
            csmp_error.notice( WARNING,
                               "SAMG_Solver::CheckConvergence",
                               "Relative solution criterion was not fulfilled" );
            return false;
        }
        cout.flush();
    }

    return true;

} // end CheckConvergence


} // end namespace csmp



#ifndef CSP_PC_WINDOWS_NT_CODE_WARRIOR
void samg_user_coo( int*, int* ndim, double*, double*, double* );
void samg_user_coo_( int*, int* ndim, double*, double*, double* );

void samg_user_coo( int*, int* ndim, double*, double*, double* ) {
    cout << "\nSAMG_Adaptor::SAMG_USER_COO: Setting *ndim to 0 " << endl;
    *ndim = 0;
}
void samg_user_coo_(int*,int* ndim, double*, double*, double* ) {
    cout << "\nSAMG_Adaptor::SAMG_USER_COO: Setting *ndim to 0 " << endl;
    *ndim = 0;
}

#endif





