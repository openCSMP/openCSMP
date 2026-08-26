#include "SAMG_Settings.h"
#include "Exception.h"

#if defined(_OPENMP)
#include "omp.h"
#endif

using std::domain_error;

namespace csmp {

SAMG_Settings::SAMG_Settings() :
    // matrix subswitches (2)
    isym_(2),
    irow_(2),
    // ifirst subswitches (1)
    itypu_(1),
    // eps switch (2)
    eps_(0.),
    rel_eps_(1.E-10),
    // nsolve subswitches (7)
    napproach_(2),
    nxtyp_(0),
    nrd_(131),
    nru_(131),
    internal_(0),
    nprim_(0),
    npr_is_dummy_(0),
    nint_weights_(0),
    nint_pat_(0),
    // ncyc subswitches (4)
    igam_(1),
    ncgrad_(1),
    nkdim_(0),
    ncycle_(30),
    // iswtch subswitches (5)
    iswit_(5),
    iextent_(1),
    ndefault_(10),
    norm_typ_(0),
#if defined(_OPENMP)
    iordered_omp_(1),
    samg_omp_num_threads_external_(0),
    irestriction_openmp_(2),
#else
    iordered_omp_(0),
#endif
    ioscratch_(0),
    // chktol switch (1)
    chktol_(-1.0),
    // idump subswitches (4)
    idmp_(1),
    igdp_(0),
    iadp_(0),
    iwdp_(0),
    // iout subswitches (2)
    iout1_(1),
    iout2_(0),
    // a_cmplx switch (1)
    a_cmplx_(1.0),
    // g_cmplx switch (1)
    g_cmplx_(1.0),
    // p_cmplx switch (1)
    p_cmplx_(0.0),
    // w_avrge switch (1)
    w_avrge_(1.0),
    // ncg subswitches (5)
    ncgtyp_(1),
    nred_(0),
    nredlev_(0),
    nxf_clean_(0),
    npcol_(0),
    // other hidden parameters
    levelx_(25),
    clsolver_finest_(0),
    nptmax_(200),
    ioform_(102),
    ioform_length_(1),
    filnam_dump_("level"),
    filnam_dump_length_(5),
    iter_pre_(0),
    // negative sign for switches
    negative_nsolve_(false),
    negative_ncyc_(false),
    negative_iout_(false),
    negative_idump_(false),
    explicit_secondary_(false),
    // output variables
    ncyc_done_(0),
    ncyc_best_(0),
    mode_mess_(0)
#ifdef SAMG_MULTIPLE_INSTANCES
    , solver_instance_(0)
#endif
{
    filnam_dump_Array_[0] = 108;  // 'l'
    filnam_dump_Array_[1] = 101;  // 'e'
    filnam_dump_Array_[2] = 118;  // 'v'
    filnam_dump_Array_[3] = 101;  // 'e'
    filnam_dump_Array_[4] = 108;  // 'l'
}

// ============================================================================
// GETTERS
// ============================================================================

int32_t SAMG_Settings::Get_matrix() const { return 10*isym_ + irow_; }
int32_t SAMG_Settings::Get_ifirst() const { return itypu_; }
double  SAMG_Settings::Get_eps()    const { return eps_;   }
double  SAMG_Settings::Get_rel_eps() const { return -rel_eps_; }

int32_t SAMG_Settings::Get_nsolve() const
{
    int32_t temp = 10000000 * napproach_ + 1000000 * nxtyp_ + 100000 * internal_
                 + 1000 * nprim_ + 100 * npr_is_dummy_ + 10 * nint_weights_ + nint_pat_;
    return ( negative_nsolve_ ? -temp : temp );
}

int32_t SAMG_Settings::Get_ncyc() const
{
    std::stringstream stream;
    if ( negative_ncyc_ ) stream << "-";
    stream << igam_ << ncgrad_ << nkdim_ << ncycle_;
    return static_cast<int32_t>( std::atoi( stream.str().c_str() ) );
}

int32_t SAMG_Settings::Get_iswit()  const { return iswit_; }

int32_t SAMG_Settings::Get_iswtch() const
{
    return 1000000 * iswit_ + 100000 * iextent_ + 1000 * ndefault_ + 100 * norm_typ_ + ioscratch_;
}

double  SAMG_Settings::Get_chktol() const { return chktol_; }

int32_t SAMG_Settings::Get_idump() const
{
    int32_t temp = 1000 * idmp_ + 100 * igdp_ + 10 * iadp_ + iwdp_;
    return ( negative_idump_ ? -temp : temp );
}

int32_t SAMG_Settings::Get_idmp() const { return idmp_; }

int32_t SAMG_Settings::Get_iout() const
{
    int32_t temp = 10 * iout1_ + iout2_;
    return ( negative_iout_ ? -temp : temp );
}

double  SAMG_Settings::Get_a_cmplx() const { return a_cmplx_; }
double  SAMG_Settings::Get_g_cmplx() const { return g_cmplx_; }
double  SAMG_Settings::Get_p_cmplx() const { return p_cmplx_; }
double  SAMG_Settings::Get_w_avrge() const { return w_avrge_; }
int32_t SAMG_Settings::Get_levelx()  const { return levelx_;  }

int32_t SAMG_Settings::Get_clsolver_finest() const { return clsolver_finest_; }
int32_t SAMG_Settings::Get_nptmax()          const { return nptmax_;          }

int32_t SAMG_Settings::Get_ioform()        const { return ioform_;        }
int32_t SAMG_Settings::Get_ioform_length() const { return ioform_length_; }

/**
    Returns the filnam_dump array, rebuilding it from the string each time
    to ensure it stays in sync with any changes made via Set_filnam_dump().
*/
int* SAMG_Settings::Get_filnam_dump()
{
    const char* ch = filnam_dump_.c_str();
    for ( int32_t i = 0; i < filnam_dump_length_; ++i )
        filnam_dump_Array_[i] = int( ch[i] );
    return this->filnam_dump_Array_;
}

std::string SAMG_Settings::Get_filnam_dump_str() const { return filnam_dump_; }

int32_t SAMG_Settings::Get_filnam_dump_length() const
{
    if ( filnam_dump_length_ > 50 ) {
        std::cout << "SAMG_Settings::Get_filnam_dump_length: "
                  << "filename may not exceed 50 characters; received '"
                  << filnam_dump_ << "'\n";
    }
    return filnam_dump_length_;
}

int32_t SAMG_Settings::Get_ncg() const
{
    return 10000 * ncgtyp_ + 1000 * nred_ + 100 * nredlev_ + 10 * nxf_clean_ + npcol_;
}

int32_t SAMG_Settings::Get_ncyc_done() const { return ncyc_done_; }
int32_t SAMG_Settings::Get_ncyc_best() const { return ncyc_best_; }
int32_t SAMG_Settings::Get_iter_pre()  const { return iter_pre_;  }

bool SAMG_Settings::ExplicitSecondary()     const { return explicit_secondary_; }
bool SAMG_Settings::UsePointBasedApproach() const { return ( napproach_ > 2 );  }

#ifdef SAMG_MULTIPLE_INSTANCES
int32_t SAMG_Settings::GetSolverInstance() const { return solver_instance_; }
#endif

int32_t SAMG_Settings::Get_mode_mess() const { return mode_mess_; }
int32_t SAMG_Settings::Get_nxtyp()     const { return nxtyp_;     }
int32_t SAMG_Settings::Get_nrd()       const { return nrd_;       }
int32_t SAMG_Settings::Get_nru()       const { return nru_;       }

#ifdef _OPENMP
int32_t SAMG_Settings::Get_icolor_omp()                    const { return icolor_omp_;                    }
int32_t SAMG_Settings::Get_iordered_omp()                  const { return iordered_omp_;                  }
int32_t SAMG_Settings::Get_irestriction_openmp()           const { return irestriction_openmp_;           }
int32_t SAMG_Settings::Get_samg_omp_num_threads_external() const { return samg_omp_num_threads_external_; }
#endif

// ============================================================================
// SETTERS
// ============================================================================

void SAMG_Settings::Set_isym( int32_t isym )
{
    if ( isym != 1 && isym != 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_isym",
                               "isym must have a value of 1 or 2" );
    isym_ = isym;
}

void SAMG_Settings::Set_irow( int32_t irow )
{
    if ( irow != 1 && irow != 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_irow",
                               "irow must have a value of 1 or 2" );
    irow_ = irow;
}

void SAMG_Settings::Set_itypu( int32_t itypu )
{
    if ( itypu != 0 && itypu != 1 && itypu != 2 && itypu != 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_itypu",
                               "itypu must have a value of 0, 1, 2 or 3" );
    itypu_ = itypu;
}

void SAMG_Settings::Set_eps( double eps )     { eps_ = -eps; }

void SAMG_Settings::Set_rel_eps( double rel_eps )
{
    if ( rel_eps <= 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_rel_eps",
                               "rel_eps must have a positive value" );
    rel_eps_ = rel_eps;
}

void SAMG_Settings::Set_napproach( int32_t napproach )
{
    if ( napproach < 1 || napproach > 5 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_napproach",
                               "napproach must have an integer value between 1 and 5" );
    napproach_ = napproach;
}

void SAMG_Settings::Set_nxtyp( int32_t nxtyp )
{
    if ( nxtyp != 0 && nxtyp != 1 && nxtyp != 2 && nxtyp != 3 && nxtyp != 5 && nxtyp != 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nxtyp",
                               "nxtyp must have a value of 0, 1, 2, 3, 5 or 9" );
    nxtyp_ = nxtyp;
}

void SAMG_Settings::Set_nrd( int32_t nrd )
{
    if ( nrd > 999 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nrd",
                               "nrd must have one of three possible values. Read the SAMG manual" );
    nrd_ = nrd;
}

void SAMG_Settings::Set_nru( int32_t nru )
{
    if ( nru > 999 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nru",
                               "nru must have one of three possible values. Read the SAMG manual" );
    nru_ = nru;
}

void SAMG_Settings::Set_internal( int32_t internal )
{
    if ( internal < 0 || internal > 4 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_internal",
                               "internal must have an integer value between 0 and 4" );
    internal_ = internal;
}

void SAMG_Settings::Set_nprim( int32_t nprim )
{
    if ( nprim < 0 || nprim > 99 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nprim",
                               "nprim must have an integer value between 0 and 99" );
    nprim_ = nprim;
}

void SAMG_Settings::Set_npr_is_dummy( int32_t npr_is_dummy )
{
    if ( npr_is_dummy != 0 && npr_is_dummy != 1 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_npr_is_dummy",
                               "npr_is_dummy must have a value of 0 or 1" );
    npr_is_dummy_ = npr_is_dummy;
}

void SAMG_Settings::Set_nint_weights( int32_t nint_weights )
{
    if ( nint_weights < 0 || nint_weights > 4 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nint_weights",
                               "nint_weights must have an integer value between 0 and 4" );
    nint_weights_ = nint_weights;
}

void SAMG_Settings::Set_nint_pat( int32_t nint_pat )
{
    if ( nint_pat < 0 || nint_pat > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nint_pat",
                               "nint_pat must have an integer value between 0 and 2" );
    nint_pat_ = nint_pat;
}

void SAMG_Settings::Set_igam( int32_t igam )
{
    if ( igam < 1 || igam > 6 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_igam",
                               "igam must have an integer value between 1 and 6" );
    igam_ = igam;
}

void SAMG_Settings::Set_ncgrad( int32_t ncgrad )
{
    if ( ncgrad < 0 || ncgrad > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ncgrad",
                               "ncgrad must have an integer value between 0 and 3" );
    ncgrad_ = ncgrad;
}

void SAMG_Settings::Set_nkdim( int32_t nkdim )
{
    if ( nkdim < 0 || nkdim > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nkdim",
                               "nkdim must have an integer value between 0 and 9" );
    nkdim_ = nkdim;
}

void SAMG_Settings::Set_ncycle( int32_t ncycle )
{
    if ( ncycle < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ncycle",
                               "ncycle must have a positive integer value" );
    ncycle_ = ncycle;
}

void SAMG_Settings::Set_ncyc_done( int32_t ncyc_done ) { ncyc_done_ = ncyc_done; }
void SAMG_Settings::Set_ncyc_best( int32_t ncyc_best ) { ncyc_best_ = ncyc_best; }

void SAMG_Settings::Set_iswit( int32_t iswit )
{
    if ( iswit < 1 || iswit > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iswit",
                               "iswit must have an integer value between 1 and 9" );
    iswit_ = iswit;
}

void SAMG_Settings::Set_iextent( int32_t iextent )
{
    if ( iextent < 0 || iextent > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iextent",
                               "iextent must have an integer value between 0 and 3" );
    iextent_ = iextent;
}

void SAMG_Settings::Set_ndefault( int32_t ndefault )
{
    // 0 is valid: means no defaults are set (ETHZ version)
    if ( !(ndefault == 0)                          &&
         !( ndefault >= 10 && ndefault <= 13 )     &&
         !( ndefault >= 15 && ndefault <= 18 )     &&
         !( ndefault >= 20 && ndefault <= 23 )     &&
         !( ndefault >= 25 && ndefault <= 28 )     &&
         !( ndefault >= 30 && ndefault <= 38 )     &&
         !( ndefault >= 40 && ndefault <= 48 ) )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ndefault",
                               "ndefault must be 0 or in one of the ranges "
                               "[10,13], [15,18], [20,23], [25,28], [30,38], [40,48]" );
    ndefault_ = ndefault;
}

void SAMG_Settings::Set_norm_typ( int32_t norm_typ )
{
    if ( norm_typ < 0 || norm_typ > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_norm_typ",
                               "norm_typ must have an integer value between 0 and 2" );
    norm_typ_ = norm_typ;
}

void SAMG_Settings::Set_ioscratch( int32_t ioscratch )
{
    if ( ioscratch < 0 || ioscratch > 99 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ioscratch",
                               "ioscratch must have an integer value between 0 and 99" );
    ioscratch_ = ioscratch;
}

void SAMG_Settings::Set_chktol( double chktol ) { chktol_ = chktol; }

void SAMG_Settings::Set_idmp( int32_t idmp )
{
    if ( idmp > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_idmp",
                               "idmp must have an integer value between 0 and 9" );
    idmp_ = idmp;
}

void SAMG_Settings::Set_igdp( int32_t igdp )
{
    if ( igdp > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_igdp",
                               "igdp must have an integer value between 0 and 2" );
    igdp_ = igdp;
}

void SAMG_Settings::Set_iadp( int32_t iadp )
{
    if ( iadp > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iadp",
                               "iadp must have an integer value between 0 and 3" );
    iadp_ = iadp;
}

void SAMG_Settings::Set_iwdp( int32_t iwdp )
{
    if ( iwdp < 0 || iwdp > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iwdp",
                               "iwdp must have an integer value between 0 and 3" );
    iwdp_ = iwdp;
}

void SAMG_Settings::Set_iout1( int32_t iout1 )
{
    if ( iout1 > 4 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iout1",
                               "iout1 must have an integer value between 0 and 4" );
    iout1_ = iout1;
}

void SAMG_Settings::Set_iout2( int32_t iout2 )
{
    if ( iout2 > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iout2",
                               "iout2 must have a value less than or equal to 3" );
    iout2_ = iout2;
}

void SAMG_Settings::Set_a_cmplx( double a_cmplx )
{
    if ( a_cmplx < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_a_cmplx",
                               "a_cmplx must have a positive value" );
    a_cmplx_ = a_cmplx;
}

void SAMG_Settings::Set_g_cmplx( double g_cmplx )
{
    if ( g_cmplx < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_g_cmplx",
                               "g_cmplx must have a positive value" );
    g_cmplx_ = g_cmplx;
}

void SAMG_Settings::Set_p_cmplx( double p_cmplx )
{
    if ( p_cmplx < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_p_cmplx",
                               "p_cmplx must have a positive value" );
    p_cmplx_ = p_cmplx;
}

void SAMG_Settings::Set_w_avrge( double w_avrge )
{
    if ( w_avrge < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_w_avrge",
                               "w_avrge must have a positive value" );
    w_avrge_ = w_avrge;
}

void SAMG_Settings::Set_ncgtyp( int32_t ncgtyp )
{
    if ( ncgtyp < 1 || ncgtyp > 5 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ncgtyp",
                               "ncgtyp must have an integer value between 1 and 5" );
    ncgtyp_ = ncgtyp;
    ExplicitSecondary( true );
}

void SAMG_Settings::Set_nred( int32_t nred )
{
    if ( nred < 0 || nred > 6 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nred",
                               "nred must have an integer value between 0 and 6" );
    nred_ = nred;
    ExplicitSecondary( true );
}

void SAMG_Settings::Set_nredlev( int32_t nredlev )
{
    if ( nredlev < 0 || nredlev > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nredlev",
                               "nredlev must have an integer value between 0 and 9" );
    nredlev_ = nredlev;
    ExplicitSecondary( true );
}

void SAMG_Settings::Set_nxf_clean( int32_t nxf_clean )
{
    if ( nxf_clean < 0 || nxf_clean > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nxf_clean",
                               "nxf_clean must have an integer value between 0 and 2" );
    nxf_clean_ = nxf_clean;
    ExplicitSecondary( true );
}

void SAMG_Settings::Set_npcol( int32_t npcol )
{
    if ( npcol < 0 || npcol > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_npcol",
                               "npcol must have an integer value between 0 and 2" );
    npcol_ = npcol;
    ExplicitSecondary( true );
}

/**
    Deals with the handling of negative values in the diagonal of the solution matrix.

    Integer. Defines the maximum number of non-positive diagonal entries allowed in
    computing the Galerkin operators before SAMG gives up its attempts to modify
    interpolation and continues without further checks.
    If neg_diag < 0, checking of the diagonal is completely de-activated.

    @warning In the latter case you should know what you are doing!
*/
void SAMG_Settings::Set_neg_diag( int neg_diag )
{
    neg_diag_ = neg_diag;
}

/**
    SAMG safety limit: maximum number of levels to be created.
    Special debug option: giving levelx a negative sign causes SAMG to
    display its convergence history also on coarser levels.
    Default = 25.
*/
void SAMG_Settings::Set_levelx( int32_t levelx )
{
    if ( levelx < -25 || levelx > 25 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_levelx",
                               "levelx must have an integer value between -25 and 25" );
    levelx_ = levelx;
    ExplicitSecondary( true );
}

/**
    RTFM SAMG manual. This parameter is rarely used, unless very small problems
    are being solved. Default value is 0.
*/
void SAMG_Settings::Set_clsolver_finest( int32_t clsolver_finest )
{
    clsolver_finest_ = clsolver_finest;
    ExplicitSecondary( true );
}

void SAMG_Settings::Set_nptmax( int32_t nptmax )
{
    nptmax_ = nptmax;
}

/**
    Sets the matrix output format for SAMG file output.
    Reading or writing can be formatted, unformatted or binary.
    With the char*1 variable ch being 'f', 'u', or 'b', issue
    "call samg_set_ioform(ch,1)" to specify the requested format,
    the default being formatted.
*/
void SAMG_Settings::Set_ioform( std::string ioform )
{
    if ( ioform != "b" && ioform != "f" && ioform != "u" )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ioform",
                               "ioform must be 'f' (formatted), 'u' (unformatted) or 'b' (binary)" );

    std::cout << "\n*** SAMG_Settings::Set_ioform: format set to '" << ioform << "' ***\n";

    if ( idmp_ < 2 || idmp_ > 10 )
        std::cout << "*** SAMG_Settings::Set_ioform: to take effect, set idmp between 2 and 10 ***\n";

    const char* ch = ioform.c_str();
    ioform_        = int( ch[0] );
    ioform_length_ = 1;
    ExplicitSecondary( true );
}

/**
    Sets the file root name which SAMG uses for matrix dump output.
    If idmp > 1, the matrices are dumped to files using this root name.
    The default name is 'level'.
*/
void SAMG_Settings::Set_filnam_dump( const std::string& filnam_dump )
{
    if ( idmp_ < 2 || idmp_ > 10 )
        std::cout << "*** SAMG_Settings::Set_filnam_dump: to take effect, set idmp between 2 and 10 ***\n";

    filnam_dump_        = filnam_dump;
    filnam_dump_length_ = static_cast<int32_t>( filnam_dump_.length() );

    const char* ch = filnam_dump_.c_str();
    for ( int32_t i = 0; i < filnam_dump_length_; ++i )
        filnam_dump_Array_[i] = int( ch[i] );

    ExplicitSecondary( true );
}

/**
    Sets the hidden parameter iter_pre. If iter_pre > 0, SAMG will be used
    as a pre-conditioner for the BI-CGStab solution.
    A maximum of 2 SAMG pre-condition calls is possible, i.e. 0 <= iter_pre <= 2.
*/
void SAMG_Settings::Set_iter_pre( int32_t iter_pre )
{
    if ( iter_pre < 0 || iter_pre > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iter_pre",
                               "iter_pre must have an integer value between 0 and 2" );
    iter_pre_ = iter_pre;
    ExplicitSecondary( true );
}

void SAMG_Settings::SetNegative_nsolve( bool negative_nsolve ) { negative_nsolve_ = negative_nsolve; }
void SAMG_Settings::SetNegative_ncyc(   bool negative_ncyc   ) { negative_ncyc_   = negative_ncyc;   }
void SAMG_Settings::SetNegative_idump(  bool negative_idump  ) { negative_idump_  = negative_idump;  }
void SAMG_Settings::SetNegative_iout(   bool negative_iout   ) { negative_iout_   = negative_iout;   }

void SAMG_Settings::ExplicitSecondary( bool explicit_secondary )
{
    explicit_secondary_ = explicit_secondary;
}

// ============================================================================
// SAMG_MULTIPLE_INSTANCES
// ============================================================================

#ifdef SAMG_MULTIPLE_INSTANCES
/**
    Sets the SAMG solver instance index (0..5).
    Only available when compiled with SAMG_MULTIPLE_INSTANCES.
    
    @param instance  solver instance index in range [0, 5]
*/
void SAMG_Settings::SetSolverInstance( int32_t instance )
{
    if ( instance < 0 || instance > 5 )
        throw csmp::Exception( ERROR, "SAMG_Settings::SetSolverInstance",
                               "solver instance must be in range [0, 5]" );
    solver_instance_ = instance;
}

#endif  // SAMG_MULTIPLE_INSTANCES

// ============================================================================
// MODE_MESS
// ============================================================================

/**
    Controls SAMG screen output verbosity.

    By setting the primary parameters iout and idump to something negative,
    all output regarding SAMG's performance (setup, cycling history, etc.)
    will be suppressed. However, general messages, warnings and error messages
    will still be printed unless mode_mess is adjusted:

      mode =  1   All kinds of messages are printed.
              0   Control messages are skipped.
             -1   In addition: general messages are skipped.
             -2   In addition: warning messages are skipped.
             -3   In addition: error messages are skipped.
*/
void SAMG_Settings::Set_mode_mess( int32_t mode )
{
    mode_mess_ = mode;
}

// ============================================================================
// OPENMP
// ============================================================================

#ifdef _OPENMP
void SAMG_Settings::Set_icolor_omp( int32_t icolor_omp )
{
    if ( icolor_omp < 0 || icolor_omp > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_icolor_omp",
                               "icolor_omp must have an integer value between 0 and 2" );
    icolor_omp_ = icolor_omp;
}

void SAMG_Settings::Set_iordered_omp( int32_t iordered_omp )
{
    if ( iordered_omp < 0 || iordered_omp > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iordered_omp",
                               "iordered_omp must have an integer value between 0 and 2" );
    iordered_omp_ = iordered_omp;
}

/**
    Switch to choose between different variants to perform the fine-to-coarse
    residual restriction in parallel.

      2   Dynamic choice between 0 and -1: for single-thread runs selects 0,
          otherwise -1. (SAMG internal default)
      1   Actual computation of residual is done in parallel, but restriction
          (via transposed interpolation) is sequential. Reasonable parallel
          efficiency, no additional memory required.
      0   Sequential computation, no OpenMP.
     -1   Transposed interpolation matrix is computed at the beginning of the
          solution phase and kept in memory. Highest parallel efficiency at
          the expense of extra memory.
*/
void SAMG_Settings::Set_irestriction_openmp( int32_t irestriction_openmp )
{
    if ( irestriction_openmp < -1 || irestriction_openmp > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_irestriction_openmp",
                               "irestriction_openmp must have an integer value between -1 and 2" );
    irestriction_openmp_ = irestriction_openmp;
}
#endif  // _OPENMP

} // end namespace csmp

