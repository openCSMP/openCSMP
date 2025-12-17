#include "SAMG_Settings.h"
#include "Exception.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using std::domain_error;

namespace csmp {

/**

Default constructor. The standard settings are according to the SAMG user
manual for version 21b1, if available. These standard settings apply well
to scalar problems.
*/
SAMG_Settings::SAMG_Settings() :
    // matrix subswitches (2)
    isym_(2),
    irow_(2),
    // ifirst subswitches (1)
    itypu_(1), // First approximation using a solution vector u=0. (Required for initialization, otherwise u must be initialized)
    // eps switch (2)
    eps_(0.), // solves until roundoff error
    rel_eps_(1.E-10),
    // nsolve subswitches (7)
    napproach_(2),
    nxtyp_(0),
    nrd_(131), //default value from manual (pre-smoothing steps setup)
    nru_(131), //default value from manual (post-smoothing steps setup)
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
#if defined(_OPENMP )
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
    // other hidden parameters (6)
    levelx_(25),
    clsolver_finest_(0),
    nptmax_(200),
    ioform_(102),  // equifalent ASCII definition for "f" returns formatted matrix output
    ioform_length_(1), // lenght of character of ioform by default can be only equal to 1
    filnam_dump_("level"),
    filnam_dump_length_(5),
    iter_pre_(0), // using SAMG as preconditioner for BICGStab
    // negative sign for switches (5)
    negative_nsolve_(false),
    negative_ncyc_(false),
    negative_iout_(false),
    negative_idump_(false),
    explicit_secondary_(false),
    // output varibles (3)
    ncyc_done_(0),
    ncyc_best_(0),
    mode_mess_(0) // mode_mess default from manual
{
    filnam_dump_Array_[0] = 108;
    filnam_dump_Array_[1] = 101;
    filnam_dump_Array_[2] = 118;
    filnam_dump_Array_[3] = 101;
    filnam_dump_Array_[4] = 108;
}

/** Returns value of matrix switch.

Involved subswitches: isym, irow
*/
int32_t SAMG_Settings::Get_matrix() const
{
    return 10*isym_ + irow_;
}

/** Returns value of ifirst switch

Involved subswitches: itypu

If ifirst == 1, (end residual) <= eps,
if ifirst == 0, (end residual) / (start residual) <= eps

rel_eps can be accessed with the function Get_rel_eps
*/
int32_t SAMG_Settings::Get_ifirst() const {
    return itypu_;
}

/** Returns value of eps switch

Involved subswitches: eps
*/
double SAMG_Settings::Get_eps() const {
    return eps_;
}

/** Returns value of eps switch in the case of an explicit first guess vector
(itypu == 0).

Involved subswitches: eps
*/
double SAMG_Settings::Get_rel_eps() const {
    return -rel_eps_;
}

/** Returns value of nsolve switch

Involved subswitches: napproach, nxtyp, internal, nprim, npr_is_dummy,
nint_weights, nint_pat

For negative values of nsolve, call SetNegative_nsolve(true) first
*/
int32_t SAMG_Settings::Get_nsolve() const {
    int32_t temp = 10000000 * napproach_ + 1000000 * nxtyp_ + 100000 * internal_
            + 1000 * nprim_ + 100 * npr_is_dummy_ + 10 * nint_weights_ + nint_pat_;
    return ( negative_nsolve_ ? -temp : temp );
}

/** Returns value of ncyc switch

Involved subswitches: igam, ncgrad, nkdim, ncycle

For negative values of ncyc, call SetNegative_ncyc( true ) first
*/
int32_t SAMG_Settings::Get_ncyc() const {
    std::stringstream stream;
    if ( negative_ncyc_ )
        stream << "-";
    stream << igam_ << ncgrad_ << nkdim_ << ncycle_;
    return static_cast<int32_t>( std::atoi( stream.str().c_str()));
}

/** Returns value of iswit

Involved subswitches: iswit, iextent, ndefault, norm_typ, ioscratch
 */
int32_t SAMG_Settings::Get_iswit() const {
    return iswit_;
}

/** Returns value of iswtch switch

Involved subswitches: iswit, iextent, ndefault, norm_typ, ioscratch
 */
int32_t SAMG_Settings::Get_iswtch() const {
    return 1000000 * iswit_ + 100000 * iextent_ + 1000 * ndefault_ + 100 * norm_typ_ + ioscratch_;
}

/** Returns value of chktol switch

Involved subswitches: chktol
 */
double SAMG_Settings::Get_chktol() const {
    return chktol_;
}

/** Returns value of idump switch

Involved subswitches: idmp, igdp, iadp, iwdp

For negative values of idump, call SetNegative_idump( true ) first
 */
int32_t SAMG_Settings::Get_idump() const {
    int32_t temp = 1000 * idmp_ + 100 * igdp_ + 10 * iadp_ + iwdp_;
    return ( negative_idump_ ? -temp : temp );
}

int32_t SAMG_Settings::Get_idmp() const {
    return idmp_;
}

/** Returns value of iout switch

Involved subswitches: iout1, iout2

For negative values of iout, call SetNegative_iout( true ) first
 */
int32_t SAMG_Settings::Get_iout() const {
    int32_t temp = 10 * iout1_ + iout2_;
    return ( negative_iout_ ? -temp : temp );
}

/** Returns value of a_cmplx switch
 */
double SAMG_Settings::Get_a_cmplx() const {
    return a_cmplx_;
}

/** Returns value of g_cmplx switch
 */
double SAMG_Settings::Get_g_cmplx() const {
    return g_cmplx_;
}

/** Returns value of p_cmplx switch
 */
double SAMG_Settings::Get_p_cmplx() const {
    return p_cmplx_;
}

/** Returns value of w_avrge switch
 */
double SAMG_Settings::Get_w_avrge() const {
    return w_avrge_;
}
/** Returns value of optional levelx switch
 */
int32_t SAMG_Settings::Get_levelx() const {
    return levelx_;
}
/** Returns value of clsolver_finest switch
 */
int32_t SAMG_Settings::Get_clsolver_finest() const {
    return clsolver_finest_;
}
/** Returns value of nptmax (maximum system size supported via direct solvers)
 */
int32_t SAMG_Settings::Get_nptmax() const {
    return nptmax_;
}
/** Returns value of optional ioform switch
 */
int32_t SAMG_Settings::Get_ioform() const{
    return ioform_;
}
/** Returns value of ioform ASCII character length
 */
int32_t SAMG_Settings::Get_ioform_length() const{
    return ioform_length_;
}

/** Returns value of optional filnam_dump switch
 */
int* SAMG_Settings::Get_filnam_dump()
{
    const char *ch = filnam_dump_.c_str();
    for( int32_t i=0; i < filnam_dump_length_; ++i ){
        filnam_dump_Array_[i] = int( ch[i] );
    }
    return this->filnam_dump_Array_;
}

std::string SAMG_Settings::Get_filnam_dump_str() const
{
    return filnam_dump_;
}

/** Returns value of filnam_dump ASCII character length
 */
int32_t SAMG_Settings::Get_filnam_dump_length() const{
    if ( filnam_dump_length_ > 50 ){
        std::cout << "SAMG_Settings::Get_filnam_dump_length 'The length of the output filename may not exceed 50 characters' \n";
        std::cout << "SAMG_Settings::Get_filnam_dump_length 'received '"+filnam_dump_+"'\n";
    }

    return filnam_dump_length_;
}
/** Returns value of optional ncg switch

Involved subswitches: ncgtyp, nred, nredlev, nxf_clean, npcol

Whether the optional switch should be set can be determined with the
function ExplicitSecondary(). SAMG_Solver uses this method to decide
whether to stick to the defaults or to use the settings stored in the
SAMG_Settings object. A modification of a subswitch of ncg triggers the
use of the ncg switch automatically.
 */
int32_t SAMG_Settings::Get_ncg() const {
    return 10000 * ncgtyp_ + 1000 * nred_ + 100 * nredlev_ + 10 * nxf_clean_ + npcol_;
}

#if defined(_OPENMP )
int32_t SAMG_Settings::Get_icolor_omp() const {
    return icolor_omp_;
}

int32_t SAMG_Settings::Get_iordered_omp() const {
    return iordered_omp_;
}

int32_t SAMG_Settings::Get_irestriction_openmp() const {
    return irestriction_openmp_;
}
#endif

/** Query to determine how many iteration cycles have been performed.
 */
int32_t SAMG_Settings::Get_ncyc_done() const {
    return ncyc_done_;
}

/** Query to determine the minimum number of iteration cycles performed.
 */
int32_t SAMG_Settings::Get_ncyc_best() const {
    return ncyc_best_;
}

/** Gets iter_pre which allows using SAMG as pre-conditioner for Bi-CGStab
 */
int32_t SAMG_Settings::Get_iter_pre() const {
    return iter_pre_;
}

/** Query to determine whether user defined ncg values should be used. This
function is called in the SAMG_Solver.
 */
bool SAMG_Settings::ExplicitSecondary() const {
    return explicit_secondary_;
}

/** Query to determine if point-based approach is used. When the function returns
false, an unknown-based approach is selected by the SAMG_Settings object.
 */
bool SAMG_Settings::UsePointBasedApproach() const {
    return ( napproach_ > 2 );
}

/** Query to determine which SAMG multiple instance is used. The returned number indicates the SAMG library instance used.
 */
#if defined( LEGACY_SAMG )
int32_t SAMG_Settings::GetSolverInstance() const{
    return solver_instance_;
}
#endif

/** Query to determine which SAMG multiple instance is used. The returned number indicates the SAMG library instance used.
 */
int32_t SAMG_Settings::Get_mode_mess() const{
    return mode_mess_;
}

/** Query the number of threads that samg uses in openMP context. If zero, it is
 * whatever is given by the system in the environment variable OMP_NUM_THREADS
 * */
#if defined(_OPENMP )
int32_t SAMG_Settings::Get_samg_omp_num_threads_external() const{
    return samg_omp_num_threads_external_;
}
#endif




/** Subswitch of matrix.

@section arguments Input Arguments

[1:1], isym         1	A is symmetric
                    2	A is not symmetric

*/
void SAMG_Settings::Set_isym(int32_t isym ) {
    isym_ = isym;
    if ( isym != 1 && isym != 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_isym",
                               "isym must have a value of 1 or 2" );
}

/** Subswitch of matrix.

@section arguments Input Arguments
[2:2], irow0        1	A is a zero rowsum matrix. (For such matrices,
                    the solution will be normalized)
                    2	A is not a zero rowsum matrix
 */
void SAMG_Settings::Set_irow( int32_t irow ) {
    irow_ = irow;
    if ( irow != 1 && irow != 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_irow",
                               "irow must have a value of 1 or 2" );
}

/** Subswitch of ifirst.

@section arguments Input Arguments

[1:1], itypu    0 Actual content of u is chosen as first approximation
                1	First approximation u==0
                2	First approximation u==1
                3	First approximation is a random function

If itypu == 0, a different convergence criterion should be set based on
the norm of the right-hand side vector:

  eps = ||b|| * rel_eps

rel_eps can be accessed with the function Get_rel_eps

*/
void SAMG_Settings::Set_itypu( int32_t itypu ) {
    itypu_ = itypu;
    if ( itypu != 0 && itypu != 1 && itypu != 2 && itypu != 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_itypu",
                               "itypo must have a value of 1, 2 or 3" );
}

/** Stopping criterion.

@section arguments Input Arguments

= 0.0	Stopping criteria based in eps is deactivated (only round-off error)
> 0.0	Iteration stops if res <= eps.res0 (res0 = starting residual)
< 0.0	Iteration stops if res <= |eps|
 */
void SAMG_Settings::Set_eps( double eps ) {
    eps_ = -eps;
}

/** Stopping criterion when first guess is set ( itypu == 0 ).

@section arguments Input Arguments

Relative convergence is used as stopping criterion "res <= eps.res0" (res0 = starting residual).
*/
void SAMG_Settings::Set_rel_eps( double rel_eps ) {
    rel_eps_ = rel_eps;
    if ( rel_eps <= 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_rel_eps",
                               "rel_eps must have a positive value" );
}

/** Subswitch of nsolve.

@section arguments Input Arguments

 [1:1],napproach  1    Scalar approach (regardless of nsys)
                  2	Unknown-based (if used in scalar system napproach will be reset to 1)
                  3-5	Point32_t based approaches - selects type of interpolation to use
                  3: interp. is separate for each unknown
                  4: interp. is same " " "
                  5: interp. is point- (block-) wise
   */
void SAMG_Settings::Set_napproach( int32_t napproach ) {
    napproach_ = napproach;
    if ( napproach < 1 || napproach > 5 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_napproach",
                               "napproach must have an integer value between 1 and 5" );
}

/** Subswitch of nsolve.

@section arguments Input Arguments
  [2:2],nxtyp       0	Gauss-Seidel relaxation
                    1	ILU(0) (substantial increase in required memory)
                    2	ILUT
                    3	Special box relaxation
                    5	Gauss-Seidel blockwise
*/
void SAMG_Settings::Set_nxtyp( int32_t nxtyp ) {
    nxtyp_ = nxtyp;
    if ( nxtyp != 0 && nxtyp != 1 && nxtyp != 2 && nxtyp != 3 && nxtyp != 5 && nxtyp != 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nxtyp",
                               "nxtyp must have a value of 0, 1, 2 or 5" );
}

int32_t SAMG_Settings::Get_nxtyp() const {
    return this->nxtyp_;
}

void SAMG_Settings::Set_nrd( int32_t nrd ) {
    nrd_ = nrd;
    if ( nrd > 999 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nrd",
                               "nrd must have one of three possible values. Read the samg manual" );
}

int32_t SAMG_Settings::Get_nrd() const {
    return this->nrd_;
}

void SAMG_Settings::Set_nru( int32_t nru ) {
    nru_ = nru;
    if ( nru > 999 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nru",
                               "nru must have one of three possible values. Read the samg manual" );
}

int32_t SAMG_Settings::Get_nru() const {
    return this->nru_;
}


/** Subswitch of nsolve.

@section arguments Input Arguments

[3:3],internal    = 0	primary matrix P is user-defined (nprim has to be !0)
                  > 0	P is defined internally

 */
void SAMG_Settings::Set_internal( int32_t internal ) {
    internal_ = internal;
    if ( internal < 0 || internal > 4 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_internal",
                               "internal must have an integer value between 0 and 4" );
}

/** Subswitch of nsolve.

@section arguments Input Arguments
[4:5],nprim
 */
void SAMG_Settings::Set_nprim( int32_t nprim ) {
    nprim_ = nprim;
    if ( nprim < 0 || nprim > 99 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nprim",
                               "nprim must have an integer value between 0 and 99" );
}

/** Subswitch of nsolve.

@section arguments Input Arguments

[6:6],npr_is_dummy    0	Primary unknown is a physical unknown
                      1	Primary unknown is a dummy

  */
void SAMG_Settings::Set_npr_is_dummy( int32_t npr_is_dummy ) {
    npr_is_dummy_ = npr_is_dummy;
    if ( npr_is_dummy != 0 && npr_is_dummy != 1 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_npr_is_dummy",
                               "npr_is_dummy must have a value of 0 or 1" );
}

/** Subswitch of nsolve.

@section arguments Input Arguments
[7:7], nint_weights
 */
void SAMG_Settings::Set_nint_weights( int32_t nint_weights ) {
    nint_weights_ = nint_weights;
    if ( nint_weights < 0 || nint_weights > 4 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nint_weights",
                               "nint_weights must have an integer value between 0 and 4" );
}

/** Subswitch of nsolve.

@section arguments Input Arguments
[8:8], nint_pat
 */
void SAMG_Settings::Set_nint_pat( int32_t nint_pat ) {
    nint_pat = nint_pat;
    if ( nint_pat < 0 || nint_pat > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nint_pat",
                               "nint_pat must have an integer value between 0 and 2" );
}

/** Subswitch of ncyc.

@section arguments Input Arguments
[1:1],igam          1   V-cycle (standard)
                    2	F-cycle
                    3	W-cycle
                    4	WW-cycle (very expensive)
                    5   Diagonally pre?conditioned Krylov method
                    6   Classical ICCG (? Sec on 7.7.2).  This method is only applicable to
                        symmetric matrices for which just the lower tridiagonal part is
                        supplied (? isym in Section 5.1).

*/
void SAMG_Settings::Set_igam( int32_t igam ) {
    igam_ = igam;
    if ( igam < 1 || igam > 6 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_igam",
                               "igam must have an integer value between 1 and 6" );
}

/** Subswitch of ncyc.

@section arguments Input Arguments
[2:2],ncgrad        0	default accelerator
(accelerator)       1	Preconditioner for CG (standard)
                    2	Precon. for BI-CGSTAB
                    3	Precon. for GMRES

 */
void SAMG_Settings::Set_ncgrad( int32_t ncgrad ) {
    ncgrad_ = ncgrad;
    if ( ncgrad < 0 || ncgrad > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ncgrad",
                               "nc_grad must have an integer value between 0 and 3" );
}

/** Subswitch for ncyc.

@section arguments Input Arguments

[3:3],nkdim       0	Select default dimension
(Krylov space     1-8	Dimension = nkdim+1
dimension)        9	Dimension = 20

*/
void SAMG_Settings::Set_nkdim( int32_t nkdim ) {
    nkdim_ = nkdim;
    if ( nkdim < 0 || nkdim > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nkdim",
                               "nkdim must have an integer value between 0 and 9" );
}

/** Subswitch of ncyc.

@section arguments Input Arguments
[4: ]			Max. number of cycles to be performed
 */
void SAMG_Settings::Set_ncycle( int32_t ncycle ) {
    ncycle_ = ncycle;
    if ( ncycle < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ncycle",
                               "ncycle must have a positive integer value" );
}

/** Storage for cycle iterations performed

@section arguments Output Arguments
[4: ]			Number of cycles performed
 */
void SAMG_Settings::Set_ncyc_done( int32_t ncyc_done ) {
    ncyc_done_ = ncyc_done;
}

/** Storage for best cycle iterations performed

@section arguments Output Arguments
[4: ]			Benchmark variable for best number of cycles performed
 */
void SAMG_Settings::Set_ncyc_best( int32_t ncyc_best ) {
    ncyc_best_ = ncyc_best;
}

/** Subswitch of iswitch.

@section arguments Input Arguments

[1:1], iswit	Controls re-use of SAMG decompositions during repeated calls.

          5   Complete SAMG run. Upon return, memory is released
          4   Same as 5 except memory not released
          3   partial setup: Re-use coarser grids and interpolation but update
              Galerkin operators. Memory not released
          2   No setup: Re-use coarser grids, interp. and Galerkin from prev. run
          1   Same as 2 except SAMG assumes matrix A to be the same as prev. run

With settings iswit(4 ,..., 1), it is the user’s responsibility to employ some kind of “outer”
convergence control. Such a control is done automatically if the primary solver control is activated.
This is done by just using the settings iswit( 6,..., 9) instead of iswit(4 ,..., 1). More precisely,
in solving a series of linear systems, the primary control process starts automatically whenever some
iswit(>5) is used, and it will terminate as soon as some iswit(≤ 5) is used.

 */
void SAMG_Settings::Set_iswit( int32_t iswit ) {
    iswit_ = iswit;
    if ( iswit < 1 || iswit > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iswit",
                               "iswit must have an integer value between 1 and 9" );
}

/** Subswitch of iswitch.

@section arguments Input Arguments
 [2:2],iextent	Memory extension switch. Selects beahaviour when limits of initial
               dimensioning have been reached

          0	SAMG returns with error code
          1	SAMG allocates extra memory and continues (if no core space,
               writes prev. allocated data to disk
          2	SAMG allocates ext. memory and continues (if no core space,
               SAMG terminates)
          3	SAMG allocates ext. memory and continues (prev. allocated data
               is written to disk)

 */
void SAMG_Settings::Set_iextent( int32_t iextent ) {
    iextent_ = iextent;
    if ( iextent < 0 || iextent > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iextent",
                               "iextent must have an integer value between 0 and 3" );
}

/**Subswitch of iswitch.

@section arguments Input Arguments
  [3:4],n_default	Selects default secondary parameters (if n_default = 0, no defaults
     are set and must be expicitly defined)
          10-13	Used if no 'critical'positive off-diagonal entries
          15-18	Same as prev but more effort in construction of interpolation
          20-23	alternative to 10-13 for 'critical' positive off-diagonal entries
          25-28	Same as prev with more effort in construction of interpolation

tested: */
void SAMG_Settings::Set_ndefault( int32_t ndefault ) {
    ndefault_ = ndefault;
    if ( !( ndefault >= 10 && ndefault <= 13 ) &&
         !( ndefault >= 15 && ndefault <= 18 ) &&
         !( ndefault >= 20 && ndefault <= 23 ) &&
         !( ndefault >= 25 && ndefault <= 28 ) &&
         !( ndefault >= 30 && ndefault <= 38 ) &&
         !( ndefault >= 40 && ndefault <= 48 ) )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ndefault",
                               "ndefault must be an integer value in one of the ranges [10,13], [15,18], [20,23], or [25,28]" );
}

/** Subswitch of iswitch.

@section arguments Input Arguments
[5:5],norm_typ	Selects type of norm to be used in computing residuals

          0	L2-norm
          1	L1-norm
          2	Maximum norm

 */
void SAMG_Settings::Set_norm_typ( int32_t norm_typ ) {
    norm_typ_ = norm_typ;
    if ( norm_typ < 0 || norm_typ > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_norm_typ",
                               "norm_typ must have an integer value between 0 and 2" );
}

#if defined(_OPENMP )
void SAMG_Settings::Set_icolor_omp( int32_t icolor_omp ) {
    icolor_omp_ = icolor_omp;
    if ( icolor_omp < 0 || icolor_omp > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_icolor_omp",
                               "icolor_omp must have an integer value between 0 and 2" );
}

void SAMG_Settings::Set_iordered_omp( int32_t iordered_omp ) {
    iordered_omp_ = iordered_omp;
    if ( iordered_omp_ < 0 || iordered_omp_ > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iordered_omp",
                               "iordered_omp must have an integer value between 0 and 2" );
}
/**
Integer. Switch to choose between different variants to perform the
irestriction_openmp  fine?to?coarse residual restriction in parallel.
(samg internal default = 2)
2  ->  Dynamic choice between 0 and ?1: For single?thread runs it selects ì0î, otherwise ì?1î.

1  ->  While the actual computation of the residual is done in parallel, its restriction
       (via the transposed interpolation) is sequential: Parallel efficiency is reasonable,
       no additional memory required.
0  ->  Sequential computation, no OpenMP

?1 ->  The transposed interpolation matrix (for restricting the residual) is computed at the beginning of the solution phase and then kept in memory. This
results in highest parallel efficiency at the expense of quite some extra memory.
*/

void SAMG_Settings::Set_irestriction_openmp( int32_t irestriction_openmp ) {
    irestriction_openmp_ = irestriction_openmp;
    if ( irestriction_openmp_ < -1 || irestriction_openmp_ > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_irestriction_openmp",
                               "irestriction_openmp must have an integer value between -1 and 2" );
}

#endif

/** Subswitch of iswitch.

@section arguments Input Arguments

[6:7],ioscratch	Unit number for scratch files used by SAMG for memory management

 */
void SAMG_Settings::Set_ioscratch( int32_t ioscratch ) {
    ioscratch_ = ioscratch;
    if ( ioscratch < 0 || ioscratch > 99 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ioscratch",
                               "ioscratch must have an integer value between 0 and 99" );
}

/** Checking of input matrix

@section arguments Input Arguments

Used to control the amount of checking of the input matrix.

  <0.0d0 No checking. This is the standard for production runs.
  =0.0d0 Standard checking for logical correctness.
  >0.0d0 Enhanced checking. The concrete value of chktol serves as a tolerance.
  Standard value: 1.0d-7.
 */
void SAMG_Settings::Set_chktol( double chktol ) {
    chktol_ = chktol;
    if ( /* DISABLES CODE */ (false) )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_chktol",
                               "chktol has an invaid value" );
}

/** Subswitch of idump.

@section arguments Input Arguments
[1:1], idmp   0 Coarsening history
              1 Standard print output (coarsening history).

 */
void SAMG_Settings::Set_idmp( int32_t idmp ) {
    idmp_ = idmp;
    if ( idmp > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_idmp",
                               "idmp must have an integer value between 0 and 9" );
}

/** Subswitch of idump.

@section arguments Input Arguments

Selects print output regarding the coarse levels.

[2:2], igdp  >1 is only relevant for coupled systems. Otherwise: ignored.
              0 No particular output.
              1 Display table on grids (full problem).
              2 Same for all submatrices (only if nsys>1).
     */
void SAMG_Settings::Set_igdp( int32_t igdp ) {
    igdp_ = igdp;
    if ( igdp > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_igdp",
                               "igdp must have an integer value between 0 and 2" );
}

/** Subswitch of idump.

@section arguments Input Arguments

Selects print output regarding the coarse-level matrices.

[3:3], iadp   >1 is only relevant for coupled systems.
               0 No particular output.
               1 Display table on coarse-level matrices (full problem).
               2 In addition: same info for all submatrices.
               3 In addition: connectivity info between unknowns.
 */
void SAMG_Settings::Set_iadp( int32_t iadp ) {
    iadp_ = iadp;
    if ( iadp > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iadp",
                               "iadp must have an integer value between 0 and 2" );
}

/** Subswitch of idump.

@section arguments Input Arguments

Selects print output regarding the interpolation matrices.

[4:4], iwdp   >1 is only relevant for coupled systems.
              0 No particular output.
              1 Display table on interpolation matrices (full problem).
              2 In addition: same info for all submatrices.
              3 In addition: connectivity info between unknowns.
 */
void SAMG_Settings::Set_iwdp( int32_t iwdp ) {
    iwdp_ = iwdp;
    if ( iwdp < 0 || iwdp > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iwdp",
                               "iwdp must have an integer value between 0 and 2" );
}

/** Subswitch of iout.

@section arguments Input Arguments

[1:1], iout1, 1 Table of input data and work statistics.
              2 Standard history of cycling process.
              3 Extended history: including all levels, full smoothing steps.
              4 Extended history: including all levels and even partial smoothing steps.

 */
void SAMG_Settings::Set_iout1( int32_t iout1 ) {
    iout1_ = iout1;
    if ( iout1 > 4 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iout1",
                               "iout1 must have an integer value between 0 and 4" );
}

/** Subswitch of iout.

@section arguments Input Arguments

[2:2], iout2  0 No action.
              1 Display most relevant SAMG hidden parameters.
              2 Display all SAMG hidden parameters.
              3 Display all SAMG hidden parameters in a single list.

 */
void SAMG_Settings::Set_iout2( int32_t iout2 ) {
    iout2_ = iout2;
    if ( iout2 > 3 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_iout2",
                               "iout2 must have a value less than or equal to 3" );
}

/** Subswitch of a_cmplx.

@section arguments Input Arguments

a_cmplx - Should be an upper limit for the operator complexity which is defined as the ratio between
the total number of matrix entries (summed over all AMG levels) and the number of entries in the given
matrix (= nna). Depending on the problem and the strategy chosen, this may be a value as low as 1.2,
but it may also be as high as 4.0, say. Typical values are 1.5-3.0.
 */
void SAMG_Settings::Set_a_cmplx( double a_cmplx ) {
    a_cmplx_ = a_cmplx;
    if ( a_cmplx < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_a_cmplx",
                               "a_cmplx must have a positive value" );
}

/** Subswitch of g_cmplx.

@section arguments Input Arguments

g_cmplx - Should be an upper limit for the grid complexity which is defined as the ratio between the
total number of variables (summed over all AMG levels) and the number of variables in the given
problem (= nnu). Depending on the problem and the strategy chosen, this may be a value as low as 1.2.
Usually, an upper limit is 2.0.
  */
void SAMG_Settings::Set_g_cmplx( double g_cmplx ) {
    g_cmplx_ = g_cmplx;
    if ( g_cmplx < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_g_cmplx",
                               "g_cmplx must have a positive value" );
}

/** Subswitch of p_cmplx.

@section arguments Input Arguments

p_cmplx - This parameter is relevant only if any of the point-based approaches is selected. p_cmplx
should then be an upper limit for the point complexity which is defined as the ratio between the total
number of points (summed over all AMG levels) and the number of points in the given problem.
Depending on the problem and the strategy chosen, this may be a value as low as 1.2. Usually, an upper
limit is 2.0.
*/
void SAMG_Settings::Set_p_cmplx( double p_cmplx ) {
    p_cmplx_ = p_cmplx;
    if ( p_cmplx < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_p_cmplx",
                               "p_cmplx must have a positive value" );
}

/** Subswitch of w_avrge.

@section arguments Input Arguments

w_avrge - Should be an upper limit for the average row length of interpolation, that is, the total
number of interpolation weights used by SAMG, summed over all levels, divided by the total number of
variables (summed over all levels). Depending on the problem and the strategy chosen, this may be a
value as low as 1.5, but it may also be as high as 6.0, say. A typical average value is 3.0.
*/
void SAMG_Settings::Set_w_avrge( double w_avrge ) {
    w_avrge_ = w_avrge;
    if ( w_avrge < 0 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_w_avrge",
                               "w_avrge must have a positive value" );
}


/** Subswitch of ncg.

@section arguments Input Arguments

Selects the process of defining strong connectivity.
Standard choices: ncgtyp=1 or ncgtyp=4.

[1:1], ncgtyp   1 Standard process. This is supposed to be used if A has
                  mostly negative off-diagonals. Positive off-diagonal
                  elements (if any) should be small. Variables with only
                  positive couplings will become C-variables.
 
                2 Standard process except that variables which have only
                  positive couplings are treated by absolute value.
 
                3 Standard process except that, for mixed-sign rows, all
                  "large" positive entries (threshold parameter ewt2) are
                  eliminated before a decision on strong connectivity is
                  made. If, for some variable i, this does not lead to a clear
                  decision, i will become a C-variable.
 
                4 Same as 3 except that, if the elimination of positive
                  couplings of variable i fails to give a clear picture, this
                  option temporarily switches to the standard process 1.
 
                5 Same as 4 except that variables which have only
                  positive couplings are treated by absolute value.

 */
void SAMG_Settings::Set_ncgtyp( int32_t ncgtyp ) {
    if ( ncgtyp < 1 || ncgtyp > 5 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ncgtyp",
                               "ncgtyp must have an integer value between 1 and 5" );
    ncgtyp_ = ncgtyp;
    ExplicitSecondary( true );
}



/** Subswitch of ncg.

@section arguments Input Argument

Specifies speed of coarsening, ie, how fast the number of variables
is reduced from one level to the next. Standard choice: nred=0.

[2:2], nred   0 Standard coarsening.
              1-4 Aggressive coarsening. 1 is most, 4 is least
                  aggressive; 2-3 are in between. Recommendation:
                  use 1 for anisotropic and 2 for isotropic problems.
              5 Cluster coarsening & piecewise constant interpolation.
              6 Cluster coarsening & multi-pass interpolation.

 */
void SAMG_Settings::Set_nred( int32_t nred ) {
    if ( nred < 0 || nred > 6 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nred",
                               "nred must have an integer value between 0 and 6" );
    nred_ = nred;
    ExplicitSecondary( true );
}

/** Subswitch of ncg.

@section arguments Input Arguments

Specifies the levels where aggressive or cluster coarsening is to be
used. On the remaining levels, standard coarsening will be applied.
Only relevant if nred>0. Standard choice: nredlev=0.

[3:3], nredlev    0-8 From level 1 up to level nredlev+1.
                  9 On all levels (not recommended).

 */
void SAMG_Settings::Set_nredlev( int32_t nredlev ) {
    if ( nredlev < 0 || nredlev > 9 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nredlev",
                               "nredlev must have an integer value between 0 and 9" );
    nredlev_ = nredlev;
    ExplicitSecondary( true );
}

/** Subswitch of ncg.

@section arguments Input Arguments

Defines how to treat exceptional F-variables (XF) left over at the end
of standard or cluster coarsening. Standard choice: nxf_clean=0.

[4:4], nxf_clean    0 Standard procedure (currently equivalent to 1).
                    1 Ensure that all XF-variables have a strong coupling to
                      (at least) one regular F-variable.
                    2 Re-set all XF-variables to C-variables.

 */
void SAMG_Settings::Set_nxf_clean( int32_t nxf_clean ) {
    if ( nxf_clean < 0 || nxf_clean > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_nxf_clean",
                               "nxf_clean must have an integer value between 0 and 2" );
    nxf_clean_ = nxf_clean;
    ExplicitSecondary( true );
}

/** Subswitch of ncg.

@section arguments Hidden Input Arguments

Enforce particular conditions. This is done a posteriori by adding
extra C-variables. Thus, npcol>0 should not be used in combination
with aggressive or cluster coarsening. Standard choice: npcol=0.

[5:5], npcol    0 No action.
                1 Enforce weak F-to-F diagonal dominance (1.0).
                2 Enforce strong F-to-F diagonal dominance (0.75).

 */
void SAMG_Settings::Set_npcol( int32_t npcol ) {
    if ( npcol < 0 || npcol > 2 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_npcol",
                               "npcol must have an integer value between 0 and 2" );
    npcol_ = npcol;
    ExplicitSecondary( true );
}


/**
    Deals with the handling of negative values in the diagonal of the solution matrix
    (SAMG manual on Galerkin coarse-;eve; matrices.
 
    Integer. Defines the maximum number of non‐positive diagonal entries allowed in computing the Galerkin operators before SAMG gives up its attempts to modify interpolation and continues without further checks.
    If neg_diag<0, checking of the diagonal is completely de‐activated. Warning: In the latter case you should know what you are doing!
*/
void SAMG_Settings::Set_neg_diag( int neg_diag )
 {
    neg_diag_ = neg_diag;
 }





/** Switch levelx

@section arguments Hidden Input Arguments

SAMG Safety limit: maximim number of levels to be created.
Special debug option: Giving levelx a negative sign causes SAMG to
display its convergence history also on coarser levels.
Default = 25

*/
void SAMG_Settings::Set_levelx( int32_t levelx ) {
    if ( levelx < -25 || levelx > 25 )
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_levelx",
                               "levelx must have an integer value between -25 and 25" );
    levelx_ = levelx;
    ExplicitSecondary( true );
}

/** RTFM Samg manual.  This parameter is rarely used, unless very small problems are being solved.
    default value is 0.
*/
void SAMG_Settings::Set_clsolver_finest( int32_t clsolver_finest ) {

    clsolver_finest_ = clsolver_finest;
    ExplicitSecondary( true );
}

void SAMG_Settings::Set_nptmax( int32_t nptmax ) {

    nptmax_=nptmax;
}

/** Switch ioform

@section arguments Hidden Input Arguments

Reading  or  writing  can  be  formatted,  unformatted  or  binary.  In  case  of  reading,  SAMG  knows  the
format  by  looking  into  a  special  “format”  file  (suffix  “.frm”, → Section  7.8.2).  In  case  of  writing
(→Sections  7.8.3‐7.8.5),  the  hidden  char*1  variable  ioform  defines  the  format  SAMG  is  supposed  to
use. That is, with the char*1 variable ch being ‘f’, ‘u’, or ‘b’, issue "call samg_set_ioform(ch,1)" to specify
the requested format, the default being formatted.

*/
void SAMG_Settings::Set_ioform( std::string ioform ){
    if ( ioform == "b" || ioform == "f" || ioform == "u" ){
        std::cout <<"\n\n*** SAMG_Settings::Set_ioform 'ioform received the matrix output format setting " << ioform << " ***\n\n";
    }
    else {
        throw csmp::Exception( ERROR, "SAMG_Settings::Set_ioform",
                               "ioform must have a string character of either ‘f’, ‘u’ or ‘b’" );
    }
    if (idmp_ < 2 || idmp_ > 10 )
        std::cout <<"\n\n*** SAMG_Settings::Set_ioform 'To take effect set idmp parameter between 2 and 10' ***\n\n";

    const char *ch = ioform.c_str();
    ioform_ = int( ch[0] );
    ioform_length_ = 1;
    ExplicitSecondary( true );
}

/** Switch filnam_dump

@section arguments Hidden Input Arguments

Character. File root name which SAMG tries to find for reading. If  idmp>1,  the  matrices  as  described  above  are  dumped  to  files  using  the  root  name  filnam_dump,
the  default  name  being  ‘level’.

*/
void SAMG_Settings::Set_filnam_dump( const std::string& filnam_dump ){
    if ( idmp_ < 2 || idmp_ > 10 )
        std::cout <<"\n\n*** SAMG_Settings::Set_filnam_dump 'To take effect set idmp parameter between 2 and 10' ***\n\n\n";

    filnam_dump_ = filnam_dump;
    filnam_dump_length_ = static_cast<int32_t>(filnam_dump_.length());
    const char *ch = filnam_dump_.c_str();
    for( auto i = 0; i < filnam_dump_length_; ++i ){
        filnam_dump_Array_[i] = int( ch[i] );
    }
    ExplicitSecondary( true );
}

/** Subswich of ncg. Sets the hidden parameter iter_pre. If iter_pre > 0, SAMG will be used
    as a pre-conditioner for the BI-CGStab solution. A maximum of 2 SAMG pre-
    condition calls is possible, i.e. 0 <= iter_pre <= 2;
 */
void SAMG_Settings::Set_iter_pre(int32_t iter_pre) {
    if (iter_pre < 0 || iter_pre > 2)
      throw csmp::Exception(ERROR, "SAMG_Settings::Set_iter_pre:",
                                   "iter_pre must have an integer value between 0 and 2");
    iter_pre_ = iter_pre;
    ExplicitSecondary(true);
}

void SAMG_Settings::SetNegative_nsolve( bool negative_nsolve ) {
    negative_nsolve_ = negative_nsolve;
}

void SAMG_Settings::SetNegative_ncyc( bool negative_ncyc ) {
    negative_ncyc_ = negative_ncyc;
}

void SAMG_Settings::SetNegative_idump( bool negative_idump ) {
    negative_idump_ = negative_idump;
}

void SAMG_Settings::SetNegative_iout( bool negative_iout ) {
    negative_iout_ = negative_iout;
}

void SAMG_Settings::ExplicitSecondary( bool explicit_secondary ) {
    explicit_secondary_ = explicit_secondary;
}

#if defined( LEGACY_SAMG )
void SAMG_Settings::SetSolverInstance( int32_t instance ) {
    solver_instance_ = instance;
    if ( instance < 0 || instance > 5 ) {
        throw csmp::Exception( ERROR, "SAMG_Settings::SetSolverInstance",
                               "SAMG Multiple Instances settings are out of range, available instances are 0, and 1" );
    }
}
#endif


void SAMG_Settings::Set_mode_mess(int32_t mode)
{
    // This is a further attempt to manipulate samg output to the screen

    // By setting the primary parameters iout and idump to
    // something negative, all output regarding SAMG’s performance
    // (setup, cycling history, etc) will be suppressed (Section 5.8).
    // HOWEVER: General messages, warnings and error messages will still
    // be printed unless the hidden parameter mode_mess is adjusted:
    //  mode = 1  All kinds of messages are printed.
    //         0  Control‐messages35 are skipped.
    //        ‐1  In addition: General messages are skipped.
    //        ‐2  In addition: Warning messages are skipped.
    //        ‐3  In addition: Error messages are skipped.
    mode_mess_=mode;
}



} // end namespace csmp
