#ifndef CSMP_ODE_STIFF_SOLVER_H
#define CSMP_ODE_STIFF_SOLVER_H

#include "SumOfProducts.h"
#include "SumOfProductsWithExponents.h"

namespace csmp {

/**

@brief Burlisch-Stoer solver for a stiff systemof ordinary differential equations.
  
@section motivation Motivation

A great variety of physical processes are best described by non-linear
systems of ordinary differential equations (ODEs) which describe the 
change of the system. Typically, these ODEs can only be integrated 
(solved) numerically as is the case for rate-dependent 
chemical reactions. If the rate coefficients differ among the 
equations by many orders of magnitude, the systems are called stiff. 
The main numerical problem imposed by the 'stiffness' is that, for instance
in a time-dependent system, some equilibria described by individual ODE
equations are obtained very fast while other equations evolve very slowly. 
Coupling such different equations together, therefore requires very small time 
steps to limit the change in each step over periods where rapid changes occur.
Large timesteps are however essential in periods of minor change to obtain a
solution in a reasonable amount of time.  

Solution methods which provide
the necessary adaptive timestepping are known as stiff ODE solvers, and
the motivation for the design of the class ODE_StiffSolver is to
make a robust solver of this kind available to CSP.   
 
 
@section Design Intent
 

The class ODE_StiffSolver encapsulates a semi-implicit form of a 
Burlisch-Stoer method as described by Bader and Deuflhard (1983) into
a code module which is also capable of the house keeping and analytical 
differentiation of the actual ordinary differential equations. Textfile
input provides an interface for the definition and description
of the system of equations in human-readable form. The equations are 
linked in a matrix scheme which uses LU decomposition and back substitution
of the results after the inversion of the matrix in each timestep. 
The integration routine 'odeint' controls the timestep size using the 
cross derivatives of all of the individual ODEs. For this purpose, the 
derivatives are assembled into a Jacobian matrix. Additional methods are
provided to change coefficients in individual equations without having
to redine them. This capability is useful, for instance, if activities
or surface area changes in an equilibration calculation.  
 
 
@section Applicability
 

The class ODE_StiffSolver is useful only for relatively small systems
of ODEs (<50 degrees of freedom). Runtime scales with number of ODEs
which is also readily explained on a physical basis. Since the
ODE coefficients will influence the condition number of the solution 
matrix, these must not differ from one-another by more than
the numerical (double64) precision. Otherwise the matrix will become singular. 
The solution speed also limits the applicability of the stiff ODE
integration algorithm in serial application in a finite-element model.
 
 
@section Implementation
 

The implementation of the ODE solver stems from the second edition of
Numerical Recipes in C by Press et al. 1996 (p. 742ff) which also provides
detailed documentation of functions included. The original names of the functions
which have been integrated as private methods of the class ODE_StiffSolver
are derivs(), jacobn(), ludcmp(), lubksb(), simpr() pzextr(), odeint(), and
stifbs(). Additional auxiliary methods manage storage in the peculiar 
numerical recipe way. 

The method implemented herein goes back to Bader, G. and Deuflhard, P. 1983,
Numerische Mathematik 41, 373-398.
 
 
@section Participants
 

The objects SumOfProductsWithExponents and SumOfProducts are used to
represent individual OD equations in the stiff solver object. 
 
 
@section Collaboration
 

Inside the CSP framework the class ODE_StiffSolver is used by the 
ReactionVisitor object. This object uses the solver to find fluid-rock
equilibria in individual finite elements.
 
 
@section Consequences
 

The ODE_StiffSolver form a sequentially applied reaction module in the 
object-oriented implementation of an decoupled incremental transport followed 
by reaction scheme for reactive transport calculations with CSP.
 
 
@section Application Examples
 

The ODE_StiffSolver can be applied using its Solve() method inside of
other CSP objects or it can be used as a stand-alone tool to solve 
non-linear ODEs. The system of equations is defined in a textfile 
starting with an underlined headline followed by a blank line, e.g.:

@code
carbonate - dissolution example
===============================

independent variables
---------------------

dependent variable definitions for which is solved for (aqueous, gaseous, and solid)
------------------------------------------------------------------------------------
OH_1m  CO2_aq  CC  Ca_2p  H_1p  H2O_aq  HCO3_1m  

equilibrium constants (in sequence of listed reactions, same number)
--------------------------------------------------------------------
KF=4.37e-2  KB1=3.4593938e-12  KB2=8.336811e+10  KB3=9727.4722

reactions (odes) A + B -> C (K):  dA/dt = -kf [A][B] + kf [C] / K
-----------------------------------------------------------------
OH_1m    = -KF H_1p OH_1m    +  KB1 H2O_aq
CO2_aq   = -KF CO2_aq H2O_aq +  KB2 HCO3_1m H_1p
CC       = -KF CC H_1p       +  KB3 HCO3_1m Ca_2p
Ca_2p    =  KF CC H_1p       + -KB3 HCO3_1m Ca_2p
H_1p     =  KF CO2_aq H2O_aq + -KB2 HCO3_1m H_1p + -KF CC H_1p + KB3 HCO3_1m Ca_2p + -KF H_1p OH_1m + KB1 H2O_aq
H2O_aq   =  KF H_1p OH_1m    + -KB1 H2O_aq
HCO3_1m  =  KF CO2_aq H2O_aq + -KB2 HCO3_1m H_1p + KF CC H_1p + -KB3 HCO3_1m Ca_2p
@endcode
 
Independent variables are specified next (such as constant coefficients). 
In this example, there are no independent variables. Again
there is an underlined headline followed by the data and a blank line. 

The next block defines the names of the dependent variables in the 
equations (i.e., those variables which shall be calculated as a function
of time). Importantly, no '+' or '-' signs must be used here since these 
symbols are reserved by the equation parser for algebraic operations.
A similar block follows, specifying the numerical values of the
equation coefficients using the assigment operator without blank space.
In the example, forward and backward reaction rates are specified in
terms of dissolution data and the equilibrium constants. This first
order kinetic approach is described by Carl Steefel in the Min. Soc.
America volume on reactive transport processes.  

A following 'reaction' block defines the actual ordinary differential
equations in terms of the notation established above. Importantly, the 
ODEs must be listed in exactly the same order in which the dependent
variables were defined before. There must be exactly one equation for each
dependent variable and the stochiometry of the aqueous species/minerals
determines the exponents and multipliers of the equation coefficients.
This derives directly from the related mass-action expressions.  

Everything which follows after the blank line below the equation block
is ignored. Here the user can insert comments and a discussion of the 
system he defined and its time-dependent behaviour. 

Of course, the described syntax can also be used to define arbitrary ODE 
systems other than sets of reaction equations.
 
 
copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts

*/
class  ODE_StiffSolver {
  public:
    explicit ODE_StiffSolver( const std::vector<SumOfProductsWithExponents>& equations );
    ~ODE_StiffSolver();
    
    // assign initial conditions and read computed values (equations 1...n)
    double64& operator()( size_t i );
    // retrieve the time derivatives from the last step before the solver quit
    double64  operator[]( size_t i ) const;
    
    void       Tolerance( double64 epsilon );
    double64  Tolerance() const;
    size_t  Equations() const;
    void       IntermediateStepsToSave( size_t s );
    size_t  IntermediateStepsToSave() const;
    void       InitialConditions( const std::vector<double64>& init );
    
    void    ReadODEsFrom( const char* file );
    
    void    ReadODEsFrom( const char* file, std::vector<std::pair<int32,std::string> >& coefs, 
                                            std::vector<std::pair<int32,std::string> >& dependent_comps,
                                            std::vector<std::pair<int32,std::string> >& independent_comps );
    
    bool    Solve( double64 time, double64 init_stepsize, 
                   size_t n_savesteps, double64 tolerance=1.0e-10, bool print_step_numbers=false );

    void    Out() const { Out(std::cout); }
    void    Out(std::ostream& os) const;  
    void    OutputResults() const;
    void    OutputResultsTo( const char* textfile ) const;
    void    OutputLastDerivatives() const;

    void    EvaluateIndependentProductTerms( const double64 dependent_vals[], 
                                             std::vector<double64>& indep_vals );

    void    ResetIndependentProductTerms( std::vector<std::pair<int32,std::string> >& ind, 
                                          std::vector<double64>& vals );

    void    AssignIndependentProductTerms( std::vector<std::pair<int32,std::string> >& ind, 
                                           std::vector<double64>& vals );
	                   
    void    Out( const std::vector<int32>& vec ) const;  
    void    Out( const std::vector<double64>& vec ) const;                
    void    Out( const std::vector<std::vector<double64> >& vec ) const;                
    void    Out( int32* vec, int32 n );                
    void    Out( double64* vec, int32 n );                
    void    Out( double64** vec, int32 m, int32 n );
    
  private:
    const double64 SAFETY;
    const double64 GROW;
    const double64 PGROW;
    const double64 SHRNK;
    const double64 PSHRNK;
    const double64 ERRCON;
    const double64 MAXTRY;
    const double64 GAM;
    const double64 TINY;
    const double64 MAXSTP;
    const double64 A21; const double64 A31; const double64 A32;
    const double64 C21; const double64 C31; const double64 C32; const double64 C41; const double64 C42; const double64 C43;
    const double64 B1;  const double64 B2;  const double64 B3;  const double64 B4;
    const double64 E1;  const double64 E2;  const double64 E3;  const double64 E4;
    const double64 C1X; const double64 C2X; const double64 C3X; const double64 C4X;
    const double64 A2X; const double64 A3X;
    const int32     NR_END;
    
    // for stiff Bulirsch-Stoer method
    const size_t  KMAXX;
    const size_t  IMAXX;
    const double64     SAFE1;
    const double64     SAFE2;
    const double64     REDMAX;
    const double64     REDMIN;
    const double64     SCALMX;
    const bool       correct_interpretation;

    double64                           eps;   // required accuracy
    size_t                           dof;   // n-equations = n-unknowns
    std::vector<SumOfProductsWithExponents>  odes;  // ODEs
    
    // used by 'pzextr'
    double64  **pzextr_d, *pzextr_x;
    // used by 'odeint'
    double64 *odeint_yscal, *odeint_y, *odeint_dydx;
    // used by 'stifbs'
    double64 *stifbs_dfdx,**stifbs_dfdy,*stifbs_err,*stifbs_yerr,*stifbs_ysav,*stifbs_yseq;
    // initial conditions
    double64* init_cond;
    // saving intermediate results
    double64                             dxsav;
    size_t                             kmax, kount;
    std::vector<double64>                result_xp; // time_increments (size=number of timesteps=kount)
    std::vector<std::vector<double64> >  result_yp; // results nvar rows x timestep cols

    // not usually used:
    std::vector<SumOfProducts>             eqns;  // if there are independent algebraic equations
    std::vector<std::pair<int32,std::string> >  independent_terms;  // independent algebraic equations labels
    std::vector<double64>                 independent_values; // independent variable values

    // auxialiary methods for solver
    int32     *ivector(size_t nl, size_t nh);
    double64  *dvector(size_t nl, size_t nh);
    int32    **imatrix(size_t nrl, size_t nrh, size_t ncl, size_t nch);
    double64 **dmatrix(size_t nrl, size_t nrh, size_t ncl, size_t nch);
    void     free_ivector(int32 *v, size_t nl);
    void     free_dvector(double64 *v, size_t nl);
    void     free_imatrix(int32 **m, size_t nrl, size_t ncl);
    void     free_dmatrix(double64 **m, size_t nrl, size_t ncl);
    
    void     derivs( double64 y[], double64 dydx[] );
    void     jacobn( double64 y[], double64 dfdx[], double64 **dfdy );
    void     LU_Decomposition( double64 **a, size_t n, int32 *indx, double64& d );
    void     LU_BackSubstitution(double64 **a, size_t n, int32 *indx, double64 b[] );
    
    void     simpr( double64 y[], double64 dydx[], double64 dfdx[], double64 **dfdy,
		            size_t n, double64 xs, double64 htot, int32 nstep, double64 yout[] );

    void     pzextr(int32 iest, double64 xest, double64 yest[], double64 yz[], double64 dy[], size_t nv);


    void     odeint( double64 ystart[], size_t nvar, double64 x1, double64 x2, 
                     double64 eps, double64 h1, double64 hmin, size_t& nok, size_t& nbad );
	             
    void     StiffBulirschStoer( double64 y[], double64 dydx[], size_t nv, double64& xx, double64 htry, double64 eps,
	                               double64 yscal[], double64& hdid, double64& hnext );
	                              
    void     AllocateSolverMemory();
    void     FreeSolverMemory();
    void     ZeroAllStorage();
    void     TokenizeString( char* s, const char* delim, std::list<std::string>& li ) const;
    
    double64   DSQR( double64 a ) const;
 };


/// assign initial conditions and read computed values (equations 1...n)
inline  double64& ODE_StiffSolver::operator()( size_t i ) 
  { 
      assert ( i>0 && i<=dof ); 
      return init_cond[i]; 
  }
  
/// retrieve the time derivatives from the last step before the solver quit
inline  double64  ODE_StiffSolver::operator[]( size_t i ) const 
  { 
      assert ( i>0 && i<=dof ); 
      return odeint_dydx[i]; 
  }


inline  void       ODE_StiffSolver::Tolerance( double64 epsilon ) { eps=epsilon; }
inline  double64  ODE_StiffSolver::Tolerance() const           { return eps; }
inline  size_t  ODE_StiffSolver::Equations() const           { return dof; }
inline  size_t  ODE_StiffSolver::IntermediateStepsToSave() const { return kmax; }

inline  double64  ODE_StiffSolver::DSQR( double64 a ) const { return (a==0.0) ? 0.0 : a*a; }

} // csmp

#endif

