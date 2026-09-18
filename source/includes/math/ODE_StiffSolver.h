// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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
the numerical (double) precision. Otherwise the matrix will become singular. 
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
    double& operator()( size_t i );
    // retrieve the time derivatives from the last step before the solver quit
    double  operator[]( size_t i ) const;
    
    void       Tolerance( double epsilon );
    double  Tolerance() const;
    size_t  Equations() const;
    void       IntermediateStepsToSave( size_t s );
    size_t  IntermediateStepsToSave() const;
    void       InitialConditions( const std::vector<double>& init );
    
    void    ReadODEsFrom( const char* file );
    
    void    ReadODEsFrom( const char* file, std::vector<std::pair<int32_t,std::string> >& coefs, 
                                            std::vector<std::pair<int32_t,std::string> >& dependent_comps,
                                            std::vector<std::pair<int32_t,std::string> >& independent_comps );
    
    bool    Solve( double time, double init_stepsize, 
                   size_t n_savesteps, double tolerance=1.0e-10, bool print_step_numbers=false );

    void    Out() const;  
    void    OutputResults() const;
    void    OutputResultsTo( const char* textfile ) const;
    void    OutputLastDerivatives() const;

    void    EvaluateIndependentProductTerms( const double dependent_vals[], 
                                             std::vector<double>& indep_vals );

    void    ResetIndependentProductTerms( std::vector<std::pair<int32_t,std::string> >& ind, 
                                          std::vector<double>& vals );

    void    AssignIndependentProductTerms( std::vector<std::pair<int32_t,std::string> >& ind, 
                                           std::vector<double>& vals );
	                   
    void    Out( const std::vector<int32_t>& vec ) const;  
    void    Out( const std::vector<double>& vec ) const;                
    void    Out( const std::vector<std::vector<double> >& vec ) const;                
    void    Out( int32_t*  vec, int32_t n );                
    void    Out( double* vec, int32_t n );                
    void    Out( double** vec, int32_t m, int32_t n );
    
  private:
    const double SAFETY;
    const double GROW;
    const double PGROW;
    const double SHRNK;
    const double PSHRNK;
    const double ERRCON;
    const double MAXTRY;
    const double GAM;
    const double TINY;
    const double MAXSTP;
    const double A21; const double A31; const double A32;
    const double C21; const double C31; const double C32; const double C41; const double C42; const double C43;
    const double B1;  const double B2;  const double B3;  const double B4;
    const double E1;  const double E2;  const double E3;  const double E4;
    const double C1X; const double C2X; const double C3X; const double C4X;
    const double A2X; const double A3X;
    const int32_t     NR_END;
    
    // for stiff Bulirsch-Stoer method
    const size_t  KMAXX;
    const size_t  IMAXX;
    const double     SAFE1;
    const double     SAFE2;
    const double     REDMAX;
    const double     REDMIN;
    const double     SCALMX;
    const bool       correct_interpretation;

    double                           eps;   // required accuracy
    size_t                           dof;   // n-equations = n-unknowns
    std::vector<SumOfProductsWithExponents>  odes;  // ODEs
    
    // used by 'pzextr'
    double  **pzextr_d, *pzextr_x;
    // used by 'odeint'
    double *odeint_yscal, *odeint_y, *odeint_dydx;
    // used by 'stifbs'
    double *stifbs_dfdx,**stifbs_dfdy,*stifbs_err,*stifbs_yerr,*stifbs_ysav,*stifbs_yseq;
    // initial conditions
    double* init_cond;
    // saving intermediate results
    double                             dxsav;
    size_t                             kmax, kount;
    std::vector<double>                result_xp; // time_increments (size=number of timesteps=kount)
    std::vector<std::vector<double> >  result_yp; // results nvar rows x timestep cols

    // not usually used:
    std::vector<SumOfProducts>             eqns;  // if there are independent algebraic equations
    std::vector<std::pair<int32_t,std::string> >  independent_terms;  // independent algebraic equations labels
    std::vector<double>                 independent_values; // independent variable values

    // auxialiary methods for solver
    int32_t     *ivector(size_t nl, size_t nh);
    double  *dvector(size_t nl, size_t nh);
    int32_t    **imatrix(size_t nrl, size_t nrh, size_t ncl, size_t nch);
    double **dmatrix(size_t nrl, size_t nrh, size_t ncl, size_t nch);
    void     free_ivector(int32_t *v, size_t nl);
    void     free_dvector(double *v, size_t nl);
    void     free_imatrix(int32_t **m, size_t nrl, size_t ncl);
    void     free_dmatrix(double **m, size_t nrl, size_t ncl);
    
    void     derivs( double y[], double dydx[] );
    void     jacobn( double y[], double dfdx[], double **dfdy );
    void     LU_Decomposition( double **a, size_t n, int32_t *indx, double& d );
    void     LU_BackSubstitution(double **a, size_t n, int32_t *indx, double b[] );
    
    void     simpr( double y[], double dydx[], double dfdx[], double **dfdy,
		            size_t n, double xs, double htot, int32_t nstep, double yout[] );

    void     pzextr(int32_t iest, double xest, double yest[], double yz[], double dy[], size_t nv);


    void     odeint( double ystart[], size_t nvar, double x1, double x2, 
                     double eps, double h1, double hmin, size_t& nok, size_t& nbad );
	             
    void     StiffBulirschStoer( double y[], double dydx[], size_t nv, double& xx, double htry, double eps,
	                               double yscal[], double& hdid, double& hnext );
	                              
    void     AllocateSolverMemory();
    void     FreeSolverMemory();
    void     ZeroAllStorage();
    void     TokenizeString( char* s, const char* delim, std::list<std::string>& li ) const;
    
    double   DSQR( double a ) const;
 };


/// assign initial conditions and read computed values (equations 1...n)
inline  double& ODE_StiffSolver::operator()( size_t i ) 
  { 
      assert ( i>0 && i<=dof ); 
      return init_cond[i]; 
  }
  
/// retrieve the time derivatives from the last step before the solver quit
inline  double  ODE_StiffSolver::operator[]( size_t i ) const 
  { 
      assert ( i>0 && i<=dof ); 
      return odeint_dydx[i]; 
  }


inline  void       ODE_StiffSolver::Tolerance( double epsilon ) { eps=epsilon; }
inline  double  ODE_StiffSolver::Tolerance() const           { return eps; }
inline  size_t  ODE_StiffSolver::Equations() const           { return dof; }
inline  size_t  ODE_StiffSolver::IntermediateStepsToSave() const { return kmax; }

inline  double  ODE_StiffSolver::DSQR( double a ) const { return (a==0.0) ? 0.0 : a*a; }

} // csmp

#endif

