#include "ODE_StiffSolver.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

/**
 
Builds the ODE solver, initializing the Kaps-Rentrop parameters for the
adaptive stepsize algorithm and a number of other constants used by the
solver. GROW and SHRNK define by how much the stepsize can increase or
decrease in one step, and TINY defines a threshold value below which 
entries in the solution matrix will be zeroed. Memory is allocated for 
the solution and jacobian matrices and vectors which store temporary 
variables. The system of equations is build. 

@param equations A vector of ordinary differential equations which define system to be
integrated by the solver. The equations are stored as STL maps inside
of the 'SumOfProductsWithExponents' class objects. 

@section implementation Implementation

The constructor calls the private method AllocateSolverMemory() to 
allocate the matrices and vectors. 

@section application Application

Standard constructor for the ODE solver. 
 */
ODE_StiffSolver::ODE_StiffSolver( const vector<SumOfProductsWithExponents>& equations )
  : NR_END( 1 ),
    SAFETY(   0.9 ),
    GROW(     1.5 ),
    PGROW(   -0.25 ),
    SHRNK(    0.5 ),
    PSHRNK(  -1.0/3.0),
    ERRCON(   0.1296), 
    MAXTRY(  40),
    GAM (    1.0/2.0),
    A21(     2.0),
    A31 (   48.0/25.0),
    A32 (    6.0/25.0),
    C21(    -8.0),
    C31 (372.0/25.0),
    C32 (12.0/5.0),
    C41 (-112.0/125.0),
    C42 (-54.0/125.0),
    C43 (-2.0/5.0),
    B1 (19.0/9.0),
    B2 (1.0/2.0),
    B3 (25.0/108.0),
    B4 (125.0/108.0),
    E1 (17.0/54.0),
    E2 (7.0/36.0),
    E3( 0.0),
    E4 (125.0/108.0),
    C1X (1.0/2.0),
    C2X (-3.0/2.0),
    C3X (121.0/50.0),
    C4X (29.0/250.0),
    A2X( 1.0),
    A3X (3.0/5.0),
    // Bulirsch-Stoer method
    KMAXX(7),
    IMAXX(7+1),
    SAFE1(0.25),
    SAFE2(0.7),
    REDMAX(1.0e-5),
    REDMIN(0.7),
    SCALMX(0.1),
    //
    TINY( 1.0e-50 ),                // values below this will be zeroed
    MAXSTP( 1000000 ),              // maximum timesteps of ODE solver that are permitted
    odes( equations ),
    dof( equations.size() ),        // number of ODEs in system
    eps(1.0e-10),                   // permitted tolerance from step to step
    kmax(0),                        // maximum number of intermediate timesteps to save
    correct_interpretation(false)   // terminates program when asserted
 {
     AllocateSolverMemory();
 }




ODE_StiffSolver::~ODE_StiffSolver()
 {
     FreeSolverMemory();
 }


/**
 
Allocates matrix and vector memory using the auxiliary routines of the
numerical recipes 'nrutil.c' source. The size of these arrays is
determined by the solverinternal parameters 'dof' and 'KMAXX'. The first
refers to the number of equations (degrees of freedom). The second 
denotes how many intermediate steps of the solution shall be saved.   

@section implementation Implementation

Method uses the Numerical Recipes functions dvector() and dmatrix().  

@section application Application

Called  by the constructor of the ODE solver.  
*/
void ODE_StiffSolver::AllocateSolverMemory()
 {
    // starting values and final values
    init_cond    = dvector(1,dof);
 
    // odeint
    odeint_yscal = dvector(1,dof);
    odeint_y     = dvector(1,dof);
    odeint_dydx  = dvector(1,dof);
	
	// stifbs
    pzextr_d=dmatrix(1,dof,1,KMAXX);
    stifbs_dfdx=dvector(1,dof);
    stifbs_dfdy=dmatrix(1,dof,1,dof);
	stifbs_err=dvector(1,KMAXX);
	pzextr_x=dvector(1,KMAXX);
	stifbs_yerr=dvector(1,dof);
	stifbs_ysav=dvector(1,dof);
	stifbs_yseq=dvector(1,dof);
 }



/** De-allocates the temporary storage vectors and matrices used by the solver.

@section implementation Implementation

Calls the Numerical Recipes functions free_dvector() and free_dmatrix(). 

@section application Application

Called by the destructor of the ODE solver. 
*/
void ODE_StiffSolver::FreeSolverMemory()
 {
    // starting values and final values
    free_dvector(init_cond,1);
 
    // odeint
	free_dvector(odeint_dydx,1);
	free_dvector(odeint_y,1);
	free_dvector(odeint_yscal,1);
	
	// stifbs
    free_dvector(stifbs_yseq,1);
	free_dvector(stifbs_ysav,1);
	free_dvector(stifbs_yerr,1);
	free_dvector(pzextr_x,1);
	free_dvector(stifbs_err,1);
	free_dmatrix(stifbs_dfdy,1,1);
	free_dvector(stifbs_dfdx,1);
	free_dmatrix(pzextr_d,1,1);
 } 


/** Initializes all the values of the temporary storage arrays to zero.
*/
void ODE_StiffSolver::ZeroAllStorage()
 {
    uint32_t i, j;

    // starting values and final values
    for ( i=1; i<=dof; i++ )
      {
// is assigned         init_cond[i]    = 0.0;
         odeint_yscal[i] = 0.0;
         odeint_y[i]     = 0.0;
         odeint_dydx[i]  = 0.0;
         stifbs_dfdx[i]  = 0.0;
	     stifbs_yerr[i]  = 0.0;
	     stifbs_ysav[i]  = 0.0;
	     stifbs_yseq[i]  = 0.0;
      }
    for ( i=1; i<=dof; i++ )
      for ( j=1; j<=dof; j++ ) stifbs_dfdy[i][j] = 0.0;

    if ( KMAXX > 0 )
      {
         for ( i=1; i<=dof; i++ )
           for ( j=1; j<=KMAXX; j++ ) pzextr_d[i][j] = 0.0;
	
         for ( j=1; j<=KMAXX; j++ ) 
           {
	      stifbs_err[j] = 0.0;
	      pzextr_x[j]   = 0.0;
           }
      }
 } 


/**
 
Defines the number of intermediate solution steps which shall be saved 
before the integration interval is covered. The storage for these steps is
allocated correspondingly. The solver attempts to distribute the save
steps evenly over the timespan of interest.  

@param s The number of steps which shall be saved.
*/
void ODE_StiffSolver::IntermediateStepsToSave( size_t s )
 {
    if ( s > result_xp.size() )
      {
         result_xp.reserve(s);
         for ( size_t i=kmax-1; i<s; i++ ) result_xp.push_back(0.0);
         for ( size_t j{0U}; j<dof; j++ ) 
           {
               result_yp[j].reserve(s);
               for ( size_t i=kmax-1; i<s; i++ ) result_yp[j].push_back(0.0); 
           }
      }
    kmax = s;
 }



/**
 
Sums the righthand sides of each ordinary differential equation in the
system given the current values of the dependent variable and returns
the results. This procedure is also illustrated in Numerical Recipes
2nd edition p. 742.   

@param y a C-array which supplies the values of the
dependent variables in the system. 

@param dydx The second argument is a C-array into which the current sums over the ODEs
will be returned. 

@section implementation Implementation

The method calls the Sum() method of each of the SumOfProductsWithExponents
objects which hold the ODEs. 

@section application Application

Called at each step which the ODE solver takes. 
*/
void ODE_StiffSolver::derivs( double y[], double dydx[] )
 {
    vector<SumOfProductsWithExponents>::const_iterator  oit;
    int32_t n;
    
    for ( n=1, oit=odes.begin(); oit!=odes.end(); oit++ )
      dydx[n++] = (*oit).Sum( y );

 } // end derivatives



/**
 
Fills the jacobian matrix with the current derivatives of the ODEs. The
jacobian matrix has one row for each ordinary differential equation. This
row contains the derivatives of this ODE with regard to each 
dependent variable in the system in the same order as the differential
equations describe these. See example in Numerical Recipes p. 742.  

@section arguments Input Arguments

The first and the second C-style array arguments supply the current 
values of the dependent variable and the time-derivatives thereof to
jacobn(). The time derivatives are simply set to an initial guess of zero.  

@param dfdy The jacobian matrix is returned into the third C-style 2D array argument.

@section implementation Implementation

Calls the DerivativeWithRespectTo() method of the SumOfProductsWithExponents
class to obtain the specific derivative entries for the jacobian matrix. 

@section application Application

Called at each step of the ODE solver. 
*/
void ODE_StiffSolver::jacobn( double y[], double dfdx[], double** dfdy )
 {
    size_t i, j;
    
    for ( i=1; i<=dof; i++ ) dfdx[i] = 0.0;
    for ( i=1; i<=odes.size(); i++ )
      for ( j=1; j<=odes.size(); j++ )
        // differentiate the i'th equation w.r.t. j'th species
        dfdy[i][j] = odes[i-1].DerivativeWithRespectTo( static_cast<int32_t>(j)-1, y );        
    
 } // end jacobn





/**
 
Performs an LU decomposition of the global solution matrix as is described
in section 2.3 of Numerical Recipes (p. 43). 

@section implementation Implementation

The method is described in detail in Numerical Recipes p. 46. 

@section application Application

Is used together with lubksb() to invert the system of ODEs at each
solution step. 

@section messages Messages

The method may print a warning if the solution matrix becomes singular
in the diagonalization process. 
*/
void ODE_StiffSolver::LU_Decomposition( double **a, size_t n, int32_t *indx, double& d )
{
	int32_t i,imax = std::numeric_limits<double>::quiet_NaN(),j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) cout <<"ODE_StiffSolver::LU_Decomposition: Singular matrix in routine ludcmp"<< endl;
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1);
}



/**
 
Returns the solution vector of the system Ax = b given an LU decomposition of 
the solution matrix A, a righthand vector b and an indicator whether an even 
or an odd number of partial pivots were performed.  

@section implementation Implementation

The method is documented in Numerical Recipes p. 47. 

@section application Application

LU_BackSubstitution() is called once for each step of the solver. 
*/
void ODE_StiffSolver::LU_BackSubstitution(double **a, size_t n, int32_t *indx, double b[])
{
	size_t  i,ii=0,j;
	int32_t      ip;
	double  sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}




/**
 
Sets the initial values of the dependent variables for the starting step
of the calculation.
In the case of a chemical rate-dependent speciation calculation with the
ODE solver, initial conditions are the starting activities of the 
aqueous species and minerals. 

@param init A constant reference to an STL vector with the initial dependent variable
values. 
*/
void ODE_StiffSolver::InitialConditions( const vector<double>& init )
 {
     assert( init.size() == dof+1 );
     for ( size_t i=1; i<=dof; i++ ) init_cond[i] = init[i];
 }





/**
 
Evolves the system of equations from the starting time (t=0) over a
user-specified period of time and returns the dependent variable values
at the end of the timespan, as well as intermediate values if this is
requested.  

@section arguments Input Arguments

The first argument defines the timespan over which the system shall be 
evolved and the second argument sets the size of the first timestep which
shall be calculated. If this initial timestep is too large and it already
bridges pronounced changes which would normally occur in the system, no 
solution may be obtained.  

The third argument gives the number of intermediate steps which shall 
be saved in the solution process. The default value is zero since one 
is usually only interested in the final solution. In this case the 
solver will perform significantly faster.  

The fourth argument defines the precision to which the dependent variable
values shall be approximated to at each timestep. This parameter strongly 
influences the number of iterations which the solver must perform to find 
a solutions. 

@return a boolean variable indicating whether a solution was
found or not. The solution can be retrieved simply with the () operator
of the ODE-solver object. 

@section implementation Implementation

The method relies on the Numerical Recipes C functions odeint(), stifbs(),
simpr(), and pzextr() as documented on p. 744ff in the book. The results
will be saved into a file with the same name as the equation systems
file but with the extension '.results'. 

@section application Application

For practical purposes the initial step size and tolerance which can
be obtained varies considerable with the type of OD equations which
shall be integrated over time. For the chemical reaction calculations 
carried out so far, an initial step size of 1.0e-17 and a tolerance of
1.0e-10 has proved useful. The solver will rapidly increase the 
stepsize if it is too small, but it will fail to converge if the first
step was too large. It seems impossible to achieve convergence at a 
tolerance smaller than 1.0e-12 in most of the examined systems. 

Before using a system of equations in a finite-element model, a series
of tests for a host of initial values is recommended to get a feel for
the time-characteristics of the system behaviour and to tune the 
Solve() arguments accordingly. This is best achieved by saving about
100 steps and examining these with a plot program to understand 
system behaviour. 

@section messages Messages

Upon completion, the solver returns the number of successful and bad 
solution steps which were performed. Other messages may include those
issued by methods which are called internally, for instance, when the 
matrix inversion routine failed to find a solution. Usually, a large
percentage of bad steps indicates a problem with the system of equations.
The most common problem may be that the system is too stiff, i.e. that 
the dependent variable coefficients diverge by more than 16 orders of
magnitude (about the limit for this double-precision implementation of
the solver). 
*/
bool ODE_StiffSolver::Solve(  double time, double init_stepsize, 
                              size_t n_savesteps, double tolerance,
                              bool print_step_numbers )
 {
//    ZeroAllStorage();

    // compute from x1 to x2
	double  hstart, x1 = 0., x2 = time;
	size_t  nbad, nok, i, j;

    eps    = tolerance;
    hstart = init_stepsize;
	kmax   = n_savesteps;
    if ( n_savesteps > 0 ) dxsav = init_stepsize;
    else dxsav = 0.0;
    
    // setting up storage for intermediate steps
    // O.K.
    if ( kmax > result_xp.size() )
      {
         result_xp.reserve(kmax); 
         for ( i=result_xp.size(); i<kmax; i++ ) result_xp.push_back(0.0);

         // building result matrix (output to elements 1...n)
         result_yp.reserve(dof);
         for ( i=0; i<dof; i++ )
            {
               result_yp.push_back( vector<double>() );
               result_yp[i].reserve(kmax);
               for ( j=0; j<kmax; j++ )  result_yp[i].push_back(0.0);
            }
      }
    // solving differential equations (hmin=0.0) 
    odeint(init_cond,dof,x1,x2,eps,hstart,0.0,nok,nbad);
    
    if ( print_step_numbers ) {
 	     printf("\n%s %13s %3lu\n","successful steps:"," ",nok);
         printf("%s %20s %3lu\n","bad steps:"," ",nbad);
      }
    return true;

 } // end Solve




/** Prints the results including intermediate steps to the screen.
*/
void ODE_StiffSolver::OutputResults() const
 {
	size_t i, j;

	cout <<"\nODE_StiffSolver::Solve: solution" << endl;
	cout <<"time   results ODE_n(t) 1...n" << endl;
	if ( kmax > 0 )
	  {
	      for ( i=0; i<kount; i++ )
	        {
	            // time
	            cout << result_xp[i] <<"  ";
	            // solutions ODE_0 to ODE_n
	            for ( j=0; j<dof; j++ ) cout << result_yp[j][i] <<"  ";
	            cout << endl;
	        } 
	  }
	else 
	  {
	     cout <<"\nNo intermediate steps were savecd. Final values:"<< endl;
	     for ( i=1; i<=dof; i++ ) cout << init_cond[i] <<"  ";
	     cout << endl;  
	  }

 } // end OutputResults





/**

Writes the results of a calculation to a user-defined textfile with the
extension '.txt'. Each line will give the timestep and the values of each
dependent variable. The output may be corrupted when not enough savesteps
were specified to accurately portray the system behaviour. 

@param textfile The name of the textfile to which the results shall be written to.
*/
void ODE_StiffSolver::OutputResultsTo( const char* textfile ) const
 {
    char name[200];
    strcpy ( name, textfile );
    strcat( name, ".txt" );

    ofstream  ofs( name );
    
    assert( ofs.is_open() );
   
    ofs.setf( ios::scientific );
    
	size_t i, j;

	ofs <<"ODE_StiffSolver::Solve: solution obtained:" << endl;
	ofs <<"time   results ODE_n(t) 1...n" << endl;
	if ( kmax > 0 )
	  {
	      for ( i=0; i<kount; i++ )
	        {
	            // time
	            ofs << result_xp[i] <<"  ";
	            // solutions ODE_0 to ODE_n
	            for ( j=0; j<dof; j++ ) ofs << result_yp[j][i] <<"  ";
	            ofs << endl;
	        } 
	  }
	else 
	  {
	     ofs <<"\nNo intermediate steps were savecd. Final values:"<< endl;
	     for ( i=1; i<=dof; i++ ) ofs << init_cond[i] <<"  ";
	     ofs  << endl;  
	  }
    ofs.close();

    cout <<"\nODE_StiffSolver::OutputResultsTo: file '"<< name;
    cout <<"' written successfully..." << endl;

 } // end OutputResultsTo






/**
 
odeint() is the integration routine which is used in most of the 
Numerical Recipe algorithms for differential equations. It is basically
a scheme that uses LU decomposition and back substitution to invert the
system of equations and to implement timestepping on the basis of criteria
which are established by sub-functions. The method is documented in
Numerical Recipes p. 721. 

 
@param ystart starting values to integrate over (starting concentrations)
@param nvar number of equations
@param input_eps tolerance
@param h1 guess of initial stepsize
@param hmin minimum permitted stepsize
@param nok good steps taken
@param nbad unecessary steps taken but fixed
 

@section implementation Implementation

See numerical recipes pages 721 and 744ff. 

@section application Application

Used by Solve() to evolve the system of ODEs. 

@section messages Messages

The method will exit and issue a warning if the limit of steps has been
reached without finding a solution. 
*/
void ODE_StiffSolver::odeint( double ystart[], size_t nvar, double x1, 
                              double x2, double input_eps, double h1,
	                            double hmin, size_t& nok, size_t& nbad )
{
	size_t  nstp, i;
  double nan(numeric_limits<double>::quiet_NaN());
	double xsav(nan), x, hnext, hdid;
	x=x1;
//	h=copysign(h1,x2-x1); this is not available on many platforms
	double h = fabs(h1) * sign(x2-x1);
	nok = nbad = kount = 0;
	for (i=1;i<=nvar;i++) odeint_y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;

	for (nstp=1;nstp<=MAXSTP;nstp++) 
          {
                // DERIVATIVES
		derivs(odeint_y,odeint_dydx);
		for (i=1;i<=nvar;i++) odeint_yscal[i]=fabs(odeint_y[i])+fabs(odeint_dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav) )
                  {
	             result_xp[kount++]=x;
		     for (i=1;i<=nvar;i++) result_yp[i-1][kount-1]=odeint_y[i];
		     xsav=x;
		  }
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

                // STIFF SOLVER
		StiffBulirschStoer(odeint_y,odeint_dydx,nvar,x,h,input_eps,odeint_yscal,hdid,hnext);

		if (hdid == h) ++nok; else ++nbad;

		if ((x-x2)*(x2-x1) >= 0.0) 
                  {
		     for (i=1;i<=nvar;i++) ystart[i]=odeint_y[i];
	             if (kmax) 
                       {
			  result_xp[kount++]=x;
			  for (i=1;i<=nvar;i++) result_yp[i-1][kount-1]=odeint_y[i];
		       }
		     return;
		  }
		if (fabs(hnext) <= hmin) 
                  {
                     cout <<"ODE_StiffSolver::odeint: Step size too small in odeint"<< endl;
                     assert( correct_interpretation );
                  }
		h=hnext;

                if ( !independent_terms.empty() ) {
                     EvaluateIndependentProductTerms( odeint_y, independent_values );
                     ResetIndependentProductTerms( independent_terms, independent_values );
                  }
	}
      cout <<"ODE_StiffSolver::odeint: Too many steps in routine odeint"<< endl;
}





/**
 
Semi-implicit extrapolation method by  Bader and Deuflhard (1986) for
stiff systems of ordinary differential equations. See Numerical Recipes
pages 742ff for documentation. 


@param y array of dependent variables
@param dydx array of (time)derivatives of dep. vars
@param nv number of equations
@param xx starting value of independent variable, e.g., t=0
@param htry initial try of stepsize
@param input_eps tolerance
@param yscal vector of values to scale error against
@param hdid accomplished stepsize
@param hnext stepsize if a new step was to follow
 
New values of x(time) and y(derivatives) are returned into respective
arguments. 

@section implementation Implementation

See Numerical Recipes p. 742ff. 

@section application Application

Called by odeint(). 

@section messages Messages

If the stepsize reaches zero, the method has failed and a stepsize 
underflow is reported. 
*/
void ODE_StiffSolver::StiffBulirschStoer( double y[], double dydx[], size_t nv, 
                                          double& xx, double htry, double input_eps,
	                                        double yscal[], double& hdid, double& hnext )
{
	int32_t i,iq,k,kk,km= std::numeric_limits<double>::quiet_NaN();
	static int32_t first=1,kmax1,kopt,nvold = -1;
	static double epsold = -1.0,xnew;
	double eps1,errmax= std::numeric_limits<double>::quiet_NaN(),fact,h,red,scale=std::numeric_limits<double>::quiet_NaN(),work,wrkmin,xest;
	static double a[9];                            // 9=IMAXX+1
	static double alf[8][8];                       // 8=KMAXX+1
	static int32_t nseq[9]={0,2,6,10,14,22,34,50,70};  // 8=IMAXX+1
	int32_t reduct,exitflag=0;

	if(input_eps != epsold || nv != nvold) {
		hnext = xnew = -1.0e29;
		eps1=SAFE1 * input_eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[1]+1.0)*(2*k+1))));
		}
		epsold=input_eps;
		nvold=static_cast<int32_t>(nv);
		a[1] += nv;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax1=kopt;
	}
	h=htry;
	for (i=1;i<=nv;i++) stifbs_ysav[i]=y[i];
	jacobn(y,stifbs_dfdx,stifbs_dfdy);
	if (xx != xnew || h != hnext) {
		first=1;
		kopt=kmax1;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax1;k++) {
			xnew=xx+h;
			if (xnew == xx) 
                          {
                             cout <<"ODE_StiffSolver::StiffBulirschStoer: step size underflow in stifbs"<< endl;
                             assert( correct_interpretation );
                          }
			simpr(stifbs_ysav,dydx,stifbs_dfdx,stifbs_dfdy,nv,xx,h,nseq[k],stifbs_yseq);
			xest=DSQR(h/nseq[k]);
			pzextr(k,xest,stifbs_yseq,y,stifbs_yerr,nv);
			if (k != 1) {
				errmax=TINY;
                                for (i=1;i<=nv;i++) errmax=std::max(errmax,fabs(stifbs_yerr[i]/yscal[i]));
				errmax /= input_eps;
				km=k-1;
				stifbs_err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax1 || k == kopt+1) {
					red=SAFE2/stifbs_err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < stifbs_err[km]) {
						red=1.0/stifbs_err[km];
						break;
					}
				else if (kopt == kmax1 && alf[km][kmax1-1] < stifbs_err[km]) {
						red=alf[km][kmax1-1]*SAFE2/stifbs_err[km];
						break;
					}
				else if (alf[km][kopt] < stifbs_err[km]) {
					red=alf[km][kopt-1]/stifbs_err[km];
					break;
				}
			}
		}
		if (exitflag) break;
                red=std::min(red,REDMIN);
                red=std::max(red,REDMAX);
		h *= red;
		reduct=1;
	}
	xx=xnew;
	hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=1;kk<=km;kk++) {
            fact=std::max(stifbs_err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	hnext=h/scale;
	if (kopt >= k && kopt != kmax1 && !reduct) {
            fact=std::max(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			hnext=h/fact;
			kopt++;
		}
	}
}



/**
 
Polynomial extrapolation method used by the stiff Burlisch-Stoer method. 
See Numerical Recipes p. 731 for documentation. 

@section application Application

See Numerical Recipes p. 731 for documentation. 

@section application Application

Used by StiffBulirschStoer(). 
*/
void ODE_StiffSolver::pzextr( int32_t iest, double xest, double yest[], 
                              double yz[], double dy[], size_t nv )
{
	int32_t k1,j;
	double q,f2,f1,delta,*c;

	c=dvector(1,nv);
	pzextr_x[iest]=xest;
	for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 1) {
		for (j=1;j<=nv;j++) pzextr_d[j][1]=yest[j];
	} else {
		for (j=1;j<=nv;j++) c[j]=yest[j];
		for (k1=1;k1<iest;k1++) {
			delta=1.0/(pzextr_x[iest-k1]-xest);
			f1=xest*delta;
			f2=pzextr_x[iest-k1]*delta;
			for (j=1;j<=nv;j++) {
				q=pzextr_d[j][k1];
				pzextr_d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=1;j<=nv;j++) pzextr_d[j][iest]=dy[j];
	}
	free_dvector(c,1);
}



/**
 
Performs semi-implicit midpoint-rule based timestepping inside of 
odeint(). See documentation p. 743 in Numerical Recipes. 

@section implementation Implementation

See documentation p. 743 in Numerical Recipes. 

@section application Application

Called by odeint(). 
*/
void ODE_StiffSolver::simpr( double y[], double dydx[], double dfdx[], 
                             double **dfdy, size_t n, double xs, double htot, 
                             int32_t nstep, double yout[] )
{
	int32_t i,j,nn,*indx;
	double d,h,x,**a,*del,*ytemp;

	indx=ivector(1,n);
	a=dmatrix(1,n,1,n);
	del=dvector(1,n);
	ytemp=dvector(1,n);
	h=htot/nstep;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
		++a[i][i];
	}
	LU_Decomposition(a,n,indx,d);
	for (i=1;i<=n;i++)
		yout[i]=h*(dydx[i]+h*dfdx[i]);
	LU_BackSubstitution(a,n,indx,yout);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	derivs(ytemp,yout);
	for (nn=2;nn<=nstep;nn++) {
		for (i=1;i<=n;i++)
			yout[i]=h*yout[i]-del[i];
		LU_BackSubstitution(a,n,indx,yout);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		derivs(ytemp,yout);
	}
	for (i=1;i<=n;i++)
		yout[i]=h*yout[i]-del[i];
	LU_BackSubstitution(a,n,indx,yout);
	for (i=1;i<=n;i++)
		yout[i] += ytemp[i];
	free_dvector(ytemp,1);
	free_dvector(del,1);
	free_dmatrix(a,1,1);
	free_ivector(indx,1);
}


// allocate an int vector with subscript range v[nl..nh] 
int32_t *ODE_StiffSolver::ivector(size_t nl, size_t nh)
{
	int32_t* v =(int32_t*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int32_t)));
	if (!v) cout<<"\nallocation failure in ivector()"<< endl;
	return v-nl+NR_END;
}


// allocate a double vector with subscript range v[nl..nh] 
double *ODE_StiffSolver::dvector(size_t nl, size_t nh)
{
	double* v =(double*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) cout<<"\nallocation failure in dvector()"<< endl;
	return v-nl+NR_END;
}


// allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] 
double **ODE_StiffSolver::dmatrix(size_t nrl, size_t nrh, size_t ncl, size_t nch)
{
	size_t i, nrow=nrh-nrl+1,ncol=nch-ncl+1;

	// allocate pointers to rows 
	double** m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) cout<<"\nallocation failure 1 in matrix()"<< endl;
	m += NR_END;
	m -= nrl;

	// allocate rows and set pointers to them 
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) cout<<"\nallocation failure 2 in matrix()"<< endl;
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	// return pointer to array of pointers to rows
	return m;
}


// allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] 
int32_t **ODE_StiffSolver::imatrix(size_t nrl, size_t nrh, size_t ncl, size_t nch)
{
	size_t i, nrow=nrh-nrl+1,ncol=nch-ncl+1;

	// allocate pointers to rows 
	int32_t** m=(int32_t **) malloc((size_t)((nrow+NR_END)*sizeof(int32_t* )));
	if (!m) cout<<"\nallocation failure 1 in matrix()"<< endl;
	m += NR_END;
	m -= nrl;


	// allocate rows and set pointers to them 
	m[nrl]=(int32_t*) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int32_t)));
	if (!m[nrl]) cout<<"\nallocation failure 2 in matrix()"<< endl;
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	// return pointer to array of pointers to rows
	return m;
}


// free an int vector allocated with ivector()
void ODE_StiffSolver::free_ivector( int32_t *v, size_t nl )
{
	free( (char*) (v+nl-NR_END));
}


// free a double vector allocated with dvector() 
void ODE_StiffSolver::free_dvector( double *v, size_t nl )
{
	free( (char*) (v+nl-NR_END));
}


// free a double matrix allocated by dmatrix() 
void ODE_StiffSolver::free_dmatrix( double **m, size_t nrl, size_t ncl )
{
	free( (char*) (m[nrl]+ncl-NR_END));
	free( (char*) (m+nrl-NR_END));
}


// free an int matrix allocated by imatrix() 
void ODE_StiffSolver::free_imatrix( int32_t **m, size_t nrl, size_t ncl )
{
	free( (char*) (m[nrl]+ncl-NR_END));
	free( (char*) (m+nrl-NR_END));
}






/**
 
Prints the defined differential and potential algebraic equations
involving independent variables to screen. The differential equation are
stored as SumOfProductsWithExponents objects and the algebraic equations
as SumOfProducts objects. The Out() methods of these classes are used
to obtain the output.   

@section implementation Implementation

Calls the Out() interfaces of the SumOfProductsWithExponents and
SumOfProducts class objects. 

@section application Application

To externalize the state of the solver to examine whether the system
of equations has been correctly parsed from the input file. 
*/
void ODE_StiffSolver::Out() const
 {
    int32_t n;
    cout <<"\nODE_StiffSolver::Out:"<< endl;
    cout <<"\nrighthand sides of ODEs: "<< endl;
    vector<SumOfProductsWithExponents>::const_iterator  it;
    vector<SumOfProducts>::const_iterator               ait;
    for ( n=0, it=odes.begin(); it!=odes.end(); it++, n++ ) 
      {
         cout <<"\nEquation: "<< n; 
         (*it).Out();
      }
    cout << endl;
    if ( !eqns.empty() )
      {
          cout <<"\nAdditional algebraic equations: "<< endl;
          for ( n=-1, ait=eqns.begin(); ait!=eqns.end(); ait++, n-- )
            {
               cout << n <<" ";
               (*ait).Out();
            }
      }
 }                  



/** Prints an stl vector, or a C 1- or 2D-array indexed from 1...n to screen.

@param vec The template object to be printed.

@section implementation Implementation

All of these output routines assume an indexing from 1...n. 
*/
void ODE_StiffSolver::Out( const vector<int32_t>& vec ) const
 {
    cout <<"\nvector n ("<< vec.size() <<")"<< endl;
    for ( size_t i{0U}; i<vec.size(); i++ ) cout << vec[i] <<"  ";
    cout << endl;
 }                


void ODE_StiffSolver::Out( int32_t*  vec, int32_t n ) 
 {
    cout <<"\nvector n ("<< n <<")"<< endl;
    for ( int32_t i=1; i<=n; i++ ) cout << vec[i] <<"  ";
    cout << endl;
 }


void ODE_StiffSolver::Out( double* vec, int32_t n ) 
 {
    cout <<"\nvector n ("<< n <<")"<< endl;
    for ( int32_t i=1; i<=n; i++ ) cout << vec[i] <<"  ";
    cout << endl;
 }
 
                
void ODE_StiffSolver::Out( double** vec, int32_t m, int32_t n )
 {
    cout <<"\nmatrix m x n ("<< m <<","<< n <<")"<< endl;
    for ( int32_t i=1; i<=m; i++ )
      {
         for ( int32_t j=1; j<=n; j++ ) cout << vec[i][j] <<"  ";
         cout << endl;
      }
    cout << endl;  
 }


void ODE_StiffSolver::Out( const vector<double>& vec ) const
 {
    cout <<"\nvector n ("<<vec.size() <<")"<< endl;
    for ( size_t i{0U}; i<vec.size(); i++ ) cout << vec[i] <<"  ";
    cout << endl;
 }                



void ODE_StiffSolver::Out( const vector<vector<double> >& vec ) const
 {
    cout <<"\nmatrix n x m ("<<vec.size()-1 <<","<< vec[0].size()-1 <<")"<< endl;
    for ( size_t i=1; i<vec.size(); i++ )
      {
         for ( size_t j=1; j<vec[i].size(); j++ ) cout << vec[i][j] <<"  ";
         cout << endl;
      }
    cout << endl;  
 }                



/**

Creates an STL list of tokens from the input string using user-defined token
delimiters. 

@param s name of the string to be tokenized (separated into segments)
@param delim a list of single character delimters ("_,+,\0" for instance).

@param li the list of tokens supplied as third argument. If this list already
contains tokens, these are erased before new tokens are added. 

*/
void ODE_StiffSolver::TokenizeString( char* s, const char* delim, 
                                      list<string>& li ) const
  {
     li.erase( li.begin(), li.end() );
     char* tok;
     tok = strtok( s, delim );
     do {
           li.push_back( string(tok) );
        }
     while ( (tok=strtok( NULL, delim )) != NULL ); 
  }



/**
 
Initializes the system of ordinary differential equations maintained by 
the ODE solver from a textfile with the format specified in the documentation
of the ODE solver object (\par ApplicationExample section). 

@param file the name of the textfile which contains the definition of the system of
differential equations. 

The textfile itself consists of 5 text blocks always starting with an
underlined heading and separated by blank lines. The blocks are: the file 
header (1), the independent variables (2), the dependent variables (3), 
the equation coefficients (4), and the equations themselves (5). Text 
after the blank line following the last equation is ignored. This gives
the user the possibility to add comments. 

@section implementation Implementation

The method uses STL functionality to parse the textfile and store the 
interpreted equations into a vector of SumOfProductsWithExponents objects. 

@section messages Messages

The method tries to analyze the reasons for parsing mistakes and reports
corresponding messages to the user. Common mistakes include missing blank
lines between definition blocks or mispelled 
independent or dependent variable names or the use of reserved '+,-'
notation in variable names. An error is reported also when the number
of differential equations does not correspond to the number of 
dependent variables. 
*/
void ODE_StiffSolver::ReadODEsFrom( const char* file )
 {
    char    c, text_line[3000];
    double  exp, K= std::numeric_limits<double>::quiet_NaN();
    int32_t coeff_idx(0), idx, i, eq_counter = 0;
    char    *token, *sub1, *sub2, 
            temp[3000], temp1[3000], temp2[3000];
    const char* const white_delims =" ,\t,\n,\r";
    const char* const lhs_delim    ="=";
    const char* const sum_delim    = "+";
    const char* const coeff_delim  =" ";
         
    map<int32_t,pair<double,double> >  product; 
	  SumOfProductsWithExponents         ode;
	  list<string>                       tokens, tokens1;
	  list<string>::iterator    tik, tik1;

    strcpy( text_line, file );
    ifstream  ifs( text_line );    
    assert( ifs.is_open() );

       
    // 1. reading headlines and species names later used in the equations O.K.
    // ------------------------------------------------------------------
    ifs.getline( text_line, 256 ); // title line
    cout <<"\nODE_StiffSolver::ReadODEsFrom: reading: "<< text_line << endl;
    ifs.getline( text_line, 256 ); // line of ========
    ifs.getline( text_line, 256 ); // empty line
    ifs.getline( text_line, 256 ); // definition of species line
    ifs.getline( text_line, 256 ); // line of --------
    ifs.getline( text_line, 256 ); // SPECIES definition
    map<string,int>  species;
    int32_t n = 0; // species names must not contain whitespace
    token = strtok( text_line, white_delims ); 
    do {
          species[ string(token) ] = n++;
       }
    while ( (token=strtok(NULL, white_delims)) != NULL );
    ifs.getline( text_line, 256 ); // get next blank line
   
    // 2. reading equilibrium constants for reactions (same sequence as species) O.K.
    // -------------------------------------------------------------------------
    // spec.idx
    map<string,pair<int,double> >  coeffs;
    string  kname;
    ifs.getline( text_line, 256 ); // title line
    ifs.getline( text_line, 256 ); // line of -------
    ifs.getline( text_line, 256 ); // COEFFICIENT definition
    TokenizeString( text_line, white_delims, tokens );
    for ( n=0, tik=tokens.begin(); tik!=tokens.end(); tik++ ) 
      {
          strcpy( temp1, (*tik).c_str() );
          // now the substrings must be broken into names and values
          // name
          sub1  = strtok( temp1, lhs_delim );
          kname = sub1;
          // constant
          sub1  = strtok( NULL, lhs_delim );
          K     = atof( sub1 );
          // recording
          coeffs[ kname ] = make_pair(n++,K);
       }
    ifs.getline( text_line, 256 ); // get next blank line
   
   
    // 3. reading ordinary differential equations using the terminology established
    //    above (if expressions do not comply errors are raised)
    // ----------------------------------------------------------------------------
    map<string,pair<int32_t,double> >::const_iterator  cit;
    map<string,int32_t>::const_iterator               sit;
    map<string,int32_t>::const_iterator               current_lhs;
    ifs.getline( text_line, 256 ); // title line
    ifs.getline( text_line, 256 ); // line of -------
    
    // removing previously defined reactions
    odes.erase( odes.begin(), odes.end() );

    // 4. reading all the ODEs in the file
    // -----------------------------------
    while ( !ifs.eof() )
      {
          // 4.1 reading ODE lhs and rhs
          // ---------------------------
          ifs.getline( text_line, 256 ); 
          token = strtok( text_line, lhs_delim ); // lhs
          if ( token == NULL ) break; // if empty line is reached
          strcpy( temp1, token );
          // terminating string where whitespace starts
          sub1=strchr( temp1,' '); *sub1 = '\0';
          // making sure that lhs is listed in order specified in species list (1)
          if ( (current_lhs=species.find(string(temp1))) != species.end() )
            {
               // checking ordering of equation list in file which must correspond to
               // listing of species above
               if ( eq_counter != (*current_lhs).second )
                 {
                     cout <<"\nODE_StiffSolver::ReadODEsFrom: ODEs are not listed in the order. ";
                     cout <<"of definitions of lefthandside terms given above." << endl;
                     ifs.close();
                     assert( correct_interpretation );
                 }
            }
          else
            {
                cout <<"\nODE_StiffSolver::ReadODEsFrom: ODEs lefthandside term must be. ";
                cout <<"ODE identifier given above. Unknown term: ";
                cout << temp1 << endl;
                ifs.close();
                assert( correct_interpretation );
            }

          // 4.2 tokenizing righthandside ODE into products (separator + )
          // -------------------------------------------------------------
          token = strtok( NULL, lhs_delim ); // rhs
          // saving rhs into a list  
          TokenizeString( token, sum_delim, tokens );
          odes.reserve( odes.size()+1 );

          for ( tik=tokens.begin(); tik!=tokens.end(); tik++ )  // for each product in 'rhs'
            {
               // 4.2.1 create subtokens corresponding to terms that make up each product
               // -----------------------------------------------------------------------
               strcpy( temp, (*tik).c_str() );
               TokenizeString( temp, coeff_delim, tokens1 );

               // 4.2.2coefficients always precede product terms
               // -----------------------------------------------
               tik1 = tokens1.begin();
               strcpy( temp1, (*tik1).c_str() );

               // test for plain coefficient
               if ( (cit=coeffs.find(string(temp1))) != coeffs.end() )
                 {
                    // getting coefficient (=equilibrium constant)
                    coeff_idx = (*cit).second.first;
                    K         = (*cit).second.second;
                 }
               else  // coefficients may be preceded by a minus which must be removed 
                 {   // -------------------------------------------------------------
                     for ( i=1; i<=strlen(temp1); i++ ) temp[i-1] = temp1[i];
                     if ( (cit=coeffs.find(string(temp))) != coeffs.end() ) 
                       {
                           coeff_idx = (*cit).second.first;
                           K         = -1.0 * (*cit).second.second;
                       }
                     else
                      {
                          cout <<"\nError when parsing coefficient in RHS of ODE. ";
                          cout <<"Unknown coefficient: "<< temp1 << endl;
                          ifs.close();
                          assert( correct_interpretation );
                      }
                 }
               tik1++;

               // 4.2.3 reading product terms
               // ---------------------------
               while ( tik1 != tokens1.end() )  // for each product term
                 {
                     strcpy( temp1, (*tik1).c_str() );
                     
                     // parsing product term
                     if ( (sit=species.find(string(temp1))) != species.end() )
                        {
                           idx = (*sit).second;
                           exp = 1.0;
                        }
                     // if product term is not found in the species list, 
                     // the term is checked for whether it 
                     // contains an exponent which needs to be separated
                     else
                       {
                           // return pointer to first occurrence of "^" in string
                           sub2 = strchr( temp1, '^' );
                           if ( sub2 == NULL ) {
                                 cout <<"\nError: Could not parse product term: "<< (*tik).c_str();
                                 cout <<" in ODE. Unknown term: "<< temp1 << endl;
                                 // checking for accidential '0' instead of 'O'
                                 if ( strchr(temp1,'0') != NULL )
                                    cout <<"Term contains zero(=0) character. Is this intended ?"<< endl;
                                 ifs.close();
                                 assert( correct_interpretation );
                             }
                           else
                             {
                                 // reading species name
                                 n=0;
                                 while ( (c=temp1[n]) != '^' ) temp2[n++] = c;
                                 temp2[n] = '\0';
                                 if ( (sit=species.find(string(temp2))) == species.end() ) {
                                      cout <<"\nError: Could not parse product term: "<< (*tik).c_str();
                                      cout <<" in ODE. Unknown term: "<< temp2 << endl;
                                      // checking for accidential '0' instead of 'O'
                                      if ( strchr(temp2,'0') != NULL )
                                        cout <<"Term contains zero(=0) character. Is this intended ?"<< endl;
                                      ifs.close();
                                      assert( correct_interpretation );
                                   }
                                 idx = (*sit).second;
                                 // parsing exponent (omitting ^ character)
                                 i=0;
                                 while ( (c=temp1[++n]) != '\0' ) temp2[i++] = c;
                                 temp2[i] = '\0';
                                 exp = atof( temp2 );
                              }
                        }
                     // adding new term to product map  
                     product[ idx ] = make_pair(1.0,exp);
                     tik1++;
                 }
               // 4.2.4 recording product terms
               // -----------------------------
               ode.AddProduct( coeff_idx, K, product );
               product.erase( product.begin(), product.end() );

            } // end for each product separated by " "

      	  // 5. Adding ODE to system of ODEs
	      // -------------------------------
	      odes.push_back( ode );
	      ode.Erase();
	      eq_counter++;
	      
      } // end while eof
      
    ifs.close();
    
    assert( odes.size() == species.size() );
    
    cout <<"\nODE_StiffSolver::ReadODEsFrom: file '"<< file;
    cout <<"' read successfully..." << endl;
    
    // 5. setting solver up for new system of equations
    // ------------------------------------------------
    FreeSolverMemory();
    dof = odes.size();
    AllocateSolverMemory();
    
} // end 



/**
 
Initializes the system of ordinary differential equations maintained by
the ODE solver from a textfile with the format specified in the documentation
of the ODE solver object (\par ApplicationExample section). The variable
coefficients, the dependent and independent variable identifiers are 
returned int STL vectors. 

@param file The name of the textfile which contains the definition of the system of
differential equations. 

The textfile itself consists of 5 text blocks always starting with an
underlined heading and separated by blank lines. The blocks are: the file 
header (1), the independent variables (2), the dependent variables (3), 
the equation coefficients (4), and the equations themselves (5). Text 
after the blank line following the last equation is ignored. This gives
the user the possibility to add comments. 

Returs three STL vectors of int / string pairs
corresponding to the coefficients of the dependent variables and their
indices giving their position in the differential equations (1), the names
and equation-position indices for the dependent variables (2) and the 
same for the independent variables (3). 

@section implementation Implementation

The method uses STL functionality to parse the textfile and store the 
interpreted equations into a vector of SumOfProductsWithExponents objects. 

Any previously defined system of equations is deleted before the new
system is initialized. 

@section messages Messages

The method tries to analyze the reasons for parsing mistakes and reports
corresponding messages to the user. Common mistakes include missing blank
lines between definition blocks or mispelled 
independent or dependent variable names or the use of reserved '+,-'
notation in variable names. An error is reported also when the number
of differential equations does not correspond to the number of 
dependent variables. 
*/
void ODE_StiffSolver::ReadODEsFrom( const char* file, 
                                    vector<pair<int32_t,string> >& coefs, 
                                    vector<pair<int32_t,string> >& dependent_comps,
                                    vector<pair<int32_t,string> >& independent_comps )
 {
    char    c, text_line[5000];
    double  exp, K=std::numeric_limits<double>::quiet_NaN();
    int32_t coeff_idx{UNSPECIFIED}, idx, i, n, eq_counter = 0;
    char    *token, *sub1, 
            temp[5000], temp1[5000], temp2[5000];
    const char* const white_delims =" ,\t,\n,\r";
    const char* const lhs_delim    ="=";
    const char* const sum_delim    = "+";
    const char* const coeff_delim  =" ";
    bool  is_constraint{false}, negative_term{false};
         
    map<int32_t,pair<double,double> >  product; 
	  SumOfProductsWithExponents         ode;
	  list<string>                       tokens, tokens1;
	  list<string>::iterator             tik, tik1;
    SumOfProducts                      algebraic_eqn;

    strcpy( text_line, file );
    ifstream  ifs( text_line );
    assert( ifs.is_open() );
       
    // --------------------------------------------------------------------------------------
    // 1. reading headline, species names, and independent values later used in the equations 
    // --------------------------------------------------------------------------------------
    ifs.getline( text_line, 1000 ); // title line
    cout <<"\nODE_StiffSolver::ReadODEsFrom: reading: "<< text_line << endl;
    ifs.getline( text_line, 1000 ); // line of ========
    ifs.getline( text_line, 1000 ); // empty line
    
    // 1.1 independent species = constraints (are also added to first list)
    //     (always have negative species indices since they are not part of
    //      the system of equations).
    // --------------------------------------------------------------------
    ifs.getline( text_line, 1000 ); // title of definition of species line
    ifs.getline( text_line, 1000 ); // line of --------
    ifs.getline( text_line, 1000 ); // SPECIES definition
    map<string,pair<int,double> >  indep_species;
    string  spec_name;
    double      constraint;
    bool        empty_line = true;
  
    n = 0;
    while ( (c=text_line[n++]) != '\0' ) if ( !isspace(c) ) { empty_line = false; break; }

    if ( !empty_line )
      {   
         TokenizeString( text_line, white_delims, tokens );
         for ( n=-1, tik=tokens.begin(); tik!=tokens.end(); tik++ ) 
           {
              strcpy( temp1, (*tik).c_str() );
              // substring is broken into name and value part
              spec_name  = strtok( temp1, lhs_delim );
              constraint = atof( strtok( NULL, lhs_delim ) );
              // listing 
              indep_species[ spec_name ] = make_pair( n--, constraint );
           }
         ifs.getline( text_line, 1000 ); // get next blank line
      }

    // 1.2 dependent species 
    // ---------------------
    ifs.getline( text_line, 1000 ); // title of definition of species line
    ifs.getline( text_line, 1000 ); // line of --------
    ifs.getline( text_line, 1000 ); // SPECIES definition
    map<string,int>  dep_species;
    n = 0; // species names must not contain whitespace
    token = strtok( text_line, white_delims ); 
    do {
          dep_species[ string(token) ] = n++;
       }
    while ( (token=strtok(NULL, white_delims)) != NULL );
    ifs.getline( text_line, 1000 ); // get next blank line
   
    
    // 1.3 generating a list of all involved species
    // ---------------------------------------------
    map<string,int32_t>::const_iterator               sit;
    map<string,pair<int32_t,double> >::const_iterator  cit;
    map<string,int32_t>                               species;
    for ( n=0, sit=dep_species.begin(); sit!=dep_species.end();   sit++ ) species[ (*sit).first ] = n++;
    for (    cit=indep_species.begin(); cit!=indep_species.end(); cit++ ) species[ (*cit).first ] = n++;
   
    // -------------------------------------------------------------------------
    // 2. reading equilibrium constants for reactions (same sequence as species) 
    // -------------------------------------------------------------------------
    // spec.idx
    map<string,pair<int,double> >  coeffs;
    string  kname;
    ifs.getline( text_line, 1000 ); // title line
    ifs.getline( text_line, 1000 ); // line of -------
    ifs.getline( text_line, 1000 ); // COEFFICIENT definition
    TokenizeString( text_line, white_delims, tokens );
    for ( n=0, tik=tokens.begin(); tik!=tokens.end(); tik++ ) 
      {
          strcpy( temp1, (*tik).c_str() );
          // now the substrings must be broken into names and values
          // name
          sub1  = strtok( temp1, lhs_delim );
          kname = sub1;
          // constant
          sub1  = strtok( NULL, lhs_delim );
          K     = atof( sub1 );
          // recording
          coeffs[ kname ] = make_pair(n++,K);
       }
    ifs.getline( text_line, 1000 ); // get next blank line
   
    // ----------------------------------------------------------------------------
    // 3. reading ordinary differential equations using the terminology established
    //    above (if expressions do not comply errors are raised)
    // ----------------------------------------------------------------------------
    map<string,pair<int32_t,double> >::const_iterator  vit; // constraint species
    map<string,int32_t>::const_iterator               current_lhs;
    ifs.getline( text_line, 1000 ); // title line
    ifs.getline( text_line, 1000 ); // line of -------
    // removing previously defined reactions
    odes.erase( odes.begin(), odes.end() );

    while ( !ifs.eof() )
      {
          // 3.1 reading ODE lhs
          // -------------------
          ifs.getline( text_line, 1000 ); 
          token = strtok( text_line, lhs_delim ); // lhs
          if ( token == NULL ) break; // if empty line is reached
          strcpy( temp1, token );
          // terminating string where whitespace starts
          sub1=strchr( temp1,' '); *sub1 = '\0';
          // making sure that lhs is listed in order specified in species list (1)
          if ( (current_lhs=dep_species.find(string(temp1))) != dep_species.end() )
            {
               // checking ordering of equation list in file which must correspond to
               // listing of species above
               if ( eq_counter != (*current_lhs).second )
                 {
                     cout <<"\nODE_StiffSolver::ReadODEsFrom: ODEs are not listed in the order. ";
                     cout <<"of definitions of lefthandside terms given above." << endl;
                     ifs.close();
                     assert( correct_interpretation );
                 }
            }
          else
            {
                cout <<"\nODE_StiffSolver::ReadODEsFrom: ODEs lefthandside term must be. ";
                cout <<"ODE identifier given above. Unknown term: ";
                cout << temp1 << endl;
                ifs.close();
                assert( correct_interpretation );
            }

          // 3.2 tokenizing righthandside ODE into products (separator + )
          // -------------------------------------------------------------
          token = strtok( NULL, lhs_delim ); // rhs
          // saving rhs into a list  
          TokenizeString( token, sum_delim, tokens );
          odes.reserve( odes.size()+1 );

          // 3.3 PRODUCTS LOOP
          // -----------------
          for ( tik=tokens.begin(); tik!=tokens.end(); tik++ )  // for each product in 'rhs'
            {
               // 3.2.1 create subtokens corresponding to terms that make up each product
               // -----------------------------------------------------------------------
               strcpy( temp, (*tik).c_str() );
               TokenizeString( temp, coeff_delim, tokens1 );

               // 3.2.2 reading product coefficients
               // ----------------------------------
               tik1 = tokens1.begin();
               strcpy( temp1, (*tik1).c_str() );

               // test for plain coefficient
               if ( (cit=coeffs.find(string(temp1))) != coeffs.end() )
                 {
                    // getting coefficient (=equilibrium constant)
                    coeff_idx = (*cit).second.first;
                    K         = (*cit).second.second;
                 }
               else  // coefficients may be preceded by a minus which must be removed 
                 {   // -------------------------------------------------------------
                     for ( i=1; i<=strlen(temp1); i++ ) temp[i-1] = temp1[i];
                     if ( (cit=coeffs.find(string(temp))) != coeffs.end() ) 
                       {
                           coeff_idx = (*cit).second.first;
                           K         = -1.0 * (*cit).second.second;
                       }
                     else
                      {
                          cout <<"\nError when parsing coefficient in RHS of ODE. ";
                          cout <<"Unknown coefficient: "<< temp1 << endl;
                          ifs.close();
                          assert( correct_interpretation );
                      }
                 }
               tik1++;

               // 3.2.3 reading products (term by term including exponents)
               // ---------------------------------------------------------
               while ( tik1 != tokens1.end() )  // for each product term
                 {
                     is_constraint = false;
                     strcpy( temp1, (*tik1).c_str() );
                     
                     // parsing product term
                     if ( (sit=species.find(string(temp1))) != species.end() )
                        {
                           // checking whether independent variable must be substituted
                           if ( (vit=indep_species.find(string(temp1))) != indep_species.end() ) is_constraint = true;
                           else
                             {
                                sit=dep_species.find(string(temp1));
                                assert( sit != dep_species.end() );
                                idx = (*sit).second;
                             }
                           exp = 1.0;
                        }
                     // if product term is not found in the species list, 
                     // term is checked for whether it is preceded by a minus sign
                     // and/or whether it contains an exponent which needs to be separated
                     else
                       {
                           // trying to parse product term without minus sign if it works the sign is removed
                           if ( temp1[0] == '-' )
                             {
                                negative_term = true;
                                for ( i=1; i<=strlen(temp1); i++ ) temp[i-1] = temp1[i];
                                temp[i-1] = '\0';
                                strcpy( temp1, temp );
                             }
                           if ( strchr( temp1,'^' ) == NULL )
                             {
                                cout <<"\nError: Could not parse product term: "<< (*tik).c_str();
                                cout <<" in ODE. Unknown term: "<< temp1 << endl;
                                // checking for accidential '0' instead of 'O'
                                if ( strchr(temp1,'0') != NULL )
                                  cout <<"Term contains zero(=0) character. Is this intended ?"<< endl;
                                ifs.close();
                                assert( correct_interpretation );
                             }
                           else 
                             {
                                n=0;
                                while ( (c=temp1[n]) != '^' ) temp2[n++] = c;
                                temp2[n] = '\0';
                                if ( (sit=species.find(string(temp2))) == species.end() ) 
                                  {
                                     cout <<"\nError: Could not parse product term: "<< (*tik).c_str();
                                     cout <<" in ODE. Unknown term: "<< temp2 << endl;
                                     // checking for accidential '0' instead of 'O'
                                     if ( strchr(temp2,'0') != NULL )
                                       cout <<"Term contains zero(=0) character. Is this intended ?"<< endl;
                                     ifs.close();
                                     assert( correct_interpretation );
                                  }
                                 // checking whether independent variable must be substituted
                                 if ( (vit=indep_species.find(string(temp2))) != indep_species.end() )
                                   is_constraint = true;
                                 else
                                   {
                                       sit=dep_species.find(string(temp2));
                                       assert( sit != dep_species.end() );
                                       // reading dependent variable information
                                       idx = (*sit).second;
                                   }
                                 // parsing exponent (omitting ^ character)
                                 i=0;
                                 while ( (c=temp1[++n]) != '\0' ) temp2[i++] = c;
                                 temp2[i] = '\0';
                                 exp = atof( temp2 );
                              }
                        }
                     // if product term involves a dependent variable, a new term is added to product map
                     if ( !is_constraint ) 
                       {
                          if ( negative_term ) product[ idx ] = make_pair(-1.0,exp);
                          else                 product[ idx ] = make_pair(1.0,exp);
                       }
                     // otherwise constraint is multiplied into coefficient
                     else 
                       {
                          assert( vit != indep_species.end() );
                          if ( negative_term ) 
                             //             negative index                              const value
                             product[ (*vit).second.first ] = make_pair( (*vit).second.second * -1.0, exp ); 
                          else
                             product[ (*vit).second.first ] = make_pair( (*vit).second.second, exp ); 
                       }
                     negative_term = false;
                     tik1++;
                 }
               // 3.2.4 recording product terms
               // -----------------------------
               ode.AddProduct( coeff_idx, K, product );
               product.erase( product.begin(), product.end() );

            } // end for each product separated by " "

	      // -------------------------------
      	  // 4. Adding ODE to system of ODEs
	      // -------------------------------
	      odes.push_back( ode );
	      ode.Erase();
	      eq_counter++;
	      
      } // end while eof

    // ----------------------------------------------------------------------------
    // 5. reading independent algebraic equations using the terminology established
    //    above. For each independent variable one algebraic equation is expected.
    //    These equations are expected in the order of the independent species listing
    //    from above.
    // ----------------------------------------------------------------------------
    ifs.getline( text_line, 1000 ); // title line
    ifs.getline( text_line, 1000 ); // line of -------
    eqns.erase( eqns.begin(), eqns.end() );
    eq_counter = -1;

    while ( !ifs.eof() )
      {
          // 5.1 reading lhs
          // -------------------
          ifs.getline( text_line, 1000 ); 
          token = strtok( text_line, lhs_delim ); // lhs
          if ( token == NULL ) break; // if empty line is reached
          strcpy( temp1, token );
          // terminating string where whitespace starts
          sub1=strchr( temp1,' '); *sub1 = '\0';
          // making sure that lhs is listed in order specified in species list (1)
          if ( (vit=indep_species.find(string(temp1))) != indep_species.end() )
            {
               // checking ordering of equation list in file which must correspond to
               // listing of species above
               if ( eq_counter != (*vit).second.first )
                 {
                     cout <<"\nODE_StiffSolver::ReadODEsFrom: mass balances are not listed in the order ";
                     cout <<"of definitions of lefthandside terms given above." << endl;
                     ifs.close();
                     assert( correct_interpretation );
                 }
            }
          else
            {
                cout <<"\nODE_StiffSolver::ReadODEsFrom: algebraic equations lefthandside term must be. ";
                cout <<"independent term identifier given above. Unknown term: ";
                cout << temp1 << endl;
                ifs.close();
                assert( correct_interpretation );
            }
          // reading algebraic equations (rhs)
          // ---------------------------------
          token = strtok( NULL, lhs_delim ); // rhs
          algebraic_eqn.InitializeFrom( token, indep_species, dep_species ); 
          eqns.reserve( eqns.size()+1 );
          eqns.push_back( algebraic_eqn );
          algebraic_eqn.Erase();

          eq_counter--;          

       } // end while !eof
      
    ifs.close();

    // --------------------------------
    // 5. Generating vectors for return
    // --------------------------------
    coefs.erase( coefs.begin(), coefs.end() );
    dependent_comps.erase( dependent_comps.begin(), dependent_comps.end() ); 
    independent_comps.erase( independent_comps.begin(), independent_comps.end() );
    coefs.reserve( coeffs.size() );
    dependent_comps.reserve( dep_species.size() );
    independent_comps.reserve( indep_species.size() );
    for ( cit=coeffs.begin(); cit!=coeffs.end(); cit++ ) 
      coefs.push_back( make_pair((*cit).second.first,(*cit).first) );
    for ( sit=dep_species.begin(); sit!=dep_species.end(); sit++ ) 
      dependent_comps.push_back( make_pair((*sit).second,(*sit).first) );
    for ( cit=indep_species.begin(); cit!=indep_species.end(); cit++ ) 
      independent_comps.push_back( make_pair((*cit).second.first,(*cit).first) );
   
    assert( odes.size() == dep_species.size() );
    
    cout <<"\nODE_StiffSolver::ReadODEsFrom: file '"<< file;
    cout <<"' read successfully..." << endl;
    
    // ------------------------------------------------
    // 6. setting solver up for new system of equations
    // ------------------------------------------------
    FreeSolverMemory();
    dof = odes.size();
    AllocateSolverMemory();
    
} // end 




/**
 
The method uses position integers and new values to change 
independent product terms in a system of differential
equations which was already defined earlier. 

@param ind An STL vector of integer - string pairs defines the name of
the independent variables/coefficient, its position index in the equation
and its numerical value. 

@section implementation Implementation

The method first calls ResetIndependentProductTerms() and then re-assigns
new values using STL vector assignment. 

@section application Application

To change independent variables in equations without redefining the 
system of equations. This may apply for instance if PT-dependent 
equilibrium constants are used in a model in which the temperature varies. 
*/
void ODE_StiffSolver::AssignIndependentProductTerms( vector<pair<int,string> >& ind, 
                                                     vector<double>&                vals )
 {
    ResetIndependentProductTerms( ind, vals );

    independent_terms  = ind;
    independent_values = vals;

 } // end AssignIndependentProductTerms






/**
 
The need to change independent equation coefficients (such as activity 
coefficients) as the values of dependent variables change in an evolution 
of a system is accounted for by this method. 

@param ind The method takes vector arguments which define the new product coefficients.
The individual products in the equations are identified by integers
(e.g. product 1, product 2, etc.) and negative integers denote 
independent product coefficients (those that are not calculated by ODEs).
The values correspond to the exponents of the terms. 

@section implementation Implementation

The SumOfProductWithExpondents is defined as pairs of c.idx/coefficient 
and e.idx/exponents. The method performs on all independent variables 
in the ODEs, on every ODE, in every product which the ODE contains, if the 
tested independent variable is part of this product. Also, for all independent 
variables in algebraic equations, in every algebraic equation, in every 
product which the equation contains, if the independent variable is part 
of this product, the coefficient will be modified. The correct index is 
abs(idx-1) since all independent var indices are set as negative (-1...-n). 
The exponents stay untouched. 

@section application Application

The method is called inside of odeint() if independent variables are part
of the system of differential equations. The capability to change 
coefficients while the ODE solver converges to a solution allows to 
modify, for instance, the activities of chemical species as their 
concentrations change over time. 
*/
void ODE_StiffSolver::ResetIndependentProductTerms( vector<pair<int32_t,string> >& ind, 
                                                    vector<double>&  vals )
 {
    vector<pair<int32_t,string> >::const_iterator  it;
    vector<SumOfProductsWithExponents>::iterator     pit;
    vector<sumofproducts_pair>::iterator             sit;
    map<int32_t,pair<double,double> >::iterator  tit;

    // SumOfProductWithExpondents is defined as 
    //    c.idx/coefficient    e.idx/exponent
    // pair<pair<int,double>,map<int,pair<double,double>,less<int> > >  

    // 1. for all independent variables in ODES
    // ----------------------------------------
    for ( it=ind.begin(); it!=ind.end(); it++ )
      // in every ODE
      //-------------
      for ( pit=odes.begin(); pit!=odes.end(); pit++ )
        // in every product which the ODE contains
        //----------------------------------------
        for ( sit=(*pit).Begin(); sit!=(*pit).End(); sit++ )
          // if the independent variable is part of this product
          //----------------------------------------------------
          if ( (tit=(*sit).second.find( (*it).first )) != (*sit).second.end() )
            // correct index is abs(idx-1) since all independent var indices 
            // are negative (-1...-n). The exponent stays untouched
            (*tit).second.first = vals[ static_cast<uint32_t>(abs((*it).first)-1) ];

    // 2. for all independent variables in algebraic equations
    // -------------------------------------------------------
    vector<SumOfProducts>::iterator         oit;
    vector<map<int32_t,double> >::iterator soit;
    map<int32_t,double>::iterator          ot;

    for ( it=ind.begin(); it!=ind.end(); it++ )
      // in every algebraic equation
      //----------------------------
      for ( oit=eqns.begin(); oit!=eqns.end(); oit++ )
        // in every product which the equation contains
        //---------------------------------------------
        for ( soit=(*oit).Begin(); soit!=(*oit).End(); soit++ )
          // if the independent variable is part of this product
          //----------------------------------------------------
          if ( (ot=(*soit).find( (*it).first )) != (*soit).end() )
            // correct index is abs(idx-1) since all independent var indices 
            // are negative (-1...-n). The exponent stays untouched
            (*ot).second = vals[ static_cast<uint32_t>(abs((*it).first)-1) ];

 } // end




/**
 
A change of the coefficients in the system of equations, for instance by
ResetIndependentProductTerms() necessitates a reavaluation of the 
product terms as done in derivs(). EvaluateIndependentProductTerms()
was written for this purpose and is called inside odeint().  

@section arguments Input Arguments

The methods requires the values of the dependent variables as first 
argument. 

@param vals the new derivatives (righthand sides of the ODEs) are returned into the
STL vector which is supplied as second argument. 

@section implementation Implementation

The righthand side values are calculated using the Sum() method of 
the SumOfProducts class object. The terms are added in the output vector 
at the locations 1...n. 

@section application Application

The method is called inside of odeint(), after ResetIndependentProductTerms(),
if and only if independent parameters were modified between solution steps. 

@section messages Messages

Currently, the method reports the calculated values for each product term. 
*/
//                                                                  from 1...n               from 0...n-1
void ODE_StiffSolver::EvaluateIndependentProductTerms( const double dependent_vals[], vector<double>& vals )
 {
    assert ( eqns.size() == vals.size() );
    for ( size_t i{0U}; i<eqns.size(); i++ ) 
      {
         vals[i] = eqns[i].Sum( dependent_vals );
         cout <<"\nindependent value: "<< vals[i];
         cout.flush();
      }

 } // end EvaluateIndependentProductTerms

} // end namespace csmp
