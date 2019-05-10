#include "Solver.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include "ModelTime.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"
#include "SAMG_Solver.h"

using namespace std;

namespace csmp {

/**
 
Default constructor chooses the algebraic multigrid method as solution
method, sets the pivot factor to 0.5, and the residual multiplication
factor to 1.0e-15. Verbose is disabled and therefore no solution
matrices or vectors will be printed.  

*/
Solver::Solver()
 : verbose_(false),
   solver_settings_(nullptr)
 {
 }

Solver::Solver( SolverSettings* settings )
 : verbose_(false) ,
   solver_settings_(settings)
 {
 }

 
Solver::~Solver()
 {
 }
 


/**
 
When verbose is switched on, the Solver will print the global solution
matrix, the solution vector, and the righthand side to stdout. 
*/ 
void Solver::Verbose( bool verb ) { verbose_ = verb; }

bool Solver::Verbose() const { return verbose_; }



/**
 
Calls solution methods from the Meschach++ and algebraic multigrid
libraries to invert the global solution matrix. If the solution process
fails, the solution matrix and the righthand vector are output to
textfiles whose name contains the runtime in seconds and which 
have the extension '.text'. 

@param A The sparse solution matrix A, 
@param b the righthand vector 'b' with boundary conditions,
@param x the solution vector.

After successful matrix inversion, the solution vector is returned into 'x'.

@section implementation Implementation 

SolveMatrixEquation() first tries to calculate a "typical" size of the
residual (the difference between the exact and the numerical solution).
Secondly, based on this residual it sets the 'tolerance' for the matrix 
inversion. This involves an initial guess with correct boundary 
conditions and assumes that the righthand vector b has been setup with 
Dirichlet boundary conditions and the solution matrix with corresponding 
1's on its diagonal. 

If the matrix inversion fails, SolveMatrixEquation() calls the method
Out() on the solution matrix and the righthand vector. 

@section application Application 

Everytime an Algorithm is passed to the Model, its computation 
involves the solution of a matrix equation. 

@section messages Messages 

If verbose is enabled the method will print out the solution matrix,
the solution vector, and the righthand vector. Also it will report
a gues of the initial residual and the actual residual of the 
obtained solution: 

@code
Typical initial Residue   = 5.21798
Residual after solution   = 3.58661e-16
@endcode

Important is that this output gives the ratio between the typical initial 
residue and the residual after the solution (here 1.0e+16). This ratio 
is essentially the signal-to-noise ratio and must (as a minimum for 
the solution to be meaningful), be greater than the ratio 
between the largest and the smallest dependent-variable value which
was specified as initial or as boundary condition.  

When the AMG solver is used, the residual relative to the initial residual
is given inside of the AMG output (see above). It is found as the L2
vector norm.  
*/ 
void  Solver::Solve( SparseMatrix& A, 
                     std::vector<double64>& b,
                     std::vector<double64>& x,
                     size_t no_unknowns )
 {
    // checking for consistent sizes 
    if ( A.Rows() != b.size() )
      throw csmp::Exception( ERROR, "Solver::Solve", "sparse matrix G cols is not equal to RHS rows");

    x.resize( b.size() );
    vector<double64>(x).swap(x);

    if ( Verbose() ) {
      if ( A.Symmetric() )
        cout <<"\nSolver::SolveMatrixEquation: Solution matrix is symmetric."<< endl;

      cout <<"\nSolver::SolveMatrixEquation: Global matrix before solution:"<< endl;
      A.Out();
      cout <<"\nSolver::SolveMatrixEquation: Righthand vector:"<< endl;
      out( b );
    }
    
    // Setup initial guess with correct BC, This assumes that b has been setup with Dirichlet
    // BC and A with corresponding 1's on diagonal
    // * Will prevent using x as a first guess
    // x = b;
    
    // Delegate solution process (this is a purely virtual function)
	
	// cout << "\nSolving Ax=b with GuessSidel\n\n"; 
	SolveMatrixEquation(A, b, x, no_unknowns);
                                     
    if ( Verbose() ) {
         cout <<"\nSolver::SolveMatrixEquation: Global matrix after solution:"<< endl;
         A.Out();
         cout <<"\nSolver::SolveMatrixEquation: Righthand vector:"<< endl;
         out( b );
         cout <<"\nSolver::SolveMatrixEquation: Solution vector:"<< endl;
         out( x );
      }
                                    
 } // end SolveMatrixEquation

/*

*/
void  Solver::Solve( CompressedRowMatrix& A,
                     std::vector<double64>& b,
                     std::vector<double64>& x,
                     size_t no_unknowns )
{
  // checking for consistent sizes 
  if ( A.Rows() != b.size() )
    throw csmp::Exception( ERROR, "Solver::Solve", "sparse matrix G cols is not equal to RHS rows" );

  x.resize( b.size() );
  vector<double64>( x ).swap( x );

  if ( Verbose() ) {    
    cout << "\nSolver::SolveMatrixEquation: Global matrix before solution:" << endl;
    A.Out();
    cout << "\nSolver::SolveMatrixEquation: Righthand vector:" << endl;
    out( b );
  }

  // Setup initial guess with correct BC, This assumes that b has been setup with Dirichlet
  // BC and A with corresponding 1's on diagonal
  // * Will prevent using x as a first guess
  // x = b;

  // Delegate solution process (this is a purely virtual function)

  // cout << "\nSolving Ax=b with GuessSidel\n\n"; 
  SolveMatrixEquation( A, b, x, no_unknowns );

  if ( Verbose() ) {
    cout << "\nSolver::SolveMatrixEquation: Global matrix after solution:" << endl;
    A.Out();
    cout << "\nSolver::SolveMatrixEquation: Righthand vector:" << endl;
    out( b );
    cout << "\nSolver::SolveMatrixEquation: Solution vector:" << endl;
    out( x );
  }

} // end SolveMatrixEquation

/**
 
Outputs the solution matrix to a textfile in SparseMatrix output
format. 

@section arguments Input Arguments

The first argument is a constant reference to the solution matrix and the
second argument is the name of the file which the matrix will be output to.
To this file name, the global runtime of the model and the extension '.text' 
will be appended.  

@section implementation Implementation 

Out() calls the Meschach function sp_foutput(). 

@section application Application 

To independently test the solution matrix or to transfer it to another
solver. 

@section messages Messages 

A successful write will be reported. 
*/
void  Solver::Out( const SparseMatrix& mat, const char* fname ) const
 {
    double64& model_time( ModelTime::Instance().modelTime );
    char   file[200], gtime[30];

    // appending the global runtime to the file name
    // ---------------------------------------------
    sprintf( gtime, "%lf", model_time );
    strcpy( file, fname );
    strcat( file, gtime );
    strcat( file,".text" );

    mat.Out( file );

    cout <<"\nSolver:Out: file '"<< file <<"' written successfully."<< endl;

 } // end Out


/**
 
Outputs the Meschach++ vector 'vec' to a textfile in Meschach++ output
format. 

@section arguments Input Arguments

The first argument is a constant reference to the Meschach vector and the
second argument is the name of the file which the vector will be output to.
To this file name, the global runtime of the model and the extension '.text' 
will be appended.  

@section implementation Implementation 

Out() calls the Meschach function v_foutput() on a standard C filestream.
 

@section application Application 

To independently examine a solution or righthand vector or to transfer it 
to another solver. 

@section messages Messages 

A failure to open the output file or a successful write will be reported. 
*/
void  Solver::Out( const vector<double64>& vec, const char* fname ) const
 {
    double64& model_time( ModelTime::Instance().modelTime );
    FILE  *fp;
    char   file[200], gtime[30];

    // appending the global runtime to the file name
    // ---------------------------------------------
    sprintf( gtime, "%lf", model_time );
    strcpy( file, fname );
    strcat( file, gtime );
    strcat( file,".text" );

    if ((fp = fopen (file,"wt")) == NULL )
      throw csmp::Exception( ERROR, "Solver::Out (vector<double64>)", "output file could not be created" );

    out( vec );

    fclose( fp );
    cout <<"\nSolver:Out: file '"<< file <<"' written successfully."<< endl;

 } // end Out vector<double64>   




/**
 
Calculates the size of the residual error by comparing the product A*x with
the righthand vector b.  

@section arguments Input Arguments

The global solution matrix, the righthand vector, and the solution vector. 
 

@section application Application 

To calculate the residual of a solution to see whether the solver reached 
some convergence and how good the solution actually is.  
*/
double64 Solver::CalculateResidual( const SparseMatrix& A, 
                                      const vector<double64>& b, 
                                      const vector<double64>& x ) const
 {
    if ( b.size() != x.size() ) {
         cout <<"\nSolver::CalculateResidual: Size of solution vector 'x' differs from ";
         cout <<" righthand 'b' vector (x vs. b): "<< x.size() <<" vs. "<< b.size() << endl;
         throw range_error("Solver::CalculateResidual");
      }
 
    const size_t   len(x.size());
    vector<double64> tvec(len,0.);
    
    // calculating residual for actual solution
    for ( size_t i=0U; i<len; i++ ) {
         for ( size_t j=0; j<len; j++ ) tvec[i] += A(i,j) * x[j];
         tvec[i] = b[i] - tvec[i];
      }
    
    // calculate the L1 norm
    vector<double64> scale_vec;
    return vector_norm2( tvec, scale_vec );
 }

SolverSettings* Solver::GetSolverSettings() {
    return solver_settings_;
}

} // end namespace csmp
