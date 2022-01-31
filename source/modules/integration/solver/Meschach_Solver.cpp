#include "Meschach_Solver.h"
#include "Solver.h"
#include "SparseMatrix.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"

#ifdef CSMP_WITH_MESCHACH
#include "iter.h"
#include "sparse2.h"


using namespace std;

namespace csmp {

Meschach_Solver::Meschach_Solver() :
    Solver(),
    residual_factor_(1e-15)    // signal to noise ratio
{}

Meschach_Solver::~Meschach_Solver(){}

void Meschach_Solver::SolveMatrixEquation(SparseMatrix& A,
                                          std::vector<double>& b,
                                          std::vector<double>& x,
                                          size_t )
{
    ErrorHandler&  skm_err( ErrorHandler::Instance() );
    double typical = GuessResidual( A, b, x );
    
    // setting the 'tolerance' for the matrix inversion
    double meschach_tolerance = fmin(typical * residual_factor_, residual_factor_);
    
    if (Verbose()) {
        cout << "\n\tUser aspired precision of solution: " << meschach_tolerance << endl;
    }
    
    SolveWithMeschach(A, b, x, meschach_tolerance);
    
    double residual = CalculateResidual( A, b, x );

    if (Verbose()) {
        cout << "\nInferred initial residual   = " << typical << endl;
        cout << "Residual after solution   = " << residual << endl;
        cout << "\n\tActually achieved precision of solution: " << residual/typical << endl;
    }

    if ( residual > typical ) {
        // outputting matrices to text files
        Out( A );
        Out( b );
        skm_err.notice( FATAL_ERROR, "Solver::SolveMatrixEquation",
                        "No convergence occurred..." );
    }
}

/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  double  Solver::GuessResidual( SpMat& A, const vector<double>& b, vector<double>& x )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Tries to calculate a "typical" size of initial residual by randomizing the
solution vector and then calculation the product A*x and comparing it with
the righthand vector b. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->

The global solution matrix, the righthand vector, and the solution vector. <p>

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->

The solution vector is used as temporary storage by this method. It is
therefore modified. <p>

<H4>Application:</H4><!--------------------------------------------------->

To get an initial guess of the problem such that it can later be identified
whether the solver reached some convergence. <p>

<!------------------------------------------------------------------------>
tested: O.K. */
double Meschach_Solver::GuessResidual(const SparseMatrix& A,
                                        const vector<double>& b,
                                        vector<double>& x)
{
    const size_t   len(b.size());
    vector<double> tvec(len,0.);

    random_generator rng;
    vector_randomize( rng, x );
    
    // calculating residual for potential solution
    for ( size_t i=0; i<len; i++ ) {
        for ( size_t j=0; j<len; j++ ) tvec[i] += A(i,j) * x[j];
        tvec[i] = b[i] - tvec[i];
    }
    
    // calculate the L1 norm
    vector<double> scale_vec;
    return vector_norm2( tvec, scale_vec );
}



/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
double  Solver::CalculateResidual( const SparseMatrix& A,
                                      const vector<double>& b,
                                      const vector<double>& x )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Calculates the size of the residual error by comparing the product A*x with
the righthand vector b. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->

The global solution matrix, the righthand vector, and the solution vector. 
<p>

<H4>Application:</H4><!--------------------------------------------------->

To calculate the residual of a solution to see whether the solver reached 
some convergence and how good the solution actually is. <p>

<!------------------------------------------------------------------------>
tested: O.K. */
double Meschach_Solver::CalculateResidual(const SparseMatrix& A,
                                            const vector<double>& b,
                                            const vector<double>& x)
{
    if ( b.size() != x.size() ) {
        cerr <<"\nSolver::CalculateResidual: Size of solution vector 'x' differs from ";
        cerr <<" righthand 'b' vector (x vs. b): "<< x.size() <<" vs. "<< b.size() << endl;
        throw range_error("Solver::CalculateResidual");
    }

    const size_t   len(x.size());
    vector<double> tvec(len,0.0);
    
    // calculating residual for actual solution
    for ( size_t i=0; i<len; i++ ) {
        for ( size_t j=0; j<len; j++ ) tvec[i] += A(i,j) * x[j];
        tvec[i] = b[i] - tvec[i];
    }
    
    // calculate the L1 norm
    vector<double> scale_vec;
    return vector_norm2( tvec, scale_vec );
}

}

#endif
