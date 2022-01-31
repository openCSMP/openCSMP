#include "Gauss_Solver.h"
#include "SparseMatrix.h"
#include "meschach.h"

#ifdef CSMP_WITH_MESCHACH
#include "iter.h"
#include "sparse2.h"

using namespace std;

namespace csmp {

Gauss_Solver::Gauss_Solver():
    Meschach_Solver(),
    pivot_factor(0.5) // empirically determined
{}

Gauss_Solver::~Gauss_Solver()
{}

void Gauss_Solver::SolveWithMeschach( SparseMatrix& A,
                                      vector<double>& b,
                                      vector<double>& x,
                                      double )
{
    if (Verbose()) {
        cout <<"\nGauss_Solver::SolveMatrixEquation: Allocating memory..." << endl;
    }
    // build solution matrix with 2 elements per row
    SPMAT*  M   = sp_get( static_cast<int32_t>(A.Rows()), static_cast<int32_t>(A.Cols()), 3 );
    convert_CSP_SparseMatrix_to_SPMAT( A, M );
    SPMAT*  LU  = sp_copy( M );
    VEC*    rhs = v_get(b.size());
    VEC*    sol = v_get(x.size());
    
    for ( uint32_t i=0; i<x.size(); i++ ) {
        rhs->ve[i] = b[i];
        sol->ve[i] = x[i];
    }
    
    // Apply Meschach's direct solver, solution returned in x
    if (Verbose()) {
        cout <<"\nGauss_Solver::SolveMatrixEquation: Linear solution in progress..." << endl;
    }
    PERM*  pivot = px_get(x.size());
    spLUfactor( LU , pivot, pivot_factor );
    spLUsolve ( LU , pivot, rhs, sol );
    
    for ( uint32_t i=0; i<x.size(); i++ ) x[i] = sol->ve[i];
    
    // giving the memory back to the system
    sp_free( M );
    sp_free( LU );
    PX_FREE( pivot );
    V_FREE( rhs );
    V_FREE( sol );
}

} // end namespace csmp

#endif
