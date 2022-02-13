#include "GaussJordan_Solver.h"
#include "Exception.h"
#include <cmath>

using namespace std;

namespace csmp {


GaussJordan_Solver::GaussJordan_Solver(){
}


GaussJordan_Solver::~GaussJordan_Solver(){
}


/** Function to solve Ax = b with single column RHS
*/
void GaussJordan_Solver::SolveMatrixEquation( SparseMatrix& A,
                                              vector<double>& b,
                                              vector<double>& x,
                                              size_t )
{
    (*this).GaussJordan( A, b );
    x = b;
}

void GaussJordan_Solver::SolveMatrixEquation( CompressedRowMatrix& A,
                                              vector<double>& b,
                                              vector<double>& x,
                                              size_t )
{
  throw csmp::Exception( ERROR, "GaussJordan_Solver::SolveMatrixEquation", "Method not implemented for CompressedSparseRowMatrix yet" );
}



/** Gauss-Jordan Solver

Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. The input matrix
is a[0..n-1][0..n-1]. b is input containing the m right-hand side.
On output, a is replaced by its matrix inverse, and b is replaced by the corresponding solution
vector.
@param A the Matrix A of the linear system Ax=b in nxn
@param b the RHS vector b of the linear system Ax=b in n
*/
void GaussJordan_Solver::GaussJordan( SparseMatrix& A, vector<double>& b )
 {
    size_t          i, icol(0), irow(0), j, k, l, ll;
    const size_t    n = A.Rows();
    double       big, dum, pivinv;
    vector<uint32_t>  indxc(n,0), indxr(n,0), ipiv(n,0);

    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
//        cout<<" i : "<<i<<endl;
        big = 0.;
        for (j=0;j<n;j++)
            if (ipiv[j] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (abs(A( j,k )) >= big) {
                            big=abs(A( j,k ));
                            irow=j;
                            icol=k;
                        }
                    }

                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=0;l<n;l++) SwapSparseMatrixElements( A, irow, l, icol, l );
            swap( b[irow], b[icol] );
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if ( A( icol, icol ) == 0.0) throw("gaussj: Singular Matrix");
        pivinv = 1.0 / A( icol, icol );
        A.Assign( icol, icol, 1. );
        for (l=0;l<n;l++) A.Assign( icol,l, A(icol,l) * pivinv);
        b[icol] *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll != icol) {
                dum=A(ll,icol);
                A.Assign(ll,icol, 0.);
                for (l=0;l<n;l++) A.Assign(ll,l, A(ll,l)-A(icol,l)*dum);
                b[ll] -= b[icol]*dum;
            }
    }

    for (long l=n-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
                SwapSparseMatrixElements(A, k, indxr[l], k, indxc[l]);
    }

} // end



/** Swaps elements of a csmp::SparseMatrix

function to swap elements of a csmp::SparseMatrix

@section motivation Motivation

Helper class to the GaussJordan Solver

@param M the Matrix of which elements are to be swapt
@param row1 the row of element 1
@param col1 the column of element 1
@param row2 the row of element 2
@param col2 the column of element 2
*/
void GaussJordan_Solver::SwapSparseMatrixElements( SparseMatrix& M,
                                                   long row1,
                                                   long col1,
                                                   long row2,
                                                   long col2 ) const
{
    double cache = M( row1, col1);
    M.Assign(row1,col1, M(row2,col2));
    M.Assign(row2,col2, cache);
}

} // csmp
