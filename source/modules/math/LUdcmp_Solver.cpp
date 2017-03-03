#include "LUdcmp_Solver.h"
#include "ErrorHandler.h"


using namespace std;

namespace csmp {


LUdcmp_Solver::LUdcmp_Solver()
  : tiny_(numeric_limits<double>::epsilon())
  {
  }


LUdcmp_Solver::~LUdcmp_Solver()
  {
  }


/// Solves Ax = b with rhs[1]
void LUdcmp_Solver::SolveMatrixEquation( SparseMatrix& A,
                                         vector<double64>& b,
                                         vector<double64>& x,
                                         size_t )
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( A.Rows() > 500U )
       csmp_error.notice( CSMP_WARNING, "LUdcmp_Solver::SolveMatrixEquation:",
                         "the matrix is rather large; try SAMG to get a result in a decent time." );


    tiny_ = numeric_limits<double>::epsilon();
    n_ = A.Rows();
    index_= std::vector<size_t>( n_ );

    //luout();
    ludcmp( A, n_, index_ );
    lubksb( A, n_, index_, x, b );
  }


/// LU decomposition
void LUdcmp_Solver::ludcmp( SparseMatrix& a,
                            long n,
                            std::vector<size_t>& indx)
  {
      vector<double64> vv( n, 0.);
      long i, j, k;
      double64 d = 1.;
      double64 big, cache;

      for(i=0; i<n; i++){
          big = 0.;
          for(j=0; j<n; j++)
              if((cache = abs(a( i,j ))) > big) big = cache;
          if( big == 0.) {
              throw csmp::Exception(CSMP_ERROR,
                      "LUdcmp_Solver::ludcmp", "LU decomposition Solver: Singular Matrix!");
          }

          vv[i] = 1./big;
      }

      long imax;

      for(k=0; k<n; k++){
          big = 0.;
          for(i=k; i<n; i++){
              cache = vv[i] * abs( a( i,k ));
              if( cache > big){
                  big = cache;
                  imax = i;
              }
          }
          if( k != imax ){
              for(size_t j=0; j<n; j++){
                  cache = a( static_cast<size_t>(imax), j );
                  a.Assign( imax,j, a( k,j ));
                  a.Assign( k,j, cache);
              }
              d = -d;
              vv[imax] = vv[k];
          }

          indx[k] = imax;

          if( a( k,k ) == 0. ) a.Assign( k,k, tiny_ );

          for(i = k+1; i<n; i++){
              a.Assign( i,k, a( i,k ) / a( k,k ) );
              cache = a ( i,k );
              for(j=k+1; j<n; j++){
                  a.Assign( i, j, a( i,j ) - cache * a( k, j ) );
              }
          }
      }
  }
  
  
  
  

/// substitution and establishing of solution vector x
void LUdcmp_Solver::lubksb( SparseMatrix& a,
                            long n,
                            std::vector<size_t>& indx,
                            std::vector<double64>& x,
                            std::vector<double64>& b )
  {
      long i, j, ip, ii(0);
      double64 sum;

      for(i=0; i<n; i++)
          x[i] = b[i];
      for(i=0; i<n; i++){
          ip = indx[i];
          sum = x[ip];
          x[ip] = x[i];
          if(ii != 0)
              for(size_t j=ii-1; j<i; j++)
                  sum -= a( i,j ) * x[j];
          else if (sum != 0.)
              ii = i+1;
          x[i] =sum;
      }
      for( long i=n-1; i>=0; i--){
          sum = x[i];
          for(j=i+1; j<n; j++)
              sum -= a( i,j ) * x[j];
          x[i] = sum / a( i,i );
      }
  }





void LUdcmp_Solver::luout()
  {
    cout << "\n\nLU decomposition, solving system of ";
    cout << n_ << " degrees of freedom\n\n";
  }

} // csmp

