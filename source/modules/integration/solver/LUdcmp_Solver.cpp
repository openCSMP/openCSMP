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
                                         vector<double>& b,
                                         vector<double>& x,
                                         size_t )
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( A.Rows() > 500U )
       csmp_error.Note( WARNING, "LUdcmp_Solver::SolveMatrixEquation:",
                         "the matrix is rather large; try SAMG to get a result in a decent time." );


    tiny_ = numeric_limits<double>::epsilon();
    n_ = A.Rows();
    index_= std::vector<size_t>( n_ );

    //luout();
    ludcmp( A, n_, index_ );
    lubksb( A, n_, index_, x, b );
  }


void LUdcmp_Solver::SolveMatrixEquation( CompressedRowMatrix& A,
                                         vector<double>& b,
                                         vector<double>& x,
                                         size_t )
  {
      throw csmp::Exception( ERROR, "LUdcmp_Solver::SolveMatrixEquation", "Method not implemented for CompressedRowMatrix yet");
  }




/** LU decomposition
      
    @note this is super slow because it operates directly on the SparseMatrix class, creating and deleting entries!
  
    @todo replace with Eigen / work of a DenseMatrix object with size limited to 10k rows =  columns (100 million entries!)
*/
void LUdcmp_Solver::ludcmp( SparseMatrix& a,
                            long n,
                            std::vector<size_t>& indx)
  {
      vector<double> vv( n, 0.);
      long i, j, k;
      double d = 1.;
      double big, cache;

      for(i=0; i<n; i++){
          big = 0.;
          for(j=0; j<n; j++)
              if((cache = abs(a( i,j ))) > big) big = cache;
          if( big == 0.) throw("LU decomposition Solver: Singular Matrix!");
          vv[i] = 1./big;
      }

  long imax{0};

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
              for( j=0; j<n; j++ ){
                  cache = a( static_cast<uint32_t>(imax), j );
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
                            std::vector<double>& x,
                            std::vector<double>& b )
  {
      long i, j, ip, ii(0);
      double sum;

      for(i=0; i<n; i++)
          x[i] = b[i];
      for(i=0; i<n; i++){
          ip = indx[i];
          sum = x[ip];
          x[ip] = x[i];
          if(ii != 0)
              for( j=ii-1; j<i; j++)
                  sum -= a( i,j ) * x[j];
          else if (sum != 0.)
              ii = i+1;
          x[i] =sum;
      }
      for( i=n-1; i>=0; i--){
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

