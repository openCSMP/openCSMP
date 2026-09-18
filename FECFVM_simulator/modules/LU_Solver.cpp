// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "LU_Solver.h"
#include <cmath>

using namespace std;

namespace csmp {


LU_Solver::LU_Solver()
{
}


LU_Solver::~LU_Solver()
{
}


/** Solves Ax = b with rhs[1]
*/
void LU_Solver::Solve(DenseMatrix<DM_MIN>& A,
                      vector<double>& b,
                      vector<double>& x,
                      uint32_t n)
{

  tiny_ = 1.e-40;
  index_.resize(n);

  ludcmp(A, n, index_);
  lubksb(A, n, index_, x, b);

}

/** LU decomposition
*/
void LU_Solver::ludcmp(DenseMatrix<DM_MIN>& a,
                       uint32_t n,
                       std::vector<uint32_t>& indx)
{
      vector<double> vv(n, 0.);
      int32_t i, j, k;
      double d(1.);
      double big, cache;

      for(i=0; i<n; i++){
          big = 0.;
          for(j=0; j<n; j++)
              if((cache = abs(a( i,j ))) > big) big = cache;
          if(big == 0.) throw("LU decomposition Solver: Singular Matrix!");
          vv[i] = 1./big;
      }

      int32_t imax = numeric_limits<int32_t>::max();

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
              for( j=0; j<n; j++){
                  cache = a(static_cast<uint32_t>(imax), j);
                  a(imax, j) = a(k, j);
                  a(k, j) =  cache;
              }
              d = -d;
              vv[imax] = vv[k];
          }

          indx[k] = imax;

          if(a(k, k) == 0.) a(k, k) = tiny_;

          for(i = k+1; i<n; i++){
              a(i, k) = a(i, k) / a(k, k);
              cache = a(i, k);
              for(j=k+1; j<n; j++){
                  a(i, j) = a(i, j) - cache * a(k, j);
              }
          }
      }
}



/** substitution and establishing of solution vector x
*/
void LU_Solver::lubksb(DenseMatrix<DM_MIN>& a,
                       uint32_t n,
                       std::vector<uint32_t>& indx,
                       std::vector<double>& x,
                       std::vector<double>& b )
{
      int32_t i, j, ip, ii(0);
      double sum;

      for(i=0; i<n; i++)
          x[i] = b[i];
          
      for(i=0; i<n; i++)
        {
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
      
      for( i = n-1; i >= 0; i--)
        {
            sum = x[i];
            for(j=i+1; j<n; j++)
                sum -= a( i,j ) * x[j];
            x[i] = sum / a( i,i );
        }
}


} // csmp

