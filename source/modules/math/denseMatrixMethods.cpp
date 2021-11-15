//
//  denseMatrixMethods.cpp
//  Coding exercise with the new Melbourne team.
//
//  Created by Stephan Matthai on 17/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "denseMatrixMethods.h"

namespace csmp {

template<size_t dim>
void LU_Decomposition( DenseMatrix<dim>& A )
 {
    assert( A.Rows() == A.Columns() /* LU decompositions makes sense only for square matrices */ );
    const size_t n(A.Rows());
    for (size_t k=0;k<n-1U;k++)
      {
        for (size_t i=k+1U;i<n;i++) {
              A(i,k)=A(i,k)/A(k,k);
              for (size_t j=k+1U;j<n;j++) {
                  A(i,j)=A(i,j)-A(i,k)*A(k,j);
                }
          }
      }
 } // end LU_Decomposition


template<size_t dim>
void LU_BackSubstitution( const DenseMatrix<dim>& A, std::vector<double>& z )
 {
    assert( A.Rows() == A.Columns() /* LU decompositions makes sense only for square matrices */ );
    z.resize(A.Rows());
    const size_t n(A.Rows());

    for (size_t k=0;k<n;k++){
        for(size_t i=0;i<n;i++){
            if (i==k){z[i]=1;} else {z[i]=0;}
            for(size_t j=0;j<i;j++){
                z[i]=z[i]-A(i,j)*z[j];
            }
            std::cout<<z[i]<< std::endl;
        }

    DenseMatrix<dim> X(A.Rows(),A.Columns(),0.);
    
        for(size_t i=n-1;i+1>0;i--){
            X(i,k)=z[i]/A(i,i);
            for(size_t j=n-1;j>i;j--){
                X(i,k)=X(i,k)-A(i,j)/A(i,i)*X(j,k);
            }
        }

    }
 }


} // end csmp
