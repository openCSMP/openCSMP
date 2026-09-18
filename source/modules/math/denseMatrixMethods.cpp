// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  denseMatrixMethods.cpp
//  Coding exercise with the new Melbourne team.
//
//  Created by Stephan Matthai on 17/08/2015.
//

#include "denseMatrixMethods.h"
#include "CSMP_definitions.h"

namespace csmp {

// =============================================================================
// validateCoordinateMatrix — not performance-critical, defined in source
// =============================================================================

template<typename MatrixType>
void validateCoordinateMatrix( const MatrixType& XY, uint32_t expected_nodes, uint32_t dim )
 {
    const uint32_t r = XY.Rows();
    const uint32_t c = XY.Cols();

    // 1. Structural check
    if ( r < expected_nodes || c < dim )
      throw std::runtime_error( std::format(
        "validateCoordinateMatrix: found ({}x{}), expected at least ({}x{})",
        r, c, expected_nodes, dim ) );

    // 2. Value integrity — NaN/Inf and all-zero check
    double absolute_sum{0.0};
    for ( uint32_t i{0U}; i < expected_nodes; ++i )
      for ( uint32_t j{0U}; j < dim; ++j )
        {
          const double val = XY(i,j);
          if ( !std::isfinite(val) )
            throw std::runtime_error( std::format(
              "validateCoordinateMatrix: non-finite value at node {} coord {}", i, j ) );
          absolute_sum += std::abs(val);
        }

    if ( absolute_sum < 1.e-18 )
      throw std::runtime_error(
        "validateCoordinateMatrix: matrix contains only zeros — CoordinateMatrix() failed to load data" );

    // 3. Collocation check — squared distance helper
    auto dist_sq = [&]( uint32_t n1, uint32_t n2 ) -> double
      {
        double ds{0.0};
        for ( uint32_t d{0U}; d < dim; ++d )
          {
            const double diff = XY(n1,d) - XY(n2,d);
            ds += diff * diff;
          }
        return ds;
      };

    const double L_sq = dist_sq(0U, 1U);
    if ( L_sq < 1.e-14 )
      throw std::runtime_error(
        "validateCoordinateMatrix: degenerate element — corner nodes 0 and 1 are collocated" );

    // 4. Quadratic mid-side check
    if ( expected_nodes == 3U )
      {
        if ( dist_sq(0U,2U) < 1.e-14 || dist_sq(1U,2U) < 1.e-14 )
          throw std::runtime_error(
            "validateCoordinateMatrix: mid-side node 2 is collocated with a corner node" );

        // dot product check — mid-side node should lie between the two corners
        double dot{0.0};
        for ( uint32_t d{0U}; d < dim; ++d )
          dot += ( XY(1U,d) - XY(0U,d) ) * ( XY(2U,d) - XY(0U,d) );

        if ( dot < 0.0 || dot > L_sq )
          std::cerr << "validateCoordinateMatrix: warning — mid-side node is outside the segment; "
                  "Jacobian will be erratic.\n";
      }
 }


template<uint32_t dim>
void LU_Decomposition( DenseMatrix<dim>& A )
 {
    assert( A.Rows() == A.Columns() /* LU decompositions makes sense only for square matrices */ );
    const size_t n(A.Rows());
    for (uint32_t k=0;k<n-1U;k++)
      {
        for (size_t i=k+1U;i<n;i++) {
              A(i,k)=A(i,k)/A(k,k);
              for (size_t j=k+1U;j<n;j++) {
                  A(i,j)=A(i,j)-A(i,k)*A(k,j);
                }
          }
      }
 } // end LU_Decomposition


template<uint32_t dim>
void LU_BackSubstitution( const DenseMatrix<dim>& A, std::vector<double>& z )
 {
    assert( A.Rows() == A.Columns() /* LU decompositions makes sense only for square matrices */ );
    z.resize(A.Rows());
    const auto n(A.Rows());

    for (uint32_t k=0;k<n;k++){
        for(uint32_t i{0U};i<n;i++){
            if (i==k){z[i]=1;} else {z[i]=0;}
            for(uint32_t j{0U};j<i;j++){
                z[i]=z[i]-A(i,j)*z[j];
            }
            std::cout<<z[i]<< std::endl;
        }

    DenseMatrix<dim> X(A.Rows(),A.Columns(),0.);
    
        for(auto i=n-1;i+1>0;i--){
            X(i,k)=z[i]/A(i,i);
            for(auto j=n-1;j>i;j--){
                X(i,k)=X(i,k)-A(i,j)/A(i,i)*X(j,k);
            }
        }

    }
 }


} // end csmp
