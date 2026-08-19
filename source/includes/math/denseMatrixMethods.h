//
//  denseMatrixMethods.h
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 17/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_DENSE_MATRIX_METHODS_H
#define CSMP_DENSE_MATRIX_METHODS_H

#include "DenseMatrix.h"

namespace csmp {

/**
    @brief Applies row-sum lumping to a square matrix in-place.

    Each diagonal entry absorbs the sum of all entries in its row.
    All off-diagonal entries are zeroed. Conservation is preserved:
    the column sums of the lumped matrix equal those of the consistent
    matrix.

    @param matrix  Square matrix to lump in-place.
    
    @attention do not use for diffusion operators as it will produce zeroes in the diagonal.
*/
template<uint32_t mn_max>
inline void lumpByRowSum( DenseMatrix<mn_max>& matrix )
{
    const uint32_t n = matrix.Rows();
    for ( uint32_t i{ 0U }; i < n; ++i ) {
        double row_sum = 0.;
        for ( uint32_t j{ 0U }; j < n; ++j ) {
            row_sum += matrix(i,j);
            if ( i != j ) matrix(i,j) = 0.;
        }
        matrix(i,i) = row_sum;
    }
}


/**
   This lumping preserves the sign structure of an element stiffness matrix —
   the diagonal is positive and the off-diagonals are negative for a standard diffusion operator.
   This ensures the lumped matrix still satisfies the discrete maximum principle.

    @param matrix  Square matrix to lump in-place.
    
 */
inline void lumpByNegativeOffDiagonalSum( DenseMatrix<DM_MIN>& matrix )
{
    const uint32_t n = matrix.Rows();
    for ( uint32_t i{ 0U }; i < n; ++i ) {
        double off_diagonal_sum = 0.;
        for ( uint32_t j{ 0U }; j < n; ++j )
            if ( i != j ) {
                off_diagonal_sum += matrix(i,j);
                matrix(i,j) = 0.;
            }
        matrix(i,i) = -off_diagonal_sum;
    }
}


/**
    @brief Applies volume lumping to a square matrix in-place.

    Each diagonal entry is set to the physical element volume divided
    by the number of nodes, distributing the element volume equally
    among its nodes. All off-diagonal entries are zeroed.

    This is equivalent to replacing the consistent mass matrix with
    a diagonal matrix whose entries are the nodal volumes, and is
    appropriate when the material property is uniform over the element.

    @param matrix       Square matrix to lump in-place.
    @param node_volume  Physical volume per node = element_volume / n_nodes.
*/
template<uint32_t mn_max>
inline void lumpByVolume( DenseMatrix<mn_max>& matrix, double node_volume )
{
    const uint32_t n = matrix.Rows();
    for ( uint32_t i{ 0U }; i < n; ++i ) {
        for ( uint32_t j{ 0U }; j < n; ++j )
            matrix(i,j) = 0.;
        matrix(i,i) = node_volume;
    }
}


/**
    @brief Dispatches to the appropriate lumping function.

    Reads the lumping strategy from the base class flag. If lumping
    by volume is requested, the physical element volume must be
    supplied so that it can be distributed equally among the nodes.

    @param matrix           Square matrix to lump in-place.
    @param is_lumped        Whether lumping is requested at all.
    @param element_volume   Physical element volume, used only for volume lumping.
    @param n_nodes          Number of nodes in the element.
*/
template<uint32_t mn_max>
inline void applyLumping( DenseMatrix<mn_max>& matrix,
                          bool                 is_lumped,
                          bool                 lump_by_volume,
                          double               element_volume,
                          uint32_t             n_nodes )
{
    if ( !is_lumped ) return;

    if ( lump_by_volume )
        lumpByVolume( matrix, element_volume / static_cast<double>(n_nodes) );
    else
        lumpByRowSum( matrix );
}


/**
 * @brief Rigorous geometric and structural validation of a finite element coordinate matrix.
 * @param XY             Coordinate matrix — layout XY(node_idx, coord_idx)
 * @param expected_nodes Number of nodes the element expects (e.g. 3 for quadratic line)
 * @param dim            Model dimension (2 or 3)
 * Declared here, defined in DenseMatrix.cpp — not performance-critical.
 */
template<typename MatrixType>
void validateCoordinateMatrix( const MatrixType& XY, uint32_t expected_nodes, uint32_t dim );


/// turns the original matrix into upper diagonal form
template<uint32_t dim> void LU_Decomposition( DenseMatrix<dim>& );

template<uint32_t dim> void LU_BackSubstitution( const DenseMatrix<dim>&, std::vector<double>& solution );

} // end csmp

#endif /* defined(CSMP_DENSE_MATRIX_METHODS_H) */
