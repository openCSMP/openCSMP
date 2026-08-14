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
