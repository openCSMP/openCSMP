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

/// turns the original matrix into upper diagonal form
template<size_t dim> void LU_Decomposition( DenseMatrix<dim>& );

template<size_t dim> void LU_BackSubstitution( const DenseMatrix<dim>&, std::vector<double64>& solution );

} // end csmp

#endif /* defined(CSMP_DENSE_MATRIX_METHODS_H) */
