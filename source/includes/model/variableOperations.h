// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef VARIABLE_OPERATIONS_H
#define VARIABLE_OPERATIONS_H

#include "Point.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

namespace csmp {

/// generic dot product for vectors
template<uint32_t dim>
double dotProduct( const std::vector<double>& v1, const std::vector<double>& v2 );

/// Point class related operators
template<uint32_t dim>
Point<dim> operator+(const Point<dim>&, const VectorVariable<dim>& );

template<uint32_t dim>
Point<dim> operator-(const Point<dim>&, const VectorVariable<dim>& );

template<uint32_t dim>
Point<dim> operator*(const Point<dim>&, const VectorVariable<dim>& );

template<uint32_t dim>
Point<dim> operator/(const Point<dim>&, const VectorVariable<dim>& );

template<uint32_t dim>
double dotProduct( const VectorVariable<dim>&, const csmp::Point<dim>& );

template<uint32_t dim>
double dotProduct( const VectorVariable<dim>&, const VectorVariable<dim>& );

VectorVariable<1U> crossProduct( const VectorVariable<1U>&, const csmp::Point<1U>& );
VectorVariable<2U> crossProduct( const VectorVariable<2U>&, const csmp::Point<2U>& );
VectorVariable<3U> crossProduct( const VectorVariable<3U>&, const csmp::Point<3U>& );

VectorVariable<1U> crossProduct( const VectorVariable<1U>&, const VectorVariable<1U>& );
VectorVariable<2U> crossProduct( const VectorVariable<2U>&, const VectorVariable<2U>& );
VectorVariable<3U> crossProduct( const VectorVariable<3U>&, const VectorVariable<3U>& );

template<uint32_t dim>
VectorVariable<dim> multiplyTensorByVector( const TensorVariable<dim>&, const VectorVariable<dim>& );

template<uint32_t dim>
VectorVariable<dim> multiplyHorizontalVectorByTensor( const VectorVariable<dim>&, const TensorVariable<dim>&  );

template<uint32_t dim>
TensorVariable<dim> multiplyTensorByTensor( const TensorVariable<dim>&, const TensorVariable<dim>& );

template<uint32_t dim>
bool isDiagonalTensor( const TensorVariable<dim>& ts );

template<uint32_t dim>
void minMaxEigenValues( const TensorVariable<dim>& ts, double& tmin, double& tmax );

template<uint32_t dim>
double  angleBetween( const VectorVariable<dim>&, const VectorVariable<dim>& );

} // end namespace csmp

#endif  //VARIABLE_OPERATIONS_H












