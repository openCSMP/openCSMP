#ifndef VARIABLE_OPERATIONS_H
#define VARIABLE_OPERATIONS_H

#include "Point.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

namespace csmp {

/// generic dot product for vectors
template<size_t dim>
double64 dotProduct( const std::vector<double64>& v1, const std::vector<double64>& v2 );

/// Point class related operators
template<size_t dim>
Point<dim> operator+(const Point<dim>&, const VectorVariable<dim>& );

template<size_t dim>
Point<dim> operator-(const Point<dim>&, const VectorVariable<dim>& );

template<size_t dim>
Point<dim> operator*(const Point<dim>&, const VectorVariable<dim>& );

template<size_t dim>
Point<dim> operator/(const Point<dim>&, const VectorVariable<dim>& );

template<size_t dim>
double64 dotProduct( const VectorVariable<dim>&, const csmp::Point<dim>& );

template<size_t dim>
double64 dotProduct( const VectorVariable<dim>&, const VectorVariable<dim>& );

/// this function treats the dot product of two scalars as if they were 1D vectors, hence the formula below.
//double64 DotProduct( const ScalarVariable& s1, const ScalarVariable& s2) {return (s1.Data()*s2.Data())*s2.Data()/(s2.Data()*s2.Data());}

VectorVariable<1U> crossProduct( const VectorVariable<1U>&, const csmp::Point<1U>& );
VectorVariable<2U> crossProduct( const VectorVariable<2U>&, const csmp::Point<2U>& );
VectorVariable<3U> crossProduct( const VectorVariable<3U>&, const csmp::Point<3U>& );

VectorVariable<1U> crossProduct( const VectorVariable<1U>&, const VectorVariable<1U>& );
VectorVariable<2U> crossProduct( const VectorVariable<2U>&, const VectorVariable<2U>& );
VectorVariable<3U> crossProduct( const VectorVariable<3U>&, const VectorVariable<3U>& );

template<size_t dim>
VectorVariable<dim> multiplyTensorByVector( const TensorVariable<dim>&, const VectorVariable<dim>& );

template<size_t dim>
VectorVariable<dim> multiplyHorizontalVectorByTensor( const VectorVariable<dim>&, const TensorVariable<dim>&  );

template<size_t dim>
TensorVariable<dim> multiplyTensorByTensor( const TensorVariable<dim>&, const TensorVariable<dim>& );

template<size_t dim>
bool isDiagonalTensor( const TensorVariable<dim>& ts );

template<size_t dim>
void minMaxEigenValues( const TensorVariable<dim>& ts, double64& tmin, double64& tmax );

template<size_t dim>
double64  angleBetween( const VectorVariable<dim>&, const VectorVariable<dim>& );

/// applies any of the standard math library functions to each element of the variable
template<size_t, template<size_t> class Var>
void applyFunction();


} // end namespace csmp

#endif  //VARIABLE_OPERATIONS_H












