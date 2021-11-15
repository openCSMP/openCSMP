//
//  tensorOperations.h
//
//  Created by Stephan Matthai on 16/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#ifndef CSMP_TENSOR_OPERATIONS_H
#define CSMP_TENSOR_OPERATIONS_H

#include "TensorVariable.h"

namespace csmp {

/// recovers scalar value, vector length, or average of tensor eigenvalues for variable with given csmp::Index
template<size_t dim,template<size_t> class PLACE>
double magnitude( const PLACE<dim>&, const csmp::Index& );

/// returns magnitude of tensor in the direction (transect) given by the supplied vector variable
template<size_t dim>
inline double tensorMagnitudeInDirection( const TensorVariable<dim>& ts, const VectorVariable<dim>& vc ) {
    VectorVariable<dim>  prod(ts * vc);
    return prod.Length();
 }


///  t = (ts . n) n + n x (ts x n)
double tensorProjectionOnPlane( const TensorVariable<3U>&,
                                  const Point<3U>&,
                                  double& component_n );

double tensorProjectionOnPlane( const TensorVariable<3U>&,
                                  const VectorVariable<3U>&,
                                  double& component_n );


///  tensor component acting tangential to a plane: n x (ts x n) (Cauchy's formula); returns component tc and its magnitude
double tangentialTensorProjection( const TensorVariable<3U>&, const Point<3U>& unit_normal, Point<3U>& tangential_component );


///  tensor component acting normal to a plane: (ts . n) n (Cauchy's formula)
inline double normalTensorProjection( const TensorVariable<3U>& ts, const Point<3U>& un )
 {
    // (t . n) n (eqn. 6.49, P&F, p.216)
    return ts(0,0) * un[0] * un[0] + ts(1,1) * un[1] * un[1] +
           ts(2,2) * un[2] * un[2] + 2. * ts(0,1) * un[0] * un[1] +
           2. * ts(1,2) * un[1] * un[2] + 2. * ts(2,0) * un[2] * un[0];
 }

inline double normalTensorProjection( const TensorVariable<3U>& ts, const VectorVariable<3U>& un )
 {
    // (t . n) n (eqn. 6.49, P&F, p.216)
    return ts(0,0) * un[0] * un[0] + ts(1,1) * un[1] * un[1] +
           ts(2,2) * un[2] * un[2] + 2. * ts(0,1) * un[0] * un[1] +
           2. * ts(1,2) * un[1] * un[2] + 2. * ts(2,0) * un[2] * un[0];
 }

/// @todo to rotate tensor T, implement:  T' = (R . T) . R^T, where R is the rotation matrix (www.continuummechanics.org/cm/rotationmatrix.html)
  
} // end csmp

#endif /* defined(CSMP_TENSOR_OPERATIONS_H) */
