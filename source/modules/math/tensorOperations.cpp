//
//  tensorOperations.cpp
//
//  Created by Stephan Matthai on 16/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#include "tensorOperations.h"
#include "TensorVariable.h"
#include "Region.h"
#include "Element.h"
#include "ErrorHandler.h"
#include "Point.h"
#include "vectorOperations.h"

namespace csmp {

/**
   recovers scalar value, vector length, 
   or the average of the eigenvalues of the tensor property 
   with the given csmp::Index
*/
template<size_t dim,template<size_t> class PLACE>
double64 magnitude( const PLACE<dim>& site, const csmp::Index& prop_key )
 {
    if (  prop_key.type == SCALAR ) return site.Read( prop_key );
   
    // calculate vector length
    if (  prop_key.type == VECTOR ) {
         VectorVariable<dim> vc;
         site.Read( prop_key, vc );
         return vc.Length();
      }
    // calculate average of Eigen values
    if (  prop_key.type == TENSOR ) {
         TensorVariable<dim> ts;
         site.Read( prop_key, ts );
         VectorVariable<dim> Eigenvalues;
         ts.EigenValues( Eigenvalues );
         return valueAverage(Eigenvalues);
      }

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( ERROR, "propertyMagnitude:", parseType(prop_key.type).c_str(),
                      "method for calculating magnitude of this variable type undefined; returning NaN.");

    return std::numeric_limits<double64>::quiet_NaN();
 }

template double64 magnitude( const Region<3>&, const csmp::Index& );
template double64 magnitude( const Element<3>&, const csmp::Index& );
template double64 magnitude( const Node<3>&, const csmp::Index& );

template double64 magnitude( const Region<2>&, const csmp::Index& );
template double64 magnitude( const Element<2>&, const csmp::Index& );
template double64 magnitude( const Node<2>&, const csmp::Index& );


/**

Resolves surface-normal- and tangential components of the supplied tensor, ts.
The surface is defined by its outward-pointing unit normal n.
The components are resolved using Cauchy's formula.

 t = (ts . n) n + n x (ts x n)
 
@return The normal component is returned into
the 3rd function argument and the and tangential component
is the return value
 
*/
double64 tensorProjectionOnPlane( const TensorVariable<3U>& ts,
                                  const Point<3U>& un,
                                  double64& component_n )
 {
    // Cauchy's formula applied to find traction vector components, P&F, p. 213
    // VectorVariable<3U> t = ts * un.Coordinates();
    const double64  tx = ts(0,0) * un[0] + ts(1,0) * un[1] + ts(2,0) * un[2];
    const double64  ty = ts(0,1) * un[0] + ts(1,1) * un[1] + ts(2,1) * un[2];
    const double64  tz = ts(0,2) * un[0] + ts(1,2) * un[1] + ts(2,2) * un[2];
 
    // (t . n) n (eqn. 6.49, P&F, p.216)
    component_n  = ts(0,0) * un[0] * un[0] + ts(1,1) * un[1] * un[1] + ts(2,2) * un[2] * un[2];
    component_n += 2. * ts(0,1) * un[0] * un[1] + 2. * ts(1,2) * un[1] * un[2] + 2. * ts(2,0) * un[2] * un[0];
    
    // n x (t x n) (eqn. 6.52, P&F, p. 216) -> vector product, vp
    const double64  vpx = ((1. - un[0] * un[0]) * tx - un[0] * un[1] * ty - un[0] * un[2] * tz);  // * ex;
    const double64  vpy = (-un[0] * un[1] * tx + (1. - un[1] * un[1]) * ty - un[1] * un[2] * tz); // * ey;
    const double64  vpz = (-un[2] * un[0] * tx - un[2] * un[1] * ty + (1. - un[2] * un[2]) * tz); // * ez;
    
    // the shear stress is the magnitude of the vector product
    return sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    
 } // end tensorProjectionOnPlane


double64 tensorProjectionOnPlane( const TensorVariable<3U>& ts,
                                  const VectorVariable<3U>& un,
                                  double64& component_n )
 {
    // Cauchy's formula applied to find traction vector components, P&F, p. 213
    // VectorVariable<3U> t = ts * un.Coordinates();
    const double64  tx = ts(0,0) * un[0] + ts(1,0) * un[1] + ts(2,0) * un[2];
    const double64  ty = ts(0,1) * un[0] + ts(1,1) * un[1] + ts(2,1) * un[2];
    const double64  tz = ts(0,2) * un[0] + ts(1,2) * un[1] + ts(2,2) * un[2];
 
    // (t . n) n (eqn. 6.49, P&F, p.216)
    component_n  = ts(0,0) * un[0] * un[0] + ts(1,1) * un[1] * un[1] + ts(2,2) * un[2] * un[2];
    component_n += 2. * ts(0,1) * un[0] * un[1] + 2. * ts(1,2) * un[1] * un[2] + 2. * ts(2,0) * un[2] * un[0];
    
    // n x (t x n) (eqn. 6.52, P&F, p. 216) -> vector product, vp
    const double64  vpx = ((1. - un[0] * un[0]) * tx - un[0] * un[1] * ty - un[0] * un[2] * tz);  // * ex;
    const double64  vpy = (-un[0] * un[1] * tx + (1. - un[1] * un[1]) * ty - un[1] * un[2] * tz); // * ey;
    const double64  vpz = (-un[2] * un[0] * tx - un[2] * un[1] * ty + (1. - un[2] * un[2]) * tz); // * ez;
    
    // the shear stress is the magnitude of the vector product
    return sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    
 } // end tensorProjectionOnPlane



///  tensor component acting tangential to a plane: n x (ts x n) (Cauchy's formula); returns component tc and its magnitude
double64 tangentialTensorProjection( const TensorVariable<3U>& ts, const Point<3U>& un, Point<3U>& tc )
 {
    // Cauchy's formula applied to find traction vector components, P&F, p. 213
    // VectorVariable<3U> t = ts * un.Coordinates();
    const double64  tx = ts(0,0) * un[0] + ts(1,0) * un[1] + ts(2,0) * un[2];
    const double64  ty = ts(0,1) * un[0] + ts(1,1) * un[1] + ts(2,1) * un[2];
    const double64  tz = ts(0,2) * un[0] + ts(1,2) * un[1] + ts(2,2) * un[2];
 
    // n x (t x n) (eqn. 6.52, P&F, p. 216) -> vector product, vp
    tc[0] = ((1. - un[0] * un[0]) * tx - un[0] * un[1] * ty - un[0] * un[2] * tz);  // * ex;
    tc[1] = (-un[0] * un[1] * tx + (1. - un[1] * un[1]) * ty - un[1] * un[2] * tz); // * ey;
    tc[2] = (-un[2] * un[0] * tx - un[2] * un[1] * ty + (1. - un[2] * un[2]) * tz); // * ez;
    
    // the shear stress is the magnitude of the vector product
    return std::sqrt( tc[0] * tc[0] + tc[1] * tc[1] + tc[2] * tc[2] );
 }
 

} // end csmp


