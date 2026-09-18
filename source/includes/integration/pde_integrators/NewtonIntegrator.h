// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NEWTON_RAPHSON_INTEGRATOR_H
#define NEWTON_RAPHSON_INTEGRATOR_H

#include "IterativeIntegrator.h"

namespace csmp {

template<uint32_t> class Element;

/** 
 
@brief Algorithm implementing a Newton iteration scheme for non-linear processes.

@author Adrian Burri
@date 2004

Based on the IterativeAlgorithm base class.

@section collaboration Collaboration

Models, Visitors and Interrelations. Calculations on single Region
objects is not supported.  

*/
template<uint32_t dim,template<uint32_t> class CELLTYPE=Element>
class NewtonIntegrator : public IterativeIntegrator<dim,CELLTYPE> {
public:
  NewtonIntegrator();
  
  virtual double Residual();
};


} // end namespace csmp

#endif
