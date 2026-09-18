// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PICARD_INTEGRATOR_H
#define PICARD_INTEGRATOR_H

#include "IterativeIntegrator.h"

namespace csmp {

template<uint32_t> class Element;

/**
 
@class PicardIntegrator
@author Adrian Burri
@date 2004

Algorithm implementing a Picard iteration scheme for non-linear processes.
Derived from the IterativeAlgorithm base class.
 
*/
template<uint32_t dim,template<uint32_t> class CELLTYPE=Element>
class PicardIntegrator : public IterativeIntegrator<dim,CELLTYPE> {
  public:
    PicardIntegrator();
    
    virtual double Residual();
    
  private:
    std::vector<double> resid_;
};

} // end namespace csmp
#endif
