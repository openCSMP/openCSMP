//
//  DiffusionLHS.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 5/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_DIFFUSION_LHS_H
#define CSMP_DIFFUSION_LHS_H

#include "MatrixOperator.h"
#include "DenseMatrix.h"

namespace csmp {

template<size_t dim> {
class DiffusionLHS {
  public:
    DiffusionLHS();
    DiffusionLHS( const DiffusionLHS& );
    virtual void ~DiffusionLHS();
  
    virtual void AccumulateFiniteVolume( Node<dim>&, SparseMatrix& );
    virtual void AccumulateStencil( Element<dim>&, SparseMatrix& );
  
  private:
    DenseMatrix<DM_MIN>  DN_, DNT_;

};

}

#endif /* CSMP_DIFFUSION_LHS_H */
