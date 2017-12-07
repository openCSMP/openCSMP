//
//  DiffusionLHS.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 5/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_DIFFUSION_LHS_H
#define CSMP_DIFFUSION_LHS_H

#include "CSMP_definitions.h"
#include "MatrixOperator.h"
#include "DenseMatrix.h"
#include "VariableSet_TracerTransferImplicit.h"

namespace csmp {

/**

@brief accumulates finite-element stencils for
diffusion into a sparse solution matrix.

@author SKM
@date 6/12/2017

*/
template<size_t dim>
class DiffusionLHS : public MatrixOperator<dim> {
  public:
    DiffusionLHS( const Model<dim>&, const char* diffusivity );
    DiffusionLHS( const DiffusionLHS& );
    virtual ~DiffusionLHS() {}
  
    /// not available for this Element-based stencil
  virtual void AccumulateFiniteVolume( const Node<dim>*, SparseMatrix& ) const { }
  
    /// element-by-element accumulation of matrix terms
    virtual void AccumulateStencil( const Element<dim>*, SparseMatrix& ) const;
  
  private:
    mutable DenseMatrix<DM_MIN>  DN_, DNT_, RESULT_;
    const csmp::Index    diff_key_;  ///< diffusion coefficient
    size_t               dof_;       ///< degrees of freedom (for scalar variables dof=1, for array variables their size)
};

} // end csmp

#endif /* CSMP_DIFFUSION_LHS_H */
