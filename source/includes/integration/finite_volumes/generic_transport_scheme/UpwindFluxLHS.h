// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_UPWIND_FLUX_LHS_H
#define CSMP_UPWIND_FLUX_LHS_H

#include "MatrixOperator.h"

namespace csmp {

/**
     First-order upwind flux.
     
       Incoming fluxes create off-diagonal couplings.
       Outgoing fluxes are added to the matrix diagonal.
       
       TODO: create version for tensor permeabilities where the flux is actually calculated from the pressure gradient
*/
template<uint32_t dim>
class UpwindFluxLHS : public MatrixOperator<dim> {
  public:
    explicit UpwindFluxLHS( const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>& vol_flow_across_facet );
  
    virtual ~UpwindFluxLHS() {}

    virtual void AccumulateFiniteVolume( const Node<dim>&, SparseMatrix& ) const;
    
    virtual void AccumulateStencil( const Element<dim>&, SparseMatrix& ) const;
//    virtual void AccumulateStencil( const Face<dim>&, SparseMatrix& ) const;
//    virtual void AccumulateStencil( const InterFace<dim>&, SparseMatrix& ) const;

  private:
    const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>& ff_key_; ///< reference to volumetric facet flux key
};

} // end csmp

#endif /* CSMP_UPWIND_FLUX_LHS_H */
