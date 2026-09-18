// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_PORE_VOLUME_LHS_H
#define CSMP_PORE_VOLUME_LHS_H

#include "MatrixOperator.h"

namespace csmp {

/**
    The pore volume is accumulated into the diagonal of the sparse (solution) matrix
*/
template <uint32_t dim>
class PoreVolumeLHS : public MatrixOperator<dim> {	
	public: 
	  PoreVolumeLHS( const csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& spv_key,
                   const csmp::INDEX<SCALAR,NODE>& fpv_key );
    
    virtual ~PoreVolumeLHS() {}

	  virtual void AccumulateFiniteVolume( const Node<dim>& nptr, SparseMatrix& A ) const;

    virtual void AccumulateStencil( const Element<dim>& eptr, SparseMatrix& A ) const;
//    virtual void AccumulateStencil( const Face<dim>& eptr, SparseMatrix& A ) const;
//    virtual void AccumulateStencil( const InterFace<dim>& eptr, SparseMatrix& A ) const;

	private:
	  const INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& spv_key_;
	  const INDEX<SCALAR,NODE>& fpv_key_; ///<  integrated pore volume of the entire finite volume
};

} // csmp

#endif // PORE_VOLUME_LHS_H
