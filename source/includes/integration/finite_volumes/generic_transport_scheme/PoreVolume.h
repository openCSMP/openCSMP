#ifndef PORE_VOLUME_H
#define PORE_VOLUME_H

#include "MatrixOperator.h"

namespace csmp {

template <size_t dim>
class PoreVolume : public MatrixOperator<dim> {	
	public: 
	PoreVolume(Model<dim>& ,const char *pv, const char* porosity);
	virtual void AccumulateFiniteVolume(Node<dim>& nptr, SparseMatrix& A) const;
  virtual void AccumulateStencil(Element<dim>& eptr, SparseMatrix& A) const;

	private:
	csmp::Index pv_key_;
  csmp::Index phi_key_;

};

} // csmp

#endif // PORE_VOLUME_H
