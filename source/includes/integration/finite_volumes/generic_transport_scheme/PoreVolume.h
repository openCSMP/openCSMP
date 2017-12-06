#ifndef PORE_VOLUME_H
#define PORE_VOLUME_H

#include "MatrixOperator.h"

namespace csmp {

template <size_t dim>
class PoreVolume : public MatrixOperator<dim> {	
	public: 
	PoreVolume(Model<dim>& ,const char *);
	virtual void AccumulateFiniteVolume(const Node<dim>* nptr, SparseMatrix& A) const;
	
	private:
	csmp::Index pv_key_;
};

} // csmp

#endif // PORE_VOLUME_H
