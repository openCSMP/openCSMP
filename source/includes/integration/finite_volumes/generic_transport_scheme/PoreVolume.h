#ifndef PORE_VOLUME_H
#define PORE_VOLUME_H

#include "MatrixOperator.h"

namespace csmp {

template <size_t dim>
class PoreVolume : public MatrixOperator<dim> {	
	public: 
	PoreVolume(Model<dim>& ,const char * ,double64 delta_t);
	virtual void AccumulateFV(Node<dim>* nptr, SparseMatrix& A, double64 dt);
	
	private:
	double64 dt_;
	csmp::Index pv_key_;
};

} // csmp

#endif // PORE_VOLUME_H