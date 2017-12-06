#include "PoreVolume.h"

using namespace std;

namespace csmp {

	template<size_t dim>
		PoreVolume<dim>::PoreVolume(Model<dim>& model, const char* porevolume, double64 delta_t)
			: MatrixOperator<dim>(ADD_ACCUMULATE), pv_key_(model.Database().StorageKey(porevolume), dt_(delta_t))
		{

		}

		template<size_t dim>
		void AccumulateFV(Node<dim> *nptr, SparseMatrix &A) {
			double64 pv = nptr->Read(pv_key_);
			size_t i = nptr->Idx();
			A(i, i) = pv / dt_;
		}
	
} // end csmp
