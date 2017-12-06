#include "PoreVolume.h"
#include "Model.h"
#include "Node.h"

using namespace std;

namespace csmp {

	template<size_t dim>
		PoreVolume<dim>::PoreVolume(Model<dim>& model, const char* porevolume)
			: MatrixOperator<dim>(ADD_ACCUMULATE), pv_key_(model.Database().StorageKey(porevolume))
		{

		}

  template<size_t dim>
  void PoreVolume<dim>::AccumulateFiniteVolume(const Node<dim> *nptr, SparseMatrix &A) const {
			double64 pv = nptr->Read(pv_key_);
			size_t i = nptr->Idx();
    if (this->multiply_with_dt_) {
      A.Assign(i, i, pv / this->dt_);
    }
    else {
      A.Assign(i, i, pv);
    }
  
		}

  template class PoreVolume<1U>;
  template class PoreVolume<2U>;
  template class PoreVolume<3U>;

} // end csmp
