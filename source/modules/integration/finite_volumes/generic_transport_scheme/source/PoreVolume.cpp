#include "PoreVolume.h"
#include "Model.h"
#include "Node.h"

using namespace std;

namespace csmp {

	template<size_t dim>
		PoreVolume<dim>::PoreVolume(Model<dim>& model, const char* porevolume, const char* porosity)
			: MatrixOperator<dim>(ADD_ACCUMULATE), pv_key_(model.Database().StorageKey(porevolume)),
  phi_key_(model.Database().StorageKey(porosity))
		{

		}
  template<size_t dim>
  void PoreVolume<dim>::AccumulateStencil(const Element<dim> *eptr, SparseMatrix &A) const {
    const size_t iNrNodes(eptr->Nodes());
    const double64 phi(eptr->Read(phi_key_));
    
    for (size_t iNode = 0; iNode < iNrNodes; ++iNode ) {
      const auto nptr = eptr->N(iNode);
      const double64 sector_volume = eptr->SectorVolume(iNode);
      const size_t i = nptr->Idx();
      
      if (this->multiply_with_dt_) {
        A.Add(i, i, phi * sector_volume * this->dt_);
      }
      else {
        A.Add(i, i, phi * sector_volume);
      }
    }
  }
  

  template<size_t dim>
  void PoreVolume<dim>::AccumulateFiniteVolume(const Node<dim> *nptr, SparseMatrix &A) const {
			double64 pv = nptr->Read(pv_key_);
			size_t i = nptr->Idx();
    if (this->multiply_with_dt_) {
      A.Assign(i, i, pv * this->dt_);
    }
    else {
      A.Assign(i, i, pv);
    }
  
		}

  template class PoreVolume<1U>;
  template class PoreVolume<2U>;
  template class PoreVolume<3U>;

} // end csmp
