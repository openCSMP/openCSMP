
#include "CompressiblePressureSourceLHS.h"
#include "SparseMatrix.h"
#include "Model.h"
#include "Element.h"
#include "Node.h"

using namespace std;

namespace csmp {
  template<size_t dim>
  CompressiblePressureSourceLHS<dim>::CompressiblePressureSourceLHS(Model<dim>& model, const char* porosity, const char* compressibility, const char* previous_pressure, const char* current_pressure)
  : MatrixOperator<dim>(0),
  key_PHI(model.Database().StorageKey(porosity)),
  key_CT(model.Database().StorageKey(compressibility)),
  key_PF0(model.Database().StorageKey(previous_pressure)),
  key_PF1(model.Database().StorageKey(current_pressure))
  {
    if ( key_PHI.place != ELEMENT || key_PHI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "CompressiblePressureSourceLHS::CompressiblePressureSourceLHS:",
                            "The porosity variable must be SCALAR and placed on ELEMENT"  );
    if ( key_CT.place != NODE || key_CT.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "CompressiblePressureSourceLHS::CompressiblePressureSourceLHS:",
                            "The compressibility variable must be SCALAR and placed on NODE"  );
    if ( key_PF0.place != NODE || key_PHI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "CompressiblePressureSourceLHS::CompressiblePressureSourceLHS:",
                            "The previous pressure variable must be SCALAR and placed on NODE"  );
    if ( key_PF1.place != NODE || key_PHI.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "CompressiblePressureSourceLHS::CompressiblePressureSourceLHS:",
                            "The current pressure variable must be SCALAR and placed on NODE"  );
  } // end constructor
  
  
  template<size_t dim>
  void CompressiblePressureSourceLHS<dim>::AccumulateStencil( Element<dim>& fe, SparseMatrix& mat ) const
  {
    const double64 phi = fe.Read( key_PHI );
    const double64 ct = fe.Read( key_CT );
    const double64 dt = this->MultiplyWithTimeIncrement() ? this->dt_ : 1.0;
    for (auto sip: fe.AllSectorIntegrationPoints()) {
      const double64 sector_volume = sip.SectorVolume();
      const double64 pore_volume = phi * sector_volume;
      const double64 pf0 = sip.Obtain(key_PF0);
      const double64 pf1 = sip.Obtain(key_PF1);
      
      auto node_idx = sip.NodeIdx();
      mat.Add(node_idx,node_idx,pore_volume*ct*(pf1-pf0)*dt);
    }
  } // end AccumulateStencil
  
  
  template<size_t dim>
  void CompressiblePressureSourceLHS<dim>::AccumulateFiniteVolume( Node<dim>& fv, SparseMatrix& mat ) const
  {
    const double64 dt = this->MultiplyWithTimeIncrement() ? this->dt_ : 1.0;
    for (auto sip: fv.AllSectorIntegrationPoints()) {
      const double64 sector_volume = sip.SectorVolume();
      const double64 pf0 = sip.Obtain(key_PF0);
      const double64 pf1 = sip.Obtain(key_PF1);
      const double64 phi = sip.Obtain( key_PHI );
      const double64 ct = sip.Obtain( key_CT );
      const double64 pore_volume = phi * sector_volume;

      auto node_idx = sip.NodeIdx();
      mat.Add(node_idx,node_idx,pore_volume*ct*(pf1-pf0)*dt);
    }
  } // end AccumulateFiniteVolume
  
  template class CompressiblePressureSourceLHS<1U>;
  template class CompressiblePressureSourceLHS<2U>;
  template class CompressiblePressureSourceLHS<3U>;
  
} // end csmp


