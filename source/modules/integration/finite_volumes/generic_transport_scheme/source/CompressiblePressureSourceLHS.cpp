
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
   phi_key_(model.Database().StorageKey(porosity)),
   CT_key_(model.Database().StorageKey(compressibility)),
   pf0_key_(model.Database().StorageKey(previous_pressure)),
   pf1_key_(model.Database().StorageKey(current_pressure))
{
} // end constructor


template<size_t dim>
void CompressiblePressureSourceLHS<dim>::AccumulateStencil( const Element<dim>* fe, SparseMatrix& mat ) const
{
    const size_t sectors(fe->Sectors());
    const double64 phi = fe->Read( phi_key_ );
    const double64 ct = fe->Read( CT_key_ );
    const double64 dt = this->MultiplyWithTimeIncrement() ? this->dt_ : 1.0;

    for ( size_t i=0U; i<sectors; ++i ) {
        const double64 sector_volume = fe->SectorVolume(i);
        const double64 pore_volume = phi * sector_volume;
        const double64 pf0 = fe->PropertyValueAtSectorIntegrationPoint(i,0,pf0_key_);
        const double64 pf1 = fe->PropertyValueAtSectorIntegrationPoint(i,0,pf1_key_);

        mat.Add(fe->N(i)->Idx(),fe->N(i)->Idx(),pore_volume*ct*(pf1-pf0)*dt);
    }
} // end AccumulateStencil


template<size_t dim>
void CompressiblePressureSourceLHS<dim>::AccumulateFiniteVolume( const Node<dim>* nd, SparseMatrix& mat ) const
{
    const size_t parent_elements(nd->Parents());
    const double64 dt = this->MultiplyWithTimeIncrement() ? this->dt_ : 1.0;
    for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<dim>* const fe = nd->Parent(i);
         assert( fe != NULL );
         const double64 phi = fe->Read( phi_key_ );
         const double64 ct = fe->Read( CT_key_ );

         const size_t pnid(nd->ParentNodeNumber(i));
         const double64 sector_volume = fe->SectorVolume(pnid);
         const double64 pore_volume = phi * sector_volume;
         const double64 pf0 = fe->PropertyValueAtSectorIntegrationPoint(i,0,pf0_key_);
         const double64 pf1 = fe->PropertyValueAtSectorIntegrationPoint(i,0,pf1_key_);
         mat.Add(nd->Idx(),nd->Idx(),pore_volume*ct*(pf1-pf0)*dt);
    }         
} // end AccumulateFiniteVolume

template class CompressiblePressureSourceLHS<1U>;
template class CompressiblePressureSourceLHS<2U>;
template class CompressiblePressureSourceLHS<3U>;

} // end csmp


