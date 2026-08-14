#include "PoreVolumeLHS.h"
#include "Element.h"
#include "Node.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
PoreVolumeLHS<dim>::PoreVolumeLHS( const INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& spv_key,
                                   const INDEX<SCALAR,NODE>& fpv_key )
  : spv_key_(spv_key),
    fpv_key_(fpv_key)
  {
  }
  
    
template<uint32_t dim>
void PoreVolumeLHS<dim>::AccumulateStencil( const Element<dim>& e, SparseMatrix& A ) const 
 {
    assert( e.IntegrationPointsPerSector() == 1U );
    const size_t sector_ipoints(e.Sectors()); 
    
    for ( auto i{0U}; i <sector_ipoints; ++i ) {
        const double sector_pore_volume = e.Read( i, 0U, spv_key_ );
        const size_t   idx = e.N(i)->Idx();
        
        A.Add( idx, idx, sector_pore_volume * this->Factor() );
    }
}
  

template<uint32_t dim>
void PoreVolumeLHS<dim>::AccumulateFiniteVolume( const Node<dim>& fv, SparseMatrix& A ) const 
 {
    const double fpv = fv.Read( fpv_key_ );
    const size_t   idx = fv.Idx();
    
    A.Assign( idx, idx, fpv * this->Factor() );
  }

//template class PoreVolumeLHS<1U>;
//template class PoreVolumeLHS<2U>;
template class PoreVolumeLHS<3U>;

} // end csmp
