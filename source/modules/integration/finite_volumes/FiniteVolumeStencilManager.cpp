#include "FiniteVolumeStencilManager.h"
#include "FiniteElementManager.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim>
FiniteVolumeStencilManager<dim>::FiniteVolumeStencilManager()
 {
 }


template<size_t dim>
FiniteVolumeStencilManager<dim>::FiniteVolumeStencilManager( const FiniteElementManager& fem_manager )
 {
    Initialize( fem_manager );
 }


template<size_t dim>
FiniteVolumeStencilManager<dim>::~FiniteVolumeStencilManager()
 {
 }



template<size_t dim>
void FiniteVolumeStencilManager<dim>::Initialize( const FiniteElementManager& fem_manager )
 {
    if ( !type_mapping.empty() ) return;
    // loop over the FEM currently connected to the FEM manager and create corresponding
    // finite volume stencils which are then stored in the stencils list
    list<CSMP_FEM_TYPE>  etypes;
    
    fem_manager.CurrentElementTypes( etypes );
    
    // 1. creating mapping between finite-element types and stencil array entries
    size_t  type_counter(0U);

    for ( typename list<CSMP_FEM_TYPE>::const_iterator it=etypes.begin(); it!=etypes.end(); it++ )
      type_mapping.insert( make_pair( *it, type_counter++ ) );
    
    // 2. Initializing the finite-volume stencils
    type_counter = 0U;
    for ( typename list<CSMP_FEM_TYPE>::const_iterator it=etypes.begin(); it!=etypes.end(); it++ )
      stencils[type_counter++].Initialize( parseFiniteElementType(*it) );
   
#ifndef NDEBUG
    cout <<"\nFiniteVolumeStencilManager<"<< dim <<">::Initialize: stencil manager initialized successfully.\n";
#endif
 }


template class FiniteVolumeStencilManager<1U>;
template class FiniteVolumeStencilManager<2U>;
template class FiniteVolumeStencilManager<3U>;

} // end namespace csmp
