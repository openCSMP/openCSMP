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


/**
    returns pointer to corresponding FV stencil type
*/
template<size_t dim>
const FiniteVolumeStencil<dim>* FiniteVolumeStencilManager<dim>::Stencil( CSMP_FEM_TYPE etype ) const
 {
    std::map<CSMP_FEM_TYPE,size_t>::const_iterator  it = type_mapping.find(etype);
    if ( it == type_mapping.end() ) {
         std::cerr <<"\nFiniteVolumeStencilManager<dim>::FiniteVolumeStencil: FV Stencil for requested ";
         std::cerr << parseFiniteElementType(etype) <<" type is not available."<< std::endl;
         return nullptr;
      }
    return &stencils[ (*it).second ];
 }


template class FiniteVolumeStencilManager<1U>;
template class FiniteVolumeStencilManager<2U>;
template class FiniteVolumeStencilManager<3U>;

} // end namespace csmp
