#ifndef CSMP_FINITE_VOLUME_STENCIL_MANAGER_H
#define CSMP_FINITE_VOLUME_STENCIL_MANAGER_H

#include "FiniteElement.h"
#include "FiniteVolumeStencil.h"

namespace csmp {

class FiniteElementManager;

template<size_t dim>
class FiniteVolumeStencilManager {
  public:
    FiniteVolumeStencilManager();
    FiniteVolumeStencilManager( const FiniteElementManager& fem_manager );
    ~FiniteVolumeStencilManager();
    void Initialize( const FiniteElementManager& fem_manager );
    
    // returns pointer to corresponding FV stencil type
    const FiniteVolumeStencil<dim>* Stencil( CSMP_FEM_TYPE etype ) const;
    
  private:
    // there is a fixed number of stencils as there is a fixed number of element types
    enum { STENCIL_POLYTYPES = 7 };
    FiniteVolumeStencil<dim>        stencils[STENCIL_POLYTYPES];
    std::map<CSMP_FEM_TYPE,size_t>  type_mapping;
};


// returns pointer to corresponding FV stencil type
template<size_t dim>
inline  const FiniteVolumeStencil<dim>* FiniteVolumeStencilManager<dim>::Stencil( CSMP_FEM_TYPE etype ) const
 {
    std::map<CSMP_FEM_TYPE,size_t>::const_iterator  it = type_mapping.find(etype);
    if ( it == type_mapping.end() ) {
         std::cout <<"\nFiniteVolumeStencilManager<dim>::FiniteVolumeStencil: FV Stencil for requested ";
         std::cout << parseFiniteElementType(etype) <<" type is not available."<< std::endl;
         return NULL;
      }
    return &stencils[ (*it).second ];
 }

} // end namespace csmp

#endif
