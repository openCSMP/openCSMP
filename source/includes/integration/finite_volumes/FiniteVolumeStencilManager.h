#ifndef CSMP_FINITE_VOLUME_STENCIL_MANAGER_H
#define CSMP_FINITE_VOLUME_STENCIL_MANAGER_H

#include "FiniteElement.h"
#include "FiniteVolumeStencil.h"

namespace csmp {

class FiniteElementManager;

/**
@brief Finite element - finite volume stencils with a parametric space representation, see Paluszny et al. (2007, Geofluids).

@author SKM
@author Adriana Paluszny
@date 2001
*/
template<size_t dim>
class FiniteVolumeStencilManager {
  public:
    /// creates finite volume stencil objects for the element types supported by the current finite-element manager
    FiniteVolumeStencilManager( const FiniteElementManager& fem_manager );
    ~FiniteVolumeStencilManager();
    
    /// returns pointer to corresponding FV stencil type
    const FiniteVolumeStencil<dim>* const Stencil( CSMP_FEM_TYPE etype ) const;
  
  private:
    FiniteVolumeStencilManager() = delete;
    void Initialize( const FiniteElementManager& fem_manager );
    // currently the number of supported stencils is fixed to 7
    enum { STENCIL_POLYTYPES = 7 };
    FiniteVolumeStencil<dim>        stencils[STENCIL_POLYTYPES];
    std::map<CSMP_FEM_TYPE,size_t>  type_mapping;
};

} // end namespace csmp

#endif
