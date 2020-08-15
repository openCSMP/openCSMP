//
//  PoreVolumeRHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "PoreVolumeRHS.h"
#include "Element.h"
#include "Node.h"

namespace csmp {

template<size_t dim>
PoreVolumeRHS<dim>::PoreVolumeRHS( const INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& spv_key,
                                   const INDEX<SCALAR,NODE>& fpv_key,
                                   const INDEX<SCALAR,NODE>& transported_variable_key )
   :  spv_key_(spv_key),
      fpv_key_(fpv_key),
      adv_key_(transported_variable_key)
  {
  }
  

/**
 Expects that the facet fluxes are uptodate as precomputed before by the FluxEvaluator.
 Furthermore this assumes that the node indexes have been updated for the computational
 domain.
 */
template<size_t dim>
void PoreVolumeRHS<dim>::AccumulateStencil( const Element<dim>& e, std::vector<double64>& rhs ) const
 {
    // assumptions
    assert( e.IntegrationPointsPerSector() == 1U );
    assert( e.Sectors() == e.Nodes() );

    const size_t sector_ipoints(e.Sectors());

    for (size_t i=0U; i < sector_ipoints; ++i ) 
      {
         const double64 sector_pore_volume = e.Read( i, 0U, spv_key_ );
         const double64 advected_var_value = e.N(i)->Read( adv_key_ );
        
         rhs[ e.N(i)->Idx() ] += sector_pore_volume * advected_var_value * this->Factor();
      }
  }
  
  


/**
    Expects that the facet fluxes are uptodate as precomputed before by the FluxEvaluator.
    Furthermore this assumes that the node indexes have been updated for the computational
    domain.
*/
template<size_t dim>
void PoreVolumeRHS<dim>::AccumulateFiniteVolume( const Node<dim>& fv, std::vector<double64>& rhs ) const
 {
     const double64 pore_volume(fv.Read(fpv_key_));
     const double64 advected_var_value(fv.Read(adv_key_));
   
     rhs[ fv.Idx() ] += pore_volume * advected_var_value * this->Factor();
 }


template class PoreVolumeRHS<1U>;
template class PoreVolumeRHS<2U>;
template class PoreVolumeRHS<3U>;

} // end csmp
