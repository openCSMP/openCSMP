//
//  PoreVolumeRHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "PoreVolumeRHS.h"
#include "Model.h"
#include "Node.h"

namespace csmp {

template<size_t dim>
PoreVolumeRHS<dim>::PoreVolumeRHS( const Model<dim>& model,
                                   const char* pv_variable, const char* advected_variable )
 :  VectorOperator<dim>(ADD_ACCUMULATE),
    pv_key_(model.Database().StorageKey(pv_variable)),
    adv_key_(model.Database().StorageKey(advected_variable))
{
 }
  


/**
    Expects that the facet fluxes are uptodate as precomputed before by the FluxEvaluator.
    Furthermore this assumes that the node indexes have been updated for the computational
    domain.
*/
template<size_t dim>
void PoreVolumeRHS<dim>::AccumulateFV( const Node<dim>* nptr, std::vector<double64>& rhs ) const
 {
     const double64 pore_volume(nptr->Read(pv_key_));
     const double64 advection_value(nptr->Read(adv_key_));
   
   if (this->multiply_with_dt_) {
     rhs[ nptr->Idx() ] += pore_volume * advection_value * this->dt_;
   }
   else {
     rhs[ nptr->Idx() ] += pore_volume * advection_value;
   }
 }


template class PoreVolumeRHS<1U>;
template class PoreVolumeRHS<2U>;
template class PoreVolumeRHS<3U>;

} // end csmp
