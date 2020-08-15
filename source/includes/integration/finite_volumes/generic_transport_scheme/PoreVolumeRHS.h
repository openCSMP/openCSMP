//
//  PoreVolumeRHS.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_PORE_VOLUME_RHS_H
#define CSMP_PORE_VOLUME_RHS_H

#include "VectorOperator.h"

namespace csmp {

/**
   @brief Accumulates the pore-volume - advected variable product into the 
   right-hand vector as divided by the time increment.
   
   This method will always be the first ADD_ACCUMULATE task for the righthandside.
   
   @author SKM
   @date 6/12/2016
*/
template<size_t dim>
class PoreVolumeRHS : public VectorOperator<dim> {
  public:
    PoreVolumeRHS( const csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& spv_key,
                   const csmp::INDEX<SCALAR,NODE>& fpv_key,
                   const csmp::INDEX<SCALAR,NODE>& transported_variable_key );
                   
    virtual ~PoreVolumeRHS() {}
  
    virtual void AccumulateFiniteVolume( const Node<dim>&, std::vector<double64>& right_hand_vector ) const;

    virtual void AccumulateStencil( const Element<dim>&, std::vector<double64>& right_hand_vector ) const;
    // virtual void AccumulateStencil( const Face<dim>&, std::vector<double64>& right_hand_vector ) const;
    // virtual void AccumulateStencil( const InterFace<dim>&, std::vector<double64>& right_hand_vector ) const;
  
  private:
    const INDEX<SCALAR,NODE>&                      fpv_key_;  ///< the pore volume of the FV as integrated over the sectors surrounding the node
    const INDEX<SCALAR,SECTOR_INTEGRATION_POINT>&  spv_key_;  ///< the pore volume of each finite volume sector
    const INDEX<SCALAR,NODE>&                      adv_key_;  ///< the transported variable
};


} // end csmp

#endif /* PORE_VOLUME_RHS_H */
