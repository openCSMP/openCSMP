//
//  PoreVolumeRHS.h
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 6/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef PORE_VOLUME_RHS_H
#define PORE_VOLUME_RHS_H

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
    PoreVolumeRHS( const Model<dim>& model,
                   const char* pv_variable, const char* advected_variable,
                   double64 time_increment );
  
    virtual void AccumulateFV( const Node<dim>*, std::vector<double64>& right_hand_vector ) const;
  
  private:
    const csmp::Index  pv_key_;   ///< the pore volume of the FV as integrated over the sectors surrounding the node
    const csmp::Index  adv_key_;  ///< the transported variable
    double64           dt_;       ///< time increment (over which the transport will be computed)
};


} // end csmp

#endif /* PORE_VOLUME_RHS_H */
