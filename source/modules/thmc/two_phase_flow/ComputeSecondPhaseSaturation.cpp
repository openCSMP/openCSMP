// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ComputeSecondPhaseSaturation.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ComputeSecondPhaseSaturation<dim>::ComputeSecondPhaseSaturation( const PropertyDatabase<dim>& p,
                                              const char* reference_saturation, const char* second_saturation ) 
      : Interrelation<dim>(p),
        ref_sat( Interrelation<dim>::GlobalProperty(reference_saturation) ),
        new_sat( Interrelation<dim>::GlobalProperty(second_saturation) )
        
 {
    Interrelation<dim>::OutputCondition( new_sat, PLAIN );
    Interrelation<dim>::ResultProperty(second_saturation);
 }


/// The saturation of the first phase is adjusted such that the
/// two saturations add up to one
template<uint32_t dim>
void ComputeSecondPhaseSaturation<dim>::Calculate()
 {
    ref_sat.AssignTo( sat );
    
    new_sat = 1.0 - sat();
         
 } // end 

template class ComputeSecondPhaseSaturation<1U>;
template class ComputeSecondPhaseSaturation<2U>;
template class ComputeSecondPhaseSaturation<3U>;

}
