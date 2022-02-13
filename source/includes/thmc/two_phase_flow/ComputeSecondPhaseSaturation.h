#ifndef ComputeSecondPhaseSaturation_h__
#define ComputeSecondPhaseSaturation_h__

#include "CSMP_definitions.h"
#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class ComputeSecondPhaseSaturation : public Interrelation<dim> {
    Operand<dim>&  ref_sat;  
    Operand<dim>&  new_sat;  
    ScalarVariable  sat;
    
  public:
    ComputeSecondPhaseSaturation( const PropertyDatabase<dim>& p,
                           const char* reference_saturation, const char* new_saturation ); 
    ~ComputeSecondPhaseSaturation() {};
    void Calculate();
};

/**
 
@class ComputeSecondPhaseSaturation  ComputeSecondPhaseSaturation "two_phase_flow/ComputeSecondPhaseSaturation.h"
@author S.K. Matthaei
@author S. Geiger
@author S. Roberts
@date 2001

 
@section motivation Motivation
 
Since only one phase is propagated in the FiniteVolumeFluidPhaseVisitor,
the saturation of the second phase must be computed afterwards, which
is done by this simple interrelation, such that all saturation dependend
properties can be computed afterwards. The user must define
in the constructor of this object the saturation that was computed
in the FiniteVolumeFluidPhaseVisitor (i.e., reference saturation) and
the saturation which shall be computed by this interrelation (i.e., new oder
second phase saturation) by subtracting the reference saturation from 1.0.
 

Of course, the saturation of the second phase can also be computed by
the FiniteVolumeFluidPhaseVisitor, which is, however, computationally
more expensive than applying this interrelation.

 
@section implementation Implementation
 
The user must provide the character strings of the reference saturation
(first character string in the constructor) and the saturation to be 
computed (second character string in the constructor) to the constructor
of this interrelation. 

As part of the Interrelation base class, this interrelation is passed
to the Model in the main() program. It should be done during the transient loop, such
that the saturations are both constantly updated. Some functions of the
FiniteVolumeFluidPhaseVisitor also require a reference to this interrelation,
since the saturation of the second phase must be updated within the time-loop
constantly. See the CSMP documentation for more detail on the Interrelation 
base class and their derived classes and the FiniteVolumeFluidPhaseVisitor 
documentation for details on the different member functions.

 
@section examples Application Examples

@code
// oil is propagated, computing the water saturation <br>
// when the interrelation is passed to the Model <br>
ComputeSecondPhaseSaturation second_phase( p_ref, "saturation oil", "saturation water" );  
@endcode
 */
}

#endif

