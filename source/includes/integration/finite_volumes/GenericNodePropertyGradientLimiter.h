#ifndef GENERIC_NODE_PROPERTY_GRADIENT_LIMITER_H
#define GENERIC_NODE_PROPERTY_GRADIENT_LIMITER_H

#include "GenericNodePropertyGradient.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "FiniteVolumeStencilManager.h"
#include "FiniteVolumeStencil.h"
#include "PropertyHandle.h"
#include "ErrorHandler.h"

namespace csmp {

/** 

@brief GenericNodePropertyGradientLimiter - limits the nodal (finite-volume) gradient 
of target variable.

If a property (saturation, solute concentration, etc...) is advected with second
order accuarcy, some sort of slope limiting is required to avoid spurious oscillations. This
objects computes the slope limiting factor phi using the MINMOD limiter for the nodal gradients
that were computed employing the object  NodePropertyGradient. This procedure
is usually carreid out automatically int the various FiniteVolume<fT, dim>Visitor objects.
The resulting gradient limiter is stored automatically in the CSP_PropertyDatabase
using the CSP_Operands.

The FiniteVolumeManager<fT, dim> and the FiniteVolume<fT, dim> classes are participants.

The GenericNodePropertyGradientLimiter object collaborates with the FiniteVolume<fT, dim>AdvectionVisitor,
the FiniteVolume<fT, dim>GradientVisitor, and the FiniteVolume<fT, dim>FluidPhaseVisitor

A single GenericNodePropertyGradientLimiter object is automatically constructed when
one of the above-mentioned visitors is constructed. The GenericNodePropertyGradientLimiter
object should not be constructed from the main() file

copyright (c) 2001 by Sebastian Geiger, Stephan K. Matthaei & Stephen G. Roberts 

*/
template<size_t dim>
class GenericNodePropertyGradientLimiter {
       
  public:
    GenericNodePropertyGradientLimiter( Model<dim>&, const char* region, const char* prop ="node property gradient" );
    GenericNodePropertyGradientLimiter( Model<dim>&, const char* region, std::vector<char*> prop_names );
  
    ~GenericNodePropertyGradientLimiter(); 
    
    // computation of property gradient
    void CalculateGenericNodalGradient();

    // compute the gradient limiter
    void CalculateSlopeLimiter( const std::vector<std::pair<double64,double64> >& MINMAX, int counter = 1 );

    // setting the property key
    void SetPropertyKey( csmp::Index& key );
    
    // returning the size of the object
    double64 SizeOf() const;

    
  private:
    Region<dim>                        gref_; ///< handle to the application domain of the limiter
    std::vector<PropertyHandle<dim>*>  limiter;
    GenericNodePropertyGradient<dim>   node_prop_grad;
    double64                           tolerance;
    csmp::Index                        u_key;
    const csmp::Index                  grad_key, lim_key, mctr_key;
 };

}

#endif
