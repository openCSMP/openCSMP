#ifndef POINT_SOURCE_RHSOP_H
#define POINT_SOURCE_RHSOP_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/** @brief Assignment of scalar nodal (Neumann) source or sink terms.

    This operator applies (absolute) nodal source terms to the RHS vector as is. 
    These values must be precomputed elsewhere. One possibility to do this
    is via the method

    @code
      NodeCenteredFiniteVolumeTranport::TransformScalarBoundaryValuesIntoNeumannConditions()
    @endcode

    It computes gradients normal to the outer model boundary from influxes (sources)
    assigned and flagged as NEUMANN elsewhere. To do this it needs material parameters
    from the corresponding flow law, for instance the hydraulic conductivity
    in the case of Darcy's law.
    
    @attention when applying this term, remember that model must have 
    at least one Dirichlet constraint assigned to make solution unique.

@date 2000 
@author Stephan K. Matthaei

*/      
template<size_t dim,class SIMPLEX=Element<dim> >
class PointSource_rhsop : public MathOperatorRHS<dim> {
  public:
    PointSource_rhsop( const PropertyDatabase<dim>&, const char* nodal_src, const char* test );
    virtual void GetOperands( SIMPLEX& );
    virtual void ComputeContribution( SIMPLEX& );
    virtual PointSource_rhsop<dim,SIMPLEX>* clone() const { return new PointSource_rhsop<dim,SIMPLEX> (*this); }
  private:
    std::vector<ScalarVariable>  SRC_;
};

} // csmp

#endif
















