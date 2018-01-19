#ifndef CSMP_TWO_PHASE_MODEL_H
#define CSMP_TWO_PHASE_MODEL_H

#include "Node.h"
#include "Element.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "FiniteVolumeStencil.h"
#include "TensorVariable.h"
#include "TensorVariable1.h"
#include "TensorVariable2.h"

namespace csmp {

/** base class for 2-phase flow models, excluding any fluid related terms

base class for 2-phase flow models (Linear, Brooks Corey, Van Genuchten, Richards)
now you can inherit BrooksCoreyWetting, BrooksCoreyNonWetting... 

@todo (3) To improve calculation speed, change relperm models from dynamic to static polymorphism (C)
 
 */
template<size_t dim, template<size_t> class USER>
class SaturationFunction {
  public:
  
    virtual ~SaturationFunction();

    /// maximum value of dpcdS
    double64 MaxCapillaryPressure( size_t phase = 1U ) const;

    /// always of the wetting phase by convention
    virtual double64  EffectiveSaturation() const;
    virtual double64  SeffToSw() const;
    virtual double64  SeffToSw( double64 seff) const;

    /// relative permeabilities
    virtual double64 krn() const = 0;
    virtual double64 krw() const = 0;

    /// derivatives of relative permeabilities
    virtual double64 dkrnds() const;
    virtual double64 dkrwds() const;

    /// capillary pressure (limit this to 4e7, the max strength of the rock)
    /// do this by computing seff for which pc=4e7, then use this seff as
    /// a limiting value @attention absolute saturation is used
    virtual double64 pc() const = 0;

    /// capillary pressure derivatives (treat seff as for previous function)
    virtual double64 dpcds() const = 0;

    /// inverse capillary pressure function
    virtual double64 Sw( double64 pc ) const;

    /// derivatives of inverse capillary pressure function
    virtual double64 dsdpc( double64 pc ) const;
  
    /// Numerical derivatives
    double64  krw_at( double64 ) const;
    double64  krn_at( double64 ) const;
    double64  dkrwds_at( double64 ) const;
    double64  dkrnds_at( double64 ) const;
    double64  pc_at( double64 se ) const;
    double64  dpcds_at( double64 se) const;
    double64  Sw_at( double64 se) const;
    double64  dsdpc_at( double64 se) const;

    virtual double64 dkrwds_numerical( double64 h = 0.001 ) const;
    virtual double64 dkrnds_numerical( double64 h = 0.001 ) const;
    virtual double64 dpcds_numerical(  double64 h = 0.00001 ) const;

    // TODO: check what these splines are?
    double64 spline_value( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const;
    double64 spline_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const;
    double64 spline_second_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const;

    /// default is the wetting phase
    virtual void Out() const;
    
  protected:
    SaturationFunction();
  
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
};

} // end namespace csmp

#endif 





















