#ifndef CSMP_SATURATION_FUNCTION_H
#define CSMP_SATURATION_FUNCTION_H

#include "CSMP_definitions.h"

namespace csmp {

/** base class for 2-phase flow models, excluding any fluid related terms

base class for 2-phase flow models (Linear, Brooks Corey, Van Genuchten, Richards)
now you can inherit BrooksCoreyWetting, BrooksCoreyNonWetting... 

@todo (3) To improve calculation speed, change relperm models from dynamic to static polymorphism (C)
 
 */
template<size_t dim, template<size_t> class USER>
class SaturationFunction {
  public:
    SaturationFunction() : pc_max_(1e7), dpcds_max_(1e6) {}
  
    /// always of the wetting phase by convention
    double64  EffectiveSaturation( double64 sw ) const;
    double64  SeffToSw( double64 seff ) const;

    /// relative permeabilities
    double64 krw( double64 sw ) const;
    double64 krn( double64 sw ) const;

    /// derivatives of relative permeabilities
    double64 dkrwds( double64 sw, bool evaluate_numerically=true ) const;
    double64 dkrnds( double64 sw, bool evaluate_numerically=true ) const;

    /// capillary pressure (limit this to 4e7, the max strength of the rock)
    /// do this by computing seff for which pc=4e7, then use this seff as
    /// a limiting value @attention absolute saturation is used
    double64 pc( double64 sw ) const;

    /// maximum value of pc
    double64 MaxCapillaryPressure() const;

     /// inverse capillary pressure function
    double64 Sw( double64 pc ) const;

    /// maximum value of dpcdS
    double64 MaxCapillaryPressureDerivative() const;

   /// capillary pressure derivatives; @attention always use sw and not seff
    double64 dpcds( double64 sw, bool evaluate_numerically=true ) const;

    /// derivatives of inverse capillary pressure function
    double64 dsdpc( double64 pc, bool evaluate_numerically=true ) const;
  
    /// default is the wetting phase
    void Out() const;
    
  protected:
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    /// Numerical derivatives
    double64 dkrwds_Numerical( double64 sw, double64 h = 0.001 ) const;
    double64 dkrnds_Numerical( double64 sw, double64 h = 0.001 ) const;
    double64 dpcds_Numerical(  double64 sw, double64 h = 0.00001 ) const;
  
  private:
    const double64 pc_max_, dpcds_max_;
  
};

} // end namespace csmp

#endif 





















