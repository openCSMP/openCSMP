#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H

#include "TwoPhaseModel.h"

namespace csmp {

/** @brief saturation function policy based on the Brooks-Corey (1964) model

    @attention relies on sw, swr, snr, variable values stored by the USER. 
    @attention phases are numbered 0..2 (water=0)

    @note for lambda=0, this implementation of Brooks-Corey model switches to linear
    @note for linear case capillary pressure is a constant value equal to entry pressure
    @note base class pm1 and pm2 are used for pd and lambda, respectively
*/
template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctions {
  public:
    BrooksCoreySaturationFunctions() = delete;
     ~BrooksCoreySaturationFunctions() = delete;
    
    /// always of the wetting phase by convention
    template<class TARGET_PLACEMENT>
    double64  EffectiveSaturation( TARGET_PLACEMENT& ) const;

    double64  SeffToSw( double64 seff ) const;

    /// relative permeabilities as a function of water saturation - parameters come from subclass FlowFunctions
    template<class TARGET_PLACEMENT>
    double64 krw( TARGET_PLACEMENT& ) const;

    template<class TARGET_PLACEMENT>
    double64 krn( TARGET_PLACEMENT& ) const;

    /// derivatives of relative permeabilities
    template<class TARGET_PLACEMENT>
    double64 dkrnds( TARGET_PLACEMENT& ) const;

    template<class TARGET_PLACEMENT>
    double64 dkrwds( TARGET_PLACEMENT& ) const;

    /// capillary pressure
    template<class TARGET_PLACEMENT>
    double64 pc( TARGET_PLACEMENT& ) const;

    /// maximum value of pc
    double64 MaxCapillaryPressure() const { return 1e7; /* Pa */ }

    /// inverse capillary pressure function
    template<class TARGET_PLACEMENT>
    double64 Sw( TARGET_PLACEMENT&, double64 pc ) const;

    /// maximum value of dpcdS
    double64 MaxCapillaryPressureDerivative() const { return 1e6; /* Pa m-1 */ }

    /// capillary pressure derivatives
    template<class TARGET_PLACEMENT>
    double64 dpcds( TARGET_PLACEMENT& ) const;

    /// inverse capillary pressure derivative
    template<class TARGET_PLACEMENT>
    double64 dsdpc( TARGET_PLACEMENT&, double64 pc ) const;

    template<class TARGET_PLACEMENT>
    void Out( TARGET_PLACEMENT& ) const;

  private:
    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
};

} // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H */
