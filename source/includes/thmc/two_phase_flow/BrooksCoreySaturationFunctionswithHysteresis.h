#ifndef CSMP_BROOKSCOREYSATURATIONFUNCTIONSWITHHYSTERESIS_H
#define CSMP_BROOKSCOREYSATURATIONFUNCTIONSWITHHYSTERESIS_H


#include "CSMP_definitions.h"

namespace csmp {

/**
@brief The Brooks-Corey capillary pressure model and the Corey-Burdine relative permeability and their first derivatives
 @attention relies on sw(sH2O), swr(srH2O), sor(srCO2) and so(sCO2) variable values stored by the USER.
 @attention needs the Brooks-Corey constants (a's and c's)
 
 @note The Brooks-Corey capillary pressure model and the Corey-Burdine relative permeability has been implemented from the Skjaeveland et al. 2000
 
  Svein M. Skjaeveland, L. M. Siqveland, A. Kjosavik, W. L. Hammervold and G. A. Virnovsky, (2000).
  Capillary Pressure Correlation for Mixed-Wet Reservoirs. SPE Reservoir Evaluation & Engineering. 3. 10.2118/39497-MS.
 
@date 2018
*/

template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctionswithHysteresis {
  public:
  
    BrooksCoreySaturationFunctionswithHysteresis();
  
    BrooksCoreySaturationFunctionswithHysteresis( double64 aw,  double64 ao,
                                  double64 cw,  double64 co) ;
  
    template<class TARGET_PLACEMENT>
    double64 Pc( TARGET_PLACEMENT& );               // The Capillary pressure Eq. (2) from Skjaeveland et al. 2000
  
    template<class TARGET_PLACEMENT>
    double64 EstimateOilEntryPressure( TARGET_PLACEMENT& ); // This is the same as "co" variable,Eq. 3
                                                // from Skjaeveland et al. 2000,
                                                // which can be calculated from the intersection
                                                // point of saturation axis with capillary pressure function.
  
    template<class TARGET_PLACEMENT>
    double64 dpcds( TARGET_PLACEMENT& );             // The first derivative of Capillary pressure.
  
                                                       // Hysteresis Loop Logic
    template<class TARGET_PLACEMENT>
    double64 OilResidualSaturation( TARGET_PLACEMENT& );   // The oil residual saturation estimated from Land's formula.
  
    // The relative Permeability and first derivatives are evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 krw( TARGET_PLACEMENT& ) ;
    
    template<class TARGET_PLACEMENT>
    double64 krn( TARGET_PLACEMENT& ) ;
    
    template<class TARGET_PLACEMENT>
    double64 dkrwds( TARGET_PLACEMENT& ) ;
    
    template<class TARGET_PLACEMENT>
    double64 dkrnds( TARGET_PLACEMENT& ) ;
  
    template<class TARGET_PLACEMENT>
    double64 WaterFractionalFlow( TARGET_PLACEMENT& , double64 S) ;
  
    template<class TARGET_PLACEMENT>
    double64 FirstDerivativeWaterFractionalFlow( TARGET_PLACEMENT& , double64 S) ;
  
    template<class TARGET_PLACEMENT>
    double64 OilFractionalFlow( TARGET_PLACEMENT& , double64 S) ;
  
    template<class TARGET_PLACEMENT>
    double64 InflectionPointSaturation( TARGET_PLACEMENT& ) ;
  
    template<class TARGET_PLACEMENT>
    double64 TangentPointSaturation( TARGET_PLACEMENT& ) ;
  
    template<class TARGET_PLACEMENT>
    double64 ShockFrontVelocity( TARGET_PLACEMENT& ) ;
  
    template<class TARGET_PLACEMENT>
    double64 BuckleyLeverettFunction( TARGET_PLACEMENT& , double64 S) ;
  
    template<class TARGET_PLACEMENT>
    double64 TheSecantMethod( TARGET_PLACEMENT& , double64 S1, double64 S2) ;
  
  
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    
    double64 aw_;
    double64 ao_;
    double64 cw_;
    double64 co_;
    double64 const C_land_ = 4.0;
 
};


}  // end namespace csmp

#endif /* BrooksCoreySaturationFunctionswithHysteresis_h */
