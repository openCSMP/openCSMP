#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H


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

enum phase {H2O, CO2} ;
enum process {Drainage, Imbibition} ;

template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctionsWithHysteresis {
  public:
  
    BrooksCoreySaturationFunctionsWithHysteresis();
  
    BrooksCoreySaturationFunctionsWithHysteresis( double64 aw,  double64 ao,
                                  double64 cw,  double64 co) ;
  
    template<class TARGET_PLACEMENT>
    void InitialiseBrooksCoreyParameters( TARGET_PLACEMENT& ) ;
  
    template<class TARGET_PLACEMENT>
    void InitialiseBrooksCoreyParameters( TARGET_PLACEMENT& , double64 aw, double64 ao,  double64 cw, double64 co) ;
    
    template<class TARGET_PLACEMENT>
    void InitialiseBrooksCoreyParameters( TARGET_PLACEMENT& , double64 a[2][2], double64 c[2][2], process ProcessPath) ;
    
    template<class TARGET_PLACEMENT>
    double64 EffectiveSaturation( TARGET_PLACEMENT&  ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 pc( TARGET_PLACEMENT& ) const ;               // The Capillary pressure Eq. (2) from Skjaeveland et al. 2000
  
    template<class TARGET_PLACEMENT>
    double64 EstimateOilEntryPressure( TARGET_PLACEMENT& ); // This is the same as "co" variable,Eq. 3
                                                // from Skjaeveland et al. 2000,
                                                // which can be calculated from the intersection
                                                // point of saturation axis with capillary pressure function.
  
    template<class TARGET_PLACEMENT>
    double64 dpcds( TARGET_PLACEMENT& ) const ;             // The first derivative of Capillary pressure.
  
                                                       // Hysteresis Loop Logic
    template<class TARGET_PLACEMENT>
    double64 OilResidualSaturation( TARGET_PLACEMENT& ) const ;   // The oil residual saturation estimated from Land's formula.
    
    
    template<class TARGET_PLACEMENT>
    double64 WaterResidualSaturation( TARGET_PLACEMENT& p, double64 a[2][2], double64 c[2][2]) const ; // The Water residual saturation estimated from the intersection of capillary pressures from Dranage and Imbibition curves in Imbibition process.
    
    template<class TARGET_PLACEMENT>
    void WaterAndOilResidualSaturationImbibitionToDrainage( TARGET_PLACEMENT&, double64 Swr_, double64 Sor_) ;  // Estimate the residual of water and oil for the target saturations
                                                                                                                  // for the Dranaige process from the Imbibitions curve.
    template<class TARGET_PLACEMENT>
    void WaterAndOilResidualSaturationDrainageToImbibition( TARGET_PLACEMENT&, double64 Swr_, double64 Sor_) ; // Estimate the residual of oil and water for the target saturations
                                                                                                                 // for the Imbibitions process from the Drainage curve.
  
    // The relative Permeability and first derivatives are evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    template<class TARGET_PLACEMENT>
    double64 krw( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 krw( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 krn( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 krn( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 dkrwds( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 dkrwds( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 dkrnds( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 dkrnds( TARGET_PLACEMENT& , double64 S) const ;

    template<class TARGET_PLACEMENT>
    double64 WaterFractionalFlow( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 FirstDerivativeWaterFractionalFlow( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 OilFractionalFlow( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 InflectionPointSaturation( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 TangentPointSaturation( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 ShockFrontVelocity( TARGET_PLACEMENT& ) const ;
  
    template<class TARGET_PLACEMENT>
    double64 BuckleyLeverettFunction( TARGET_PLACEMENT& , double64 S) const ;
  
    template<class TARGET_PLACEMENT>
    double64 TheSecantMethod( TARGET_PLACEMENT& , double64 S1, double64 S2) const ;
  
  // Mahyar: These are extra prameters has to be defined in the Brooks Corey with Hysteresis in the two phase model.
    double64  aw_ ;  //Brooks Corey parameter of water
    double64  ao_ ;  //Brooks Corey parameter of oil(CO2)
    double64  cw_ ;  //entry pressure of water
    double64  co_ ;  //entry pressure of oil(CO2)
    double64 const C_land_ = 2.0; //Land's parameter
  
  
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
  double64 Swmin ;
  double64 Swmax ;
 
};


}  // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H */
