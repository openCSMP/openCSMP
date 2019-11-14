#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"


namespace csmp {

template<size_t> class Element;
template<size_t> class PropertyDatabase;

  /**
   @brief The Brooks-Corey capillary pressure model and the Corey-Burdine relative permeability and their first derivatives.
   
   @section Introduction
   
   This class contains modelling the hysteresis relative permeability curves for mixed wet 2 phase model.
   
   Firstly, we model relative permeability curves suggested by Skjaeveland et al.(2000). They introduce a capillary pressure model for
   mixed- wet reservoir rock based on the Brooks-Corey capillary pressure model.
   By using these new correlations, they show the ability of modelling the Hysteresis loop logic. Also, using the Land's empirical equations
   one can estimate the residual's saturation of Non-wet phase (oil, or CO2), which can create first Imbibition curve etc..
   We follow the same logic presented in the paper (page 61, and 62 in the paper), and implement the algorithm to do Drainage and Imbibition process.
   
   In this class, we also calculate the relative permeability of the wet, and non-wet phase.
   
   In the modelling, we check the capillary pressure of any process (Imbibitions, or Drainage) to be bound in the main Drainage and Imbibitions.
   Moreover, also all other functions, also calculated from the proper state of phases, either they are Imbibing or Draining.
   
   The input parameters are, "aw, ao, cw, co" parameters, and residual saturation of wet and non-wet (oil, CO2) phase. We have these parameters
   as the property of elements...They can be read from configuration file as:
   
   with key(key_kri_param) as described as "Hysteresis relative permeability model Parameters".
   
   They are ordered as: HYSTERIC_MODEL_PARAMETERS { AWD=0, AOD=1, CWD=2, COD=3, AWI=4, AOI=5, CWI=6, COI=7} ;
   
   The results are relative permeabilities, derivatives... which they will be used in the FlowFunctions class.
   
   @attention: To have Primary Drainage, we need just use non-autamatic and use the forced intialize by set the :
   
   c_[H2O][DRAINAGE] = entery pressure
   c_[CO2][DRAINAGE] = 0
   
   and set the ProcessPath= DRAINAGE
   
   @attention relies on sw(sH2O), swr(srH2O), sor(srCO2) and so(sCO2) variable values stored by the USER.
   @attention needs the Brooks-Corey constants (a's and c's)
   
   @note The Brooks-Corey capillary pressure model and the Corey-Burdine relative permeability has been implemented from the Skjaeveland et al. 2000.
   
   @attention current implementation is suitable only for slightly compressible flow computations, but could be extended to compositional modelling.

   @section  Reference
   
   Svein M. Skjaeveland, L. M. Siqveland, A. Kjosavik, W. L. Hammervold and G. A. Virnovsky, (2000).
   Capillary Pressure Correlation for Mixed-Wet Reservoirs. SPE Reservoir Evaluation & Engineering. 3. 10.2118/39497-MS.
   
   @author Mahyar Madadi
   @author Stephan Matthai
   @date May, 2019
   
*/
template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctionsWithHysteresis {
  public:
    enum FLUID_PHASE {H2O, CO2};
    enum TWO_PHASE_FLOW_PROCESS { DRAINAGE, IMBIBITION };
    enum HYSTERIC_MODEL_PARAMETERS { AWD=0, AOD=1, CWD=2, COD=3, AWI=4, AOI=5, CWI=6, COI=7 };
  
  public:
    /// accesses the model to create or attach to associated variable storage
    BrooksCoreySaturationFunctionsWithHysteresis( PropertyDatabase<dim>& );
  
    /// Effective saturation function
    double64 EffectiveSaturation( Element<dim>* const ) const;
  
    /// Effective saturation function for saturation S
    double64 EffectiveSaturation_at( Element<dim>* const, double64 s1 ) const;
  
    /// Capillary pressure Eq. (2) from Skjaeveland et al. 2000
    double64 pc( Element<dim>* const ) const;
  
    /// Capillary pressure Eq. (2) from Skjaeveland et al. 2000 for saturation S
    double64 pc_at( Element<dim>* const , double64 s1 ) const;
  
    /// First derivative of capillary pressure as a function of saturation
    double64 dpcds( Element<dim>* const ) const;
  
    double64 dpcds_at( Element<dim>* const, double64 ) const;
  
    /// The water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 krw( Element<dim>* const ) const ;

    double64 krw_at( Element<dim>* const , double64 S) const ;
  
    /// The CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 krn( Element<dim>* const ) const ;

    double64 krn_at( Element<dim>* const, double64 S) const ;
  
    /// The first relative of water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 dkrwds( Element<dim>* const ) const ;

    double64 dkrwds_at( Element<dim>* const, double64 S) const ;
  
    /// The first relative of CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
    double64 dkrnds( Element<dim>* const ) const ;

    double64 dkrnds_at( Element<dim>* const, double64 S) const ;
  
    /// Numerical derivatives of first derivatives of relative permeability of water and CO2
    double64 dkrwds_Numerical( Element<dim>* const, double64 delta_s ) const ;

    double64 dkrwds_at_Numerical( Element<dim>* const, double64 sw, double64 delta_s ) const ;

    double64 dkrnds_Numerical( Element<dim>* const, double64 delta_s ) const ;

    double64 dkrnds_at_Numerical( Element<dim>* const, double64 sw, double64 delta_s ) const ;
  
    // -----------------------
    // NON-STANDARD INTERFACES
    // -----------------------
  
    /// Oil residual saturation estimated from Land's formula.
    double64 OilResidualSaturation( Element<dim>* const ) const;
  
    /// Water residual saturation estimated from the intersection of capillary pressures from Dranage and Imbibition curves in Imbibition process.
    double64 WaterResidualSaturation( Element<dim>* const, double64 Sor ) const;

    /// Estimate the residual of water and oil for the target saturations for the Dranaige process from the Imbibitions curve.
    void WaterAndOilResidualSaturationImbibitionToDrainage( Element<dim>* const e, double64& Swr, double64& Sor )  const;

  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    void SetBrooksCoreyCurvesParameters( const Element<dim>* const e, std::array<double64, 2>& a,std::array<double64, 2>& c ) const ;
  
    void UpdatePseudoResidualAndEndpointSaturations( Element<dim>* ) const;
  
    /// print out the two phase flow process state based on the current and new CO2 Saturation functions for each element based on the  interpolation of saturations at the bary center.
    TWO_PHASE_FLOW_PROCESS FlowProcess( const Element<dim>* const ) const;
  
    /// print out the capillary function bounds upper and lowwer limits..
    std::pair<double64,double64> CapillaryPressureLimits( Element<dim>* const ) const;
  
    /// Check the capillary pressure is in the limits and reset the pseduo residual water (newSrH2O) and CO2 saturations (newSrCO2)
    void CheckPcLimitsAndResetResiduals( Element<dim>* const, double64& newSrH2O, double64& newSrCO2, std::array<double64, 2>& a,std::array<double64, 2>& c) const;
  
    /// Land's parameter for Air and CO2 "Prather Bray Seymour Codd 2016" paper
    const double64  C_land_ = 0.89;
  
    /// Maximum Capillary pressure for the Sw < Srw
    const double64  MaxCapillaryPressure = 1e7; // Pa
    const double64  MaxCapillaryPressureDerivative = 1e6; //Pa
  
    /// print out the water (wet phase) relative permeability function bounds upper and lowwer limits..
    std::pair<double64,double64> krwLimits_at( Element<dim>* const , double64 const ) const;
  
    /// print out the co2 (non-wet phase) relative permeability function bounds upper and lowwer limits..
    std::pair<double64,double64> krnLimits_at( Element<dim>* const , double64 const ) const;
  
  private:
    csmp::Index key_SwImbToDr_, key_SwDrToImb_, key_prsH2O_, key_prsCO2_ ;
    mutable ArrayVariable  ac_params_; ///< extra parameters that are used only in this relative permeability model
};

  
}  // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H */
