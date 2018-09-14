#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"

namespace csmp {

template<size_t> class Element;
  
  
  enum FLUID_PHASE {H2O, CO2} ;
  enum TWO_PHASE_FLOW_PROCESS { DRAINAGE, IMBIBITION } ;
  enum HYSTERIC_MODEL_PARAMETERS { AWD=0, AOD=1, CWD=2, COD=3, AWI=4, AOI=5, CWI=6, COI=7} ;  /// suggested by Stephan 17/Aug/2018
  
  /**
   @brief The Brooks-Corey capillary pressure model and the Corey-Burdine relative permeability and their first derivatives.
   
   @section Introduction
   
   This class contains two main sections, firstly modelling the hysteresis relative permeability curves for mixed wet 2 phase model, and
   secondly, use the Leverett-Buckley approximation to estimate the shock wave velocity and saturation (Shock hight).
   
   Firstly, we model relative permeability curves suggested by Skjaeveland et al.(2000). They introduce a capillary pressure model for
   mixed- wet reservoir rock based on the Brooks-Corey capillary pressure model.
   By using these new correlations, they show the ability of modelling the Hysteresis loop logic. Also, using the Land's empirical equations
   one can estimate the residual's saturation of Non-wet phase (oil, or CO2), which this can create first Imbibition curve etc..
   We follow the same logic presented in the paper (page 61, and 62 in the paper), and implement the algorithm to do Drainage and Imbibition process.
   
   In this class, we also calculate the relative permeability of the wet, and non-wet phase.
   
   In the modelling, we check the capillary pressure of any process (Imbibitions, or Drainage) to be bound in the main Drainage and Imbibitions.
   Moreover, also all other functions, also calculated from the proper state of phases, either they are Imbibing or Draining.
   
   The input parameters are, "aw, ao, cw, co" parameters, and residual saturation of wet and non-wet (oil, CO2) phase. We have these parameters
   as the property of elements...They can be read from configuration file as:
   
   with key(key_HisBCparam) as discribed as "Hysteresis relative permeability model Parameters".
   
   They are ordered as: HYSTERIC_MODEL_PARAMETERS { AWD=0, AOD=1, CWD=2, COD=3, AWI=4, AOI=5, CWI=6, COI=7} ;
   
   They will be read automatically by using the "InitialiseBrooksCoreyParameters( const Element<dim>* const )". There are two other method to read input parameters:
   
   a) Force to have a process as either Drainage or Imbibitions:                                "InitialiseBrooksCoreyParameters( const Element<dim>* const, a,  c , ProcessPath)"
   b) Alternatively, Unforced and the results from modelling will be taking into account:       "InitialiseBrooksCoreyParameters( const Element<dim>* const, a, c )"
   
   The results are relative permeabilities, derivatives... which they will be used in the FlowFunctions class.
   
   Moreover, finally, we use these relative permeability functions to model shock velocity and shock hight.
   
   
   @attention: To have Primary Drainage, we need just use non-autamatic and use the forced intialize by set the :
   
   
   c_[H2O][DRAINAGE] = entery pressure
   c_[CO2][DRAINAGE] = 0
   
   and set the ProcessPath= DRAINAGE
   
   
   @attention relies on sw(sH2O), swr(srH2O), sor(srCO2) and so(sCO2) variable values stored by the USER.
   @attention needs the Brooks-Corey constants (a's and c's)
   
   @note The Brooks-Corey capillary pressure model and the Corey-Burdine relative permeability has been implemented from the Skjaeveland et al. 2000
 
   @note  Initialising the input parameter in three alternative methods :
     
     a) directly from input files.
     b) forced Process and in the fly, this can be used for Primary Drainage and other examples which we force the process to be either Drainage or Imbibition.
     c) unforced Process and in the fly... This can help to change the parameters if needed in the process.  For example for any change of rock type during modelling
        or input from rocktype parameter file.

   
   @section  Reference
   
   Svein M. Skjaeveland, L. M. Siqveland, A. Kjosavik, W. L. Hammervold and G. A. Virnovsky, (2000).
   Capillary Pressure Correlation for Mixed-Wet Reservoirs. SPE Reservoir Evaluation & Engineering. 3. 10.2118/39497-MS.
   
   @author Mahyar Madadi
   @date 2018
   
*/
template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctionsWithHysteresis {
  public:
    BrooksCoreySaturationFunctionsWithHysteresis();
  
    /// Initialision  method (a)
        void InitialiseBrooksCoreyParameters( const Element<dim>* const );
    
    /// Initialision  method (b)
        void InitialiseBrooksCoreyParameters( const Element<dim>* const, std::array<std::array<double64, 2>,2> a,
                                              std::array<std::array<double64, 2>,2> c , TWO_PHASE_FLOW_PROCESS ProcessPath );
  
    /// Initialision  method (c)
        void InitialiseBrooksCoreyParameters( const Element<dim>* const, std::array<std::array<double64, 2>,2> a, std::array<std::array<double64, 2>,2> c );
    
    /// Effective saturation function
        double64 EffectiveSaturation( const Element<dim>* const ) const;
  
    /// Effective saturation function fro saturation S
        double64 EffectiveSaturation_at( const Element<dim>* const, double64 s1 ) const;
  
    /// Capillary pressure Eq. (2) from Skjaeveland et al. 2000
        double64 pc( const Element<dim>* const ) const;
  
    /// First derivative of capillary pressure.
        double64 dpcds( const Element<dim>* const ) const;
  
    /// Oil residual saturation estimated from Land's formula.
        double64 OilResidualSaturation( const Element<dim>* const ) const;
  
    /// Water residual saturation estimated from the intersection of capillary pressures from Dranage and Imbibition curves in Imbibition process.
        double64 WaterResidualSaturation( const Element<dim>* const, double64 Sor_ ) const;

    /// Estimate the residual of water and oil for the target saturations for the Dranaige process from the Imbibitions curve.
        void WaterAndOilResidualSaturationImbibitionToDrainage( const Element<dim>* const e, double64& Swr_, double64& Sor_)  const;
  
  
    /// The water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
        double64 krw( const Element<dim>* const ) const ;
    
        double64 krw_at( const Element<dim>* const , double64 S) const ;
  
    /// The CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
        double64 krn( const Element<dim>* const ) const ;
    
        double64 krn_at( const Element<dim>* const, double64 S) const ;
  
    /// The first relative of water relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
        double64 dkrwds( const Element<dim>* const ) const ;
    
        double64 dkrwds_at( const Element<dim>* const, double64 S) const ;
  
    /// The first relative of CO2 relative Permeability is evaluated from the Brooks Corey Capillary Pressure model. See the Skaevland et al. 2000 at page 65
        double64 dkrnds( const Element<dim>* const ) const ;
    
        double64 dkrnds_at( const Element<dim>* const, double64 S) const ;
  
    /// print out the two phase flow process state based on the current and new CO2 Saturation functions for each element based on the  interpolation of saturations at the bary center.
        TWO_PHASE_FLOW_PROCESS FlowProcess( const Element<dim>* const );
  
    /// print out the capillary function bounds upper and lowwer limits..
        std::pair<double64,double64> CapillaryPressureLimits( const Element<dim>* const ) const;
  
    /// Check the capillary pressure is in the limits and reset the pseduo residual water (newSrH2O) and CO2 saturations (newSrCO2)
        void CheckPcLimitsAndResetResiduals( const Element<dim>* const , double64& newSrH2O , double64& newSrCO2 ) const;

    /// Land's parameter for Air and CO2 "Prather Bray Seymour Codd 2016" paper
    const double64  C_land_ = 0.89;
  
    /// Maximum Capillary pressure for the Sw < Srw
    const double64  MaxCapillaryPressure = 1e7; // Pa
    const double64  MaxCapillaryPressureDerivative = 1e6; //Pa
  
  
    /// Numerical derivatives of first derivatives of relative permeability of water and CO2
        double64 dkrwds_Numerical( const Element<dim>* const, double64 delta_s ) const ;

        double64 dkrwds_at_Numerical( const Element<dim>* const, double64 sw, double64 delta_s ) const ;

        double64 dkrnds_Numerical( const Element<dim>* const, double64 delta_s ) const ;
  
        double64 dkrnds_at_Numerical( const Element<dim>* const, double64 sw, double64 delta_s ) const ;

  
    
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
  
    mutable double64  aw_ ; ///< TODO: document what these guys are !
    mutable double64  ao_ ;
    mutable double64  cw_ ;
    mutable double64  co_ ;
    
    double64  awd_ ;
    double64  aod_ ;
    double64  cwd_ ;
    double64  cod_ ;
    
    double64  awi_ ;
    double64  aoi_ ;
    double64  cwi_ ;
    double64  coi_ ;
    
    std::array<std::array<double64, 2>,2> a_;
    std::array<std::array<double64, 2>,2> c_;
    
    ArrayVariable  ac_params_;
};

  
}  // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_WITH_HYSTERESIS_H */
