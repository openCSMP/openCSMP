#ifndef HETEROGENEITY_AND_RATE_AWARE_MODEL_H
#define HETEROGENEITY_AND_RATE_AWARE_MODEL_H

#include "TwoPhaseModel.h"

namespace csmp {

/**
    First attempt at implementing Maartje Boons & Sally Benson (Stanford visit May 16-23, 2019),
    see Maple worksheet 'BHP/Stanford-visit/IMPLEMENTATION/HeterogeneityAwareSaturationFunctions'
 
    @note thus far, this is a drainage only model
 
    @note use alternative capillary number: k ||grad p|| / sigma   to get around viscosity problem

    @note for lambda=0, this implementation of Brooks-Corey model switches to linear
    @note for linear case capillary pressure is a constant value equal to entry pressure
*/
template<size_t dim>
class HeterogeneityAndRateAwareModel : public TwoPhaseModel<dim> {
  public:

    HeterogeneityAndRateAwareModel( const PropertyDatabase<dim>& database,
                                    const char* rocktype, const char* total_velocity,
                                    const char* pc_entry,
                                    const bool sw_ro_mu_placement = true ); // NODE=true ELEMENT=false
  
    virtual ~HeterogeneityAndRateAwareModel();
  
    // TODO: implement
    virtual void Initialize( const Element<dim>& e );

    // relative permeabilities
    /// heterogeneity-aware, rate-dependent version, Nc is calculated in the background
    virtual double64 krw_Phase() const;
    virtual double64 krn_Phase() const;
  
    // TODO: capillary pressure
    virtual double64 pc_Phase( ) const;

    // TODO: capillary pressure derivatives (numeric)
    virtual double64 dpcds_Phase( ) const;

    // TODO: derivative of fractional flow (advection multipliers) - (numeric)
    virtual double64 dfds() const;
    
    // derivatives of gravitational flow term (advection multipliers)
    virtual double64 dGds() const;

    // maximum absolute value returned by dfdS
    virtual double64 MaxFractionalFlowDerivative() const;

    // TODO:
    virtual void Out( size_t phase ) const;
  
  private:
    // internal support functions
    /// model specific effective saturation
    double64  SwStar( double64 sw ) const { return (sw - swr_) / (1. - swr_); }
 
    /// returns CL water saturation value below which water will enter the fine-grained layers
    double64  SwStarKinkCO2( double64 sw, double64 Nc ) const;

    double64  SwStarKinkH2O( double64 sw, double64 Nc ) const;
  
    /// standard form: Nc = vt mu_CO2 / sigma
    double64  Nc_vtmuCO2_Version() const;
  
    /// pressure gradient form: Nc = k ||grad p|| / sigma
    double64  Nc_kgradP_Version() const;

  private:

    const csmp::Index  RRT_key_, bcp_key_, pd_key_, vt_key_;
    // Brooks-Corey model
    double64 k_, vt_magnitude_;
    double64 lambda_, entry_pressure_;  ///< Brooks-Corey parameter and drainage capillary threshold pressure (Brooks-Corey), CL curve

    const double64 IFT_ = 0.0035; ///< water - CO2 (N/m)
  
    // saturations
    double64  sw_, swr_, snr_;               ///< water saturations
    // specific points
    double64  Nc_, Nc_VL_, Nc_CL_;
    // model parameters (read by Initialise function)
  
};

} // end namespace csmp

#endif
