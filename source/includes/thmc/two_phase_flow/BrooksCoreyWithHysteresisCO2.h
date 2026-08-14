#ifndef BROOKS_COREY_WITH_HYSTERESIS_CO2_H
#define BROOKS_COREY_WITH_HYSTERESIS_CO2_H

#include "TwoPhaseModel.h"
#include "LinearTwoPhaseModel.h"

namespace csmp {

// for lambda=0, this implementation of Brooks-Corey model switches to linear
// base class pm1 and pm2 are used for pd and lambda, respectively
template<uint32_t dim>
class BrooksCoreyWithHysteresisCO2 : public TwoPhaseModel<dim> {
  public:
    BrooksCoreyWithHysteresisCO2( const PropertyDatabase<dim>& database,
                               const char* permeability, 
                               double viscosity_nw, double viscosity_w,
                               double density_nw, double density_w,
                               const char* lamda, const char* pc_entry );

    BrooksCoreyWithHysteresisCO2( const PropertyDatabase<dim>& database,
                 const char* lamda, const char* pc_entry );
  
    // relative permeabilities
    double krn_Phase() const override final;
    double krw_Phase() const override final;

    // derivative of fractional flow (advection multipliers)
    double dfds() const override final;
    
    // maximum absolute value returned by dfdS                              
    double MaxFractionalFlowDerivative() const override final;

    // derivatives of gravitational flow (advection multipliers)                                      
    double dGds() const override final;
    
    // capillary pressure
    double pc_Phase() const override final;
    
    // capillary pressure derivatives
    double dpcds_Phase() const override final;
      
    // linearized fractional flow derivative
    double ShockSpeed() const override final;
    double ShockHeight() const override final;
    
    void Out( uint32_t phase ) const override final;

  private:
      /// not constant as it sets (stores) the saturation inflection point on the element
    void Initialize( Element<dim>& );

  private:
    BrooksCoreyWithHysteresisCO2();
    csmp::Index   pd_key, lamda_key, sormax_key,
                  sat_previous_key, sat_inflection_key;
    const double  acc_gravity_;
    double        sat_previous,     ///< previous timestep water saturation
                  sat_inflection,   ///< brancing point to scaning curve
                  sormax;           ///< maximum residual oil saturation
    bool          imbibing,         ///< True when wetting phase increases
                  on_scaning_curve; ///< True when it is on scaning curve
    double CC, CK, Sot, Snorm, kroDswinflection, Sof_1_Snorm, kroI_1_Snorm, pm1, pm2; //temporary values for calculating krn()
};

} // end namespace csmp

#endif
