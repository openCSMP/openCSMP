#ifndef BROOKS_COREY_WITH_HYSTERESIS_H
#define BROOKS_COREY_WITH_HYSTERESIS_H

#include "TwoPhaseModel.h"

namespace csmp {

// for lambda=0, this implementation of Brooks-Corey model switches to linear
// base class pm1 and pm2 are used for pd and lambda, respectively
template<uint32_t dim>
class BrooksCoreyWithHysteresis : public TwoPhaseModel<dim> {
  public:
    BrooksCoreyWithHysteresis( const PropertyDatabase<dim>& database,
                               const char* permeability, 
                               double viscosity_nw, double viscosity_w,
                               double density_nw, double density_w,
                               const char* lamda, const char* pc_entry );

    BrooksCoreyWithHysteresis( const PropertyDatabase<dim>& database,
                               const char* lamda, const char* pc_entry );
                           
    virtual ~BrooksCoreyWithHysteresis() {}
  
    /// not constant as it sets the saturation inflection point
    virtual void InitializeAndStore( Element<dim>& );
    
    // relative permeabilities
    virtual double krn_Phase() const;
    virtual double krw_Phase() const;

    // derivative of fractional flow (advection multipliers)
    virtual double dfds() const;
    
    // maximum absolute value returned by dfdS                              
    virtual double MaxFractionalFlowDerivative() const;

    // derivatives of gravitational flow (advection multipliers)                                      
    virtual double dGds() const;
    
    // capillary pressure
    virtual double pc_Phase( size_t phase ) const;
    
    // capillary pressure derivatives
    virtual double dpcds_Phase( size_t phase ) const;
      
    // linearized fractional flow derivative
    virtual double ShockSpeed() const;
    virtual double ShockHeight() const;
    
    virtual void Out( uint32_t phase ) const;

  private:
    BrooksCoreyWithHysteresis();
    csmp::Index          pd_key, lamda_key, sormax_key, 
                         sat_previous_key, sat_inflection_key;
    const double       acc_gravity_; 
    double             sat_previous, // previous timestep water saturation
                         sat_inflection, // brancing point to scaning curve
                         sormax; //maximum residual oil saturation
    bool                 imbibing, //True when wetting phase increases
                         on_scaning_curve; //True when it is on scaning curve                     
    double CC, CK, Sot, Snorm, kroDswinflection, Sof_1_Snorm, kroI_1_Snorm, pm1, pm2; //temporary values for calculating krn()
};

} // end namespace csmp

#endif
