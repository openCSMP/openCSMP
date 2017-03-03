#ifndef BROOKS_COREY_WITH_HYSTERESIS_H
#define BROOKS_COREY_WITH_HYSTERESIS_H

#include "TwoPhaseModel.h"
#include "LinearTwoPhaseModel.h"

namespace csmp {

// for lambda=0, this implementation of Brooks-Corey model switches to linear
// base class pm1 and pm2 are used for pd and lambda, respectively
template<size_t dim>
class BrooksCoreyWithHysteresis : public TwoPhaseModel<dim> {
  public:
    BrooksCoreyWithHysteresis( const PropertyDatabase<dim>& database,
                               const char* permeability, 
                               double64 viscosity_nw, double64 viscosity_w,
                               double64 density_nw, double64 density_w,
                               const char* lamda, const char* pc_entry );

    BrooksCoreyWithHysteresis( const PropertyDatabase<dim>& database,
                 const char* lamda, const char* pc_entry );
                           
    virtual ~BrooksCoreyWithHysteresis();
    
    virtual void Initialize( Element<dim>& e );
    
    // relative permeabilities
    virtual double64 krn_Phase() const;
    virtual double64 krw_Phase() const;

    // derivative of fractional flow (advection multipliers)
    virtual double64 dfds() const;
    
    // maximum absolute value returned by dfdS                              
    virtual double64 MaxFractionalFlowDerivative() const;

    // derivatives of gravitational flow (advection multipliers)                                      
    virtual double64 dGds() const;
    
    // capillary pressure
    virtual double64 pc_Phase( size_t phase ) const;
    
    // capillary pressure derivatives
    virtual double64 dpcds_Phase( size_t phase ) const;
      
    // linearized fractional flow derivative
    virtual double64 ShockSpeed() const;
    virtual double64 ShockHeight() const;
    
    void Out(size_t phase) const { Out(std::cout, phase); }
    virtual void Out( std::ostream& os, size_t phase ) const;

  private:
    BrooksCoreyWithHysteresis();
    csmp::Index          pd_key, lamda_key, sormax_key, 
                         sat_previous_key, sat_inflection_key;
    const double64       acc_gravity_; 
    double64             sat_previous, // previous timestep water saturation
                         sat_inflection, // brancing point to scaning curve
                         sormax; //maximum residual oil saturation
    bool                 imbibing, //True when wetting phase increases
                         on_scaning_curve; //True when it is on scaning curve                     
    double64 CC, CK, Sot, Snorm, kroDswinflection, Sof_1_Snorm, kroI_1_Snorm, pm1, pm2; //temporary values for calculating krn()
};

} // end namespace csmp

#endif
