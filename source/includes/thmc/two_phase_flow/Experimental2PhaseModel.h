#ifndef EXPERIMENTAL_2PHASE_MODEL_H
#define EXPERIMENTAL_2PHASE_MODEL_H

#include "TwoPhaseModel.h"
#include "CubicSpline.h"

namespace csmp {

/// relperm- and pc-s relations from data files are interpolated with cubic splines
/// @author SKM @date 14/1/2007
template<size_t dim>
class Experimental2PhaseModel : public TwoPhaseModel<dim> {
  public:
    Experimental2PhaseModel( const PropertyDatabase<dim>& database,
                             const char* permeability, 
                             double64 viscosity_nw, double64 viscosity_w,
                             double64 density_nw, double64 density_w,
                             const char* kr1_data_file,
                             const char* kr2_data_file,
                             const char* pc_data_file,
                             const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    Experimental2PhaseModel( const PropertyDatabase<dim>& database,
                             const char* kr1_data_file,
                             const char* kr2_data_file,
                             const char* pc_data_filey,
                             const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false
                           
    virtual ~Experimental2PhaseModel();
    
    // reads model_key variable to switch between experimental and linear model
    virtual void Initialize( const Element<dim>& e );
    
    // relative permeabilities
    virtual double64 krn_Phase() const;
    virtual double64 krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double64 dkrnds_Phase() const;
    virtual double64 dkrwds_Phase() const;

    // derivative of fractional flow (advection multipliers)
    virtual double64 dfds() const;
    
    // maximum absolute value returned by dfdS                              
    virtual double64 MaxFractionalFlowDerivative() const;

    // derivatives of gravitational flow (advection multipliers)                                      
    virtual double64 dGds( ) const;
    
    // capillary pressure
    virtual double64 pc_Phase( ) const;
    
    // capillary pressure derivatives
    virtual double64 dpcds_Phase( ) const;
      
    // linearized fractional flow derivative
    virtual double64 ShockSpeed() const;
    virtual double64 ShockHeight() const;
    
    virtual void Out( size_t phase ) const;

  private:
    Experimental2PhaseModel();
    const double64       acc_gravity_; 
    CubicSpline          kr1_, kr2_, pc_;
    mutable double64     dfds_max_,   // max frac.flow derivative
                         dfds_shock_, // shock speed multiplier for vt
                         s1_shock_;   // shock height (saturation 1)
    double64             entry_pressure_,
                         lambda_;     // to distinguish linear from exp. relperm model
    csmp::Index          model_key_;  // to distinguish linear from exp. relperm model
    csmp::Index          pd_key_;     // reading in entry pressure for linear relperm model
    
    void ComputeShockSpeedAndHeight() const;
};

} // end namespace csmp

#endif
