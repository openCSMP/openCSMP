#ifndef BROOKS_COREY_H
#define BROOKS_COREY_H

#include "TwoPhaseModel.h"
#include "CubicSpline.h"

namespace csmp {

/// @note for lambda=0, this implementation of Brooks-Corey model switches to linear
/// @note for linear case capillary pressure is a constant value equal to entry pressure
/// @note base class pm1 and pm2 are used for pd and lambda, respectively
template<uint32_t dim>
class BrooksCorey : public TwoPhaseModel<dim> {
  public:

    BrooksCorey( const PropertyDatabase<dim>& database,
                 const char* lamda, const char* pc_entry,
                 const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    BrooksCorey( const PropertyDatabase<dim>& database,
                 const char* permeability, 
                 double viscosity_nw, double viscosity_w,
                 double density_nw, double density_w,
                 const char* lamda, const char* pc_entry,
                 const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false


    /// constructor for constant capillary pressure equal to entry pressure  in case if lambda = 0
    BrooksCorey( const PropertyDatabase<dim>& database,
                 const char* permeability,
                 const char* viscosity_nw, const char* viscosity_w,
                 const char* density_nw, const char* density_w,
                 const char* lamda,
                 const char* pc_entry,
                 const char* sat_w,
                 const char* res_sat_nw, const char* res_sat_w,
                 const bool sw_ro_mu_placement = true ); // NODE=true ELEMENT=fals

    /// constructor for linear capillary pressure in case if lambda = 0
    BrooksCorey( const PropertyDatabase<dim>& database,
                 const char* permeability,
                 const char* viscosity_nw, const char* viscosity_w,
                 const char* density_nw, const char* density_w,
                 const char* lamda,
                 const char* pc_entry,
                 const char* maximum_pc,
                 const char* sat_w,
                 const char* res_sat_nw, const char* res_sat_w,
                 const bool sw_ro_mu_placement = true ); // NODE=true ELEMENT=fals

    // manual assignment of two phase properties
    BrooksCorey();
                           
    virtual ~BrooksCorey();
    
    virtual void Initialize( const Element<dim>& e );

    // relative permeabilities
    virtual double krn_Phase() const;
    virtual double krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double dkrnds_Phase() const;
    virtual double dkrwds_Phase() const;

    // capillary pressure
    virtual double pc_Phase( ) const;

    // capillary pressure derivatives
    virtual double dpcds_Phase( ) const;

    // inverse capillary pressure function
    virtual double Sw_Phase( double pc_Phase ) const;

    // inverse capillary pressure derivative
    virtual double dsdpc_Phase( double pc_Phase ) const;

    // derivative of fractional flow (advection multipliers)
    virtual double dfds() const;
    
    // derivatives of gravitational flow (advection multipliers)                                      
    virtual double dGds( ) const;

    // maximum absolute value returned by dfdS
    virtual double MaxFractionalFlowDerivative() const;

    // linearized fractional flow derivative
    virtual double ShockSpeed() const;
    virtual double ShockHeight() const;
    
    virtual void Out( size_t phase ) const;


    // brooks corey parameters
    void      Lambda( double lambda );
    double  Lambda() const;
    void      Pd( double pd );
    double  Pd() const;

  private:

    csmp::Index          pd_key_, pc_max_key_, lamda_key_;
    double             lambda_, entry_pressure_,pc_max_;
    bool                 default_capillary_pressure_max_;

};

} // end namespace csmp

#endif
