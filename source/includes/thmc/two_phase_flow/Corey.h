#ifndef COREY_H
#define COREY_H

#include "TwoPhaseModel.h"


namespace csmp {

/** 
   General Corey relationships:
   - for lambda=0, this implementation of Brooks-Corey model switches to linear
     base class pm1 and pm2 are used for pd and lambda, respectively
 */
template<uint32_t dim>
class Corey : public TwoPhaseModel<dim> {
  public:

    Corey( const PropertyDatabase<dim>& database );
    
    /// constructor for constant capillary pressure equal to entry pressure  in case if lambda = 0
    Corey( const PropertyDatabase<dim>& database,
           const char* kkk,   // permeability
           const char* viscosity_nw, const char* viscosity_w,
           const char* density_nw, const char* density_w,
           const char* lambda,// corey capillary pressure exponent
           const char* pc_entry, // capillary entry pressure
           const char* exp_nw,// corey non-wetting phase exponent
           const char* exp_w, // corey wetting phase exponent
           const char* krnw,  // residual non-wetting phase end point
           const char* krw,   // residual wetting phase end point
           const char* sat,   // saturation wetting phase
           const char* snr,   // residual saturation non-wetting phase
           const char* swr,   // residual saturation wetting phase
           const bool  sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    /// constructor for linear capillary pressure in case if lambda = 0
    Corey( const PropertyDatabase<dim>& database,
           const char* kkk,   // permeability
           const char* viscosity_nw, const char* viscosity_w,
           const char* density_nw, const char* density_w,
           const char* lambda,// corey capillary pressure exponent
           const char* pc_entry, // capillary entry pressure
           const char* maximum_pc, // maximum permitted capillary pressure
           const char* exp_nw,// corey non-wetting phase exponent
           const char* exp_w, // corey wetting phase exponent
           const char* krnw,  // residual non-wetting phase end point
           const char* krw,   // residual wetting phase end point
           const char* sat,   // saturation wetting phase
           const char* snr,   // residual saturation non-wetting phase
           const char* swr,   // residual saturation wetting phase
           const bool  sw_ro_mu_placement = true); // NODE=true ELEMENT=false


    virtual ~Corey();
    
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

    // maximum absolute value returned by dfdS
    virtual double MaxFractionalFlowDerivative() const;

    virtual void Out( size_t phase ) const;

  private:

    Corey();
    Index           pd_key_, pc_max_key_,expw_key_, expn_key_, lambda_key_, krw_key_, krn_key_;
    double        krw_, krn_, expw_, expn_, lambda_, entry_pressure_,pc_max_;
    bool            default_capillary_pressure_max_;

};

} // end namespace csp

#endif
