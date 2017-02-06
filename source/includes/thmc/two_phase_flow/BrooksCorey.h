#ifndef BROOKS_COREY_H
#define BROOKS_COREY_H

#include "TwoPhaseModel.h"
#include "CubicSpline.h"

namespace csmp {

/// @note for lambda=0, this implementation of Brooks-Corey model switches to linear
/// @note for linear case capillary pressure is a constant value equal to entry pressure
/// @note base class pm1 and pm2 are used for pd and lambda, respectively
template<size_t dim>
class BrooksCorey : public TwoPhaseModel<dim> {
  public:

    BrooksCorey( const PropertyDatabase<dim>& database,
                 const char* lamda, const char* pc_entry,
                 const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    BrooksCorey( const PropertyDatabase<dim>& database,
                 const char* permeability, 
                 double64 viscosity_nw, double64 viscosity_w,
                 double64 density_nw, double64 density_w,
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
    virtual double64 krn_Phase() const;
    virtual double64 krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double64 dkrnds_Phase() const;
    virtual double64 dkrwds_Phase() const;

    // capillary pressure
    virtual double64 pc_Phase( ) const;

    // capillary pressure derivatives
    virtual double64 dpcds_Phase( ) const;

    // inverse capillary pressure function
    virtual double64 Sw_Phase( double64 pc_Phase ) const;

    // inverse capillary pressure derivative
    virtual double64 dsdpc_Phase( double64 pc_Phase ) const;

    // derivative of fractional flow (advection multipliers)
    virtual double64 dfds() const;
    
    // derivatives of gravitational flow (advection multipliers)                                      
    virtual double64 dGds( ) const;

    // maximum absolute value returned by dfdS
    virtual double64 MaxFractionalFlowDerivative() const;

    // linearized fractional flow derivative
    virtual double64 ShockSpeed() const;
    virtual double64 ShockHeight() const;
    
    virtual void Out( std::ostream& os, size_t phase ) const;


    // brooks corey parameters
    void      Lambda( double64 lambda );
    double64  Lambda() const;
    void      Pd( double64 pd );
    double64  Pd() const;

  private:

    csmp::Index          pd_key_, pc_max_key_, lamda_key_;
    double64             lambda_, entry_pressure_,pc_max_;
    bool                 default_capillary_pressure_max_;

};

template<size_t dim>
inline double64 csmp::BrooksCorey<dim>::Pd() const
  {
    return entry_pressure_;
  }


template<size_t dim>
inline void csmp::BrooksCorey<dim>::Pd( double64 pd )
  {
    entry_pressure_ = pd;
  }


template<size_t dim>
inline void csmp::BrooksCorey<dim>::Lambda( double64 lambda )
  {
    lambda_ = lambda;
  }

template<size_t dim>
inline double64 csmp::BrooksCorey<dim>::Lambda() const
  {
    return lambda_;
  }


} // end namespace csmp

#endif
