#ifndef CSMP_LINEAR_TWO_PHASE_MODEL_H
#define CSMP_LINEAR_TWO_PHASE_MODEL_H

#include "CSMP_definitions.h"
#include "TwoPhaseModel.h"

namespace csmp {

/// relative permeability and capillary pressure, assuming a linear model
template<size_t dim>
class LinearTwoPhaseModel : public TwoPhaseModel<dim>{
  public:
    /// constructor for constant capillary pressure equal to entry pressure
    LinearTwoPhaseModel(const PropertyDatabase<dim>& database,
                        const char* permeability, 
                        double64 viscosity_nw, double64 viscosity_w,
                        double64 density_nw, double64 density_w,
                        const char* pc_entry,
                        const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    /// constructor for constant capillary pressure equal to entry pressure
    LinearTwoPhaseModel(const PropertyDatabase<dim>& database,
                        const char* permeability, 
                        double64 viscosity_nw, double64 viscosity_w,
                        double64 density_nw, double64 density_w,
                        const char* pc_entry, const char* sw, const char* rsnw, const char* rsw,
                        const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false

    /// constructor for constant capillary pressure equal to entry pressure
    LinearTwoPhaseModel( const PropertyDatabase<dim>& database,
                         const char* permeability,
                         const char* viscosity_nw, const char* viscosity_w,
                         const char* density_nw, const char* density_w,
                         const char* pc_entry,
                         const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                         const bool sw_ro_mu_placement = true ); // NODE=true ELEMENT=fals

    /// constructor for linear capillary pressure
    LinearTwoPhaseModel(const PropertyDatabase<dim>& database,
                        const char* permeability,
                        double64 viscosity_nw, double64 viscosity_w,
                        double64 density_nw, double64 density_w,
                        const char* pc_entry,
                        const char* maximum_pc,
                        const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false


    /// constructor for linear capillary pressure
    LinearTwoPhaseModel( const PropertyDatabase<dim>& database,
                         const char* permeability,
                         const char* viscosity_nw, const char* viscosity_w,
                         const char* density_nw, const char* density_w,
                         const char* pc_entry,
                         const char* maximum_pc,
                         const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                         const bool sw_ro_mu_placement = true ); // NODE=true ELEMENT=fals

    virtual ~LinearTwoPhaseModel();
    
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

  private:
    LinearTwoPhaseModel();
    csmp::Index pd_key_, pc_max_key_;
    double64    entry_pressure_,pc_max_;
    bool        default_capillary_pressure_max_;
};
 
 } // end namespace csmp
 
#endif

