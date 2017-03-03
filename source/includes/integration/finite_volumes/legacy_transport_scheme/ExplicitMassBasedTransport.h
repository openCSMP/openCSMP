#ifndef EXPLICIT_MASS_BASED_TRANSPORT_H
#define EXPLICIT_MASS_BASED_TRANSPORT_H

#include "ExplicitNodeCenteredFiniteVolumeTransport.h"

namespace csmp {

template<size_t dim,template<size_t> class STP>
class ExplicitMassBasedTransport : public ExplicitNodeCenteredFiniteVolumeTransport<dim, STP> {
  public:
    
    // single phase solute advection-only constructor for a subregion identified as a group, 
    // with different lhs and rhs
    ExplicitMassBasedTransport( const char* group,
                                Model<dim>& sg, // not constant since FV's are created
                                const char* porosity,
                                const char* advected_prop_lhs,
                                const char* advected_prop_rhs,
                                const char* transp_velocity,
                                const char* nodal_source,
                                bool second_order_accuracy,
                                bool second_order_in_time=false,
                                const char* thickness=NULL);
    virtual double64  AnisotropicCourantIncrement();
    virtual void AdvectVariable( double64 time_interval);

    /// single-phase passive advection, does NOT return courant increment, single timestep calculation
    /// no checks are made for courant condition.  Assumes external checks.
    virtual void AdvectVariableSingleStep( double64 time_increment,
                                      bool apply_flux_balance_correction,
                                      bool update_pore_volumes);

    void AdvectVariable1stOrder( double64 time_interval, bool output_result_range );

  protected:
    virtual void AccumulateFluxUpwindProducts();
    void AccumulateFluxUpwindProductsOMP(std::vector<double64>& RESULT);
    virtual void AssignFluxBoundaryConditions(const size_t var_comp_nr=0);
    virtual void ComposeSolution( double64 time_interval,const size_t var_comp_nr=0);

  protected:
    const csmp::Index                   ad_rhs_key_;       ///< transported variable

};
}
#endif //EXPLICIT_MASS_BASED_TRANSPORT_H
