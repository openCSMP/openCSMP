#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H

#include "TwoPhaseModel.h"

namespace csmp {

/** @brief saturation function policy based on the Brooks-Corey (1964) model

    @attention relies on sw, swr, snr, variable values stored by the USER. 

    @note for lambda=0, this implementation of Brooks-Corey model switches to linear
    @note for linear case capillary pressure is a constant value equal to entry pressure
    @note base class pm1 and pm2 are used for pd and lambda, respectively
*/
template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctions {
  public:
    BrooksCorey( const Model<dim>&, const char* bc_param, const char* pd_param );
     ~BrooksCorey();
    
    /// relative permeabilities - parameters come from subclass FlowFunctions
    double64 krn_Phase() const;
    double64 krw_Phase() const;

    /// derivatives of relative permeabilities
    double64 dkrnds_Phase() const;
    double64 dkrwds_Phase() const;

    /// capillary pressure
    double64 pc_Phase( ) const;

    /// capillary pressure derivatives
    double64 dpcds_Phase( ) const;

    /// inverse capillary pressure function
    double64 Sw_Phase( double64 pc_Phase ) const;

    /// inverse capillary pressure derivative
    double64 dsdpc_Phase( double64 pc_Phase ) const;

    void Out( size_t phase ) const;

  private:
    /// reads Brooks-Corey  model parameters (always element properties)
    void UpdateModelParameters();

    /// shorthand for accessing the class that FacetFlux_TracerTransferExplicit is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    const csmp::Index  pd_key_, pc_max_key_, lamda_key_;
    double64           lambda_, entry_pressure_, pc_max_;
};

} // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H */
