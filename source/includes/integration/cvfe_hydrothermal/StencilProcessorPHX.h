#ifndef STENCIL_PROCESSOR_PHX_H
#define STENCIL_PROCESSOR_PHX_H

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

#include "FV_Parameter.h"
#include "DenseMatrix.h"

namespace csmp {

template<uint32_t> class Element;

  /**
     @class StencilProcessorPHX StencilProcessorPHX.h

     Stencil processor for finitive volume calculations within CVFEM scheme (Weis et al., Geofluids, 2014).

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes
  
     @section motivation Motivation
      Specialized stencil processor for PHX scheme.

     @section usage Usage
      To be used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Performs finite volume fluxes across facets.
          
     @endcode
     
     @section dependencies Dependencies
	 The functions are tailored for use in ExplicitFiniteVolumeTransportPHX
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
template<uint32_t dim>
class StencilProcessorPHX {
  public:
    StencilProcessorPHX();

    void DetermineFluxOutWithGravity( const FV_Parameter& param,
                                      const Element<dim>& e,
                                      DenseMatrix<DM_MIN>& upwind, //matrix indicating pre-defined upwind nodes
                                      std::vector<std::vector<double> >& facet_flux,  //vector storing flux across each facet
                                      std::vector<double>& flux_out, // vector storing the flux out of the control volumes
                                      csmp::Index rhs_key, // key for the right-hand side variable of fintite volume calculations
                                      csmp::Index rho_key, // density of the transported phase
                                      csmp::Index k_key);  // element permeability

	void DetermineFluxOutWithoutGravity( const FV_Parameter& param,
                                       const Element<dim>& e,
                                       DenseMatrix<DM_MIN>& upwind, //matrix indicating pre-defined upwind nodes
                                       std::vector<std::vector<double> >& facet_flux, //vector storing flux across each facet
                                       std::vector<double>& flux_out, // vector storing the flux out of the control volumes
                                       csmp::Index rhs_key);  // density of the transported phase

    void DetermineFluxIn( const Element<dim>& e,
                          std::vector<std::vector<double> >& facet_flux, //vector storing flux across each facet
                          std::vector<double>& flux_in, // vector storing the flux into of the control volumes
                          std::vector<double>& mass_balance ); // vector for mass balance corrections

	void DetermineFluxOut( const FV_Parameter& param,
                         const Element<dim>& e,
                         std::vector<std::vector<double> >& facet_flux, //vector storing flux across each facet
                         std::vector<double>& flux_out, // vector storing the flux out of the control volumes
                         csmp::Index rhs_key); // density of the transported phase

    ~StencilProcessorPHX();

  private:
    
    VectorVariable<dim> gravity_;

    size_t            eidx_;
    mutable uint32_t  inside_node_, outside_node_;

    ScalarVariable    rhs_property_, density_, perm_;
    double            grav_, g_component_, vel_;

};

} // namespace csmp

#endif




























