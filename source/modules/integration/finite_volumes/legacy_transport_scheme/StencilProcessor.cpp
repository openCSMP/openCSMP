#include "StencilProcessor.h"
#include "FV_Parameter.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "Element.h"
#include "CSMP_mathUtilities.h"
#include "TwoPhaseModel.h"

#include <limits>

using namespace std;

#define DEBUG_ImplicitStencilProcessor

namespace csmp {

template<uint32_t dim>
StencilProcessor<dim>::StencilProcessor( const csmp::Index& adv_key,   
                                         const csmp::Index& velo_key,
                                         const csmp::Index& src_key,
                                         const csmp::Index& velo_mult_key)
  : adv1_key_(adv_key),
    vel_key_(velo_key),
    src_key_(src_key),
    velo_mult_key_(velo_mult_key),
    dsdn_(dim),
    dpcdsn_(dim)
 {
 }


template<uint32_t dim>
StencilProcessor<dim>::StencilProcessor( const csmp::Index& adv_key,   
                                         const csmp::Index& velo_key,
                                         const csmp::Index& diffusivity_key,
                                         const csmp::Index& src_key,
                                         const csmp::Index& velo_mult_key)
  : adv1_key_(adv_key),
    vel_key_(velo_key),
    diff_key_(diffusivity_key),
    src_key_(src_key),
    velo_mult_key_(velo_mult_key),
    dsdn_(dim),
    dpcdsn_(dim)
 {
 }



/** Standard assignment operator.
*/
/*
template<uint32_t dim>
StencilProcessor<dim>& StencilProcessor<dim>::operator=( const StencilProcessor& tfs )
 {
    if ( &tfs != this ) {
         eidx_               = tfs.eidx_;
         diff_coeff_         = tfs.diff_coeff_;
  	     sector_pore_volume_ = tfs.sector_pore_volume_;
  	     facet_flux_         = tfs.facet_flux_;
  	     psi1_               = tfs.psi1_;
  	     ipsi1_              = tfs.ipsi1_;
  	     src_                = tfs.src_;
  	     adv1key_            = tfs.adv1key_;
  	     vel0key_            = tfs.vel0key_;
  	     diffkey_            = tfs.diffkey_;  	         
      }
    return *this;
 }
*/

 
 
template<uint32_t dim>
StencilProcessor<dim>::StencilProcessor( const StencilProcessor<dim>& tfs )
  :  eidx_(tfs.eidx_),
     diff_coeff_(tfs.diff_coeff_),
     sector_pore_volume_ (tfs.sector_pore_volume_),
     facet_flux_(tfs.facet_flux_),
     psi1_(tfs.psi1_),
     ipsi1_(tfs.ipsi1_),
     src_(tfs.src_),
     adv1_key_(tfs.adv1_key_),
     vel_key_(tfs.vel_key_),
     diff_key_(tfs.diff_key_),
     dsdn_(tfs.dsdn_),
     dpcdsn_(tfs.dpcdsn_)
 {
 }
 

 
/**
Destructor is used to write theta limiter values (if accumulated) to
a text file named "spatial-limiter-values" if the flag write_out_limiter
has been set.  
*/
template<uint32_t dim>
StencilProcessor<dim>::~StencilProcessor()
 {
 } // end destructor



/**
Initializes sector pore volume and facet flux vectors for
current finite element = finite volume stencil.  
It is implied that advected variable and nodal source are always of the same type 
(both SCALAR, both ARRAY or both FLAGGEDARRAY )
*/
template<uint32_t dim>
void StencilProcessor<dim>::InitializeFirstOrder( const FV_Parameter& param,
                                                  const Element<dim>& e,
                                                  const VARIABLE_TYPE& vt,
                                                  uint32_t var_comp_nr )
 {
     eidx_ = e.Idx();
     sector_pore_volume_.resize(e.Nodes());
     psi1_.resize(e.Nodes());
     src_.resize(e.Nodes());
     nodal_src_.resize(e.Nodes());
     if ( !(diff_key_.index == ULONG_MAX) )
       diff_coeff_ = e.Read( diff_key_ );
     bool velo_mult_flag;
     if ( velo_mult_key_.index == ULONG_MAX )
        velo_mult_flag = false;
     else
        velo_mult_flag = true;

     // getting all the node-related information
     // ----------------------------------------
     if (vt==SCALAR){
         for ( auto i{0U}; i<e.Nodes(); i++ ) {
             // sector pore volumes
             sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable
             psi1_[i]       = e.N(i)->Read( adv1_key_ );
             nodal_src_[i]  = e.N(i)->Read( src_key_ );
         }
     }
     else if (vt==ARRAY){
         for ( auto i{0U}; i<e.Nodes(); i++ ) {
             // sector pore volumes
             sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable and source
             ArrayVariable var, src;
             
             e.N(i)->Read( adv1_key_ ,var);
             psi1_[i] = var(var_comp_nr);
             
             e.N(i)->Read( src_key_ ,src);
             nodal_src_[i] = src(var_comp_nr);
         }
     }
     else if (vt==FLAGGEDARRAY){
         for ( auto i{0U}; i<e.Nodes(); i++ ) {
             // sector pore volumes
             sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable and source
             FlaggedArrayVariable var, src;
             
             e.N(i)->Read( adv1_key_ ,var);
             psi1_[i] = var(var_comp_nr);
             
             e.N(i)->Read( src_key_ ,src);
             nodal_src_[i] = src(var_comp_nr);
         }
     }

     // getting all the finite-volume facet related information
     // -------------------------------------------------------
     facet_flux_.resize(e.FV()->Facets());
     
     for ( auto i{0U}; i<e.FV()->Facets(); i++ ) {
         // project the velocities onto the facet normals to get
         // fluxes once the projections have been multiplied with
         // the surface areas
         facet_flux_[i]  = param.FacetNormalVelocity( i );
         facet_flux_[i] *= param.FacetArea( i );
         
         e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
         if (velo_mult_flag)
         {
            if (param.FacetNormalVelocity(i) > 0)
              velo_mult_ = e.N(inside_node_)->Read(velo_mult_key_);
            else
              velo_mult_ = e.N(outside_node_)->Read(velo_mult_key_);
            facet_flux_[i] *= velo_mult_;
         }
     }
     
 } // end InitializeFirstOrder

/** Computes Flux Mismatch for FVs at the Boundaries
*/
template<uint32_t dim>
void StencilProcessor<dim>::ComputeBoundaryFluxMismatch(const FV_Parameter& param, const Element<dim>& e )
  {
     eidx_ = e.Idx();
     src_.resize(e.Nodes());

     bool velo_mult_flag;
     if ( velo_mult_key_.index == ULONG_MAX )
        velo_mult_flag = false;
     else
        velo_mult_flag = true;
     
     for ( auto i{0U}; i<e.Nodes(); i++ ) {

         src_[i] = 0.;
         
         // compute flux mismatch
         if (e.N(i)->AtBoundary() != NOT)
           {
             for ( auto k{0U}; k<e.FV()->FacetsPerSector(i); k++ ) {
                 const auto n(e.FV()->FacetSurroundingSector(i,k));
                 e.FV()->FacetEdgeNodes( n, inside_node_, outside_node_ );
                 // if the sector node is the inside node then an incoming flux will create a positive source term
                 double flux( (i == inside_node_) ? -param.FacetNormalVelocity(n) * param.FacetArea(n)
                                                          :  param.FacetNormalVelocity(n) * param.FacetArea(n) );

                 if (velo_mult_flag)
                 {
                    if (param.FacetNormalVelocity(n) > 0)
                      velo_mult_ = e.N(inside_node_)->Read(velo_mult_key_);
                    else
                      velo_mult_ = e.N(outside_node_)->Read(velo_mult_key_);
                    flux *= velo_mult_;
                 }

                 src_[i] += flux;
             }
           }
     }
  
  } // end ComputeBoundaryFluxMismatch


/**
Initializes sector pore volume, facet flux, and facet advected variable
vectors for current finite element = finite volume stencil.  
*/
template<uint32_t dim>
void StencilProcessor<dim>::InitializeSecondOrder( const FV_Parameter& param,
                                                   const Element<dim>& e )
 {
     eidx_ = e.Idx();

     const auto facets(e.FV()->Facets());
     sector_pore_volume_.resize(e.Nodes());
     psi1_.resize(e.Nodes());
     src_.resize(e.Nodes());
     if ( !(diff_key_.index == ULONG_MAX) )
       diff_coeff_ = e.Read( diff_key_ );
     bool velo_mult_flag;
     if ( velo_mult_key_.index == ULONG_MAX )
        velo_mult_flag = false;
     else
        velo_mult_flag = true;

     // getting all the node-property related information
     // -------------------------------------------------
     for ( auto i{0U}; i<e.FV()->Sectors(); i++ ) {
          // sector pore volumes
          sector_pore_volume_[i] = param.SectorVolume( i );
          // advected variable (transformed by bijective mapping)
          psi1_[i] = e.N(i)->Read( adv1_key_ );
       }

     ipsi1_.resize(facets);
     facet_flux_.resize(facets);

     // interpolating values of advected quantity to facet integration points
     // --------------------------------------------------------------------
     for ( auto i{0U}; i<facets; i++ ) {
          // interpolating the property value to the (first) facet integration point
          ipsi1_[i] = ( e ).PropertyValueAtFacetIntegrationPoint( i, 0U, adv1_key_ );

          // integration of facet fluxes
          facet_flux_[i]  = param.FacetNormalVelocity( i );
          facet_flux_[i] *= param.FacetArea( i );
         
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          if (velo_mult_flag)
          {
            if (param.FacetNormalVelocity(i) > 0)
              velo_mult_ = e.N(inside_node_)->Read(velo_mult_key_);
            else
              velo_mult_ = e.N(outside_node_)->Read(velo_mult_key_);
            facet_flux_[i] *= velo_mult_;
          }
       }
     
 } // end InitializeSecondOrder



/**
Like 'Initialize()' but without the calculation of the facet
fluxes. Depending on whether bijective mapping is used or not either 
variable is interpolated to the facet integration points.  

@section arguments Input Arguments 

The FACETFLUXES at time-level t+dt are used to find the upstream 
saturation values.
@section application Application

This method is used in the non-linear iteration loops to avoid 
unnecessary floating point operations.  
 */
template<uint32_t dim>
void StencilProcessor<dim>::InitializeAdvectedVariableValues( const Element<dim>& e )
 {
     eidx_ = e.Idx();

     const auto facets(e.FV()->Facets());
     psi1_.resize(e.Nodes());

     if ( !(diff_key_.index == ULONG_MAX) ) diff_coeff_ = e.Read( diff_key_ );

     // getting all the node-property related information
     // -------------------------------------------------
     for ( auto i{0U}; i<e.Nodes(); i++ )
       // advected variable & advected variable transformed by bijective mapping
       psi1_[i] = e.N(i)->Read( adv1_key_ ); // for upstream saturations

     ipsi1_.resize(facets);

     // interpolating advected variable to finite-volume facet integration points
     // -------------------------------------------------------------------------
     for ( auto i{0U}; i<facets; i++ )
       ipsi1_[i] = (e).PropertyValueAtFacetIntegrationPoint( i, 0U, adv1_key_ );
     
 } // end InitializeAdvectedPropertyValues







/** Limiting the value of the interpolated variable using the NVD approach
(see Matthai et al. 2009, TIPM).

@section arguments Input Arguments

SMINMAX stores the min/max values of the transport variable in the v
vicinity of and including the upstream node.

@section implementation Implementation

If the facet is located at the model boundary the limiter automatically
sets the flux to first order.

@section application Application

To avoid numerical oscillations in higher-order transport schemes.
*/
template<uint32_t dim>
void  StencilProcessor<dim>::IsotropicallyLimitTransportProperties( const Element<dim>& e,
                                                                    const vector< pair<double,double> >& SMINMAX )
 {
    eidx_ = e.Idx();
    uint32_t     n_upstr, n_dnstr;
    const double zero(0.), xi(2.);

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
      {
         e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

         // ----------   finding out which connected element is upstream and which is downstream
         if ( facet_flux_[i] > zero ) { n_upstr=inside_node_;  n_dnstr=outside_node_; }
         else                         { n_upstr=outside_node_; n_dnstr=inside_node_; }  // n_upstr = n_current

         // if either node on the facet is at the boundary of the Model
         if ( e.N(inside_node_)->AtBoundary()  != NOT ||
              e.N(outside_node_)->AtBoundary() != NOT ) {
              ipsi1_[i] = psi1_[n_upstr];
           }
         else
             ipsi1_[i] = limitProperty( psi1_[n_upstr], psi1_[n_dnstr], ipsi1_[i], SMINMAX[ e.N(n_upstr)->Idx() ], xi);
     }

 } // end IsotropicallyLimitTransportProperties



template<uint32_t dim>
void  StencilProcessor<dim>::ApplyLeastSquareMethodToLimitTransportProperties( const Element<dim>& e,
                                                                              const std::vector<std::pair<double,double> >& SMINMAX,
                                                                              const csmp::Index& mass_center_key,
                                                                              const csmp::Index& grad_psi_key,
                                                                              const csmp::Index& grad_psi_limiter_key)
 {

    eidx_ = e.Idx();
    uint32_t n_upstr, n_dnstr;
    const double   zero(0.);

    for ( auto i{0U}; i<e.FV()->Facets(); i++ ){

          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          // ----------   finding out which connected element is upstream and which is downstream
          if ( facet_flux_[i] > zero ) { n_upstr=inside_node_;  n_dnstr=outside_node_; }
          else                         { n_upstr=outside_node_; n_dnstr=inside_node_; }  // n_upstr = n_current

          // if either node on the facet is at the boundary of the Model
          if ( e.N(inside_node_)->AtBoundary()  != NOT ||
               e.N(outside_node_)->AtBoundary() != NOT ) {
               ipsi1_[i] = psi1_[n_upstr];
            }
          else
              ipsi1_[i]  = limitProperty_LSMGRAD<dim>( e, mass_center_key, grad_psi_key, grad_psi_limiter_key,
                                                       n_upstr, i, psi1_[n_upstr]);

    }

 } // end ApplyLeastSquareMethodToLimitTransportProperties

/**

Limits advected variable values averaged onto the current facet so that
it can be used to compute a non-oscillatory source term.

@section arguments Input Arguments

Methods expects the StencilProcessor to be initialized.
It needs the private variables:

inside and outside node and a flux that is positive if aligned with
the facet normal and negative if flow is into the facet.

@return The method returns the limited value of the variable at the facet
integration point.

@section application Application

Limiting of source terms for explicit (for instance as in IMPES)
calculations.
 */
template<uint32_t dim>
double  StencilProcessor<dim>::LimitExplicitSourceTerm(
                                       const Element<dim>& e,
                                       const vector<pair<double,double> >& SMINMAX,
                                       double inside_var_value,
                                       double outside_var_value,
                                       double psi_dash_f,        // average value at facet
                                       double flux )
 {
     if ( flux > 0. )
         return limitProperty( inside_var_value, outside_var_value, psi_dash_f, SMINMAX[ e.N(inside_node_)->Idx() ]);

     return limitProperty( outside_var_value, inside_var_value, psi_dash_f, SMINMAX[ e.N(outside_node_)->Idx() ]);

 } // end LimitExplicitSourceTerm



/**

See Matthai et al. 2009, TIPM for the explanation of the theta-
time-level balancing method.
 */
template<uint32_t dim>
void StencilProcessor<dim>::EvaluateThetaValues( const Element<dim>& e,
                                                    const vector<double>&  FVPOREVOL, // vols of FV's
                                                    const vector<double>&  SAT0, // saturation at initial time-level
                                                    const vector<vector<double> >& FACETFLUXES0, // old timestep
                                                    const vector<vector<double> >& LTDSATS0, // old timestep
                                                    double time_increment,
                                                    bool  use_max_theta )
 {
    uint32_t  n_upstr, n_dnstr; // upstream & downstream elements
    const double zero(0.);
    theta_.resize(e.FV()->Facets());

    // for each sector divider
    for ( auto i{0U}; i<theta_.size(); i++ )
      {
         e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

         // ----------   finding out which connected element is upstream and which is downstream
         if ( facet_flux_[i] > zero ) { n_upstr=inside_node_;  n_dnstr=outside_node_; }
         else                         { n_upstr=outside_node_; n_dnstr=inside_node_; }  // n_upstr = n_current

         // computing solute fluxes (eqn. 44, Pain et al.) across face 1 at time levels t and t + dt
         double hf_t0 = FACETFLUXES0[eidx_][i] * LTDSATS0[eidx_][i]; // from last timestep
         double hf_t1 = facet_flux_[i]         * ipsi1_[i];          // for current timestep

         if ( !use_max_theta )
           theta_[i] = thetaLimiter( hf_t0, hf_t1,
                                    FVPOREVOL[ e.N(n_upstr)->Idx() ],
                                    FVPOREVOL[ e.N(n_dnstr)->Idx() ],
                                    psi1_[n_upstr] - SAT0[ e.N(n_upstr)->Idx() ],
                                    psi1_[n_dnstr] - SAT0[ e.N(n_dnstr)->Idx() ],
                                    time_increment );
         else
           // the theta value is allowed to increase but not decrease
           theta_[i] = std::max( thetaLimiter( hf_t0, hf_t1,
                                             FVPOREVOL[ e.N(n_upstr)->Idx() ],
                                             FVPOREVOL[ e.N(n_dnstr)->Idx() ],
                                             psi1_[n_upstr] - SAT0[ e.N(n_upstr)->Idx() ],
                                             psi1_[n_dnstr] - SAT0[ e.N(n_dnstr)->Idx() ],
                                             time_increment ), theta_[i] );
//         cout << theta[i] <<" ";
      }

} // end EvaluateThetaValues






/** Computes all the terms necessary to solve in the saturation equation.

The gravitational segregation of the phases is taken care via  
the second term in:

dfds vt - d lambda_overbar / ds k (rho1 - rho2) g

This requires the analytically found derivative of the lambda overbar 
function.  

@section implementation Implementation

All flow properties are computed at the element barycenter.  

Capillary spreading is modelled with the FEM, so it is not included
here.   

MAYBE FACET INTEGRATION POINT INTERPOLATION WORKS BETTER in multidimensional case:

@code
   relperm.InitializeForFacetIntegrationPoint( i, 0U, pmem_, e );  
   relperm.EffectiveSaturation(); 
   facet_flux_[i] = param.FacetNormalVelocity(i) * param.FacetArea(i) * relperm.dfds();
@endcode
*/
template<uint32_t dim>
void  StencilProcessor<dim>::AccumulateImplicitTwoPhaseSolution1( 
                                          const FV_Parameter& param,
                                          const Element<dim>& e,
                                          TwoPhaseModel<dim>& relperm )
 {
     eidx_ = e.Idx();
     // --------------------------------------------------------------------------  
     // initialization of relative permeability model for saturation interpolated
     // to element barycenter (this is the only internally consistent approach 
     // when total mobility is computed at the barycenter by the FEM computations)
     // --------------------------------------------------------------------------  
     relperm.Initialize( e );
     relperm.InitializeForBaryCenter( e );  
     relperm.EffectiveSaturation(); 
 
     assert( !isnan(relperm.AdvectionMultiplier()) );
     assert( !isnan(relperm.CapillaryDiffusionMultiplier( )) );
     assert( !isnan(relperm.GravityMultiplier_dGds()) );

     const double dfds = relperm.AdvectionMultiplier();
     const double dGds = relperm.GravityMultiplier_dGds();
     if ( !(diff_key_.index == ULONG_MAX) )
       diff_coeff_ = -relperm.CapillaryDiffusionMultiplier( );
       
     // some provisions for taking account of gravitational forces
     const size_t v( (dim==1u) ? 0u : 1u );

     // ---------------------------------------------------------------------
     // 1. keeping track of the divergence of the non-linear facet flux
     // ---------------------------------------------------------------------
     sector_pore_volume_.resize(e.Nodes());
     psi1_.resize(e.Nodes());
     src_.resize(e.Nodes());

     // for all FV sectors = nodes = matrix columns
     for ( auto j{0U}; j<e.Nodes(); j++ )
       {
          sector_pore_volume_[j] = param.SectorVolume(j);
          psi1_[j]               = e.N(j)->Read( adv1_key_ );
          // for all sectors contained in this element
          // accumulate fluid sources (+) or sinks (-) due to deviations from 
          // potentially non-conservative fluxes
          src_[j] = 0.;  
          for ( auto k{0U}; k<e.FV()->FacetsPerSector(j); k++ ) {
               const auto n(e.FV()->FacetSurroundingSector(j,k));
               // if the sector node is the inside node then an incoming flux will create a positive source term
               const double fsign( (j == e.FV()->InsideNode(n)) ? -1. : 1. ); 
               src_[j] += fsign * param.FacetArea(n) * 
                         (param.FacetNormalVelocity(n) * dfds + param.FacetNormalComponent(n,v) * dGds);
            }
       }

      // -----------------------------------------------------------------------------
      // 2. accumulation of first-order facet fluxes (pressure- & gravity driven flow)
      //    (compare with Helmig, 97, p. 108, eqn. 3.74)
      // -----------------------------------------------------------------------------
      facet_flux_.resize( e.FV()->Facets() );
      for ( auto i{0U}; i<e.FV()->Facets(); i++ )
      {
        facet_flux_[i] = param.FacetNormalVelocity(i) * dfds;

        // ---------------------------------------------------------------------------
        // 3. gravitational segregation of wetting and non-wetting phases in zones of mixed
        //    saturation
        // ---------------------------------------------------------------------------
        // only if there is a density difference
        if ( fabs( relperm.DensityWettingPhase() - relperm.DensityNonWettingPhase() ) > numeric_limits<double>::epsilon() and
             fabs( param.FacetNormalComponent(i, v) ) > numeric_limits<double>::epsilon() )
        {
          // del_lambda_overbar / del_sn * k * -g(rho_w - rho_n)
          facet_flux_[i] += param.FacetNormalComponent(i, v) * dGds;
        }
            
        facet_flux_[i] *= param.FacetArea(i);
      }

  } // end AccumulateImplicitTwoPhaseSolution1 (generic version)


template<uint32_t dim>
void  StencilProcessor<dim>::AccumulateImplicitTwoPhaseSolution2(
                                          const FV_Parameter& param,
                                          const Element<dim>& e,
                                          TwoPhaseModel<dim>& relperm )
 {
     eidx_ = e.Idx();

     // --------------------------------------------------------------------------
     // initialization of relative permeability model for saturation interpolated
     // to element barycenter (this is the only internally consistent approach
     // when total mobility is computed at the barycenter by the FEM computations)
     // --------------------------------------------------------------------------
     relperm.Initialize( e );
     relperm.InitializeForBaryCenter( e );
     relperm.EffectiveSaturation();

     assert( !isnan(relperm.AdvectionMultiplier()) );
     assert( !isnan(relperm.CapillaryDiffusionMultiplier( )) );
     assert( !isnan(relperm.GravityMultiplier_dGds()) );

     const double dfds = relperm.AdvectionMultiplier();
     const double dGds = relperm.GravityMultiplier_dGds();
     if ( !(diff_key_.index == ULONG_MAX) )
       diff_coeff_ = -relperm.CapillaryDiffusionMultiplier();

     // some provisions for taking account of gravitational forces
     const auto v( (dim==1u) ? 0u : 1u );

     // ---------------------------------------------------------------------
     // 1. keeping track of the divergence of the non-linear facet flux
     // ---------------------------------------------------------------------
     sector_pore_volume_.resize(e.Nodes());
     psi1_.resize(e.Nodes());
     src_.resize(e.Nodes());

     // getting all the node-property related information
     // for all FV sectors = nodes = matrix columns
     // -------------------------------------------------
     for ( auto j{0U}; j<e.FV()->Sectors(); j++ ) //for ( size_t j{0U}; j<e.Nodes(); j++ )
     {
       // sector pore volumes
       sector_pore_volume_[j] = param.SectorVolume(j);

       // advected variable (transformed by bijective mapping)
       psi1_[j]               = e.N(j)->Read( adv1_key_ );

       // compute flux mismatch
       // for all sectors contained in this element
       // accumulate fluid sources (+) or sinks (-) due to deviations from
       // potentially non-conservative fluxes
       src_[j] = 0.;
       for ( auto k{0U}; k<e.FV()->FacetsPerSector(j); k++ ) {
         const auto n(e.FV()->FacetSurroundingSector(j,k));
         // if the sector node is the inside node then an incoming flux will create a positive source term
         const double fsign( (j == e.FV()->InsideNode(n) ) ? -1. : 1. );
         src_[j] += fsign * param.FacetArea(n) *
                    (param.FacetNormalVelocity(n) * dfds + param.FacetNormalComponent(n,v) * dGds);
       }
     }

     // -----------------------------------------------------------------------------
     // 2. accumulation of first-order facet fluxes (pressure- & gravity driven flow)
     //    (compare with Helmig, 97, p. 108, eqn. 3.74)
     // -----------------------------------------------------------------------------
     ipsi1_.resize( e.FV()->Facets() );
     facet_flux_.resize( e.FV()->Facets() );

     // interpolating values of advected quantity to facet integration points
     // --------------------------------------------------------------------
     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
     {
       // interpolating the property value to the (first) facet integration point
       ipsi1_[i] = (e).PropertyValueAtFacetIntegrationPoint( i, 0U, adv1_key_ );

       // integration of facet fluxes
       facet_flux_[i] = param.FacetNormalVelocity(i) * dfds;

       // ---------------------------------------------------------------------------
       // 3. gravitational segregation of wetting and non-wetting phases in zones of
       //    mixed saturation
       // ---------------------------------------------------------------------------
       // only if there is a density difference
       if ( fabs( relperm.DensityWettingPhase() - relperm.DensityNonWettingPhase() ) > numeric_limits<double>::epsilon() and
            fabs( param.FacetNormalComponent(i, v) ) > numeric_limits<double>::epsilon() )
       {
         // del_lambda_overbar / del_sn * k * -g(rho_w - rho_n)
         facet_flux_[i] += param.FacetNormalComponent(i, v) * dGds;
       }

       facet_flux_[i] *= param.FacetArea(i);

     }

  } // end AccumulateImplicitTwoPhaseSolution2 (generic version)







template<uint32_t dim>
void  StencilProcessor<dim>::AccumulateImplicitTwoPhaseSolution1_NonlinearNewtonRaphson(
                                          const FV_Parameter& param,
                                          const Element<dim>& e,
                                          TwoPhaseModel<dim>& relperm,
                                          bool with_gravity_forces,
                                          bool with_capillary_spreading)
 {

    eidx_ = e.Idx();

    bool DebugOutput(false);

    sector_pore_volume_.resize(e.Nodes());
    psi1_.resize(e.Nodes());
    // for all FV sectors = nodes = matrix columns
    for ( uint32_t j{0U}; j<e.Nodes(); j++ )
    {
        sector_pore_volume_[j] = param.SectorVolume(j);
        psi1_[j]               = e.N(j)->Read( adv1_key_ );
    }

    facet_flux_.resize( e.FV()->Facets());
    fill(facet_flux_.begin(),facet_flux_.end(),0.0);

    facet_flux_rhs_.resize( e.Nodes());
    fill(facet_flux_rhs_.begin(),facet_flux_rhs_.end(),0.0);

    upstream_node_.resize( e.FV()->Facets());
    fill(upstream_node_.begin(),upstream_node_.end(),0U);


    double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
    double upstream_fn(0.0),upstream_lambda_overbar(0.0);
    double upstream_sn(0.0),upstream_sw(0.0);
    double vn_at_facet_int_point(0.0),vw_at_facet_int_point(0.0);

    double linear_flux(0.0),nonlinear_flux(0.0);
    double viscous_velocity_component(0.0),gravity_velocity_component(0.0),capillary_velocity_component(0.0);
    double dsdn(0.0),dpcdn(0.0),dpcdsn(0.0);

    const double zero(0.);
    const size_t v( (dim==1u) ? 0u : 1u );  // gravity direction

    relperm.Initialize( e );
    relperm.InitializeForBaryCenter(e);
    relperm.EffectiveSaturation();

    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );


        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( uint32_t k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
        }

        fill( dpcdsn_.begin(), dpcdsn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             relperm.InitializeForNode( e, j );
             relperm.EffectiveSaturation();
             const double dpcds = -relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative ove sn is needed
             for ( auto k{0U}; k<dim; k++ ) dpcdsn_[k] += DN_(k,j) * dpcds;
        }

    }

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
      {

        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        // Average values and effective saturation
        // ------------------------------------------
        relperm.InitializeForNode( e, inside_node_ );
        relperm.EffectiveSaturation();
        const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
        const double ln_inside_node  = relperm.MobilityPhase(NONWETTING_PHASE);
        const double lw_inside_node  = relperm.MobilityPhase(WETTING_PHASE);

        relperm.InitializeForNode( e, outside_node_ );
        relperm.EffectiveSaturation();
        const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );
        const double ln_outside_node = relperm.MobilityPhase(NONWETTING_PHASE);
        const double lw_outside_node = relperm.MobilityPhase(WETTING_PHASE);


        relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
        relperm.EffectiveSaturation();

        // RHS: Nonlinear flux calculations

        viscous_velocity_component = 0.0;
        gravity_velocity_component = 0.0;
        capillary_velocity_component = 0.0;

        if( with_gravity_forces || with_capillary_spreading ){

            double vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );

            if( with_gravity_forces ){

                vn_gravity_component_of_velocity = relperm.MobilityPhase(NONWETTING_PHASE)*relperm.GravityTerm() * param.FacetNormalComponent(i,v);
                vw_gravity_component_of_velocity = relperm.MobilityPhase(WETTING_PHASE)*relperm.GravityTerm() * param.FacetNormalComponent(i,v);

            }

            if(with_capillary_spreading){

                dsdn  = param.FacetNormalProjection( i, dsdn_ );
                dpcdn = -dsdn*relperm.dpcds_Phase( );

                vn_capillary_component_of_velocity = relperm.MobilityPhase(NONWETTING_PHASE)*relperm.Permeability() * dpcdn;
                vw_capillary_component_of_velocity = relperm.MobilityPhase(WETTING_PHASE)*relperm.Permeability() * dpcdn;
            }

            vn_at_facet_int_point = param.FacetNormalVelocity(i) - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            vw_at_facet_int_point = param.FacetNormalVelocity(i) + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity;


            // mixture moving inside the CV        mixture moving outside the CV
            //
            //       |  vt                                 /|\ vt
            //       |                                      |
            //      \|/              /|\                    |
            //     -----              |  N_up             -----
            //   /       \                              /       \
            //   \       /                              \       /
            //     -----              |                   -----
            //      /|\              \|/  N_down            |
            //       |                                      |
            //       |  vt                                 \|/ vt

            if((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point>zero)){

                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_inside_node;

                upstream_sn = sn_inside_node;
                upstream_sw = 1.0 - sn_inside_node;

            }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero)){

                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_outside_node;

                upstream_sn = sn_inside_node;
                upstream_sw = 1.0 - sn_outside_node;

            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point>zero)){

                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_inside_node;

                upstream_sn = sn_outside_node;
                upstream_sw = 1.0 - sn_inside_node;

            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point<zero)){

                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_outside_node;

                upstream_sn = sn_outside_node;
                upstream_sw = 1.0 - sn_outside_node;

            }else{

                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);

                upstream_sn = 0.5*(sn_inside_node + sn_outside_node);
                upstream_sw = 1.0 - upstream_sn;

            }

            total_mobility=upstream_mobility_n+upstream_mobility_w;
            upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);

            viscous_velocity_component = upstream_fn * param.FacetNormalVelocity(i);

            if( with_gravity_forces )
                gravity_velocity_component = upstream_lambda_overbar * relperm.GravityTerm() * param.FacetNormalComponent(i,v);

            if( with_capillary_spreading )
                // Luat fix (2/9/2019) capillary_velocity_component = upstream_fn*vn_capillary_component_of_velocity;
                capillary_velocity_component = upstream_lambda_overbar * relperm.Permeability() * dpcdn;


        }else{

            if ( param.FacetNormalVelocity(i) != zero ) {

                 // identifying the upstream saturation (to get psi_hat_c, n_upstream = n_current )
                 upstream_sn  = ( (param.FacetNormalVelocity(i) < zero) ? sn_outside_node : sn_inside_node );
                 // fsat is sn interpolated to facet integration point by relperm model
                 // NB: Improve this value using up and downstream values together with Kurganov's scheme?
                 relperm.SaturationWettingPhase( 1. - upstream_sn );
                 relperm.EffectiveSaturation();

                 viscous_velocity_component = relperm.f_Phase(NONWETTING_PHASE/*non-wetting phase*/)* param.FacetNormalVelocity(i);

            }

        }

        nonlinear_flux = ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component)* param.FacetArea(i);

        // RHS fluxes:
        facet_flux_rhs_[inside_node_]  += nonlinear_flux;
        facet_flux_rhs_[outside_node_] -= nonlinear_flux;

        // LHS: Jacobian matrix calculations ( Linearized flux calculations )

        relperm.SaturationWettingPhase( upstream_sw );
        relperm.EffectiveSaturation();

        viscous_velocity_component = relperm.dfds() * param.FacetNormalVelocity(i);
        gravity_velocity_component = 0.0;
        capillary_velocity_component = 0.0;

        if( with_capillary_spreading){
            dpcdsn = param.FacetNormalProjection( i, dpcdsn_ );
            capillary_velocity_component = relperm.Permeability() * ( -relperm.dGds( ) * dpcdn  + relperm.G( )*dpcdsn );
        }

        if( with_gravity_forces )
            gravity_velocity_component = - relperm.GravityMultiplier_dGds( )*param.FacetNormalComponent(i,v);

        linear_flux = ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component ) * param.FacetArea(i);

        facet_flux_[i] += linear_flux;

        if( facet_flux_[i] > 0.0 ){

            upstream_node_[i] = inside_node_;

        }else{

            upstream_node_[i] = outside_node_;

        }


        if(DebugOutput){
            cout<<"So="<<relperm.Saturation(NONWETTING_PHASE)<<endl;
            cout<<"Sw="<<relperm.Saturation(WETTING_PHASE)<<endl;
            cout<<"So+Sw="<<relperm.Saturation(WETTING_PHASE)+relperm.Saturation(NONWETTING_PHASE)<<endl;
            cout<<"So_outside="<<sn_outside_node<<endl;
            cout<<"Sn_inside="<<sn_inside_node<<endl;
        }

     }

  } // end AccumulateImplicitTwoPhaseSolution1NonlinearNewtonRaphson





template<uint32_t dim>
void  StencilProcessor<dim>::AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson(
                                          const vector<pair<double,double> >& SMINMAX,
                                          const FV_Parameter& param,
                                          const Element<dim>& e,
                                          TwoPhaseModel<dim>& relperm,
                                          bool with_gravity_forces,
                                          bool with_capillary_spreading)
{

    eidx_ = e.Idx();

    bool DebugOutput(false);

    sector_pore_volume_.resize(e.Nodes());
    psi1_.resize(e.Nodes());
    // for all FV sectors = nodes = matrix columns
    for ( auto j{0U}; j<e.Nodes(); j++ )
    {
        sector_pore_volume_[j] = param.SectorVolume(j);
        psi1_[j]               = e.N(j)->Read( adv1_key_ );
    }

    facet_flux_.resize( e.FV()->Facets());
    fill(facet_flux_.begin(),facet_flux_.end(),0.0);

    facet_flux_rhs_.resize( e.Nodes());
    fill(facet_flux_rhs_.begin(),facet_flux_rhs_.end(),0.0);

    upstream_node_.resize( e.FV()->Facets());
    fill(upstream_node_.begin(),upstream_node_.end(),0U);


    double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
    double upstream_fn(0.0),upstream_lambda_overbar(0.0);
    double upstream_sn(0.0),upstream_sw(0.0);
    double vn_at_facet_int_point(0.0),vw_at_facet_int_point(0.0);

    double linear_flux(0.0),nonlinear_flux(0.0);
    double viscous_velocity_component(0.0),gravity_velocity_component(0.0),capillary_velocity_component(0.0);
    double dsdn(0.0),dpcdn(0.0),dpcdsn(0.0);

    const double zero(0.);
    const size_t v( (dim==1u) ? 0u : 1u );  // gravity direction

    relperm.Initialize( e );
    relperm.InitializeForBaryCenter(e);
    relperm.EffectiveSaturation();

    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );

        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
        }

        fill( dpcdsn_.begin(), dpcdsn_.end(), 0. );
        for ( uint32_t j{0U}; j<e.Nodes(); j++ ) {
             relperm.InitializeForNode( e, j );
             relperm.EffectiveSaturation();
             const double dpcds = -relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative over sn is needed
             for ( auto k{0U}; k<dim; k++ ) dpcdsn_[k] += DN_(k,j) * dpcds;
        }

    }

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
      {

        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        // Average values and effective saturation
        // ------------------------------------------
        // Saturation at the facet
        relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
        relperm.EffectiveSaturation();
        double sn_facet        = relperm.Saturation(NONWETTING_PHASE);

        // Saturation at the nodes
        const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
        const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );


        // #######################################################################
        // Limited saturation at the inside node
        // #######################################################################


        relperm.SaturationWettingPhase( 1. - limitProperty( sn_inside_node, sn_outside_node, sn_facet,
                                             SMINMAX[ e.N( inside_node_ )->Idx() ] ) );
        relperm.EffectiveSaturation();
        const double ln_inside_node  = relperm.MobilityPhase(NONWETTING_PHASE);
        const double lw_inside_node  = relperm.MobilityPhase(WETTING_PHASE);

        // #######################################################################
        // Limited saturation at the outside node
        // #######################################################################

        relperm.SaturationWettingPhase( 1. - limitProperty( sn_outside_node, sn_inside_node, sn_facet,
                                             SMINMAX[ e.N( outside_node_ )->Idx() ] ) );
        relperm.EffectiveSaturation();
        const double ln_outside_node = relperm.MobilityPhase(NONWETTING_PHASE);
        const double lw_outside_node = relperm.MobilityPhase(WETTING_PHASE);


        relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
        relperm.EffectiveSaturation();


        // RHS: Nonlinear flux calculations

        viscous_velocity_component = 0.0;
        gravity_velocity_component = 0.0;
        capillary_velocity_component = 0.0;

        if( with_gravity_forces || with_capillary_spreading ){

            double vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );

            if( with_gravity_forces ){

                vn_gravity_component_of_velocity = relperm.MobilityPhase(NONWETTING_PHASE)*relperm.GravityTerm() * param.FacetNormalComponent(i,v);
                vw_gravity_component_of_velocity = relperm.MobilityPhase(WETTING_PHASE)*relperm.GravityTerm() * param.FacetNormalComponent(i,v);

            }

            if(with_capillary_spreading){

                dsdn  = param.FacetNormalProjection( i, dsdn_ );
                dpcdn = -dsdn*relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative over sn is needed

                vn_capillary_component_of_velocity = relperm.MobilityPhase(NONWETTING_PHASE)*relperm.Permeability() * dpcdn;
                vw_capillary_component_of_velocity = relperm.MobilityPhase(WETTING_PHASE)*relperm.Permeability() * dpcdn;
            }

            vn_at_facet_int_point = param.FacetNormalVelocity(i) - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            vw_at_facet_int_point = param.FacetNormalVelocity(i) + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity;


            // mixture moving inside the CV        mixture moving outside the CV
            //
            //       |  vt                                 /|\ vt
            //       |                                      |
            //      \|/              /|\                    |
            //     -----              |  N_up             -----
            //   /       \                              /       \
            //   \       /                              \       /
            //     -----              |                   -----
            //      /|\              \|/  N_down            |
            //       |                                      |
            //       |  vt                                 \|/ vt

            if((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point>zero)){

                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_inside_node;

                upstream_sn = sn_inside_node;
                upstream_sw = 1.0 - sn_inside_node;

            }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero)){

                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_outside_node;

                upstream_sn = sn_inside_node;
                upstream_sw = 1.0 - sn_outside_node;

            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point>zero)){

                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_inside_node;

                upstream_sn = sn_outside_node;
                upstream_sw = 1.0 - sn_inside_node;

            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point<zero)){

                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_outside_node;

                upstream_sn = sn_outside_node;
                upstream_sw = 1.0 - sn_outside_node;

            }else{

                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);

                upstream_sn = 0.5*(sn_inside_node + sn_outside_node);
                upstream_sw = 1.0 - upstream_sn;

            }

            total_mobility=upstream_mobility_n+upstream_mobility_w;
            upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);

            viscous_velocity_component = upstream_fn * param.FacetNormalVelocity(i);

            if ( with_gravity_forces )
                gravity_velocity_component = upstream_lambda_overbar * relperm.GravityTerm() * param.FacetNormalComponent(i,v);

            if ( with_capillary_spreading )
                // Luat fix (2/9/2019) capillary_velocity_component = upstream_fn*vn_capillary_component_of_velocity;
                capillary_velocity_component = upstream_lambda_overbar * relperm.Permeability() * dpcdn;


        }else{

            if ( param.FacetNormalVelocity(i) != zero ) {

                 // identifying the upstream saturation (to get psi_hat_c, n_upstream = n_current )
                 upstream_sn  = ( (param.FacetNormalVelocity(i) < zero) ? sn_outside_node : sn_inside_node );
                 // fsat is sn interpolated to facet integration point by relperm model
                 // NB: Improve this value using up and downstream values together with Kurganov's scheme?
                 relperm.SaturationWettingPhase( 1. - upstream_sn );
                 relperm.EffectiveSaturation();

                 viscous_velocity_component = relperm.f_Phase(NONWETTING_PHASE/*non-wetting phase*/)* param.FacetNormalVelocity(i);

            }

        }

        nonlinear_flux = ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component)* param.FacetArea(i);

        // RHS fluxes:
        facet_flux_rhs_[inside_node_]  += nonlinear_flux;
        facet_flux_rhs_[outside_node_] -= nonlinear_flux;

        // LHS: Jacobian matrix calculations ( Linearized flux calculations )

        relperm.SaturationWettingPhase( upstream_sw );
        relperm.EffectiveSaturation();

        viscous_velocity_component = relperm.dfds() * param.FacetNormalVelocity(i);
        gravity_velocity_component = 0.0;
        capillary_velocity_component = 0.0;

        if( with_capillary_spreading){
            dpcdsn = param.FacetNormalProjection( i, dpcdsn_ );
            capillary_velocity_component = relperm.Permeability() * ( -relperm.dGds( ) * dpcdn  + relperm.G()*dpcdsn );
        }

        if( with_gravity_forces )
            gravity_velocity_component = - relperm.GravityMultiplier_dGds( )*param.FacetNormalComponent(i,v);

        linear_flux = ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component ) * param.FacetArea(i);

        facet_flux_[i] += linear_flux;

        if( facet_flux_[i] > 0.0 ){

            upstream_node_[i] = inside_node_;

        }else{

            upstream_node_[i] = outside_node_;

        }


        if(DebugOutput){
            cout<<"So="<<relperm.Saturation(NONWETTING_PHASE)<<endl;
            cout<<"Sw="<<relperm.Saturation(WETTING_PHASE)<<endl;
            cout<<"So+Sw="<<relperm.Saturation(WETTING_PHASE)+relperm.Saturation(NONWETTING_PHASE)<<endl;
            cout<<"So_outside="<<sn_outside_node<<endl;
            cout<<"Sn_inside="<<sn_inside_node<<endl;
        }

     }

  } // end AccumulateImplicitTwoPhaseSolution2NonlinearNewtonRaphson


template<uint32_t dim>
void  StencilProcessor<dim>::AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson(
                                          const vector<pair<double,double> >& SMINMAX,
                                          const FV_Parameter& param,
                                          const Element<dim>& e,
                                          TwoPhaseModel<dim>& relperm,
                                          bool with_gravity_forces,
                                          bool with_capillary_spreading,
                                          const csmp::Index& mass_center_key,
                                          const csmp::Index& grad_sn_key,
                                          const csmp::Index& grad_sn_limiter_key)
{

    eidx_ = e.Idx();

    bool DebugOutput(false);

    sector_pore_volume_.resize(e.Nodes());
    psi1_.resize(e.Nodes());
    // for all FV sectors = nodes = matrix columns
    for ( auto j{0U}; j<e.Nodes(); j++ )
    {
        sector_pore_volume_[j] = param.SectorVolume(j);
        psi1_[j]               = e.N(j)->Read( adv1_key_ );
    }

    facet_flux_.resize( e.FV()->Facets());
    fill(facet_flux_.begin(),facet_flux_.end(),0.0);

    facet_flux_rhs_.resize( e.Nodes());
    fill(facet_flux_rhs_.begin(),facet_flux_rhs_.end(),0.0);

    upstream_node_.resize( e.FV()->Facets());
    fill(upstream_node_.begin(),upstream_node_.end(),0U);


    double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
    double upstream_fn(0.0),upstream_lambda_overbar(0.0);
    double upstream_sn(0.0),upstream_sw(0.0);
    double vn_at_facet_int_point(0.0),vw_at_facet_int_point(0.0);

    double linear_flux(0.0),nonlinear_flux(0.0);
    double viscous_velocity_component(0.0),gravity_velocity_component(0.0),capillary_velocity_component(0.0);
    double dsdn(0.0),dpcdn(0.0),dpcdsn(0.0);

    double limited_sn_inside_node(1.0), limited_sn_outside_node(1.0);

    const double zero(0.);
    const size_t v( (dim==1u) ? 0u : 1u );  // gravity direction

    relperm.Initialize( e );
    relperm.InitializeForBaryCenter(e);
    relperm.EffectiveSaturation();

    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );

        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
        }

        fill( dpcdsn_.begin(), dpcdsn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             relperm.InitializeForNode( e, j );
             relperm.EffectiveSaturation();
             const double dpcds = -relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative over sn is needed
             for ( auto k{0U}; k<dim; k++ ) dpcdsn_[k] += DN_(k,j) * dpcds;
        }

    }

    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
      {

        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        // Average values and effective saturation
        // ------------------------------------------
        // Saturation at the facet
        relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
        relperm.EffectiveSaturation();

        // Saturation at the nodes
        const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
        const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );


        limitProperty_LSMGRAD<dim>( e, mass_center_key, grad_sn_key, grad_sn_limiter_key,
                                    inside_node_,outside_node_, i,
                                    sn_inside_node, sn_outside_node,
                                    limited_sn_inside_node, limited_sn_outside_node);


        // #######################################################################
        // Limited saturation at the inside node
        // #######################################################################


        relperm.SaturationWettingPhase( 1. - limited_sn_inside_node );
        relperm.EffectiveSaturation();
        const double ln_inside_node  = relperm.MobilityPhase(NONWETTING_PHASE);
        const double lw_inside_node  = relperm.MobilityPhase(WETTING_PHASE);

        // #######################################################################
        // Limited saturation at the outside node
        // #######################################################################

        relperm.SaturationWettingPhase( 1. - limited_sn_outside_node );
        relperm.EffectiveSaturation();
        const double ln_outside_node = relperm.MobilityPhase(NONWETTING_PHASE);
        const double lw_outside_node = relperm.MobilityPhase(WETTING_PHASE);


        // LHS: Jacobian matrix calculations ( Linearized flux calculations )

        relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
        relperm.EffectiveSaturation();

        // RHS: Nonlinear flux calculations

        viscous_velocity_component = 0.0;
        gravity_velocity_component = 0.0;
        capillary_velocity_component = 0.0;

        if( with_gravity_forces || with_capillary_spreading ){

            double vn_gravity_component_of_velocity( 0.0 ),vw_gravity_component_of_velocity( 0.0 );
            double vn_capillary_component_of_velocity ( 0.0 ), vw_capillary_component_of_velocity ( 0.0 );

            if( with_gravity_forces ){

                vn_gravity_component_of_velocity = relperm.MobilityPhase(NONWETTING_PHASE)*relperm.GravityTerm() * param.FacetNormalComponent(i,v);
                vw_gravity_component_of_velocity = relperm.MobilityPhase(WETTING_PHASE)*relperm.GravityTerm() * param.FacetNormalComponent(i,v);

            }

            if(with_capillary_spreading){

                dsdn  = param.FacetNormalProjection( i, dsdn_ );
                dpcdn = -dsdn*relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative over sn is needed

                vn_capillary_component_of_velocity = relperm.MobilityPhase(NONWETTING_PHASE)*relperm.Permeability() * dpcdn;
                vw_capillary_component_of_velocity = relperm.MobilityPhase(WETTING_PHASE)*relperm.Permeability() * dpcdn;
            }

            vn_at_facet_int_point = param.FacetNormalVelocity(i) - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
            vw_at_facet_int_point = param.FacetNormalVelocity(i) + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity;


            // mixture moving inside the CV        mixture moving outside the CV
            //
            //       |  vt                                 /|\ vt
            //       |                                      |
            //      \|/              /|\                    |
            //     -----              |  N_up             -----
            //   /       \                              /       \
            //   \       /                              \       /
            //     -----              |                   -----
            //      /|\              \|/  N_down            |
            //       |                                      |
            //       |  vt                                 \|/ vt

            if((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point>zero)){

                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_inside_node;

                upstream_sn = sn_inside_node;
                upstream_sw = 1.0 - sn_inside_node;

            }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero)){

                upstream_mobility_n=ln_inside_node;
                upstream_mobility_w=lw_outside_node;

                upstream_sn = sn_inside_node;
                upstream_sw = 1.0 - sn_outside_node;

            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point>zero)){

                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_inside_node;

                upstream_sn = sn_outside_node;
                upstream_sw = 1.0 - sn_inside_node;

            }else if ((vn_at_facet_int_point<zero)&&(vw_at_facet_int_point<zero)){

                upstream_mobility_n=ln_outside_node;
                upstream_mobility_w=lw_outside_node;

                upstream_sn = sn_outside_node;
                upstream_sw = 1.0 - sn_outside_node;

            }else{

                upstream_mobility_n=0.5*(ln_inside_node+ln_outside_node);
                upstream_mobility_w=0.5*(lw_inside_node+lw_outside_node);

                upstream_sn = 0.5*(sn_inside_node + sn_outside_node);
                upstream_sw = 1.0 - upstream_sn;

            }

            total_mobility=upstream_mobility_n+upstream_mobility_w;
            upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
            upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);

            viscous_velocity_component = upstream_fn * param.FacetNormalVelocity(i);

            if( with_gravity_forces )
                gravity_velocity_component = upstream_lambda_overbar * relperm.GravityTerm() * param.FacetNormalComponent(i,v);

            if ( with_capillary_spreading )
              // Luat fix (2/9/2019) capillary_velocity_component = upstream_fn*vn_capillary_component_of_velocity;
              capillary_velocity_component = upstream_lambda_overbar * relperm.Permeability() * dpcdn;


        }else{

            if ( param.FacetNormalVelocity(i) != zero ) {

                 // identifying the upstream saturation (to get psi_hat_c, n_upstream = n_current )
                 upstream_sn  = ( (param.FacetNormalVelocity(i) < zero) ? sn_outside_node : sn_inside_node );
                 // fsat is sn interpolated to facet integration point by relperm model
                 // NB: Improve this value using up and downstream values together with Kurganov's scheme?
                 relperm.SaturationWettingPhase( 1. - upstream_sn );
                 relperm.EffectiveSaturation();

                 viscous_velocity_component = relperm.f_Phase(NONWETTING_PHASE/*non-wetting phase*/)* param.FacetNormalVelocity(i);

            }

        }

        nonlinear_flux = ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component)* param.FacetArea(i);

        // RHS fluxes:
        facet_flux_rhs_[inside_node_]  += nonlinear_flux;
        facet_flux_rhs_[outside_node_] -= nonlinear_flux;

        // LHS: Jacobian matrix calculations ( Linearized flux calculations )

        relperm.SaturationWettingPhase( upstream_sw );
        relperm.EffectiveSaturation();

        viscous_velocity_component = relperm.dfds() * param.FacetNormalVelocity(i);
        gravity_velocity_component = 0.0;
        capillary_velocity_component = 0.0;

        if( with_capillary_spreading){
            dpcdsn = param.FacetNormalProjection( i, dpcdsn_ );
            capillary_velocity_component = relperm.Permeability() * ( -relperm.dGds( ) * dpcdn  + relperm.G()*dpcdsn );
        }

        if( with_gravity_forces )
            gravity_velocity_component = - relperm.GravityMultiplier_dGds( )*param.FacetNormalComponent(i,v);

        linear_flux = ( viscous_velocity_component - gravity_velocity_component - capillary_velocity_component ) * param.FacetArea(i);

        facet_flux_[i] += linear_flux;

        if( facet_flux_[i] > 0.0 ){

            upstream_node_[i] = inside_node_;

        }else{

            upstream_node_[i] = outside_node_;

        }


        if(DebugOutput){
            cout<<"So="<<relperm.Saturation(NONWETTING_PHASE)<<endl;
            cout<<"Sw="<<relperm.Saturation(WETTING_PHASE)<<endl;
            cout<<"So+Sw="<<relperm.Saturation(WETTING_PHASE)+relperm.Saturation(NONWETTING_PHASE)<<endl;
            cout<<"So_outside="<<sn_outside_node<<endl;
            cout<<"Sn_inside="<<sn_inside_node<<endl;
        }

     }

  } // end AccumulateImplicitTwoPhaseSolution2NonlinearNewtonRaphson


template<uint32_t dim>
void  StencilProcessor<dim>::CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson(
                                                        const FV_Parameter& param,
                                                        const Element<dim>& e,
                                                        TwoPhaseModel<dim>& relperm,
                                                        uint32_t pnid,
                                                        double& flux,
                                                        bool with_gravity_forces,
                                                        bool with_capillary_spreading)
 {
    eidx_ = e.Idx();

    relperm.Initialize( e );
    relperm.InitializeForNode( e, pnid );
    relperm.EffectiveSaturation();

    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> velo;
    double facetArea;

    /*
    // Accumulate sources for correction at the boundary for LHS
    src_[pnid] = 0.;
    rhs_src_[pnid] = 0.;
    relperm.InitializeForNode( e, pnid );
    relperm.EffectiveSaturation();
    for ( size_t k{0U}; k<e.FV()->FacetsPerSector(pnid); k++ ) {
          const size_t i(e.FV()->FacetSurroundingSector(pnid,k));
          // if the sector node is the inside node then an incoming flux will create a positive source term
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          const double sn_inside_node  = e.N(inside_node_)->Read( adv1key_ );
          const double sn_outside_node = e.N(outside_node_)->Read( adv1key_ );
          relperm.InitializeForFacetIntegrationPoint( i, 0U, e ); // costly (interpolates rho and mu as well)
          relperm.EffectiveSaturation();

          const double lhs_sign( (pnid == e.FV()->InsideNode(i)) ? -1. : 1. );
          const double rhs_sign( (pnid == e.FV()->InsideNode(i)) ? -1. : 1. );
          src_[pnid] += lhs_sign * relperm.dfds() * (param.FacetNormalVelocity(i)) * param.FacetArea(i);
          rhs_src_[pnid] += rhs_sign * (relperm.f_Phase(NONWETTING_PHASE) * param.FacetNormalVelocity(i)) * param.FacetArea(i);
    }
    */


    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );
        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
        }

    }

    // now the saturation dependent properties are computed
    // for all FACETS per SECTOR surrounding the finite volume at the boundary
    for ( auto k{0U}; k<e.FV()->FacetsPerSector(pnid); k++ )
      {
         auto i( e.FV()->FacetSurroundingSector(pnid,k) );

         e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
         // get velocity across FV facet

         Point<dim>  n =e.FacetNormal(i);
         e.Read( vel_key_, velo );
         facetArea = e.FacetArea(i);

         const double viscous_vel_component = velo.DotProduct(n);

         if( with_gravity_forces){

             const double gravity_vel_component = n[v];

             if( with_capillary_spreading ){

                 const double capillary_vel_component( -param.FacetNormalProjection( i, dsdn_ )*relperm.dpcds_Phase( ) );

                 if ( pnid == inside_node_ )
                     flux += (relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component  - relperm.GravityMultiplier_G()*gravity_vel_component - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;
                 else
                     flux -= (relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component  - relperm.GravityMultiplier_G()*gravity_vel_component - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;

             }else{

                 if ( pnid == inside_node_ )
                     flux += ( relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component - relperm.GravityMultiplier_G()*gravity_vel_component )*facetArea;
                 else
                     flux -= ( relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component - relperm.GravityMultiplier_G()*gravity_vel_component )*facetArea;

             }

         }else{

             if(with_capillary_spreading){

                 const double capillary_vel_component( -param.FacetNormalProjection( i, dsdn_ )*relperm.dpcds_Phase( ) );

                 if ( pnid == inside_node_ )
                     flux += (relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component  - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;
                 else
                     flux -= (relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component  - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;

             }else{

                 if ( pnid == inside_node_ )
                     flux += relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component*facetArea;
                 else
                     flux -= relperm.f_Phase(NONWETTING_PHASE) * viscous_vel_component*facetArea;

             }

         }
      }

 } // end CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson





/** Outputs the current internal variables of the FV sector to screen.
*/
template<uint32_t dim>
void  StencilProcessor<dim>::Out() const
 {
    cout <<"\n\n\nImplicitStencilProcessor<"<< dim <<">::Out: "<< endl;
    // outputting data of object
    cout <<"Finite volume stencil data."<< endl;
    cout <<"\n\tassociated sectors, pore volumes, and advected property values: ";
    for ( size_t i{0U}; i<sector_pore_volume_.size(); i++ )
      cout <<"\n\tsector "<< i+1 <<": "<< sector_pore_volume_[i] <<", "<< psi1_[i];

    cout <<"\n\nFacets between nodes i, i+1, and associated properties:";
    cout <<"\n\tintegrated volume flux and advected variable interpolated to integration points:";
    for ( size_t i{0U}; i<facet_flux_.size(); i++ )
    cout <<"\n\tfacet "<< i+1 <<": "<< facet_flux_[i] <<", "<< ipsi1_[i];

    if ( !src_.empty() ) {
         cout <<"\n\nFacet integrated source terms:";
         for ( size_t i{0U}; i<src_.size(); i++ )
         cout <<"\n\tfacet "<< i+1 <<": "<< src_[i];
      }

    if ( !ipsi1_.empty() ) {
         cout <<"\n\nadvected variable values at facet integration points: ";
         for ( size_t i{0U}; i<ipsi1_.size(); i++ )
         cout <<"\n\tfacet "<< i+1 <<": "<< ipsi1_[i];
      }

    cout << endl << endl;

 } // end Out


template struct StencilProcessor<1U>;
template struct StencilProcessor<2U>;
template struct StencilProcessor<3U>;


} // namespace csmp

