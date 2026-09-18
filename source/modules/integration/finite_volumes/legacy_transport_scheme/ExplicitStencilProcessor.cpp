// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ExplicitStencilProcessor.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "Element.h"
#include "FV_Parameter.h"
#include "TwoPhaseModel.h"

using namespace std;

#define DEBUG_ExplicitStencilProcessor

namespace csmp {

template<uint32_t dim>
ExplicitStencilProcessor<dim>::ExplicitStencilProcessor( const csmp::Index& adv_key,   
                                                         const csmp::Index& velo_key )
 : adv1_key_(adv_key),
   vel_key_(velo_key),
   diff_key_(csmp::Index()),
   grad_(dim),
   pc_grad_(dim),
   dsdn_(dim)
 {
 }


template<uint32_t dim>
ExplicitStencilProcessor<dim>::ExplicitStencilProcessor( const csmp::Index& adv_key,  
                                                         const csmp::Index& velo_key,
                                                         const csmp::Index& diffusivity_key )
 : adv1_key_(adv_key),
   vel_key_(velo_key),
   diff_key_(diffusivity_key),
   grad_(dim),
   pc_grad_(dim),
   dsdn_(dim)
 {
 }
 
template<uint32_t dim>
ExplicitStencilProcessor<dim>::ExplicitStencilProcessor()
 : adv1_key_(),
   vel_key_(),
   diff_key_(),
   grad_(dim),
   pc_grad_(dim),
   dsdn_(dim)
 {
 }
 



///  standard assignment operator  
/*
template<uint32_t dim>
ExplicitStencilProcessor<dim>& ExplicitStencilProcessor<dim>::operator=( const ExplicitStencilProcessor& tfs )
 {
    if ( &tfs != this ) {
         sector_pore_volume_ = tfs.sector_pore_volume_;
         facet_flux_         = tfs.facet_flux_;
         psi1_               = tfs.psi1_;
         src_                = tfs.src_;
         grad_               = tfs.grad_;
         pc_grad_            = tfs.pc_grad_;
         dsdn_               = tfs.dsdn_;

         // property keys
  	     adv1key_            = tfs.adv1key_;
  	     vel0key_            = tfs.vel0key_;
  	     diffkey_            = tfs.diffkey_;  	         
      }
    return *this;
 }


 
 
template<uint32_t dim>
ExplicitStencilProcessor<dim>::ExplicitStencilProcessor( const ExplicitStencilProcessor<dim>& tfs )
 {
    *this = tfs;
 }
 
*/

 
///  destructor
template<uint32_t dim>
ExplicitStencilProcessor<dim>::~ExplicitStencilProcessor()
 {
 } 

    

/**

Initializes sector pore volume and facet flux vectors for
current finite element = finite volume stencil.
*/
template<uint32_t dim>
void ExplicitStencilProcessor<dim>::InitializeFirstOrder( const FV_Parameter& param,
                                                          const Element<dim>& e , const VARIABLE_TYPE &vt, uint32_t var_comp_nr )
 {
     sector_pore_volume_.resize(e.Nodes());
     psi1_.resize(e.Nodes());

     if (vt==SCALAR){
         // getting all the node-related information
         // ----------------------------------------
         for ( auto i{0U}; i<e.Nodes(); i++ ) {
             // sector pore volumes
             sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable
             psi1_[i] = e.N(i)->Read( adv1_key_ );
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
         }
     }
     else if (vt==ARRAY){
         ArrayVariable av;
         for ( auto i{0U}; i<e.Nodes(); i++ ) {
             // sector pore volumes
             sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable
             e.N(i)->Read( adv1_key_,av );
             psi1_[i] = av(var_comp_nr);
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
         }
     }
     if (vt==FLAGGEDARRAY){
         FlaggedArrayVariable fav;
         for ( uint32_t i{0U}; i<e.Nodes(); i++ ) {
             // sector pore volumes
             sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable
             e.N(i)->Read( adv1_key_, fav );
             psi1_[i] = fav(var_comp_nr);
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
         }
     }

 } // end InitializeFirstOrder




/**

Accumulates the second-order-in-space accurate explicit solution into the result vector
using the MINMOD limiter to prevent oscillations.

*/
template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitAdvectionSolution2( 
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     vector<double>& res ) 
 {
     eidx_ = e.Idx();
            
     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          const double psi_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double psi_outside_node = e.N(outside_node_)->Read( adv1_key_ );
    
          // ------------------------------------------
          // 1. accumulation of linear facet flux, O.K.
          // ------------------------------------------
          if ( fabs(param.FacetNormalVelocity(i)) > numeric_limits<double>::epsilon() ) {

              double psi_facet   = ( e ).PropertyValueAtFacetIntegrationPoint( i, 0U, adv1_key_ );

              double linear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i);
	           
              if ( param.FacetNormalVelocity(i) > 0. ) {
	                linear_flux *= limitProperty( psi_inside_node, psi_outside_node, psi_facet, 
                                                SMINMAX[ e.N( e.FV()->InsideNode(i) )->Idx() ] );
	             } 
              else {
	                linear_flux *= limitProperty( psi_outside_node, psi_inside_node, psi_facet, 
                                                SMINMAX[ e.N( e.FV()->OutsideNode(i) )->Idx() ] );
	             }  
          
             res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += linear_flux;
             res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= linear_flux;
	       }
       }
  
} // AccumulateExplicitAdvectionDiffusionSolution2



template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitAdvectionSolution2(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     vector<double>& res,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_psi_key,
                                     const csmp::Index& grad_psi_limiter_key)
 {
     eidx_ = e.Idx();

     double  limited_psi_inside_node, limited_psi_outside_node;

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          const double psi_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double psi_outside_node = e.N(outside_node_)->Read( adv1_key_ );

          // ------------------------------------------
          // 1. accumulation of linear facet flux, O.K.
          // ------------------------------------------
          if ( fabs(param.FacetNormalVelocity(i)) > numeric_limits<double>::epsilon() ) {

              double linear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i);

              limitProperty_LSMGRAD<dim>( e, mass_center_key, grad_psi_key, grad_psi_limiter_key,
                                          inside_node_,outside_node_, i,
                                          psi_inside_node, psi_outside_node,
                                          limited_psi_inside_node, limited_psi_outside_node);

               if ( param.FacetNormalVelocity(i)> 0. )
                   linear_flux*=limited_psi_inside_node;
               else
                   linear_flux*=limited_psi_outside_node;


             res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += linear_flux;
             res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= linear_flux;
           }
       }

} // AccumulateExplicitAdvectionDiffusionSolution2 with lsm grad limiter









/**  Accumulates first-order-in-space accurate solution into res vector
*/
template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitAdvectionDiffusionSolution1(
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     vector<double>& res,
                                     const VARIABLE_TYPE &vt, 
                                     uint32_t var_comp_nr )
 {
     const double  zero(0.);

     double  linear_flux;

     eidx_ = e.Idx();
     // ---------------------------------------------------------
     // 0. computing advected variable gradient from psi at nodes
     //    and interpolating psi to facet integration point
     // ---------------------------------------------------------

     e.dN_AtBaryCenter( DN_ );
     double psi_node = std::numeric_limits<double>::quiet_NaN();
     double psi_inside_node  = std::numeric_limits<double>::quiet_NaN();
     double psi_outside_node = std::numeric_limits<double>::quiet_NaN();

     fill( grad_.begin(), grad_.end(), zero );
     for ( auto j{0U}; j<DN_.Cols(); j++ ) {
          if (vt == SCALAR)  
            psi_node = e.N(j)->Read( adv1_key_ );
          else if (vt == ARRAY)
          {
            ArrayVariable av;
            e.N(j)->Read( adv1_key_, av );
            psi_node = av[var_comp_nr];
          }
          else if (vt == FLAGGEDARRAY)
          {
            FlaggedArrayVariable fav;
            e.N(j)->Read( adv1_key_, fav );
            psi_node = fav[var_comp_nr];
          }

          for ( auto k{0U}; k<dim; k++ ) grad_[k] += DN_(k,j) * psi_node;
     }

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          if (vt == SCALAR)
            {
                psi_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
                psi_outside_node = e.N(outside_node_)->Read( adv1_key_ );
            }
          else if (vt == ARRAY)
            {
                ArrayVariable av;
                e.N(inside_node_)->Read( adv1_key_, av );
                psi_inside_node  = av[var_comp_nr]; 
                e.N(outside_node_)->Read( adv1_key_, av );                
                psi_outside_node = av[var_comp_nr];
            }
          else if (vt == FLAGGEDARRAY)
            {
                FlaggedArrayVariable fav;
                e.N(inside_node_)->Read( adv1_key_, fav );
                psi_inside_node  = fav[var_comp_nr]; 
                e.N(outside_node_)->Read( adv1_key_, fav );                
                psi_outside_node = fav[var_comp_nr];
            }

          // ------------------------------------------
          // 1. accumulation of linear facet flux, O.K.
          // ------------------------------------------
          if ( param.FacetNormalVelocity(i) != zero ) {

               linear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i);

               if ( param.FacetNormalVelocity(i) > zero )
                    linear_flux *= psi_inside_node;
               else
                    linear_flux *= psi_outside_node;
          }
          else linear_flux = zero;

          // -----------------------------------------
          // 2. diffusion (1st order flux)
          // -----------------------------------------
          // NB: a special boundary treatment seems to be required to avoid oscillations which arise at small gradients
          const double diffusive_flux = e.Read( diff_key_ ) * -param.FacetNormalProjection( i, grad_ )*param.FacetArea(i);

          res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += (linear_flux + diffusive_flux);
          res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= (linear_flux + diffusive_flux);
       }

} // AccumulateExplicitAdvectionDiffusionSolution1



/**  Accumulates second-order-in-space accurate solution into res vector
*/
template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitAdvectionDiffusionSolution2( 
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     vector<double>& res)
 {
     const double  zero(0.);

     double  linear_flux;

     eidx_ = e.Idx();
     // ---------------------------------------------------------
     // 0. computing advected variable gradient from psi at nodes
     //    and interpolating psi to facet integration point
     // ---------------------------------------------------------

     e.dN_AtBaryCenter( DN_ );

     fill( grad_.begin(), grad_.end(), zero );
     for ( auto j{0U}; j<DN_.Cols(); j++ ) {
          double psi_node = e.N(j)->Read( adv1_key_ );
          for ( auto k{0U}; k<dim; k++ ) grad_[k] += DN_(k,j) * psi_node;
     }
            
     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          const double psi_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double psi_outside_node = e.N(outside_node_)->Read( adv1_key_ );
    
          // ------------------------------------------
          // 1. accumulation of linear facet flux, O.K.
          // ------------------------------------------
          if ( param.FacetNormalVelocity(i) != zero ) {

               double psi_facet = (e).PropertyValueAtFacetIntegrationPoint( i, 0U, adv1_key_ );

               linear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i);
	           
               if ( param.FacetNormalVelocity(i) > zero )
	                linear_flux *= limitProperty( psi_inside_node, psi_outside_node, psi_facet, 
                                                SMINMAX[ e.N( e.FV()->InsideNode(i) )->Idx() ] );
               else
                    linear_flux *= limitProperty( psi_outside_node, psi_inside_node, psi_facet,
                                                SMINMAX[ e.N( e.FV()->OutsideNode(i) )->Idx() ] );
          }
          else linear_flux = zero;

          // -----------------------------------------
          // 2. diffusion (1st order flux)
          // -----------------------------------------
          // NB: a special boundary treatment seems to be required to avoid oscillations which arise at small gradients
          const double diffusive_flux = e.Read( diff_key_ ) * -param.FacetNormalProjection( i, grad_ )*param.FacetArea(i);

          res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += (linear_flux + diffusive_flux);
          res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= (linear_flux + diffusive_flux);
       }
  
} // AccumulateExplicitAdvectionDiffusionSolution2


template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitAdvectionDiffusionSolution2(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     vector<double>& res,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_psi_key,
                                     const csmp::Index& grad_psi_limiter_key)
 {
     const double  zero(0.);
     double  linear_flux;
     double  limited_psi_inside_node, limited_psi_outside_node;

     eidx_ = e.Idx();
     // ---------------------------------------------------------
     // 0. computing advected variable gradient from psi at nodes
     //    and interpolating psi to facet integration point
     // ---------------------------------------------------------

     e.dN_AtBaryCenter( DN_ );

     fill( grad_.begin(), grad_.end(), zero );
     for ( auto j{0U}; j<DN_.Cols(); j++ ) {
          double psi_node = e.N(j)->Read( adv1_key_ );
          for ( uint32_t k{0U}; k<dim; k++ ) grad_[k] += DN_(k,j) * psi_node;
     }

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          const double psi_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double psi_outside_node = e.N(outside_node_)->Read( adv1_key_ );

          // ------------------------------------------
          // 1. accumulation of linear facet flux, O.K.
          // ------------------------------------------
          if ( param.FacetNormalVelocity(i) != zero ) {

             linear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i);

             limitProperty_LSMGRAD<dim>( e, mass_center_key, grad_psi_key, grad_psi_limiter_key,
                                         inside_node_,outside_node_, i,
                                         psi_inside_node, psi_outside_node,
                                         limited_psi_inside_node, limited_psi_outside_node);

             if ( param.FacetNormalVelocity(i)> 0. )
                 linear_flux*=limited_psi_inside_node;
             else
                 linear_flux*=limited_psi_outside_node;

          }
          else linear_flux = zero;

          // -----------------------------------------
          // 2. diffusion (1st order flux)
          // -----------------------------------------
          // NB: a special boundary treatment seems to be required to avoid oscillations which arise at small gradients
          const double diffusive_flux = e.Read( diff_key_ ) * param.FacetArea(i) * -param.FacetNormalProjection( i, grad_ );
          res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += (linear_flux + diffusive_flux);
          res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= (linear_flux + diffusive_flux);
       }

} // AccumulateExplicitAdvectionDiffusionSolution2 with lsm grad limiter













/**

1st-order upwind-based advection scheme including gravitational
and capillary forces, see Helmig Buch, p. 108, eqn. 3.74.  

@section arguments Input Arguments 

See argument list of method. The base class of a two-phase 
relative permeability model is used to compute the non-linear flow 
multipliers from Helmig's formulation.  

Entries into the results and source term vectors for the finite volumes
that are part of the current finite-element finite-volume stencil.  

@section @section implementation Implementation

Key Issue: What saturation value is used to calculate the nonlinear fluxes 
(see Kurganov) is crucial for the quality of the results and the mesh
resolution dependency of the scheme.

If fluxes or facet integration point saturations are estimated well, the 
first-order scheme will outperform a 2nd-order one which uses linear 
estimates of the fluxes or facet saturations.  

@section @section application Application

In IMPES formulations of transport problems.  

@attention Comments: 

To decrease the mesh sensitivity of the scheme one would need to improve 
estimate of fsat in the presence of shocks:  

- if shock is between upstream and downstream node (see their saturations)
  speed must be used to assign sn=sn_shock or sn_d.  
  
- if sat range is in rarefaction fan, the normal interpolation should be fine
   NB: Kurganov model might be suitable as well or it may be selectively applied
       only when the shock is between the upstream and downstream nodes.           
*/
template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution1_Visc(
                                                        const FV_Parameter& param,
                                                        const Element<dim>& e,
                                                        TwoPhaseModel<dim>& relperm,
                                                        vector<double>& res)
 {
     const double zero(0.);
     // irreducible saturation water
     // irreducible saturation oil
     // Brooks Corey parameter
     relperm.Initialize( e );
     
     // computing capillary pressure gradient from pc at nodes
     eidx_ = e.Idx();

     e.dN_AtBaryCenter( DN_ );
     fill( dsdn_.begin(), dsdn_.end(), 0. );
     for ( auto j{0U}; j<e.Nodes(); j++ ) {
          const double sn = e.N(j)->Read( adv1_key_);
          for ( auto k=0; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
     }

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );
          
          // 0. average values and effective saturation
          // ------------------------------------------
          relperm.InitializeForFacetIntegrationPoint( i, 0U, e ); // costly (interpolates rho and mu as well)
          relperm.EffectiveSaturation();        


#ifdef DEBUG_ExplicitStencilProcessor

          assert( !isnan(relperm.AdvectionMultiplier()) );
          assert( !isnan(relperm.CapillaryDiffusionMultiplier( )) );
          assert( !isnan(relperm.GravityMultiplier_G()) );
#endif          

          // ---------------------------------------------
          // 1. (non-linear) (capillary) diffusion (term1)
          // ---------------------------------------------

          /// Note: this upwind technique doesn't work for capillary term, doesn't take into account co-current flow
          ///       there better upwind strategy can be found in function AccumulateExplicitTwoPhaseSolution1 with grav. and cap. boolean flags
          /// (Roman,2013)

          // compute first order diffusive flux
          // get length of edge between in and outside node

          /*
          if( with_capillary_spreading ){

              const double ds = (sn_outside_node - sn_inside_node);

              if ( ds != zero ) {

                   const double dsdn(param.FacetNormalProjection( i, dsdn_ ));

                   // identifying the upstream node
                   // flow will be out of a facet if proj. of sat-gradient onto normal is negative
                   const double  psi_hat_c( (dsdn < zero) ? sn_inside_node : sn_outside_node );

                   // computing flow properties at upstream node
                   relperm.SaturationWettingPhase( 1. - psi_hat_c );
                   relperm.EffectiveSaturation();

                   const double capillary_flux = -dsdn * relperm.CapillaryDiffusionMultiplier(NONWETTING_PHASE) * param.FacetArea(i);

                   res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += capillary_flux;
                   res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= capillary_flux;
                }
          }
          */

          // ------------------------------------------------------------------------------
          // 2. gravitational segregation in zones of mixed saturation (2nd part of term 2)
          // ------------------------------------------------------------------------------

          /// Note: this upwind technique doesn't work for gravity term, doesn't take into account co-current flow as well as capillary flow
          ///       there better upwind strategy can be found in function AccumulateExplicitTwoPhaseSolution1 with grav. and cap. boolean flags
          /// (Roman, 2013)
          /*
          const size_t v( (dim==1u) ? 0u : 1u );
          // check whether there is gravity driven flow
          if ( param.FacetNormalComponent(i,v) != zero and
               relperm.DensityWettingPhase()   != relperm.DensityNonWettingPhase() ) 
            {
 	           // identifying the upstream node
 	           // flow of the lighter phase will always be out of a facet whose normal
 	           // is pointing upward
               const double  psi_hat_c( (param.FacetNormalComponent(i,v) < zero) ? sn_outside_node : sn_inside_node );

               // computing flow properties at upstream node
               relperm.SaturationWettingPhase( 1. - psi_hat_c );
               relperm.EffectiveSaturation();

               const double gravity_induced_flux = param.FacetNormalComponent(i,v) * param.FacetArea(i) * 
                                                     relperm.GravityMultiplier_G();

               // 2nd-order fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += gravity_induced_flux; 
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= gravity_induced_flux;
            }

          */

          // --------------------------------------------------------------
          // 3. accumulation of nonlinear facet fluxes (1st part of term 2)
          // --------------------------------------------------------------
          if ( param.FacetNormalVelocity(i) != zero ) {
 	           // identifying the upstream saturation (to get psi_hat_c, n_upstream = n_current )
               const double  psi_hat_c( (param.FacetNormalVelocity(i) < zero) ? sn_outside_node : sn_inside_node );
               // fsat is sn interpolated to facet integration point by relperm model
               // NB: Improve this value using up and downstream values together with Kurganov's scheme?
	           relperm.SaturationWettingPhase( 1. - psi_hat_c );
	           relperm.EffectiveSaturation();        

               // try whether upcasting does not work
               const double nonlinear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i) * relperm.f_Phase(2u);

	           // 2nd-order fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
	           res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += nonlinear_flux; 
	           res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= nonlinear_flux;
            }
      }
  
 } // end AccumulateExplicitTwoPhaseSolution1 ( generic fully upwind version 1st-order )




/**

Second-order accurate solution with advection and source term limiting
(Geiger, 2004).

@section arguments Input Arguments


@section implementation Implementation

The following capillary flow implementation also works fine but
needs extra source term vector<double> for accumulation:

@code
{
   capillary_flux *= param.FacetArea(i) * relperm.CapillaryDiffusionMultiplier(NONWETTING_PHASE);
   const double sn_ltd = LimitExplicitSourceTerm( mapping, SMINMAX,
                                              sn_inside_node, sn_outside_node,
                                              sn_facet, capillary_flux );

   src[ mapping.FiniteVolume(e.N(inside_node_)->ID()) ]  += capillary_flux * sn_ltd;
   src[ mapping.FiniteVolume(e.N(outside_node_)->ID()) ] -= capillary_flux * sn_ltd;
}
@endcode

@attention Comments:

Method is not 100% stable where flow occurs across porosity jumps.

AssignFractionalFlowBoundaryConditions() shows how to compensate these fluxes
but it does it only for the divergence recorded by the FLUX_BALANCE
*/


template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution2_Visc(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     TwoPhaseModel<dim>& relperm,
                                     vector<double>& res)
 {
     eidx_ = e.Idx();

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          // --------------------------------------------------------
          // 0. facet saturation and relperm
          // --------------------------------------------------------
          relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
          relperm.EffectiveSaturation();
          double sn_facet        = relperm.Saturation(NONWETTING_PHASE);
          // get non-wetting-phase saturations for the inside and the outside nodes
          const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );

          // -----------------------------------------
          // 1. accumulation of nonlinear facet fluxes
          // -----------------------------------------
          if ( fabs(param.FacetNormalVelocity(i)) > numeric_limits<double>::epsilon() )
            {
               // identifying the upstream node (to get psi_hat_c, n_upstream = n_current )
               if ( param.FacetNormalVelocity(i)> 0. )
                 relperm.SaturationWettingPhase( 1. - limitProperty( sn_inside_node, sn_outside_node, sn_facet,
                                                                     SMINMAX[ e.N( e.FV()->InsideNode(i) )->Idx() ] ) );
               else
                 relperm.SaturationWettingPhase( 1. - limitProperty( sn_outside_node, sn_inside_node, sn_facet,
                                                                     SMINMAX[ e.N( e.FV()->OutsideNode(i) )->Idx() ] ) );
               // computing flow properties at upstream node
               relperm.EffectiveSaturation();
               // + div . [fn vt]
               const double nonlinear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i) * relperm.f_Phase(NONWETTING_PHASE);

               // 2nd-order fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += nonlinear_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= nonlinear_flux;
             }

      }

 } // end AccumulateTwoPhaseSolution2 ( generic fully upwind version, but without capillary displacement )



template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution2_Visc(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     TwoPhaseModel<dim>& relperm,
                                     vector<double>& res,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_sn_key,
                                     const csmp::Index& grad_sn_limiter_key)
 {
     eidx_ = e.Idx();

     double limited_sn_inside_node(1.0), limited_sn_outside_node(1.0);

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          // --------------------------------------------------------
          // 0. facet saturation and relperm
          // --------------------------------------------------------
          relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
          relperm.EffectiveSaturation();

          // get non-wetting-phase saturations for the inside and the outside nodes
          const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );


          // -----------------------------------------
          // 1. accumulation of nonlinear facet fluxes
          // -----------------------------------------
          if ( fabs(param.FacetNormalVelocity(i)) > numeric_limits<double>::epsilon() )
            {

              limitProperty_LSMGRAD<dim>( e, mass_center_key, grad_sn_key, grad_sn_limiter_key,
                                          inside_node_,outside_node_, i,
                                          sn_inside_node, sn_outside_node,
                                          limited_sn_inside_node, limited_sn_outside_node);


               if ( param.FacetNormalVelocity(i)> 0. )
                   relperm.SaturationWettingPhase( 1. - limited_sn_inside_node);
               else
                   relperm.SaturationWettingPhase( 1. - limited_sn_outside_node);

               // computing flow properties at upstream node
               relperm.EffectiveSaturation();
               // + div . [fn vt]
               const double nonlinear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i) * relperm.f_Phase(NONWETTING_PHASE);

               // 2nd-order fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += nonlinear_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= nonlinear_flux;
             }

      }

 } // end AccumulateTwoPhaseSolution2 ( generic fully upwind version, but without capillary displacement ) with lsm grad limiter



template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution2(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     TwoPhaseModel<dim>& relperm,
                                     vector<double>& res)
 {
     eidx_ = e.Idx();

     // 0. computing capillary pressure gradient from pc at nodes
     // ---------------------------------------------------------

     /*
     // what we are after is the derivative of this gradient to get an approximation
     // like dt/x2 diffusivity (psi_i-1 - 2 psi_i + psi_i+1), see Leveque, p. 60

     // note: doesn't give the right result
     relperm.Initialize( e );
     e.dN_AtBaryCenter( DN_ );

     fill( pc_grad_.begin(), pc_grad_.end(), 0. );
     for ( size_t j{0U}; j<e.Nodes(); j++ ) {
          relperm.SaturationWettingPhase( 1. - e.N(j)->Read( adv1key_ ) );
          const double pc_at_node = relperm.pc_Phase(NONWETTING_PHASE);
          for ( size_t k{0U}; k<dim; k++ ) pc_grad_[k] += DN_(k,j) * pc_at_node;
     }*/

     relperm.Initialize( e );
     e.dN_AtBaryCenter( DN_ );

     fill( dsdn_.begin(), dsdn_.end(), 0. );
     for ( auto j{0U}; j<e.Nodes(); j++ ) {
          const double sn = e.N(j)->Read( adv1_key_);
          for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
     }

     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          // --------------------------------------------------------
          // 0. interpolate property values and effective saturation
          // --------------------------------------------------------
          relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
          relperm.EffectiveSaturation();
          double sn_facet = relperm.Saturation(NONWETTING_PHASE);
          // get non-wetting-phase saturations for the inside and the outside nodes
          const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );

          // -----------------------------------------
          // 1. accumulation of nonlinear facet fluxes
          // -----------------------------------------
          if ( fabs(param.FacetNormalVelocity(i)) > numeric_limits<double>::epsilon() ) // economising on the computations
            {

               // first-order unconditionally stable scheme if at model boundary
               if ( e.N(inside_node_)->AtBoundary() or e.N(outside_node_)->AtBoundary() ) {
                 }
               // identifying the upstream node (to get psi_hat_c, n_upstream = n_current )
               if ( param.FacetNormalVelocity(i) > 0. ) {
                 if ( e.N(inside_node_)->AtBoundary() ) relperm.SaturationWettingPhase( 1. - sn_inside_node );
                 else relperm.SaturationWettingPhase( 1. - limitProperty( sn_inside_node, sn_outside_node, sn_facet,
                                                                          SMINMAX[ e.N( e.FV()->InsideNode(i) )->Idx() ] ) );
                 }
               else {
                 if ( e.N(outside_node_)->AtBoundary() ) relperm.SaturationWettingPhase( 1. - sn_outside_node );
                 else relperm.SaturationWettingPhase( 1. - limitProperty( sn_outside_node, sn_inside_node, sn_facet,
                                                                          SMINMAX[ e.N( e.FV()->OutsideNode(i) )->Idx() ] ) );
                 }


               // computing flow properties at upstream node
               relperm.EffectiveSaturation();

               // + div . [fn vt]
               const double nonlinear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i) * relperm.f_Phase(NONWETTING_PHASE);

               // 2nd-order fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += nonlinear_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= nonlinear_flux;
             }


          // -----------------------------------------
          // 2. (non-linear) (capillary) diffusion
          // -----------------------------------------
          /* comments:
             - this is not the correct approximation of the second derivative of the saturation
               gradient which is needed by the explicit scheme (Leveque, p. 64, eqn. 4.12)
             - however, it appears to give the right results using this high resolution method?!

          // note: doesn't give the right results!

          const double dpcdn(param.FacetNormalProjection( i, pc_grad_ ));

          if ( fabs(dpcdn) > numeric_limits<double>::epsilon() )
          {
               relperm.SaturationWettingPhase( 1. - sn_facet );
               relperm.EffectiveSaturation();
               // =  - div . [k * lambda_overbar * dpc/dn]
               const double capillary_flux = param.FacetArea(i) * relperm.Permeability() * relperm.G() * -dpcdn;
               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += capillary_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= capillary_flux;
          }
          */

          /// note: doesn't give the right results!
          /// ( find better approximation in function AccumulateExplicitTwoPhaseSolution2 with grav. and cap. bool flags )
          /// (Roman, 2013 )
          /*

          const double ds = (sn_outside_node - sn_inside_node);

          if ( ds != 0.0 ) {

               const double dsdn(param.FacetNormalProjection( i, dsdn_ ));

               // identifying the upstream node
               // flow will be out of a facet if proj. of sat-gradient onto normal is negative
               const double  psi_hat_c( (dsdn < zero) ? sn_inside_node : sn_outside_node );

               // computing flow properties at upstream node
               relperm.SaturationWettingPhase( 1. - psi_hat_c );
               relperm.EffectiveSaturation();

               const double capillary_flux = -dsdn * relperm.CapillaryDiffusionMultiplier(NONWETTING_PHASE) * param.FacetArea(i);

               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += capillary_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= capillary_flux;
            }
            */

      }

 } // end AccumulateTwoPhaseSolution2 ( generic fully upwind version )


template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution2(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     TwoPhaseModel<dim>& relperm,
                                     vector<double>& res,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_sn_key,
                                     const csmp::Index& grad_sn_limiter_key)
 {
     eidx_ = e.Idx();

     double limited_sn_inside_node(1.0), limited_sn_outside_node(1.0);

     // 0. computing capillary pressure gradient from pc at nodes
     // ---------------------------------------------------------

     relperm.Initialize( e );
     e.dN_AtBaryCenter( DN_ );

     fill( dsdn_.begin(), dsdn_.end(), 0. );
     for ( auto j{0U}; j<e.Nodes(); j++ ) {
          const double sn = e.N(j)->Read( adv1_key_);
          for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
     }


     for ( auto i{0U}; i<e.FV()->Facets(); i++ )
       {
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          // --------------------------------------------------------
          // 0. interpolate property values and effective saturation
          // --------------------------------------------------------
          relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
          relperm.EffectiveSaturation();

          // get non-wetting-phase saturations for the inside and the outside nodes
          const double sn_inside_node  = e.N(inside_node_)->Read( adv1_key_ );
          const double sn_outside_node = e.N(outside_node_)->Read( adv1_key_ );

          // -----------------------------------------
          // 1. accumulation of nonlinear facet fluxes
          // -----------------------------------------
          if ( fabs(param.FacetNormalVelocity(i)) > numeric_limits<double>::epsilon() ) // economising on the computations
            {

              limitProperty_LSMGRAD<dim>( e, mass_center_key, grad_sn_key, grad_sn_limiter_key,
                                          inside_node_,outside_node_, i,
                                          sn_inside_node, sn_outside_node,
                                          limited_sn_inside_node, limited_sn_outside_node);

               if ( param.FacetNormalVelocity(i)> 0. )
                   relperm.SaturationWettingPhase( 1. - limited_sn_inside_node);
               else
                   relperm.SaturationWettingPhase( 1. - limited_sn_outside_node);


               // computing flow properties at upstream node
               relperm.EffectiveSaturation();

               // + div . [fn vt]
               const double nonlinear_flux = param.FacetNormalVelocity(i) * param.FacetArea(i) * relperm.f_Phase(NONWETTING_PHASE);

               // 2nd-order fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += nonlinear_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= nonlinear_flux;
             }


          // -----------------------------------------
          // 2. (non-linear) (capillary) diffusion
          // -----------------------------------------

          /// note: doesn't give the right results!
          /// ( find better approximation in function AccumulateExplicitTwoPhaseSolution2 with grav. and cap. bool flags )
          /// (Roman, 2013 )

          /*
          const double ds = (sn_outside_node - sn_inside_node);

          if ( ds != 0.0 ) {

               const double dsdn(param.FacetNormalProjection( i, dsdn_ ));

               // identifying the upstream node
               // flow will be out of a facet if proj. of sat-gradient onto normal is negative
               const double  psi_hat_c( (dsdn < zero) ? sn_inside_node : sn_outside_node );

               // computing flow properties at upstream node
               relperm.SaturationWettingPhase( 1. - psi_hat_c );
               relperm.EffectiveSaturation();

               const double capillary_flux = -dsdn * relperm.CapillaryDiffusionMultiplier(NONWETTING_PHASE) * param.FacetArea(i);

               res[ e.N( e.FV()->InsideNode(i) )->Idx() ]  += capillary_flux;
               res[ e.N( e.FV()->OutsideNode(i) )->Idx() ] -= capillary_flux;
            }

          */
      }

 } // end AccumulateTwoPhaseSolution2 ( generic fully upwind version ) with lsm grad limiter




/// New upwind strategy for capillary and gravitational terms
/// Roman

template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution1(
                                                        const FV_Parameter& param,
                                                        const Element<dim>& e,
                                                        TwoPhaseModel<dim>& relperm,
                                                        vector<double>& res,
                                                        bool with_gravity_forces,
                                                        bool with_capillary_spreading)
 {

    eidx_ = e.Idx();
    bool DebugOutput(false);

    relperm.Initialize( e );

    const size_t v( (dim==1u) ? 0u : 1u );
    const double zero(0.);

    double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
    double upstream_fn(0.0),upstream_lambda_overbar(0.0);
    double upstream_sn(0.0),upstream_sw(0.0);
    double vn_at_facet_int_point(0.0),vw_at_facet_int_point(0.0);

    double nonlinear_flux(0.0);
    double viscous_velocity_component(0.0),gravity_velocity_component(0.0),capillary_velocity_component(0.0);
    double dsdn(0.0),dpcdn(0.0);


    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );

        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
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
                  dpcdn = -dsdn*relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative ove sn is needed

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

          res[ e.N( inside_node_)->Idx()  ] += nonlinear_flux;
          res[ e.N( outside_node_)->Idx() ] -= nonlinear_flux;


          if(DebugOutput){
              cout<<"So="<<relperm.Saturation(NONWETTING_PHASE)<<endl;
              cout<<"Sw="<<relperm.Saturation(WETTING_PHASE)<<endl;
              cout<<"So+Sw="<<relperm.Saturation(WETTING_PHASE)+relperm.Saturation(NONWETTING_PHASE)<<endl;
              cout<<"So_outside="<<sn_outside_node<<endl;
              cout<<"Sn_inside="<<sn_inside_node<<endl;
          }

    }

 } // end AccumulateExplicitTwoPhaseSolution1  ( generic fully upwind version 1st-order with other upwind technique)




template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution2(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     TwoPhaseModel<dim>& relperm,
                                     vector<double>& res,
                                     bool with_gravity_forces,
                                     bool with_capillary_spreading)
 {

    eidx_ = e.Idx();

    bool DebugOutput(false);

    relperm.Initialize( e );

    const size_t v( (dim==1u) ? 0u : 1u );
    const double zero(0.);

    double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
    double upstream_fn(0.0),upstream_lambda_overbar(0.0);
    double upstream_sn(0.0),upstream_sw(0.0);
    double vn_at_facet_int_point(0.0),vw_at_facet_int_point(0.0);

    double nonlinear_flux(0.0);
    double viscous_velocity_component(0.0),gravity_velocity_component(0.0),capillary_velocity_component(0.0);
    double dsdn(0.0),dpcdn(0.0);

    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );

        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
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
                  dpcdn = -dsdn*relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative ove sn is needed

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


          res[ e.N( inside_node_)->Idx()  ] += nonlinear_flux;
          res[ e.N( outside_node_)->Idx() ] -= nonlinear_flux;


          if(DebugOutput){
              cout<<"So="<<relperm.Saturation(NONWETTING_PHASE)<<endl;
              cout<<"Sw="<<relperm.Saturation(WETTING_PHASE)<<endl;
              cout<<"So+Sw="<<relperm.Saturation(WETTING_PHASE)+relperm.Saturation(NONWETTING_PHASE)<<endl;
              cout<<"So_outside="<<sn_outside_node<<endl;
              cout<<"Sn_inside="<<sn_inside_node<<endl;
          }
    }

} // end AccumulateTwoPhaseSolution2

template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolution2(
                                     const vector<pair<double,double> >& SMINMAX,
                                     const FV_Parameter& param,
                                     const Element<dim>& e,
                                     TwoPhaseModel<dim>& relperm,
                                     vector<double>& res,
                                     bool with_gravity_forces,
                                     bool with_capillary_spreading,
                                     const csmp::Index& mass_center_key,
                                     const csmp::Index& grad_sn_key,
                                     const csmp::Index& grad_sn_limiter_key)
 {

    eidx_ = e.Idx();

    bool DebugOutput(false);

    relperm.Initialize( e );

    const size_t v( (dim==1u) ? 0u : 1u );
    const double zero(0.);

    double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
    double upstream_fn(0.0),upstream_lambda_overbar(0.0);
    double upstream_sn(0.0),upstream_sw(0.0);
    double vn_at_facet_int_point(0.0),vw_at_facet_int_point(0.0);
    double limited_sn_inside_node(0.0),limited_sn_outside_node(0.0);

    double nonlinear_flux(0.0);
    double viscous_velocity_component(0.0),gravity_velocity_component(0.0),capillary_velocity_component(0.0);
    double dsdn(0.0),dpcdn(0.0);

    if(with_capillary_spreading){

        e.dN_AtBaryCenter( DN_ );

        fill( dsdn_.begin(), dsdn_.end(), 0. );
        for ( auto j{0U}; j<e.Nodes(); j++ ) {
             const double sn = e.N(j)->Read( adv1_key_);
             for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
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


          relperm.InitializeForFacetIntegrationPoint( i, 0U, e );
          relperm.EffectiveSaturation();

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
                  dpcdn = -dsdn*relperm.dpcds_Phase( ); // minus, because one calculate derivative of pc over sw, when the derivative ove sn is needed

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


          res[ e.N( inside_node_)->Idx()  ] += nonlinear_flux;
          res[ e.N( outside_node_)->Idx() ] -= nonlinear_flux;


          if(DebugOutput){
              cout<<"So="<<relperm.Saturation(NONWETTING_PHASE)<<endl;
              cout<<"Sw="<<relperm.Saturation(WETTING_PHASE)<<endl;
              cout<<"So+Sw="<<relperm.Saturation(WETTING_PHASE)+relperm.Saturation(NONWETTING_PHASE)<<endl;
              cout<<"So_outside="<<sn_outside_node<<endl;
              cout<<"Sn_inside="<<sn_inside_node<<endl;
          }
    }

} // end AccumulateTwoPhaseSolution2 with lsm grad limiter


template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::AccumulateExplicitTwoPhaseSolutionAtBoundary(  const FV_Parameter& param,
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

    const uint32_t v = (dim==1u) ? 0u : 1u;
    VectorVariable<dim> velo;
    double facetArea;

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

             if(with_capillary_spreading){

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

 } // end AccumulateExplicitTwoPhaseSolutionAtBoundary_WithGravity


/** print stencil state to screen
 
Outputs the current internal variables of the FV sector to screen.  

*/
template<uint32_t dim>
void  ExplicitStencilProcessor<dim>::Out() const
 {
    cout <<"\n\n\nExplicitStencilProcessor<"<< dim <<">::Out: "<< endl;
    // outputting data of object
    cout <<"Finite volume stencil data."<< endl;
    cout <<"\n\tassociated sectors, pore volumes, and advected property values: ";
    for ( auto i{0U}; i<sector_pore_volume_.size(); i++ )
      cout <<"\n\tsector "<< i+1 <<": "<< sector_pore_volume_[i] <<", "<< psi1_[i];
          
    if ( !src_.empty() ) {
         cout <<"\n\nFacet integrated source terms:";
         for ( size_t i{0U}; i<src_.size(); i++ )
         cout <<"\n\tfacet "<< i+1 <<": "<< src_[i];
      }
        
    cout << endl << endl;
    
 } // end Out


template struct ExplicitStencilProcessor<1U>;
template struct ExplicitStencilProcessor<2U>;
template struct ExplicitStencilProcessor<3U>;


} // namespace csmp



















