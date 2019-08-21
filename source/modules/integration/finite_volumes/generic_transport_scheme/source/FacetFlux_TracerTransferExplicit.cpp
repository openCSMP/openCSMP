//
//  FacetFlux_TracerTransferExplicit.cpp
//
//  Created by Stephan Matthai on 2/21/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "FacetFlux_TracerTransferExplicit.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "TwoPhaseModel.h"
#include "ExplicitTransport.h"

using namespace std;

namespace csmp {


/** 
    Computes A_i vD . n_i for all facets and its product with the upstream concentrations
    of element stencil, storing it there.
    Stores 'facet flux'  and 'facet flux concentration' product variables.
    
    @attention that vD can be used directly as interstitial velocity since the FV is scaled by porosity
 
    @attention the thickness of lower-dimensional elements is taken into account
 
    @tested OK SKM 17/8/19 
*/
template<size_t dim, template<size_t> class USER>
void FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O1_FluxesInterior( Element<dim>* eptr ) const
 {
   assert( eptr != nullptr );

   // element-based Darcy velocity
   eptr->Read( User()->key_V, vD_ );
   const double64 thickness = eptr->Read( User()->key_thi );

   // computing total facet fluxes by projecting vt onto facet normals
   const size_t facets(eptr->FV()->Facets());
   for ( size_t j=0U; j<facets; ++j ) {
        eptr->Read( j, 0U, User()->key_fn, nrml_ );
        // projection (dot product)
        double64 facet_flux = nrml_[0] * vD_[0];
        if ( dim != 1U ) facet_flux += nrml_[1] * vD_[1];
        if ( dim >  2U ) facet_flux += nrml_[2] * vD_[2];
        // multiplication with element thickness
        facet_flux *= thickness;
        // multiplication with facet area
        facet_flux *= eptr->Read( j, 0U, User()->key_fA );

        // storing the volumetric facet flux without altering the variables flag
        const VARIABLE_FLAG flag1 = eptr->Status( j, 0U, User()->key_ff );
        eptr->Store( j, 0U, User()->key_ff, makeScalar(flag1,facet_flux) );

        // multiplying the volumetric flux with the upstream concentration
        const size_t inside_node  = eptr->FV()->InsideNode( j );
        const size_t outside_node = eptr->FV()->OutsideNode( j );
        // fluxes are multiplied with upstream concentrations
        if ( facet_flux < 0. ) facet_flux *= eptr->N(outside_node)->Read( User()->key_C0 );
        else                   facet_flux *= eptr->N(inside_node)->Read( User()->key_C0 );

        // storing facet flux concentration product without altering the variables flag
        const VARIABLE_FLAG flag2 = eptr->Status( j, 0U, User()->key_ffC );
        eptr->Store( j, 0U, User()->key_ffC, makeScalar(flag2,facet_flux) );
     }
   
 } // end AdvectiveFluxes



/* TESTING CODE
cerr <<"\n"<< eptr->Idx() <<":"<< j <<": facet flux: "<< facet_flux;
     
     
bool error(false);

// SKM TEST
if ( fabs(1. - nrml_.Length()) > numeric_limits<double64>::epsilon() ) {
     cerr.precision(15);
     cerr <<"\nfacet "<< j <<", normal length: "<< std::scientific << nrml_.Length() <<"\n";
     nrml_.Out();
     Point<dim> pnrml = eptr->FacetNormal(j);
     cerr <<"\nnewly computed normal:";
     pnrml.Out();
     error = true;
  }

if (error ) {
    eptr->Out();
    cerr <<"\nvolume: "<< eptr->Volume();
    cerr << "\n";
 }
*/








/**
    AdvectiveFluxesAtBoundary() computes:
    
     1. flux balance of perimeter finite volume if it exists (not the case at model boundary)

     2. accumulation of facet flux-concentration products
 
     3. inflow if the the perimeter finite volume is truncated by the model boundary.
        This is determined from the AtBoundary() flag (flag!=NOT).
        
    In summary, this method computes flux balances and concentration-facet flux products where possible,
    at sliced boundaries inflow (+) concentration products are stored to C1.
    In this case the influx is found from the flux balance and FV cell's concentration.

    @attention This method is only for perimeter nodes
    
    TODO: USE NO-FLOW BOUNDARY CONDITION TO GET EXACT FLUXES AND FLUX-BALANCES through facets at such boundaries
*/
template<size_t dim, template<size_t> class USER>
double64 FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O1_FluxesAtBoundary( Node<dim>* nd_ptr ) const
  {
     assert( nd_ptr != NULL );
    
     // 1. Dirichlet FV
     // ---------------
     // nothing needs to be done for FVs the saturation of which is flagged as Dirichlet
     if ( nd_ptr->Status( User()->key_C0 ) == DIRICH ) {
          // the new value is initialised to the old one and the flag is kept
          nd_ptr->Store( User()->key_C1, makeScalar( nd_ptr->Status( User()->key_C1 ), nd_ptr->Read( User()->key_C0 ) ) );
// TODO: is this the correct treatment of the flux balance
          // the flux balance is set to zero
          nd_ptr->Store( User()->key_FB, makeScalar( nd_ptr->Status( User()->key_FB ), 0. ) );
          // and returned
          return 0.;
       }
    
     // 2. a full finite volume is available so that influxes and outfluxes can be balanced
     // -----------------------------------------------------------------------------------
     if ( nd_ptr->AtBoundary() == NOT ) {
          double64 flux_balance(0.), accumulation(0.);
          const size_t node_parent_elements(nd_ptr->Parents());
          for ( size_t t=0U; t<node_parent_elements; t++ )
            {
              Element<dim>* const eptr(nd_ptr->Parent(t));
              assert( eptr != NULL );
              const size_t pnid(nd_ptr->ParentNodeNumber(t));
              eptr->Read( User()->key_V, vD_ );
              const double64 thickness = eptr->Read( User()->key_thi );

              const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
              for ( size_t i=0U; i<sector_facets; i++ )
                {
                   const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                   const size_t inside_node(eptr->FV()->InsideNode(iFacet));
                   const size_t outside_node(eptr->FV()->OutsideNode(iFacet));

                   eptr->Read( iFacet, 0U, User()->key_fn, nrml_ );
                   const double64  vD_n = vD_.DotProduct(nrml_) * thickness;
                   const double64  facetArea = eptr->Read( iFacet, 0U, User()->key_fA );
                  
                   // finding the upstream concentration
                   const double64 C_upstream = (vD_n < 0.) ? eptr->N(outside_node)->Read( User()->key_C0 ) :
                                                             eptr->N(inside_node)->Read( User()->key_C0 );
                  
                   const double64 sign = ( pnid == inside_node ) ? 1. : -1.;
                   // storing the facet flux
                   eptr->Store( iFacet, 0U, User()->key_ff, makeScalar(eptr->Status(iFacet, 0U, User()->key_ff), sign * vD_n * facetArea * C_upstream) );
                   flux_balance += sign * vD_n * facetArea;
                   accumulation += sign * vD_n * facetArea * C_upstream;
                }
            } // end for loop for parent elements
     
          nd_ptr->Store( User()->key_C1, makeScalar( nd_ptr->Status( User()->key_C1 ), accumulation ) );
          nd_ptr->Store( User()->key_FB, makeScalar( nd_ptr->Status( User()->key_FB ), flux_balance ) );
          return flux_balance;
       }
    
     // 3. The FV is at the model boundary so that only an in- or outflux can be obtained
     // ----------------------------------------------------------------------------------
     // inflow and outflow are measured using the stencils in the interior of the computational region
     // thus inflows to model originate as positive and outflows as negative
     //
     // 3 cases for AtBoundary != NOT:
     //   3.1: prescribed C value at inflow boundary
     //   3.2: prescribed flux (which has an effect only at inflow boundary)
     //   3.3: free outflow (where flux balance missing the outflow facets is negative)
     //
     double64  inflow(0.); // (+) at an inflow boundary and negative at an outflow one
     double64  influx(0.); // the inflow upstream concentration product

     // for all FV SECTORS of FE_FV-stencils of this boundary finite volume
     const size_t node_parent_elements(nd_ptr->Parents());
     for ( size_t t=0U; t<node_parent_elements; t++ )
       {
          const Element<dim>* const eptr(nd_ptr->Parent(t));
          assert( eptr != NULL );
          const size_t pnid(nd_ptr->ParentNodeNumber(t));
          eptr->Read( User()->key_V, vD_ );
          const double64 thickness = eptr->Read( User()->key_thi );

          // for all FACETS per SECTOR surrounding the finite volume at the boundary
          // getting the volumetric fluxes only (upstream concentrations are found later)
          const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
          for ( size_t i=0U; i<sector_facets; i++ )
            {
               const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
               const size_t inside_node(eptr->FV()->InsideNode(iFacet));
               const size_t outside_node(eptr->FV()->OutsideNode(iFacet));

               eptr->Read( iFacet, 0U, User()->key_fn, nrml_ );
               const double64  normal_vel_component = vD_.DotProduct(nrml_) * thickness;
               const double64  facetArea = eptr->Read( iFacet, 0U, User()->key_fA );
              
               const double64 C_upstream = (normal_vel_component < 0.) ? eptr->N(outside_node)->Read( User()->key_C0 ) :
                                                                         eptr->N(inside_node)->Read( User()->key_C0 );
               if ( pnid == inside_node ) {
                    inflow += normal_vel_component * facetArea;
                    influx += normal_vel_component * facetArea * C_upstream;
                 }
               else {
                    inflow -= normal_vel_component * facetArea;
                    influx -= normal_vel_component * facetArea * C_upstream;
                 }
            }
       } // end for loop for parent elements

      // the volume flux balance is stored (source terms are not subtracted if such were applied)
      nd_ptr->Store( User()->key_FB, makeScalar( nd_ptr->Status( User()->key_FB ), inflow ) );

      // computing the flux - concentration products for the sliced FV
      // -------------------------------------------------------------
      // inflow: cases 3.1 and 3.2
      // -------------------------
      if (inflow < 0.) {
           // if there is a nodal fluid volume source causing inflow, it is assumed that a concentration value wat set for the FV
           if ( fabs(nd_ptr->Read(User()->key_NQV)) > numeric_limits<double64>::epsilon() ) {
               VARIABLE_FLAG influx_FV_concentration_status(nd_ptr->Status(User()->key_C0));
               assert( influx_FV_concentration_status == DIRICH or influx_FV_concentration_status == CONSTANT_FLUX );
               // source terms are accumulated later
             }
           // the influx concentration product is assigned
           nd_ptr->Store( User()->key_C1, makeScalar( nd_ptr->Status( User()->key_C1 ), inflow * nd_ptr->Read( User()->key_C0 ) ) );
        }
      // case 3.3: free outflow
      // ----------------------
      // (the upstream values of the concentration from within the model domain are used)
      else nd_ptr->Store( User()->key_C1, makeScalar( nd_ptr->Status( User()->key_C1 ), influx ) );
   
      // the influx is returned
      return inflow; // positive when outgoing
   
 } // end AdvectiveFluxesAtBoundary






/**
    stores and returns FV flux balance computed from current thickness-weighted facet fluxes.
*/
template<size_t dim, template<size_t> class USER>
double64 FacetFlux_TracerTransferExplicit<dim,USER>::FluxBalance( Node<dim>* const nptr ) const
 {
     const size_t parent_elements(nptr->Parents());
     double64 flux_balance(0.);
     for ( size_t i=0U; i<parent_elements; ++i ) {
          const Element<3U>* const eptr = nptr->Parent(i);
          const size_t sector_node      = nptr->ParentNodeNumber(i);
          for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
               const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
               const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
               flux_balance += facet_flux;
            }
       }
     nptr->Store( User()->key_FB, makeScalar(nptr->Status(User()->key_FB),flux_balance) );
   
     return flux_balance;

 } // end FluxBalance

  
  
  




template class FacetFlux_TracerTransferExplicit<3U,ExplicitTransport>;

} // end csmp
