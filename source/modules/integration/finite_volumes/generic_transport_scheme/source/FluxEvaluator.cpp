//
//  FluxEvaluator.cpp
//
//  Created by Stephan Matthai on 2/21/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "FluxEvaluator.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "TwoPhaseModel.h"
#include "ExplicitTransport.h"
#include "ImplicitTransport.h"

using namespace std;

namespace csmp {


/** 
    Computes A_i vD . n_i for all facets and its product with the upstream concentrations.
    Stores 'facet flux'  and 'facet flux concentration' product variables.
    Accumulates the flux balance and flux concentration balances in the corresponding
    node variables.
    
    @attention vD can be used directly as interstitial velocity since the FV is scaled by porosity
 
    @attention the thickness of lower-dimensional elements is taken into account
 
    @attention, assumes that the FV variables 'flux balance' and 'new concentration' were zeroed
    out before.
 
    @tested OK SKM 17/8/19 
*/
template<size_t dim, template<size_t> class USER>
void FluxEvaluator<dim,USER>::Advective_O1_FluxesAndBalances( Element<dim>* const eptr ) const
 {
   assert( eptr != nullptr );

   // element-based Darcy velocity
   eptr->Read( User()->key_V, vD_ );
   const double64 thickness = eptr->Read( User()->key_THI );

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

        // accumulating the flux balances into FV's associated with the element's
        const size_t inside_node  = eptr->FV()->InsideNode( j );
        const size_t outside_node = eptr->FV()->OutsideNode( j );
     
         // facet edge node pair by node pair
        eptr->N(inside_node)->Read( User()->key_FB, sc_ );
        eptr->N(inside_node)->Store( User()->key_FB, (sc_ += facet_flux) );
        eptr->N(outside_node)->Read( User()->key_FB, sc_ );
        eptr->N(outside_node)->Store( User()->key_FB, (sc_ -= facet_flux) );

        // multiplying the volumetric facet flux with the upstream concentration
         // fluxes get multiplied with upstream concentrations
        if ( facet_flux < 0. ) facet_flux *= eptr->N(outside_node)->Read( User()->key_C0 );
        else                   facet_flux *= eptr->N(inside_node)->Read( User()->key_C0 );

        // storing facet flux concentration product without altering the variables flag
        const VARIABLE_FLAG flag2 = eptr->Status( j, 0U, User()->key_ffC );
        eptr->Store( j, 0U, User()->key_ffC, makeScalar(flag2,facet_flux) );
     
        // accumulating the volume-flux concentration product balances into FV's associated with the element's
        // facet edge node pair by node pair
        eptr->N(inside_node)->Read( User()->key_acc_, sc_ );
        eptr->N(inside_node)->Store( User()->key_acc_, (sc_ += facet_flux) );
        eptr->N(outside_node)->Read( User()->key_acc_, sc_ );
        eptr->N(outside_node)->Store( User()->key_acc_, (sc_ -= facet_flux) );
     }
   
 } // end Advective_O1_FluxesAndBalances



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
     COMPUTING FLUX BALANCE FROM VELOCITY AND NORMAL AND FACET AREA.
 
     The divergence of the transport velocity on the current FV = node,
     and its product with the transport variable is stored.
*/
template<size_t dim, template<size_t> class USER>
void FluxEvaluator<dim,USER>::Advective_O1_FluxesAndBalances( Node<dim>* const nd_ptr ) const
  {
    assert( nd_ptr != nullptr );
  
    double64 flux_balance(0.), accumulation(0.);
  
    const size_t node_parent_elements(nd_ptr->Parents());
    for ( size_t t=0U; t<node_parent_elements; ++t )
      {
        Element<dim>* const eptr(nd_ptr->Parent(t));
        assert( eptr != NULL );
        const size_t pnid(nd_ptr->ParentNodeNumber(t));
        eptr->Read( User()->key_V, vD_ );
        const double64 thickness = eptr->Read( User()->key_THI );

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
             // computing and storing the facet fluxes and their concentration products
             const double64 sign = ( pnid == inside_node ) ? 1. : -1.;
             const double64 flux = sign * vD_n * facetArea;
             eptr->Store( iFacet, 0U, User()->key_ff, makeScalar(eptr->Status(iFacet, 0U, User()->key_ff), flux) );
             eptr->Store( iFacet, 0U, User()->key_ffC, makeScalar(eptr->Status(iFacet, 0U, User()->key_ffC), flux * C_upstream) );
             flux_balance += flux;
             accumulation += flux * C_upstream;
          }
      } // end for loop for parent elements

    // storing the balances of the FVs = nodes
    nd_ptr->Store( User()->key_FB, makeScalar( nd_ptr->Status( User()->key_FB ), flux_balance ) );
    nd_ptr->Store( User()->key_acc_, makeScalar( ANY, accumulation ) );
  
} // end Advective_O1_FluxesAndBalances





/**
    From precomputed facet fluxes, AdvectiveFluxesAtBoundary() computes:
 
     0. Nothing for FVs on the model boundary with Dirichlet concentration constraints.
    
     1. facet fluxes and flux balances of intact perimeter finite volumes that are not on the model boundary

     2. accumulation of facet flux-concentration products
 
     3. inflow if the the perimeter finite volume is truncated by the model boundary.
        This is determined from the AtBoundary() flag (flag!=NOT).
        
    In summary, this method computes flux balances and concentration-facet flux products where possible,
    at sliced boundaries inflow (+) concentration products are stored to C1.
    At such boundaries the influx is found from the flux balance and the FV cell's concentration.

    @attention This method is only for perimeter nodes
 
*/
template<size_t dim, template<size_t> class USER>
void FluxEvaluator<dim,USER>::Advective_O1_FluxBalancesAtBoundary( Node<dim>* const nd_ptr ) const
  {
     assert( nd_ptr != NULL );
  
     // INTACT FINITE VOLUMES
  
     // 1. computation of flux balances
     // -------------------------------
     // The divergence of vD is stored and returned
     // the transport variable flux balance is computed from current facet fluxes and stored in 'accumulation'
     // TODO: probably no need to recompute this
//     FluxBalancesFromFacetFluxes( nd_ptr );
  
  
     // 2. where the boundary FVs are intact influxes and outfluxes are balanced
     //    and no special treatment is needed
     // ------------------------------------------------------------------------
     // any flux balance here indeed gets counter balanced in Assemble_O1_Solution()
     if ( nd_ptr->AtBoundary() == NOT ) return;
  
  
     // TRUNCATED FINITE VOLUMES AT THE MODEL BOUNDARY

     // 3. for truncated FVs at no-flow boundaries, in- or outflux can be obtained,
     //    therefore no special treatment is needed.
     // ----------------------------------------------------------------------------------------------------
     //   (if there is a nodal fluid volume source causing a flux divergence it is accounted for)
     const VARIABLE_FLAG  pf_status(nd_ptr->Status( User()->key_PF ));
     assert( nd_ptr->AtBoundary() != NOT );
     assert( pf_status != NEUMANN );
     if (  pf_status != DIRICH  ) return;
  
  
     // 4. for truncated FV at inflow boundary, the flux balance must be corrected using their transport variable value
     // ---------------------------------------------------------------------------------------------------------------
     // truncated FV, but no Dirichlet constraint
     const double64 in_flow = -1. * this->InFlow( nd_ptr );
     if (  in_flow < 0. ) {
          double64 accumulation = nd_ptr->Read( User()->key_acc_ );
          // correction:  incoming transport variable     inflow
          accumulation += nd_ptr->Read(User()->key_C0) * -in_flow;
          // since the flux balance has been dealt with, it is set to zero
          nd_ptr->Store( User()->key_FB, makeScalar( nd_ptr->Status( User()->key_FB ), 0. ) );
          nd_ptr->Store( User()->key_acc_, makeScalar( ANY, accumulation ) );
          return;
      }


     // 3. for truncated FV's at outflow boundaries a volume correction is needed
     // -------------------------------------------------------------------------
     // concentration value is copied to next time level and flux balance set to zero
     const double64 out_flow = this->OutFlow( nd_ptr );
     if ( out_flow > 0. )  {
          double64 accumulation = nd_ptr->Read( User()->key_acc_ );
          // correction:  incoming transport variable     outflow
          accumulation += nd_ptr->Read(User()->key_C0) * -out_flow;
          // since the flux balance has been dealt with, it is set to zero
          nd_ptr->Store( User()->key_FB, makeScalar( nd_ptr->Status( User()->key_FB ), 0. ) );
          nd_ptr->Store( User()->key_acc_, makeScalar( ANY, nd_ptr->Read( User()->key_C0 ) ) );
       }
  
 } // end Advective_O1_FluxBalancesAtBoundary








/**
    Computes first-order facet fluxes, F = ff * C0 = (pre-computed) volumetric flow upstream C0 products on all facets
 
    @apply this method to the interior and halo stencils.
*/
template<size_t dim, template<size_t> class USER>
void FluxEvaluator<dim,USER>::TransportVariableFluxes( Element<dim>* const eptr ) const
 {
   assert( eptr != nullptr );

   // computing total facet fluxes by projecting vt onto facet normals
   const size_t facets(eptr->FV()->Facets());
   for ( size_t j=0U; j<facets; ++j ) {
         // reading the facet flux
        double64 facet_flux = eptr->Read( j, 0U, User()->key_ff );
        // establishing the upstream direction
        const size_t inside_node  = eptr->FV()->InsideNode( j );
        const size_t outside_node = eptr->FV()->OutsideNode( j );
        if ( facet_flux < 0. ) facet_flux *= eptr->N(outside_node)->Read( User()->key_C0 );
        else                   facet_flux *= eptr->N(inside_node)->Read( User()->key_C0 );

        // storing facet flux concentration product without altering the variables flag
        const VARIABLE_FLAG flag2 = eptr->Status( j, 0U, User()->key_ffC );
        eptr->Store( j, 0U, User()->key_ffC, makeScalar(flag2,facet_flux) );
     }
 
 } // end TransportVariableFluxes


template<size_t dim, template<size_t> class USER>
void FluxEvaluator<dim,USER>::TransportVariableFluxesAndBalances( Element<dim>* const eptr ) const
 {
   assert( eptr != nullptr );

   // computing total facet fluxes by projecting vt onto facet normals
   const size_t facets(eptr->FV()->Facets());
   for ( size_t j=0U; j<facets; ++j ) {
         // reading the facet flux
        double64 facet_flux = eptr->Read( j, 0U, User()->key_ff );
        // establishing the upstream direction
        const size_t inside_node  = eptr->FV()->InsideNode( j );
        const size_t outside_node = eptr->FV()->OutsideNode( j );
        if ( facet_flux < 0. ) facet_flux *= eptr->N(outside_node)->Read( User()->key_C0 );
        else                   facet_flux *= eptr->N(inside_node)->Read( User()->key_C0 );

        // storing facet flux concentration product without altering the variables flag
        const VARIABLE_FLAG flag2 = eptr->Status( j, 0U, User()->key_ffC );
        eptr->Store( j, 0U, User()->key_ffC, makeScalar(flag2,facet_flux) );
     
        // accumulating the volume-flux concentration product balances into FV's associated with the element's
        // facet edge node pair by node pair
        eptr->N(inside_node)->Read( User()->key_acc_, sc_ );
        eptr->N(inside_node)->Store( User()->key_acc_, (sc_ += facet_flux) );
        eptr->N(outside_node)->Read( User()->key_acc_, sc_ );
        eptr->N(outside_node)->Store( User()->key_acc_, (sc_ -= facet_flux) );
     }
 
 } // end TransportVariableFluxesAndBalances




/**
    Computes and returns FV flux balance computed from pre-calculated thickness-weighted facet fluxes.
    The divergence of the volumetric flow is stored into the variable 'flux balance'.
    The volumetric flow * upstream concentration product is stored into  'new concentration'.
    Both, the flux balance and its product with concentration are stored on the node.
 
    @attention this method relies on precomputed facet fluxes.
*/
// TODO: turn off the volumetric flow balance computation because it already exists
template<size_t dim, template<size_t> class USER>
void FluxEvaluator<dim,USER>::FluxBalancesFromFacetFluxes( Node<dim>* const nptr ) const
 {
     double64 flux_C0_products(0.);
 
     const size_t parent_elements(nptr->Parents());
     for ( size_t i=0U; i<parent_elements; ++i ) {
           const Element<3U>* const eptr = nptr->Parent(i);
           const size_t sector_node      = nptr->ParentNodeNumber(i);
           const size_t sector_facets(eptr->FV()->FacetsPerSector(sector_node));
           for ( size_t j=0U; j<sector_facets; ++j ) {
               const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
               const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
               // facet flux already is the product with element thickness
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
               // temporary accumulation of flux-concentration products into the variable 'accumulation'
               flux_C0_products += sign * facet_flux * eptr->Read( facet, 0U, User()->key_ffC );
            }
       }
 
     // variable 'accumulation' which gets subsequently overwritten
     nptr->Store( User()->key_acc_, makeScalar(ANY,flux_C0_products) );

 } // end FluxBalancesFromFacetFluxes


/* volumetric flow- and variable flux balances

template<size_t dim, template<size_t> class USER>
double64 FluxEvaluator<dim,USER>::FluxBalancesFromFacetFluxes( Node<dim>* const nptr ) const
 {
     double64 flux_balance(0.), flux_C0_products(0.);
 
     const size_t parent_elements(nptr->Parents());
     for ( size_t i=0U; i<parent_elements; ++i ) {
           const Element<3U>* const eptr = nptr->Parent(i);
           const size_t sector_node      = nptr->ParentNodeNumber(i);
           const size_t sector_facets(eptr->FV()->FacetsPerSector(sector_node));
           for ( size_t j=0U; j<sector_facets; ++j ) {
               const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
               const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
               // facet flux already is the product with element thickness
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
               flux_balance += facet_flux;
               // temporary accumulation of flux-concentration products into the variable 'accumulation'
               flux_C0_products += sign * facet_flux * eptr->Read( facet, 0U, User()->key_ffC );
            }
       }
 
     // variable 'flux balance'
     nptr->Store( User()->key_FB, makeScalar(nptr->Status(User()->key_FB),flux_balance) );
     // variable 'accumulation' which gets subsequently overwritten
     nptr->Store( User()->key_acc_, makeScalar(ANY,flux_C0_products) );

     return flux_balance;

 } // end FluxBalancesFromFacetFluxes
*/




/**
    Like FluxBalancesFromFacetFluxes(), but returns the current outflow from the cell.
 
    The latter is valid only for intact finite volumes (not at model boundary).
    Else, a ficticious outflow is detected at inflow boundaries.

    @attention this method relies on precomputed facet fluxes.
*/
template<size_t dim, template<size_t> class USER>
double64 FluxEvaluator<dim,USER>::FluxBalanceAndOutFlow( Node<dim>* const nptr ) const
 {
     double64 flux_balance(0.), out_flow(0.), flux_C0_products(0.);
 
     const size_t parent_elements(nptr->Parents());
     for ( size_t i=0U; i<parent_elements; ++i ) {
          const Element<3U>* const eptr = nptr->Parent(i);
          const size_t sector_node      = nptr->ParentNodeNumber(i);
           for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
               const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
               const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
               // facet flux already is the product with element thickness
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
               flux_balance     += facet_flux;
               // flow out the FV which is positive
               if ( facet_flux > 0. ) out_flow += facet_flux;
               // temporary accumulation of flux-concentration products into the variable 'new concentration'
               flux_C0_products += sign * facet_flux * eptr->Read( facet, 0U, User()->key_ffC );
            }
       }
 
     // variable 'flux balance'
     nptr->Store( User()->key_FB, makeScalar(nptr->Status(User()->key_FB),flux_balance) );
     // variable 'accumulation' which gets subsequently overwritten
     nptr->Store( User()->key_acc_, makeScalar(ANY,flux_balance) );

     return out_flow;

 } // end FluxBalanceAndOutFlow



/**
    Like FluxBalanceAndOutFlow(), but computes only the volumetric inflow into the cell.
 
    The latter is valid only for intact finite volumes (not at model boundary).
    Else, a ficticious inflow is detected at inflow boundaries.

    @attention this method relies on precomputed facet fluxes.
*/
template<size_t dim, template<size_t> class USER>
double64 FluxEvaluator<dim,USER>::InFlow( const Node<dim>* const nptr ) const
 {
     double64 in_flow(0.);
 
     const size_t parent_elements(nptr->Parents());
     for ( size_t i=0U; i<parent_elements; ++i ) {
          const Element<3U>* const eptr = nptr->Parent(i);
          const size_t sector_node      = nptr->ParentNodeNumber(i);
           for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
               const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
               const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
               // facet flux already is the product with element thickness
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
               // flow out the FV which is positive
               if ( facet_flux < 0. ) in_flow += facet_flux;
            }
       }

     return fabs(in_flow);

 } // end InFlow


/**
    Like FluxBalanceAndOutFlow(), but computes only the volumetric outflow from the cell.
 
    The latter is valid only for intact finite volumes (not at model boundary).
    Else, not the complete outflow may be detected at out-flow boundaries.

    @attention this method relies on precomputed facet fluxes.
*/
template<size_t dim, template<size_t> class USER>
double64 FluxEvaluator<dim,USER>::OutFlow( const Node<dim>* const nptr ) const
 {
     double64 out_flow(0.);
 
     const size_t parent_elements(nptr->Parents());
     for ( size_t i=0U; i<parent_elements; ++i ) {
          const Element<3U>* const eptr = nptr->Parent(i);
          const size_t sector_node      = nptr->ParentNodeNumber(i);
           for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
               const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
               const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
               // facet flux already is the product with element thickness
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
               // flow out the FV which is positive
               if ( facet_flux > 0. ) out_flow += facet_flux;
            }
       }

     return out_flow;

 } // end OutFlow

  
template class FluxEvaluator<3U,ExplicitTransport>;
template class FluxEvaluator<3U,ImplicitTransport>;

} // end csmp
