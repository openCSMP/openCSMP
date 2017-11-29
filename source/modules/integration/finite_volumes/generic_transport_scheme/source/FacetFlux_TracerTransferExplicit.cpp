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
#include "finiteVolumeFunctions.h"
#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp {



/** 
    Computes A_i vD . n_i for all facets and its product with the upstream concentrations
    of element stencil, storing it there.
    
    @attention that vD can be used directly as interstitial velocity since the FV is scaled by porosity
 
    @tested 
 
*/
template<size_t dim, template<size_t> class USER>
void FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O1_FluxesInterior( bool reuse_previous_velocity, Element<dim>* eptr ) const
 {
   assert( eptr != NULL );

   FiniteVolumeHelper<dim> helper(eptr);
   Point<dim> vD;
   
    if (!reuse_previous_velocity) {
      // Compute Darcy velocity
      const double64 K(eptr->Read(User()->k_key));
      const double64 mu = User()->GetModel().Read(User()->mu_key);

      vD = -K/mu * helper.GradientOfScalarNodeProperty(eptr->FV()->Barycenter(), User()->pf_key);
      VectorVariable<dim> vD_var(vD);
      eptr->Store(User()->propdb_.StorageKey("velocity"), vD_var);
    }

   // computing total facet fluxes by projecting vt onto facet normals
   const size_t iNrOfFacets(eptr->FV()->Facets());
   for ( size_t iFacet=0U; iFacet<iNrOfFacets; ++iFacet ) {
     double64 facet_flux = 0;
     if (reuse_previous_velocity) {
       facet_flux = eptr->Read( iFacet, 0U, User()->ff_key );
     }
     else {
        // Calulate local normal
        Point<dim> facetNormal(helper.NormalOfFacet(iFacet));

        // Projection. Note that facetNormal is scaled by the area.
        facet_flux = dotProduct(facetNormal, vD);

        // storing the volumetric facet flux without altering the variables flag
        const auto ff_flag = eptr->Status( iFacet, 0U, User()->ff_key );
        eptr->Store( iFacet, 0U, User()->ff_key, makeScalar(ff_flag,facet_flux) );
     }
     
     // Get concentration from upwind node
     auto upwind_node = (facet_flux < 0) ? eptr->FV()->OutsideNode(iFacet) :  eptr->FV()->InsideNode(iFacet);
     const double64 c = eptr->N(upwind_node)->Read( User()->C0_key );

     // Store facet flux concentration
     const auto ffc_flag = eptr->Status(iFacet, 0u, User()->ffC_key);
     eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
  }
   
 } // end Advective_O1_FluxesInterior

/** 
    Computes A_i vD . n_i for all facets and its product with the upstream concentrations
    of element stencil, storing it there.
    
    @attention that vD can be used directly as interstitial velocity since the FV is scaled by porosity
 
    @tested 
 
*/
template<size_t dim, template<size_t> class USER>
void FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O2_FluxesInterior( bool reuse_previous_velocity, Element<dim>* eptr ) const
 {
   assert( eptr != NULL );

   FiniteVolumeHelper<dim> helper(eptr);
   Point<dim> vD;
   auto c0_key = User()->C0_key;
   
    if (!reuse_previous_velocity) {
      // Compute Darcy velocity
      const double64 K(eptr->Read(User()->k_key));
      const double64 mu = User()->GetModel().Read(User()->mu_key);

      vD = -K/mu * helper.GradientOfScalarNodeProperty(eptr->FV()->Barycenter(), User()->pf_key);
      VectorVariable<dim> vD_var(vD);
      eptr->Store(User()->propdb_.StorageKey("velocity"), vD_var);
    }

   // computing total facet fluxes by projecting vt onto facet normals
   const size_t iNrOfFacets(eptr->FV()->Facets());
   for ( size_t iFacet=0U; iFacet<iNrOfFacets; ++iFacet ) {
     double64 facet_flux = 0;
     if (reuse_previous_velocity) {
       facet_flux = eptr->Read( iFacet, 0U, User()->ff_key );
     }
     else {
        // Calulate local normal
        Point<dim> facetNormal(helper.NormalOfFacet(iFacet));

        // Projection. Note that facetNormal is scaled by the area.
        facet_flux = dotProduct(facetNormal, vD);

        // storing the volumetric facet flux without altering the variables flag
        const auto ff_flag = eptr->Status( iFacet, 0U, User()->ff_key );
        eptr->Store( iFacet, 0U, User()->ff_key, makeScalar(ff_flag,facet_flux) );
     }

     auto outside_node = eptr->N(eptr->FV()->OutsideNode(iFacet));
     const double64 c_outside = outside_node->Read( c0_key );
     auto inside_node = eptr->N(eptr->FV()->InsideNode(iFacet));
     const double64 c_inside = inside_node->Read( c0_key );
     const double64 c_fip = helper.InterpolateScalarNodeProperty( eptr->FV()->FacetIntegrationPoint( iFacet, 0u ), c0_key );
     const auto ffc_flag = eptr->Status(iFacet, 0u, User()->ffC_key);

     if (fabs(facet_flux) > numeric_limits<double64>::epsilon()) {
       if (facet_flux > 0) {
         // XXX This is inefficient!
         std::pair<double64,double64> cminmax(c_inside, c_inside);
         for ( size_t iNeighbour=0U; iNeighbour < inside_node->Neighbors(); ++iNeighbour ) {
           const double64 c_neighbour(inside_node->Neighbor(iNeighbour)->Read( c0_key ));
           cminmax.first = std::min(cminmax.first, c_neighbour);
           cminmax.second = std::max(cminmax.second, c_neighbour);
         }

	       const double64 c = limitProperty( c_inside, c_outside, c_fip, cminmax );
         eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
       }
       else {
         // XXX This is inefficient!
         std::pair<double64,double64> cminmax(c_outside, c_outside);
         for ( size_t iNeighbour=0U; iNeighbour < outside_node->Neighbors(); ++iNeighbour ) {
           const double64 c_neighbour(outside_node->Neighbor(iNeighbour)->Read( c0_key ));
           cminmax.first = std::min(cminmax.first, c_neighbour);
           cminmax.second = std::max(cminmax.second, c_neighbour);
         }

	     const double64 c = limitProperty( c_outside, c_inside, c_fip, cminmax );
         eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
       }
     }
     else {
      const double64 c = facet_flux < 0 ? c_outside : c_inside;
      eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
     }
  }
   
 } // end Advective_O2_FluxesInterior



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
    auto c0_key = User()->C0_key;

     // 1. Dirichlet FV
     // ---------------
     // nothing needs to be done for FVs the saturation of which is flagged as Dirichlet
     if ( nd_ptr->Status( User()->C0_key ) == DIRICH ) {
          // the new value is initialised to the old one and the flag is kept
          nd_ptr->Store( User()->C1_key, makeScalar( nd_ptr->Status( User()->C1_key ), nd_ptr->Read( User()->C0_key ) ) );
          // TODO: is this the correct treatment of the flux balance
          // the flux balance is set to zero
          nd_ptr->Store( User()->fb_key, makeScalar( nd_ptr->Status( User()->fb_key ), 0. ) );
          // and returned
          return 0.;
       }
    
     // 2. a full finite volume is available so that influxes and outfluxes can be balanced
     // -----------------------------------------------------------------------------------
     if ( nd_ptr->AtBoundary() == NOT ) {
          double64 flux_balance(0.);
          const size_t node_parent_elements(nd_ptr->Parents());
          for ( size_t t=0U; t<node_parent_elements; t++ )
            {
              Element<dim>* const eptr(nd_ptr->Parent(t));
              assert( eptr != NULL );
              const size_t pnid(nd_ptr->ParentNodeNumber(t));

              const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
              for ( size_t i=0U; i<sector_facets; i++ )
                {
                   const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                   const auto ffc_flag = eptr->Status( iFacet, 0U, User()->ff_key );
                   const double64 facet_flux = eptr->Read( iFacet, 0U, User()->ff_key );
                   const size_t inside_node(eptr->FV()->InsideNode(iFacet));
                   const size_t outside_node(eptr->FV()->OutsideNode(iFacet));
                   const double64 c = (facet_flux < 0.) ? eptr->N(outside_node)->Read( c0_key ) :
                                                          eptr->N(inside_node)->Read( c0_key );

                   eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
                   const double64 sign = ( pnid == inside_node ) ? 1. : -1.;
                   flux_balance += sign * facet_flux;
                }
            } // end for loop for parent elements
     
          nd_ptr->Store( User()->fb_key, makeScalar( nd_ptr->Status( User()->fb_key ), flux_balance ) );
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
          Element<dim>* const eptr(nd_ptr->Parent(t));
          assert( eptr != NULL );
          const size_t pnid(nd_ptr->ParentNodeNumber(t));
           
          const double64 K(eptr->Read(User()->k_key));
           if (isnan(K)) {
             // XXX AJB HACK
             // Deleting boundaries during model creation means you can't set
             // properties during configuration. Retaining the boundaries means
             // that they are still included in the parents of a node.
             //
             // For now, we skip over any element which doesn't have a
             // permeability. They are not the flow domain.
             continue;
           }

          // for all FACETS per SECTOR surrounding the finite volume at the boundary
          // getting the volumetric fluxes only (upstream concentrations are found later)
          const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));

          for ( size_t i=0U; i<sector_facets; i++ )
            {
               const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
               const size_t inside_node(eptr->FV()->InsideNode(iFacet));
               const size_t outside_node(eptr->FV()->OutsideNode(iFacet));
               const double64 ff = eptr->Read( iFacet, 0u, User()->ff_key );

               const double64 C_upstream = (ff < 0.) ? eptr->N(outside_node)->Read( User()->C0_key ) :
                                                       eptr->N(inside_node)->Read( User()->C0_key );
               if ( pnid == inside_node ) {
                    inflow += ff;
                    influx += ff * C_upstream;
                 }
               else {
                    inflow -= ff;
                    influx -= ff * C_upstream;
                 }
            }
       } // end for loop for parent elements

      // the volume flux balance is stored (source terms are not subtracted if such were applied)
      nd_ptr->Store( User()->fb_key, makeScalar( nd_ptr->Status( User()->fb_key ), inflow ) );

      //   3.1: prescribed C value at inflow boundary
      //        in this case, inflow > 0 and influx > 0. The amount of concentration flux
      //        that enters through the boundary is inflow * C.
      //
      //   3.2: prescribed flux (which has an effect only at inflow boundary)
      //        this is the same as case 3.1 because nsrc should already have been added by this point.
      //
      //   3.3: free outflow (where flux balance missing the outflow facets is negative)
      //        in this case, inflow < 0.

      const double64 c0 = nd_ptr->Read( User()->C0_key );

      nd_ptr->Store( User()->C1_key, makeScalar( nd_ptr->Status( User()->C1_key ), influx - inflow * c0 ) );

      // the influx is returned
      return inflow; // positive when outgoing
   
 } // end AdvectiveFluxesAtBoundary



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
double64 FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O2_FluxesAtBoundary( Node<dim>* nd_ptr ) const
  {
     assert( nd_ptr != NULL );
     const auto c0_key = User()->C0_key;

     // 1. Dirichlet FV
     // ---------------
     // nothing needs to be done for FVs the saturation of which is flagged as Dirichlet
     if ( nd_ptr->Status( c0_key ) == DIRICH ) {
          // the new value is initialised to the old one and the flag is kept
          nd_ptr->Store( User()->C1_key, makeScalar( nd_ptr->Status( User()->C1_key ), nd_ptr->Read( c0_key ) ) );
          // TODO: is this the correct treatment of the flux balance
          // the flux balance is set to zero
          nd_ptr->Store( User()->fb_key, makeScalar( nd_ptr->Status( User()->fb_key ), 0. ) );
          // and returned
          return 0.;
       }
    
     // 2. a full finite volume is available so that influxes and outfluxes can be balanced
     // -----------------------------------------------------------------------------------
     if ( nd_ptr->AtBoundary() == NOT ) {
          double64 flux_balance(0.);
          const size_t node_parent_elements(nd_ptr->Parents());
          for ( size_t t=0U; t<node_parent_elements; t++ )
            {
              Element<dim>* const eptr(nd_ptr->Parent(t));
              assert( eptr != NULL );
              FiniteVolumeHelper<dim> helper(eptr);

              const size_t pnid(nd_ptr->ParentNodeNumber(t));

              const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
              for ( size_t i=0U; i<sector_facets; i++ )
                {
                   const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
                   const auto ffc_flag = eptr->Status( iFacet, 0U, User()->ff_key );

                   const size_t inside(eptr->FV()->InsideNode(iFacet));
                   auto inside_node = eptr->N(inside);
                   const double64 c_inside = inside_node->Read( c0_key );
                   const size_t outside(eptr->FV()->OutsideNode(iFacet));
                   auto outside_node = eptr->N(outside);
                   const double64 c_outside = outside_node->Read( c0_key );
                   const double64 c_fip = helper.InterpolateScalarNodeProperty( eptr->FV()->FacetIntegrationPoint( iFacet, 0u ), c0_key );
                   
                   const double64 facet_flux = eptr->Read( iFacet, 0U, User()->ff_key );

                   if (fabs(facet_flux) > numeric_limits<double64>::epsilon()) {
                     if (facet_flux > 0) {
                       // XXX This is inefficient!
                       std::pair<double64,double64> cminmax(c_inside, c_inside);
                       for ( size_t iNeighbour=0U; iNeighbour < inside_node->Neighbors(); ++iNeighbour ) {
                         const double64 c_neighbour(inside_node->Neighbor(iNeighbour)->Read( c0_key ));
                         cminmax.first = std::min(cminmax.first, c_neighbour);
                         cminmax.second = std::max(cminmax.second, c_neighbour);
                       }

                       const double64 c = limitProperty( c_inside, c_outside, c_fip, cminmax );
                       eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
                     }
                     else {
                       // XXX This is inefficient!
                       std::pair<double64,double64> cminmax(c_outside, c_outside);
                       for ( size_t iNeighbour=0U; iNeighbour < outside_node->Neighbors(); ++iNeighbour ) {
                         const double64 c_neighbour(outside_node->Neighbor(iNeighbour)->Read( c0_key ));
                         cminmax.first = std::min(cminmax.first, c_neighbour);
                         cminmax.second = std::max(cminmax.second, c_neighbour);
                       }

                       const double64 c = limitProperty( c_outside, c_inside, c_fip, cminmax );
                       eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
                     }
                   }
                   else {
                     // Just use the upstream concentration
                     const double64 c = facet_flux < 0 ? c_outside : c_inside;
                     eptr->Store( iFacet, 0u, User()->ffC_key, makeScalar(ffc_flag, c * facet_flux) );
                   }
                   const double64 sign = ( pnid == inside ) ? 1. : -1.;
                   flux_balance += sign * facet_flux;
                }
            } // end for loop for parent elements
     
          nd_ptr->Store( User()->fb_key, makeScalar( nd_ptr->Status( User()->fb_key ), flux_balance ) );
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
          Element<dim>* const eptr(nd_ptr->Parent(t));
          assert( eptr != NULL );
          FiniteVolumeHelper<dim> helper(eptr);
         
          const size_t pnid(nd_ptr->ParentNodeNumber(t));
           
          const double64 K(eptr->Read(User()->k_key));
           if (isnan(K)) {
             // XXX AJB HACK
             // Deleting boundaries during model creation means you can't set
             // properties during configuration. Retaining the boundaries means
             // that they are still included in the parents of a node.
             //
             // For now, we skip over any element which doesn't have a
             // permeability. They are not the flow domain.
             continue;
           }

          // for all FACETS per SECTOR surrounding the finite volume at the boundary
          // getting the volumetric fluxes only (upstream concentrations are found later)
          const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));

          for ( size_t i=0U; i<sector_facets; i++ ) {
            const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
            const size_t inside(eptr->FV()->InsideNode(iFacet));
            auto inside_node = eptr->N(inside);
            const double64 c_inside = inside_node->Read( c0_key );
            const size_t outside(eptr->FV()->OutsideNode(iFacet));
            auto outside_node = eptr->N(outside);
            const double64 c_outside = outside_node->Read( c0_key );
            const double64 c_fip = helper.InterpolateScalarNodeProperty( eptr->FV()->FacetIntegrationPoint( iFacet, 0u ), c0_key );
            
            const double64 facet_flux = eptr->Read( iFacet, 0U, User()->ff_key );
            double64 c(0.0);

            if (fabs(facet_flux) > numeric_limits<double64>::epsilon()) {
              if (facet_flux > 0) {
                // XXX This is inefficient!
                std::pair<double64,double64> cminmax(c_inside, c_inside);
                for ( size_t iNeighbour=0U; iNeighbour < inside_node->Neighbors(); ++iNeighbour ) {
                  const double64 c_neighbour(inside_node->Neighbor(iNeighbour)->Read( c0_key ));
                  cminmax.first = std::min(cminmax.first, c_neighbour);
                  cminmax.second = std::max(cminmax.second, c_neighbour);
                }

                c = limitProperty( c_inside, c_outside, c_fip, cminmax );
              }
              else {
                // XXX This is inefficient!
                std::pair<double64,double64> cminmax(c_outside, c_outside);
                for ( size_t iNeighbour=0U; iNeighbour < outside_node->Neighbors(); ++iNeighbour ) {
                  const double64 c_neighbour(outside_node->Neighbor(iNeighbour)->Read( c0_key ));
                  cminmax.first = std::min(cminmax.first, c_neighbour);
                  cminmax.second = std::max(cminmax.second, c_neighbour);
                }

                c = limitProperty( c_outside, c_inside, c_fip, cminmax );
              }
            }
            else {
              // essentially no flux means essentially no flow
              c = 0;
            }

            if ( pnid == inside ) {
              inflow += facet_flux;
              influx += facet_flux * c;
            }
            else {
              inflow -= facet_flux;
              influx -= facet_flux * c;
            }
          }
          
       } // end for loop for parent elements

      // the volume flux balance is stored (source terms are not subtracted if such were applied)
      nd_ptr->Store( User()->fb_key, makeScalar( nd_ptr->Status( User()->fb_key ), inflow ) );

      //   3.1: prescribed C value at inflow boundary
      //        in this case, inflow > 0 and influx > 0. The amount of concentration flux
      //        that enters through the boundary is inflow * C.
      //
      //   3.2: prescribed flux (which has an effect only at inflow boundary)
      //        this is the same as case 3.1 because nsrc should already have been added by this point.
      //
      //   3.3: free outflow (where flux balance missing the outflow facets is negative)
      //        in this case, inflow < 0.

      const double64 c0 = nd_ptr->Read( c0_key );

      nd_ptr->Store( User()->C1_key, makeScalar( nd_ptr->Status( User()->C1_key ), influx - inflow * c0 ) );

      // the influx is returned
      return inflow; // positive when outgoing
   
 } // end AdvectiveFluxesAtBoundary






/**
    stores and returns FV flux balance computed from current facet fluxes
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
               const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->ff_key );
               flux_balance += facet_flux;
            }
       }
     nptr->Store( User()->fb_key, makeScalar(nptr->Status(User()->fb_key),flux_balance) );
  
     //bmin = std::min( bmin, (*nit)->Read( User()->fb_key ) );
     //bmax = std::max( bmax, (*nit)->Read( User()->fb_key ) );
   
     return flux_balance;

 } // end FluxBalance

  
  
  




template class FacetFlux_TracerTransferExplicit<3U,ExplicitTransport>;

} // end csmp
