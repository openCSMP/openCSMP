//  FacetFlux_TracerTransferExplicit.cpp
//
//  Created by Stephan Matthai on 2/21/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "FacetFlux_TracerTransferExplicit.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "TwoPhaseModel.h"
#include "ExplicitTransport.h"
#include "ImplicitTransport.h"
#include "finiteVolumeFunctions.h"
#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp {


  template<size_t dim, template<size_t> class USER>
  FacetFlux_TracerTransferExplicit<dim,USER>::FacetFlux_TracerTransferExplicit()
  {
  }


  template<size_t dim, template<size_t> class USER>
  void
  FacetFlux_TracerTransferExplicit<dim,USER>::FacetFluxes( Region<dim>& gref, bool reuse_velocity, bool second_order )
  {
    if (!second_order) {
      // 1. element-by-element processing of the facet fluxes
      auto elements_end(gref.ElementsEnd());
      for ( auto eit=gref.ElementsBegin(); eit!=elements_end; ++eit )
      {
        // 1.1 computation of facet fluxes (including upstream concentrations, but no-time increment yet)
        this->Advective_O1_FluxesInterior( reuse_velocity, **eit );
      }

      // 2. processing fluxes through the FVs on regions perimeter computing outside facet fluxes as necessary
      const typename vector<Node<dim>*>::iterator pnodes_end(gref.PerimeterNodesEnd());
      for ( typename vector<Node<dim>*>::iterator
           nit=gref.PerimeterNodesBegin(); nit!=pnodes_end; ++nit )
      {
        // computing flux balances and concentration-facet flux products where possible,
        // at sliced boundaries 3-typed of conditions are applied: 1) prescribed value (only at inflow),
        // 2) prescribed flux (has consequence only where there is inflow), 3) free outflow (outflow)
        // in this case the influx is found from the flux balance and FV cell's concentration
        this->Advective_O1_FluxesAtBoundary( **nit );
      }
    }
    else {
      // 1. element-by-element processing of the facet fluxes
      auto elements_end(gref.ElementsEnd());
      for ( auto eit=gref.ElementsBegin(); eit!=elements_end; ++eit )
      {
        // 1.1 computation of facet fluxes (including upstream concentrations, but no-time increment yet)
        this->Advective_O2_FluxesInterior( reuse_velocity, **eit );
      }

      // 2. processing fluxes through the FVs on regions perimeter computing outside facet fluxes as necessary
      const typename vector<Node<dim>*>::iterator pnodes_end(gref.PerimeterNodesEnd());
      for ( typename vector<Node<dim>*>::iterator
           nit=gref.PerimeterNodesBegin(); nit!=pnodes_end; ++nit )
      {
        // computing flux balances and concentration-facet flux products where possible,
        // at sliced boundaries 3-typed of conditions are applied: 1) prescribed value (only at inflow),
        // 2) prescribed flux (has consequence only where there is inflow), 3) free outflow (outflow)
        // in this case the influx is found from the flux balance and FV cell's concentration
        this->Advective_O2_FluxesAtBoundary( **nit );
      }
    }

  }

  /**
   Computes A_i vD . n_i for all facets and its product with the upstream concentrations
   of element stencil, storing it there.

   @attention that vD can be used directly as interstitial velocity since the FV is scaled by porosity

   @tested

   */
  
  bool approxeq(double64 x, double64 y) {
    return std::abs(x-y) < 1e-3;
  }
  template<size_t dim, template<size_t> class USER>
  void FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O1_FluxesInterior( bool reuse_previous_velocity, Element<dim>& e ) const
  {
    Point<dim> vDvec;
    if (!reuse_previous_velocity) {
      // Compute Darcy velocity
      auto bctr = e.AtBarycenter();
      const auto k = bctr.Interpolate(User()->key_k);
      const auto mu = User()->GetModel().Read(User()->key_MU); // XXX Should be able to Interpolate
      const auto grad_p = bctr.Gradient(User()->key_PF);

      VectorVariable<dim> vD(bctr.Read(User()->key_V));
      vD = -k/mu * grad_p;
      vDvec = vD.P();
      bctr.Store( User()->key_V, vD );
    }
    else {
      auto bctr = e.AtBarycenter();
      VectorVariable<dim> vD(bctr.Read(User()->key_V));
      vDvec = vD.P();
    }

    // computing total facet fluxes by projecting vt onto facet normals
    size_t i = 0;
    for (auto fip : e.AllFacetIntegrationPoints()) {
      double64 facet_flux = 0;

      if (reuse_previous_velocity) {
        facet_flux = fip.Read( User()->key_ff );
      }
      else {
        // Projection of vD onto the directed area of the facet.
        facet_flux = fip.ProjectOntoDirectedArea(User()->key_V);

        // storing the volumetric facet flux without altering the variables flag
        const auto ff_flag = fip.Status( User()->key_ff );
        fip.Store( User()->key_ff, makeScalar(ff_flag,facet_flux) );
      }

      // Get concentration from upwind node
      auto upstream_node = fip.UpstreamNode( User()->key_V );
      auto c = upstream_node.Read( User()->key_C );
      auto ffc = facet_flux * c;

      // Store facet flux concentration
      const auto ffc_flag = fip.Status(User()->key_ffC);
      fip.Store( User()->key_ffC, makeScalar(ffc_flag, ffc) );
    }

  } // end Advective_O1_FluxesInterior

  /**
   Computes A_i vD . n_i for all facets and its product with the upstream concentrations
   of element stencil, storing it there.

   @attention that vD can be used directly as interstitial velocity since the FV is scaled by porosity

   @tested

   */
  template<size_t dim, template<size_t> class USER>
  void FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O2_FluxesInterior( bool reuse_previous_velocity, Element<dim>& e ) const
  {
    Point<dim> vD;
    auto& key_C = User()->key_C;

    if (!reuse_previous_velocity) {
      // Compute Darcy velocity
      auto bctr = e.AtBarycenter();
      const auto k = bctr.Interpolate(User()->key_k);
      const auto mu = User()->GetModel().Read(User()->key_MU); // XXX Should be able to interpolate
      const auto grad_p = bctr.Gradient(User()->key_PF);

      VectorVariable<dim> vD;
      bctr.Read(User()->key_V, vD);
      vD = -k/mu * grad_p;
      bctr.Store( User()->key_V, vD );
    }

    // computing total facet fluxes by projecting vt onto facet normals
    for (auto fip : e.AllFacetIntegrationPoints()) {
      double64 facet_flux = 0;
      if (reuse_previous_velocity) {
        facet_flux = fip.Interpolate( User()->key_ff );
      }
      else {
        // Projection of vD onto the directed area of the facet.
        facet_flux = fip.ProjectOntoDirectedArea(User()->key_V);

        // storing the volumetric facet flux without altering the variables flag
        const auto ff_flag = fip.Status( User()->key_ff );
        fip.Store( User()->key_ff, makeScalar(ff_flag,facet_flux) );
      }

      auto outside_node = fip.OutsideNode();
      const double64 c_outside = outside_node.Read( key_C );
      auto inside_node = fip.InsideNode();
      const double64 c_inside = inside_node.Read( key_C );
      ScalarVariable cvar_fip;
      const double64 c_fip = fip.Interpolate( key_C );
      const auto ffc_flag = fip.Status(User()->key_ffC);

      if (fabs(facet_flux) > numeric_limits<double64>::epsilon()) {
        if (facet_flux > 0) {
          // XXX Is this too inefficient?
          std::pair<double64,double64> cminmax(+std::numeric_limits<double64>::max(),
                                               -std::numeric_limits<double64>::max());
          for ( auto neighbour : inside_node.AllNeighbourNodes() ) {
            auto c_neighbour = neighbour.Read( key_C );
            cminmax.first = std::min(cminmax.first, c_neighbour);
            cminmax.second = std::max(cminmax.second, c_neighbour);
          }

          const double64 c = limitProperty( c_inside, c_outside, c_fip, cminmax );
          fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
        }
        else {
          // XXX Is this too inefficient?
          std::pair<double64,double64> cminmax(+std::numeric_limits<double64>::max(),
                                               -std::numeric_limits<double64>::max());
          for ( auto neighbour : outside_node.AllNeighbourNodes() ) {
            const double64 c_neighbour(neighbour.Read( key_C ));
            cminmax.first = std::min(cminmax.first, c_neighbour);
            cminmax.second = std::max(cminmax.second, c_neighbour);
          }

          const double64 c = limitProperty( c_outside, c_inside, c_fip, cminmax );
          fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
        }
      }
      else {
        const double64 c = facet_flux < 0 ? c_outside : c_inside;
        fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
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
  double64 FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O1_FluxesAtBoundary( Node<dim>&n ) const
  {
    auto key_C = User()->key_C;

    // 1. Dirichlet FV
    // ---------------
    // nothing needs to be done for FVs the saturation of which is flagged as Dirichlet
    if ( n.Status( User()->key_C ) == DIRICH ) {
      // the new value is initialised to the old one and the flag is kept
      n.Store( User()->key_NC, makeScalar( n.Status( key_C ), n.Read( key_C ) ) );
      // TODO: is this the correct treatment of the flux balance
      // the flux balance is set to zero
      n.Store( User()->key_FB, makeScalar( n.Status( User()->key_FB ), 0. ) );
      // and returned
      return 0.;
    }

    // 2. a full finite volume is available so that influxes and outfluxes can be balanced
    // -----------------------------------------------------------------------------------
    if ( n.AtBoundary() == NOT ) {
      double64 flux_balance(0.);
      for (auto fip : n.AllFacetIntegrationPoints()) {
        const auto ffc_flag = fip.Status( User()->key_ffC );
        const double64 facet_flux = fip.Read( User()->key_ff );
        auto upstream_node = fip.UpstreamNode(facet_flux);
        const double64 c = upstream_node.Read( key_C );

        fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
        const double64 sign = ( fip.FromInside() ) ? 1. : -1.;
        auto w = fip.IntegrationWeight();
        flux_balance += sign * w * facet_flux;
      }

      n.Store( User()->key_FB, makeScalar( n.Status( User()->key_FB ), flux_balance ) );
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
    for (auto fip : n.AllFacetIntegrationPoints()) {
      auto inside_node = fip.InsideNode();
      auto outside_node = fip.OutsideNode();
      const double64 ff = fip.Read( User()->key_ff );

      const double64 C_upstream = (ff < 0.) ? outside_node.Read( User()->key_C ) :
      inside_node.Read( User()->key_C );
      const double64 sign = ( fip.FromInside() ) ? 1. : -1.;

      auto w = fip.IntegrationWeight();

      inflow += sign * w * ff;
      influx += sign * w * ff * C_upstream;
    }

    // the volume flux balance is stored (source terms are not subtracted if such were applied)
    n.Store( User()->key_FB, makeScalar( n.Status( User()->key_FB ), inflow ) );

    //   3.1: prescribed C value at inflow boundary
    //        in this case, inflow > 0 and influx > 0. The amount of concentration flux
    //        that enters through the boundary is inflow * C.
    //
    //   3.2: prescribed flux (which has an effect only at inflow boundary)
    //        this is the same as case 3.1 because nsrc should already have been added by this point.
    //
    //   3.3: free outflow (where flux balance missing the outflow facets is negative)
    //        in this case, inflow < 0.

    const double64 c0 = n.Read( User()->key_C );

    n.Store( User()->key_NC, makeScalar( n.Status( User()->key_NC ), influx - inflow * c0 ) );

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
  double64 FacetFlux_TracerTransferExplicit<dim,USER>::Advective_O2_FluxesAtBoundary( Node<dim>& n ) const
  {
    const auto key_C = User()->key_C;

    // 1. Dirichlet FV
    // ---------------
    // nothing needs to be done for FVs the saturation of which is flagged as Dirichlet
    if ( n.Status( key_C ) == DIRICH ) {
      // the new value is initialised to the old one and the flag is kept
      n.Store( User()->key_NC, makeScalar( n.Status( User()->key_NC ), n.Read( key_C ) ) );
      // TODO: is this the correct treatment of the flux balance
      // the flux balance is set to zero
      n.Store( User()->key_FB, makeScalar( n.Status( User()->key_FB ), 0. ) );
      // and returned
      return 0.;
    }

    // 2. a full finite volume is available so that influxes and outfluxes can be balanced
    // -----------------------------------------------------------------------------------
    if ( n.AtBoundary() == NOT ) {
      double64 flux_balance(0.);
      for (auto fip : n.AllFacetIntegrationPoints()) {
        const auto ffc_flag = fip.Status( User()->key_ffC );
        auto inside_node = fip.InsideNode();
        const double64 c_inside = inside_node.Read( key_C );
        auto outside_node = fip.OutsideNode();
        const double64 c_outside = outside_node.Read( key_C );
        const double64 c_fip = fip.Interpolate( key_C );
        const double64 facet_flux = fip.Read( User()->key_ff );

        if (fabs(facet_flux) > std::numeric_limits<double64>::epsilon()) {
          if (facet_flux > 0) {
            // XXX Is this too inefficient?
            std::pair<double64,double64> cminmax(+std::numeric_limits<double64>::max(),
                                                 -std::numeric_limits<double64>::max());
            for (auto neighbour : inside_node.AllNeighbourNodes()) {
              const double64 c_neighbour = neighbour.Read( key_C );
              cminmax.first = std::min(cminmax.first, c_neighbour);
              cminmax.second = std::max(cminmax.second, c_neighbour);
            }

            const double64 c = limitProperty( c_inside, c_outside, c_fip, cminmax );
            fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
          }
          else {
            // XXX Is this too inefficient?
            std::pair<double64,double64> cminmax(+std::numeric_limits<double64>::max(),
                                                 -std::numeric_limits<double64>::max());
            for (auto neighbour : outside_node.AllNeighbourNodes()) {
              const double64 c_neighbour = neighbour.Read( key_C );
              cminmax.first = std::min(cminmax.first, c_neighbour);
              cminmax.second = std::max(cminmax.second, c_neighbour);
            }

            const double64 c = limitProperty( c_outside, c_inside, c_fip, cminmax );
            fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
          }
        }
        else {
          // Just use the upstream concentration
          const double64 c = facet_flux < 0 ? c_outside : c_inside;
          fip.Store( User()->key_ffC, makeScalar(ffc_flag, c * facet_flux) );
        }
        const double64 sign = fip.FromInside() ? 1. : -1.;
        const auto w = fip.IntegrationWeight();
        flux_balance += sign * w * facet_flux;
      }

      n.Store( User()->key_FB, makeScalar( n.Status( User()->key_FB ), flux_balance ) );
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
    for (auto fip : n.AllFacetIntegrationPoints()) {
      auto inside_node = fip.InsideNode();
      const double64 c_inside = inside_node.Read( key_C );
      auto outside_node = fip.OutsideNode();
      const double64 c_outside = outside_node.Read( key_C );
      const double64 c_fip = fip.Interpolate(key_C);

      const double64 facet_flux = fip.Read( User()->key_ff );
      double64 c(0.0);

      if (fabs(facet_flux) > numeric_limits<double64>::epsilon()) {
        if (facet_flux > 0) {
          // XXX Is this too inefficient?
          std::pair<double64,double64> cminmax(+std::numeric_limits<double64>::max(),
                                               -std::numeric_limits<double64>::max());
          for ( auto neighbour : inside_node.AllNeighbourNodes() ) {
            const double64 c_neighbour = neighbour.Read( key_C );
            cminmax.first = std::min(cminmax.first, c_neighbour);
            cminmax.second = std::max(cminmax.second, c_neighbour);
          }

          c = limitProperty( c_inside, c_outside, c_fip, cminmax );
        }
        else {
          // XXX Is this too inefficient?
          std::pair<double64,double64> cminmax(+std::numeric_limits<double64>::max(),
                                               -std::numeric_limits<double64>::max());
          for ( auto neighbour : outside_node.AllNeighbourNodes() ) {
            const double64 c_neighbour = neighbour.Read( key_C );
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

      const auto w = fip.IntegrationWeight();
      if ( fip.FromInside() ) {
        inflow += w * facet_flux;
        influx += w * facet_flux * c;
      }
      else {
        inflow -= w * facet_flux;
        influx -= w * facet_flux * c;
      }
    }

    // the volume flux balance is stored (source terms are not subtracted if such were applied)
    n.Store( User()->key_FB, makeScalar( n.Status( User()->key_FB ), inflow ) );

    //   3.1: prescribed C value at inflow boundary
    //        in this case, inflow > 0 and influx > 0. The amount of concentration flux
    //        that enters through the boundary is inflow * C.
    //
    //   3.2: prescribed flux (which has an effect only at inflow boundary)
    //        this is the same as case 3.1 because nsrc should already have been added by this point.
    //
    //   3.3: free outflow (where flux balance missing the outflow facets is negative)
    //        in this case, inflow < 0.

    const double64 c0 = n.Read( key_C );

    n.Store( User()->key_NC, makeScalar( n.Status( User()->key_NC ), influx - inflow * c0 ) );

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
      const Element<dim>* const eptr = nptr->Parent(i);
      const size_t sector_node      = nptr->ParentNodeNumber(i);
      for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
        const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
        const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
        const double64 facet_flux = sign * eptr->Read( facet, 0U, User()->key_ff );
        flux_balance += facet_flux;
      }
    }
    nptr->Store( User()->key_FB, makeScalar(nptr->Status(User()->key_FB),flux_balance) );

    //bmin = std::min( bmin, (*nit)->Read( User()->key_FB ) );
    //bmax = std::max( bmax, (*nit)->Read( User()->key_FB ) );

    return flux_balance;

  } // end FluxBalance


  template class FacetFlux_TracerTransferExplicit<1U,ImplicitTransport>;
  template class FacetFlux_TracerTransferExplicit<2U,ImplicitTransport>;
  template class FacetFlux_TracerTransferExplicit<3U,ImplicitTransport>;
  template class FacetFlux_TracerTransferExplicit<1U,ExplicitTransport>;
  template class FacetFlux_TracerTransferExplicit<2U,ExplicitTransport>;
  template class FacetFlux_TracerTransferExplicit<3U,ExplicitTransport>;

} // end csmp
