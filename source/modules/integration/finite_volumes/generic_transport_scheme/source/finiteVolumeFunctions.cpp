
//
//  finiteVolumeFunctions.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 9/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "finiteVolumeFunctions.h"

#include "Index.h"
#include "Model.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "CSMP_highLevelUtilities.h"
#include "Variables_TracerTransfer.h"

using namespace std;


namespace csmp {
  
  /**
   We store the volume of the finite volumes and their pore volumes
   
   - finite volume
   - finite volume (effective) pore volume
   - sector volume (stored at sector integration point)
   - sector weight: sector pore volume / finite volume pore volume = weighting factor
   
   Those parameters are also computed for the "halo elements" i.e., those elements
   that contribute sectors to the finite volumes that belong to the target region
   but are themselves not part of it.
   */
  template<size_t dim>
  void initializeFiniteVolumeProperties( Model<dim>& model, Region<dim>& gref )
  {
    const bool initialize_flux(true);
    const bool initialize_facet_area_perm(true);
    
    variables::Variables_TracerTransfer vars(model.Database());
    
    Point<dim> vD;

    // 0. zeroing sector pore volumes for accumulation in element loop
    // ---------------------------------------------------------------
    gref.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );

    // For the interior elements of the region compute relevant variable values
    const auto it_end(gref.ElementsEnd());
    for ( auto it=gref.ElementsBegin(); it!=it_end; ++it )
    {
      auto& e = **it;
      auto bctr = e.AtBarycenter();

      TensorVariable<dim> k;

      if ( initialize_facet_area_perm || initialize_flux ) {
        bctr.Obtain(vars.key_k, k);
      }
      
      if (initialize_facet_area_perm) {
        for ( auto fip : e.AllFacetIntegrationPoints() ) {
          VectorVariable<dim> fAk( fip.DirectedArea() * k );
          fip.Store( vars.key_fAk, fAk );
        }
      }

      // element-based total velocity
      if ( initialize_flux ) {
        const auto mu = model.Read(vars.key_MU); // XXX Should be able to interpolate
        const auto grad_p = bctr.Gradient(vars.key_PF);
        vD = -1.0/mu * (k * grad_p);
      }

      // 1. computing sector pore volumes
      // --------------------------------
      // (scaled by the cell thickness attribute=1 for volumetric elements)
      const auto phi = bctr.Read( vars.key_PHI ) * bctr.Read( vars.key_THI );
      
      for (auto n : e.AllNodes()) {
        // sector pore volume
        const double64 sector_volume = n.SectorVolume();
        // (*it)->Store( i, 0U, spv_key, makeScalar( PLAIN, phi * sector_volume ) );
        // sector volume is added to  pore volume of FV's containing this sector
        double64 pore_volume  = n.Read( vars.key_FVPV );
        // sector volume from FV traits
        pore_volume   += phi * sector_volume;
        n.Store( vars.key_FVPV, makeScalar(PLAIN,pore_volume) );
      }
      if ( initialize_flux ) {
        for (auto fip : e.AllFacetIntegrationPoints()) {
          const double64  facet_flux = fip.ProjectOntoDirectedArea(vD);
          fip.Store( vars.key_ff, makeScalar( PLAIN, facet_flux ) );
        }
      }
    }
    
    // 4. initialising sector pore volume in the elements surrounding perimeter nodes
    // -------------------------------------------------------------------------------------------------------------------
    // (here the pore volumes do not include the sectors outside the region)
#if 0
    const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
    for ( typename vector<Node<dim>*>::iterator nit=gref.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
      const size_t parent_elements((*nit)->Parents());
      for ( size_t i=0U; i<parent_elements; ++i ) {
        Element<dim>* const eptr = (*nit)->Parent(i);
        
        // ---------------------------------------
        // computing sector volumes & pore volumes
        // ---------------------------------------
        const double64 porosity = eptr->Read( phi_key );
        const size_t sectors(eptr->Sectors());
        for ( size_t j=0U; j<sectors; ++j ) {
          const double64 sector_volume = eptr->SectorVolume(j);
          eptr->Store( j, 0U, spv_key, makeScalar( PLAIN, sector_volume * porosity ) );
        }
      }
    }
#endif
    
    // 5. computing FV flux balances over the complete stencils
    // --------------------------------------------------------
    if ( initialize_flux ) {        // loop over FV stencils, computing the relevant variable values
      const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
      double64 bmin(1e30), bmax(-1e30);
      
      for ( typename vector<Node<dim>*>::iterator nit=gref.NodesBegin(); nit!=nit_end; ++nit ) {
        if ( (*nit)->AtBoundary() != NOT )
        {
          const size_t parent_elements((*nit)->Parents());
          double64 flux_balance(0.);
          for ( size_t i=0U; i<parent_elements; ++i ) {
            const Element<dim>* const eptr = (*nit)->Parent(i);
            const size_t sector_node      = (*nit)->ParentNodeNumber(i);
            for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
              const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
              const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
              const double64 facet_flux = sign * eptr->Read( facet, 0U, vars.key_ff );
              flux_balance += facet_flux;
            }
          }
          (*nit)->Store( vars.key_FB, makeScalar((*nit)->Status(vars.key_FB),flux_balance) );
          
          bmin = std::min( bmin, flux_balance );
          bmax = std::max( bmax, flux_balance );
        }
      }
      cout <<"\ninitializeFiniteVolumeProperties: initial flux balance: "<< std::max(fabs(bmin), fabs(bmax)) << endl;
    }

  } // end initializeFiniteVolumeProperties
  
  // explicit instantiation of function template in 2 and 3D
  template void initializeFiniteVolumeProperties( Model<1U>&, Region<1U>& );
  template void initializeFiniteVolumeProperties( Model<2U>&, Region<2U>& );
  template void initializeFiniteVolumeProperties( Model<3U>&, Region<3U>& );
  
  
  
  
  
  /* TESTING SECTOR INTEGRATION POINT STORAGE
   
   // sector storage: writing global node numbers to sector IP's and reading them out
   for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
   for ( size_t i=0U; i<(*it)->Sectors(); ++i )
   (*it)->Store( i, 0U, swt_key, makeScalar(PLAIN,(*it)->N(i)->Idx()) );
   
   // reading out node numbers and their double equivalents stored at the sector integration points
   for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it ) {
   cerr <<"\nelement: "<< (*it)->Idx() << endl;
   for ( size_t i=0U; i<(*it)->Sectors(); ++i ) {
   cerr << (*it)->N(i)->Idx() <<":";
   cerr << (*it)->Read( i, 0U, swt_key ) <<" ";
   }
   }
   */
  
  
  
  /*  TESTING FV volume calculations
   double volume(0.);
   for ( vector<Node<3U>*>::iterator it=ref.NodesBegin(); it!=ref.NodesEnd(); ++it )
   volume += (*it)->Read( fv_key );
   cerr <<"\nFV total volume: "<< volume;
   
   volume = 0.;
   for ( vector<Element<3U>*>::iterator it=ref.ElementsBegin(); it!=ref.ElementsEnd(); ++it )
   for ( size_t i=0U; i<(*it)->Sectors(); ++i )
   volume += (*it)->Read( i, 0U, sv_key );
   cerr <<"\nFV total volume: "<< volume;
   
   double pvolume(0.);
   for ( vector<Node<3U>*>::iterator it=ref.NodesBegin(); it!=ref.NodesEnd(); ++it )
   pvolume += (*it)->Read( fvphi_key );
   cerr <<"\nFV total volume: "<< pvolume;
   */
  
} // end csmp
