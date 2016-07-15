//
//  finiteVolumeUniversalFunctions.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 9/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "finiteVolumeUniversalFunctions.h"

#include "Index.h"
#include "Model.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "CSMP_highLevelUtilities.h"

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
   
    const csmp::Index phi_key = model.Database().StorageKey("porosity");
    const csmp::Index thi_key = model.Database().StorageKey("thickness");
    const csmp::Index vt_key  = model.Database().StorageKey("velocity");

    const csmp::Index fv_key  = model.Database().StorageKey("finite volume");
    const csmp::Index pv_key  = model.Database().StorageKey("FV pore volume");
    const csmp::Index sv_key  = model.Database().StorageKey("sector volume");
    const csmp::Index spv_key = model.Database().StorageKey("sector pore volume");
    const csmp::Index fa_key  = model.Database().StorageKey("facet area");
    const csmp::Index ff_key  = model.Database().StorageKey("facet flux");
    const csmp::Index fn_key  = model.Database().StorageKey("facet normal");
    const csmp::Index fb_key  = model.Database().StorageKey("flux balance");
   
    Point<dim>           nrml;
    VectorVariable<dim>  fnrml, vt;
   
    // 0. zeroing sector pore volumes for accumulation in element loop
    // ---------------------------------------------------------------
    gref.InputPropertyValue( "finite volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );
   
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref.ElementsBegin(); it!=it_end; ++it )
      {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // element-based total velocity
         if ( initialize_flux ) (*it)->Read( vt_key, vt );

         // 1. computing sector pore volumes
         // --------------------------------
         // (scaled by the cell thickness attribute=1 for volumetric elements)
         const double64 phi = (*it)->Read( phi_key ) * (*it)->Read( thi_key );
         for ( size_t i=0U; i<sectors; ++i ) {
              // sector pore volume
              const double64 sector_volume = (*it)->SectorVolume(i);
              (*it)->Store( i, 0U, sv_key, makeScalar( PLAIN, sector_volume ) );
              (*it)->Store( i, 0U, spv_key, makeScalar( PLAIN, phi * sector_volume ) );
              // sector volume is added to  pore volume of FV's containing this sector
              double64 finite_volume = (*it)->N(i)->Read( fv_key );
              double64 pore_volume   = (*it)->N(i)->Read( pv_key );
              // sector volume from FV traits
              finite_volume += sector_volume;
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( fv_key, makeScalar(PLAIN,finite_volume) );
              (*it)->N(i)->Store( pv_key, makeScalar(PLAIN,pore_volume) );
           }

         // 2. computing facet normals and areas
         // ------------------------------------
         for ( size_t j=0U; j<facets; ++j ) {
              // computing facet areas
              const double64 facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, fa_key, makeScalar( PLAIN, facet_area ) );
              // computing facet normals
              nrml = (*it)->FacetNormal(j);
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, fn_key, fnrml );
           
              // 3. computing total facet fluxes and flux balance
              // ------------------------------------------------
              if ( initialize_flux ) {
                   double64 facet_flux(nrml[0] * vt[0]);
                   if ( dim != 1U ) facet_flux += nrml[1] * vt[1];
                   if ( dim == 3U ) facet_flux += nrml[2] * vt[2];
                   facet_flux *= facet_area;
                   (*it)->Store( j, 0U, ff_key, makeScalar((*it)->Status( j, 0U, ff_key),facet_flux) );
                }
           }
      }

   // 4. initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // -------------------------------------------------------------------------------------------------------------------
   // (here the pore volumes do not include the sectors outside the region)
   const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
   
   for ( typename vector<Node<dim>*>::iterator nit=gref.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const size_t parent_elements((*nit)->Parents());
        for ( size_t i=0U; i<parent_elements; ++i ) {
             Element<dim>* const eptr = (*nit)->Parent(i);
             // ---------------------------------
             // computing facet normals and areas
             // ---------------------------------
             const size_t facets(eptr->Facets());
             for ( size_t j=0U; j<facets; ++j ) {
                  // computing facet areas
                  const double64 facet_area = eptr->FacetArea(j);
                  eptr->Store( j, 0U, fa_key, makeScalar( PLAIN, facet_area ) );
                  // computing facet normals
                  nrml = eptr->FacetNormal(j);
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, fn_key, fnrml );
               }
             // ---------------------------------------
             // computing sector volumes & pore volumes
             // ---------------------------------------
             const double64 porosity = eptr->Read( phi_key );
             const size_t sectors(eptr->Sectors());
             for ( size_t j=0U; j<sectors; ++j ) {
                  const double64 sector_volume = eptr->SectorVolume(j);
                  eptr->Store( j, 0U, sv_key, makeScalar( PLAIN, sector_volume ) );
                  eptr->Store( j, 0U, spv_key, makeScalar( PLAIN, sector_volume * porosity ) );
               }
          }
     }
         
   // 5. computing FV flux balances over the complete stencils
   // --------------------------------------------------------
   if ( initialize_flux ) {
        // loop over FV stencils, computing the relevant variable values
     
        const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
        double64 bmin(1e30), bmax(-1e30);
     
        for ( typename vector<Node<dim>*>::iterator nit=gref.NodesBegin(); nit!=nit_end; ++nit )
          if ( (*nit)->AtBoundary() != NOT )
            {
               const size_t parent_elements((*nit)->Parents());
               double64 flux_balance(0.);
               for ( size_t i=0U; i<parent_elements; ++i ) {
                    const Element<dim>* const eptr = (*nit)->Parent(i);
                    const size_t sector_node      = (*nit)->ParentNodeNumber(i);
                    for ( size_t j=0U; j<eptr->FV_Stencil()->FacetsPerSector(sector_node); ++j ) {
                         const size_t facet = eptr->FV_Stencil()->FacetSurroundingSector( sector_node, j );
                         const double64 sign = (sector_node==eptr->FV_Stencil()->InsideNode(facet)) ? 1. : -1.;
                         const double64 facet_flux = sign * eptr->Read( facet, 0U, ff_key );
                         flux_balance += facet_flux;
                      }
                 }
               (*nit)->Store( fb_key, makeScalar((*nit)->Status(fb_key),flux_balance) );
            
               bmin = std::min( bmin, (*nit)->Read( fb_key ) );
               bmax = std::max( bmax, (*nit)->Read( fb_key ) );
            }
        cout <<"\ninitializeFiniteVolumeProperties: initial flux balance: "<< std::max(fabs(bmin), fabs(bmax)) << endl;
     }
 
 } // end initializeFiniteVolumeProperties

// explicit instantiation of function template in 2 and 3D
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
