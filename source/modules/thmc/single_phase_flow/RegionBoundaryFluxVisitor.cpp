// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  RegionBoundaryFluxVisitor.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/30/14.
//

#include "RegionBoundaryFluxVisitor.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "PropertyDatabase.h"
#include "ErrorHandler.h"
#include "Node.h"

namespace csmp {

/** @attention EXPERIMENTAL CODE

    After construction this visitor should be passed to the region of interest
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
RegionBoundaryFluxVisitor<dim,COMPUTATION_DOMAIN>::RegionBoundaryFluxVisitor( Model<dim>& model,
                                                                              const char*  target_region,
                                                                              const char*  Darcy_velocity )
 : Visitor<dim>(MODEL,NODE),
   flux_key_(model.Database().StorageKey(Darcy_velocity)),
   FVinflux_(0.),
   delta_t_(0.),
   domain_ptr_(&model.Region(target_region))
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( domain_ptr_ == NULL )
      csmp_error.Note( ERROR,
                        "RegionBoundaryFluxVisitor<dim>(custom constructor):",
                        "the target region remains undefined.");

    if ( flux_key_.place != ELEMENT )
      csmp_error.Note( ERROR,
                        "RegionBoundaryFluxVisitor<dim>(custom constructor):",
                        Darcy_velocity, "must be an element property.");

    if ( flux_key_.type != VECTOR )
      csmp_error.Note( ERROR,
                        "RegionBoundaryFluxVisitor<dim>(custom constructor):",
                        Darcy_velocity, "must be a vector property.");
 } // end
 
  


/**
    Sitting on a boundary node of the target Region
    We look at the surrounding finite volume sectors that do not 
    form part of this region and calculate the incoming fluxes
    into the region using the user-supplied flux variable.
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void RegionBoundaryFluxVisitor<dim,COMPUTATION_DOMAIN>::Visit( Node<dim>* nd )
{
  // only if the node lies at the boundary of the region, a flux needs to be computed
  if ( !domain_ptr_->IsPerimeterNode(nd) ) return;
  
  // looping over the parent elements
  for ( auto i{0U}; i<nd->Parents(); ++i )
    {
       // ascertain that the pointer is valid
       if ( Element<dim>* const eptr = nd->Parent(i) )
         if ( domain_ptr_->Contains(eptr) == false )
           {
              assert( eptr != NULL );
              // getting the flux variable
              eptr->Read( flux_key_, vt_ );
              // find out which node were are on from a view-point of the parent element
              const auto parent_nd(nd->ParentNodeNumber(i));
              assert( eptr->FV() != NULL /* if so the FV stencils may not have been initialized yet */ );
              const auto n_sector_facets(eptr->FV()->FacetsPerSector(parent_nd));
              for ( auto j{0U}; j<n_sector_facets; j++ )
                {
                   const auto facet = eptr->FV()->FacetSurroundingSector( parent_nd, j );
                   // if the sector node is the inside node then an incoming flux will create a positive source term
                   const double  fsign = (parent_nd == nd->Parent(i)->FV()->InsideNode(facet)) ? -1. : 1.;
                   FVinflux_ += fsign * eptr->FacetArea(facet) * eptr->ProjectionOnFacetNormal(facet,vt_) * delta_t_;
                }
           }
    } //
  
} // end


/// returns the flux through the region boundary
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
double RegionBoundaryFluxVisitor<dim,COMPUTATION_DOMAIN>::InFlux() const
{
    return FVinflux_;
} // end


/// resets the internal flux counter to zero
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void RegionBoundaryFluxVisitor<dim,COMPUTATION_DOMAIN>::ResetFlux()
{
   FVinflux_ = 0.;
} // end


/// use this to set the time interval over which the flux shall be integrated
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void RegionBoundaryFluxVisitor<dim,COMPUTATION_DOMAIN>::TimeIncrement( double dt )
{
   delta_t_ = dt;
} // end



// explicit instantiation of template for 2D models
template class RegionBoundaryFluxVisitor<2U,Region>;

} // end csmp
