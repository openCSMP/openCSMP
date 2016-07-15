#include "MassBasedStencilProcessor.h"
#include "FV_Parameter.h"
#include "Element.h"

namespace csmp{



template<size_t dim>
MassBasedStencilProcessor<dim>::MassBasedStencilProcessor( const csmp::Index& adv_lhs_key,
                                                           const csmp::Index& adv_rhs_key,
                                                           const csmp::Index& velo_key )
 : ExplicitStencilProcessor<dim>(),
   adv_rhs_key_(adv_rhs_key)
 {
  this->adv1_key_ = adv_lhs_key;
  this->vel_key_ = velo_key;
 }

/**

Initializes sector pore volume and facet flux vectors for
current finite element = finite volume stencil.
*/

template<size_t dim>
void MassBasedStencilProcessor<dim>::InitializeFirstOrder(const FV_Parameter& param,
                                                          const Element<dim>& e ,
                                                          const VARIABLE_TYPE &vt,
                                                          const size_t var_comp_nr)
 {
     this->sector_pore_volume_.resize(e.Nodes());
     this->psi1_.resize(e.Nodes());
     this->psi2_.resize(e.Nodes());

     if (vt==SCALAR){
         // getting all the node-related information
         // ----------------------------------------
         for ( size_t i=0; i<e.Nodes(); i++ ) {
             // sector pore volumes
             this->sector_pore_volume_[i] = param.SectorVolume( i );
             // advected variable
             this->psi1_[i] = e.N(i)->Read( this->adv1_key_ );
             this->psi2_[i] = e.N(i)->Read( this->adv_rhs_key_ );
         }

         // getting all the finite-volume facet related information
         // -------------------------------------------------------
         this->facet_flux_.resize(e.FV_Stencil()->Facets());

         for ( size_t i=0; i<e.FV_Stencil()->Facets(); i++ ) {
             // project the velocities onto the facet normals to get
             // fluxes once the projections have been multiplied with
             // the surface areas
             this->facet_flux_[i]  = param.FacetNormalVelocity( i );
             this->facet_flux_[i] *= param.FacetArea( i );
         }
     }
     if (vt==ARRAY){
         if (adv_rhs_key_.type==FLAGGEDARRAY){
             FlaggedArrayVariable fav;
             ArrayVariable av;
             for ( size_t i=0; i<e.Nodes(); i++ ) {
                 // sector pore volumes
                 this->sector_pore_volume_[i] = param.SectorVolume( i );

                 // advected variable
                 e.N(i)->Read( this->adv1_key_,av );
                 this->psi1_[i] = av(var_comp_nr);

                 e.N(i)->Read( this->adv_rhs_key_,fav );
                 this->psi2_[i] = fav(var_comp_nr);
             }
         }
         else{
             ArrayVariable av;
             for ( size_t i=0; i<e.Nodes(); i++ ) {
                 // sector pore volumes
                 this->sector_pore_volume_[i] = param.SectorVolume( i );

                 // advected variable
                 e.N(i)->Read( this->adv1_key_,av );
                 this->psi1_[i] = av(var_comp_nr);

                 e.N(i)->Read( this->adv_rhs_key_,av );
                 this->psi2_[i] = av(var_comp_nr);
             }
         }

         this->facet_flux_.resize(e.FV_Stencil()->Facets());

         for ( size_t i=0; i<e.FV_Stencil()->Facets(); i++ ) {
             // project the velocities onto the facet normals to get
             // fluxes once the projections have been multiplied with
             // the surface areas
             this->facet_flux_[i]  = param.FacetNormalVelocity( i );
             this->facet_flux_[i] *= param.FacetArea( i );
         }

     }
     if (vt==FLAGGEDARRAY){
         if (adv_rhs_key_.type==FLAGGEDARRAY){
             FlaggedArrayVariable fav;
             for ( size_t i=0; i<e.Nodes(); i++ ) {
                 // sector pore volumes
                 this->sector_pore_volume_[i] = param.SectorVolume( i );
                 // advected variable
                 e.N(i)->Read( this->adv1_key_, fav );
                 this->psi1_[i] = fav(var_comp_nr);

                 e.N(i)->Read( this->adv_rhs_key_, fav );
                 this->psi2_[i] = fav(var_comp_nr);
             }
         }
         else
         {
             FlaggedArrayVariable fav;
             ArrayVariable av;
             for ( size_t i=0; i<e.Nodes(); i++ )
             {
                 // sector pore volumes
                 this->sector_pore_volume_[i] = param.SectorVolume( i );
                 // advected variable
                 e.N(i)->Read( this->adv1_key_, fav );
                 this->psi1_[i] = fav(var_comp_nr);

                 e.N(i)->Read( this->adv_rhs_key_, av );
                 this->psi2_[i] = av(var_comp_nr);
             }
         }

         this->facet_flux_.resize(e.FV_Stencil()->Facets());

         for ( size_t i=0; i<e.FV_Stencil()->Facets(); i++ ) {
             // project the velocities onto the facet normals to get
             // fluxes once the projections have been multiplied with
             // the surface areas
             this->facet_flux_[i]  = param.FacetNormalVelocity( i );
             this->facet_flux_[i] *= param.FacetArea( i );
         }
     }

 } // end InitializeFirstOrder

template struct MassBasedStencilProcessor<1U>;
template struct MassBasedStencilProcessor<2U>;
template struct MassBasedStencilProcessor<3U>;
 
}// end csmp

