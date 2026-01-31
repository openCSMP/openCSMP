#include "UpwindFluxLHS.h"
#include "SparseMatrix.h"
#include "Element.h"
#include "Node.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
UpwindFluxLHS<dim>::UpwindFluxLHS( const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>& vol_flow_across_facet )
 : ff_key_(vol_flow_across_facet)
{
} // end constructor


/**
       Upwinded first-order in space and time.
*/
template<uint32_t dim>
void UpwindFluxLHS<dim>::AccumulateStencil( const Element<dim>& fe, SparseMatrix& mat ) const
{
    const auto facets(fe.Facets());
    for (uint32_t i(0); i < facets; ++i)
      {
         const auto inside_node(fe.FV()->InsideNode(i));
         const auto outside_node(fe.FV()->OutsideNode(i));
         auto w = fe.FV()->FacetIntegrationWeight(i, 0U);
         double facet_flux = w * fe.Read( i, 0U, ff_key_ );

         if ( facet_flux < 0. ) {
             facet_flux *= this->Factor();
             
             // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
             //    (outside node = upstream)
             // -------------------------------------------------------------------------------------------------
             mat.Add( fe.N(inside_node)->Idx(), fe.N(outside_node)->Idx(), facet_flux );  // incoming flux
            
             // 2. outgoing fluxes are added to the matrix diagonal
             // ---------------------------------------------------
             mat.Add( fe.N(outside_node)->Idx(), fe.N(outside_node)->Idx(), -facet_flux ); // outgoing flux
        } else {
             mat.Add( fe.N(outside_node)->Idx(), fe.N(inside_node)->Idx(), -facet_flux );  // incoming flux
             mat.Add( fe.N(inside_node)->Idx(), fe.N(inside_node)->Idx(), facet_flux );  // outgoing flux
        }   
    }     
     
} // end AccumulateStencil




template<uint32_t dim>
void UpwindFluxLHS<dim>::AccumulateFiniteVolume( const Node<dim>& fv, SparseMatrix& mat ) const
 {
    const auto node_parent_elements(fv.Parents());
    for ( auto t=0U; t<node_parent_elements; t++ )
    {
      const Element<dim>* const eptr(fv.Parent(t));   
      const auto pnid(fv.ParentNodeNumber(t));
      const auto sector_facets(eptr->FV()->FacetsPerSector(pnid));
      for ( auto i{0U}; i<sector_facets; i++ )
      {
          const auto iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
          const auto inside_node(eptr->FV()->InsideNode(iFacet));
          const auto outside_node(eptr->FV()->OutsideNode(iFacet));
          auto w = eptr->FV()->FacetIntegrationWeight(iFacet, 0U);   
          double facet_flux = w * eptr->Read( iFacet, 0U, ff_key_ );
          if (pnid != inside_node) facet_flux *= -1.;
          
          if ( facet_flux < 0. ) {
              facet_flux *= this->Factor();
              
              // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
              //    (outside node = upstream)
              // -------------------------------------------------------------------------------------------------
              mat.Add( eptr->N(inside_node)->Idx(), eptr->N(outside_node)->Idx(), facet_flux );  // incoming flux
            
              // 2. outgoing fluxes are added to the matrix diagonal
              // ---------------------------------------------------
              mat.Add( eptr->N(outside_node)->Idx(), eptr->N(outside_node)->Idx(), -facet_flux ); // outgoing flux
          } else {
              mat.Add( eptr->N(outside_node)->Idx(), eptr->N(inside_node)->Idx(), -facet_flux );  // incoming flux
              mat.Add( eptr->N(inside_node)->Idx(), eptr->N(inside_node)->Idx(), facet_flux );  // outgoing flux
          }   
      }
    } 

} // end AccumulateFiniteVolume

template class UpwindFluxLHS<1U>;
template class UpwindFluxLHS<2U>;
template class UpwindFluxLHS<3U>;

} // end csmp


