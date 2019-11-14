
#include "FluxLHS.h"
#include "SparseMatrix.h"
#include "Model.h"
#include "Element.h"
#include "Node.h"

using namespace std;

namespace csmp {
template<size_t dim>
FluxLHS<dim>::FluxLHS(Model<dim>& model, const char* facet_flux)
 : MatrixOperator<dim>(0),
   ff_key_(model.Database().StorageKey(facet_flux))
{
  if ( ff_key_.place != NODE || ff_key_.type != SCALAR )
    throw csmp::Exception( FATAL_ERROR, "FluxLHS<dim>::FluxLHS:",
                          "The 'facet flux' variable must be SCALAR and placed on NODE"  );

} // end constructor


template<size_t dim>
void FluxLHS<dim>::AccumulateStencil( Element<dim>& fe, SparseMatrix& mat ) const
{
    /*
    for (auto fip : fe.AllFacetIntegrationPoints()) {
        auto inside_node = fip.InsideNode();
        auto outside_node = fip.OutsideNode();
        auto w = fip.IntegrationWeight();

        double64 facet_flux = w * fip.Read( ff_key_ );
        if (this->MultiplyWithTimeIncrement()) {
          facet_flux *= this->dt_;
        }
        if ( facet_flux < 0. ) {
            // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
            //    (outside node = upstream)
            // -------------------------------------------------------------------------------------------------
            mat.Add( inside_node.NodeIdx(), outside_node.NodeIdx(), facet_flux );  // incoming flux
            
            // 2. outgoing fluxes are added to the matrix diagonal
            // ---------------------------------------------------
            mat.Add( outside_node.NodeIdx(), outside_node.NodeIdx(), -facet_flux ); // outgoing flux
        } else {
            mat.Add( outside_node.NodeIdx(), inside_node.NodeIdx(), -facet_flux );  // incoming flux
            mat.Add( inside_node.NodeIdx(), inside_node.NodeIdx(), facet_flux );  // outgoing flux
        }
    } 
    */  

    const size_t facets(fe.Facets());
    for (size_t i(0); i < facets; ++i) {
        const size_t inside_node(fe.FV()->InsideNode(i));
        const size_t outside_node(fe.FV()->OutsideNode(i));
        auto w = fe.FV()->FacetIntegrationWeight(i, 0U);
        double64 facet_flux = w * fe.Read( i, 0U, ff_key_ );
        if (this->MultiplyWithTimeIncrement()) {
          facet_flux *= this->dt_;
        }
        if ( facet_flux < 0. ) {
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


template<size_t dim>
void FluxLHS<dim>::AccumulateFiniteVolume( Node<dim>& fv, SparseMatrix& mat ) const
{
    /*
    for (auto fip : fv.AllFacetIntegrationPoints()) {
      auto inside_node = fip.InsideNode();
      auto outside_node = fip.OutsideNode();
      auto w = fip.IntegrationWeight();

      double64 facet_flux = w * fip.Read( this->ff_key_ );
      if (!fip.FromInside()) facet_flux *= -1.;
      if (this->MultiplyWithTimeIncrement()) {
        facet_flux *= this->dt_;
      }

      if ( facet_flux < 0. ) {
          // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
          //    (outside node = upstream)
          // -------------------------------------------------------------------------------------------------
          mat.Add( inside_node.NodeIdx(), outside_node.NodeIdx(), facet_flux );  // incoming flux
      
          // 2. outgoing fluxes are added to the matrix diagonal
          // ---------------------------------------------------
          mat.Add( outside_node.NodeIdx(), outside_node.NodeIdx(), -facet_flux ); // outgoing flux
      } else {
          mat.Add( outside_node.NodeIdx(), inside_node.NodeIdx(), -facet_flux );  // incoming flux
          mat.Add( inside_node.NodeIdx(), inside_node.NodeIdx(), facet_flux );  // outgoing flux
      }
    } 
    */ 

    const size_t node_parent_elements(fv.Parents());
    for ( size_t t=0U; t<node_parent_elements; t++ )
    {
      Element<dim>* const eptr(fv.Parent(t));   
      const size_t pnid(fv.ParentNodeNumber(t)); 
      const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
      for ( size_t i=0U; i<sector_facets; i++ )
      {
          const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,i) );
          const size_t inside_node(eptr->FV()->InsideNode(iFacet));
          const size_t outside_node(eptr->FV()->OutsideNode(iFacet)); 
          auto w = eptr->FV()->FacetIntegrationWeight(iFacet, 0U);   
          double64 facet_flux = w * eptr->Read( iFacet, 0U, ff_key_ );
          if (pnid != inside_node) facet_flux *= -1.;
          if (this->MultiplyWithTimeIncrement()) {
              facet_flux *= this->dt_;
          }

          if ( facet_flux < 0. ) {
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

template class FluxLHS<1U>;
template class FluxLHS<2U>;
template class FluxLHS<3U>;

} // end csmp


