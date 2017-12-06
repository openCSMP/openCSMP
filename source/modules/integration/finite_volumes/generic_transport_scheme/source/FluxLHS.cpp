
#include "FluxLHS.h"
#include "SparseMatrix.h"
#include "Model.h"
#include "Element.h"
#include "Node.h"

using namespace std;

namespace csmp {
template<size_t dim>
FluxLHS<dim>::FluxLHS(Model<dim>& model, const char* facet_flux)
 : MatrixOperator<dim>(ADD_ACCUMULATE),
   ff_key_(model.Database().StorageKey(facet_flux))
{
} // end constructor


template<size_t dim>
void FluxLHS<dim>::AccumulateStencil( const Element<dim>* fe, SparseMatrix& mat ) const
{
    const size_t iNrOfFacets(fe->FV()->Facets());
    for ( size_t iFacet=0U; iFacet<iNrOfFacets; ++iFacet ) {
        const size_t inside(fe->FV()->InsideNode(iFacet));
        auto inside_node = fe->N(inside);
        const size_t outside(fe->FV()->OutsideNode(iFacet));
        auto outside_node = fe->N(outside);

        const double64 facet_flux = fe->Read( iFacet, 0U, this->ff_key_ );
        if ( facet_flux < 0. ) {
            // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
            //    (outside node = upstream)
            // -------------------------------------------------------------------------------------------------
            mat.Add( inside_node->Idx(), outside_node->Idx(), facet_flux );  // incoming flux
            
            // 2. outgoing fluxes are added to the matrix diagonal
            // ---------------------------------------------------
            mat.Add( outside_node->Idx(), outside_node->Idx(), -facet_flux ); // outgoing flux
        } else {
            mat.Add( outside_node->Idx(), inside_node->Idx(), -facet_flux );  // incoming flux
            mat.Add( inside_node->Idx(), inside_node->Idx(), facet_flux );  // outgoing flux
        }
    }        
} // end AccumulateStencil


template<size_t dim>
void FluxLHS<dim>::AccumulateFiniteVolume( const Node<dim>* nd, SparseMatrix& mat ) const
{
    const size_t parent_elements(nd->Parents());
    for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<dim>* const eptr = nd->Parent(i);
         assert( eptr != NULL );
         const size_t pnid(nd->ParentNodeNumber(i));
         const size_t sector_facets(eptr->FV()->FacetsPerSector(pnid));
         for ( size_t j=0U; j<sector_facets;++j ) {
              const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,j) );
              const size_t inside(eptr->FV()->InsideNode(iFacet));
              auto inside_node = eptr->N(inside);
              const size_t outside(eptr->FV()->OutsideNode(iFacet));
              auto outside_node = eptr->N(outside);

              double64 facet_flux = eptr->Read( iFacet, 0U, this->ff_key_ );
              if (pnid != inside) facet_flux *= -1.;

              if ( facet_flux < 0. ) {
                  // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
                  //    (outside node = upstream)
                  // -------------------------------------------------------------------------------------------------
                  mat.Add( inside_node->Idx(), outside_node->Idx(), facet_flux );  // incoming flux
            
                  // 2. outgoing fluxes are added to the matrix diagonal
                  // ---------------------------------------------------
                  mat.Add( outside_node->Idx(), outside_node->Idx(), -facet_flux ); // outgoing flux
              } else {
                  mat.Add( outside_node->Idx(), inside_node->Idx(), -facet_flux );  // incoming flux
                  mat.Add( inside_node->Idx(), inside_node->Idx(), facet_flux );  // outgoing flux
              }        
          }   
    }         
} // end AccumulateFiniteVolume

template class FluxLHS<1U>;
template class FluxLHS<2U>;
template class FluxLHS<3U>;

} // end csmp


