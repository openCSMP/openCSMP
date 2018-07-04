
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
} // end AccumulateStencil


template<size_t dim>
void FluxLHS<dim>::AccumulateFiniteVolume( Node<dim>& fv, SparseMatrix& mat ) const
{
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
} // end AccumulateFiniteVolume

template class FluxLHS<1U>;
template class FluxLHS<2U>;
template class FluxLHS<3U>;

} // end csmp


