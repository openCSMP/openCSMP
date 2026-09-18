// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FluidCompressibilitySourceLHS.h"
#include "SparseMatrix.h"
#include "Element.h"
#include "Node.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
FluidCompressibilitySourceLHS<dim>::FluidCompressibilitySourceLHS( const csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& key_SPV,
                                                                   const csmp::INDEX<SCALAR,ELEMENT>& key_PHI,
                                                                   const csmp::INDEX<SCALAR,ELEMENT>& key_CT,
                                                                   const csmp::INDEX<SCALAR,NODE>& key_PF0,
                                                                   const csmp::INDEX<SCALAR,NODE>& key_PF1,
                                                                   bool interpolate_pf_to_sector_ip  )
  : key_SPV_(key_SPV),
    key_PHI_(key_PHI),
    key_CT_(key_CT),
    key_PF0_(key_PF0),
    key_PF1_(key_PF1),
    interpolate_pf_to_sector_ip_(interpolate_pf_to_sector_ip)
  {
  } // end constructor


/**
    Element-wise accumulation of source term/
*/  
template<uint32_t dim>
void FluidCompressibilitySourceLHS<dim>::AccumulateStencil( const Element<dim>& e, SparseMatrix& mat ) const
  {
    const double phi = e.Read( key_PHI_ );
    const double ct  = e.Read( key_CT_ );

    const auto sectors{ e.Sectors() };
    for ( auto i{0U}; i<sectors; ++i )
     {
        const double sector_PV = (interpolate_pf_to_sector_ip_) ? e.Read( i, 0U, key_SPV_ ) : e.SectorVolume(i) * phi;  
        double pf0(0.), pf1(0.);
        // getting pressures p0 and p1 at sector integration points 
        if ( interpolate_pf_to_sector_ip_ == true ) {
             const uint32_t s_ip(0U), n_nodes(e.Nodes()); // sector integration point and element nodes
             e.N_AtSectorIntegrationPoint( i, s_ip );
             // interpolate pressures to sector integration point
             for ( uint32_t j{0U}; j<n_nodes; ++j ) {
                  pf0 += e.FE()->NRST[j] * e.N(j)->Read( key_PF0_ );
                  pf1 += e.FE()->NRST[j] * e.N(j)->Read( key_PF1_ );
               }
          }
        else { // the node=FV pressure values are used directly
             pf0 = e.N(i)->Read( key_PF0_ );
             pf1 = e.N(i)->Read( key_PF1_ );
          }
        const size_t node_idx = e.N(i)->Idx();
        mat.Add( node_idx, node_idx, sector_PV * ct * (pf1-pf0) * this->Factor() );
     }

  } // end AccumulateStencil
  
  
  /**
    Since each sector of the finite volume may consist of a different rock type, the total systems compressibility
    must be summed up sector-by-sector, honouring the different contributions.
        
        @note using a single pressure value is perhaps not rigorous enough. One should rather interpolate
        the pressure to the sector integration in order to capture the source term more accurately.
  */ 
  template<uint32_t dim>
  void FluidCompressibilitySourceLHS<dim>::AccumulateFiniteVolume( const Node<dim>& fv, SparseMatrix& mat ) const
  {
    // getting properties from parent elements
    const auto parent_elements{fv.Parents()};
    for ( auto i{0U}; i<parent_elements; ++i ) {
        const Element<dim>* const eptr = fv.Parent(i);
        const auto n_node = fv.ParentNodeNumber(i);

        const double sector_PV = (interpolate_pf_to_sector_ip_) ? eptr->Read( n_node, 0U, key_SPV_ ) : 
                                                                  eptr->SectorVolume(n_node) * eptr->Read( key_PHI_ );
        const double ct = eptr->Read( key_CT_ ); 
        double pf0(0.), pf1(0.);
        // getting pressures p0 and p1 at sector integration points 
        if ( interpolate_pf_to_sector_ip_ == true ) {
             const auto s_ip(0U), n_nodes(eptr->Nodes()); // sector integration point and element nodes
             eptr->N_AtSectorIntegrationPoint( n_node, s_ip );
             // interpolate pressures to sector integration point
             for ( auto j{0U}; j<n_nodes; ++j ) {
                  pf0 += eptr->FE()->NRST[j] * eptr->N(j)->Read( key_PF0_ );
                  pf1 += eptr->FE()->NRST[j] * eptr->N(j)->Read( key_PF1_ );
               }
          }
        else { // the node=FV pressure values are used directly
             pf0 = fv.Read( key_PF0_ );
             pf1 = fv.Read( key_PF1_ );
          }
        // perform sector-by-sector accumulation
        mat.Add( fv.Idx(), fv.Idx(), sector_PV * ct * (pf1 - pf0) * this->Factor()  );
     }

  } // end AccumulateFiniteVolume
  
//  template class FluidCompressibilitySourceLHS<1U>;
//  template class FluidCompressibilitySourceLHS<2U>;
  template class FluidCompressibilitySourceLHS<3U>;
  
} // end csmp


