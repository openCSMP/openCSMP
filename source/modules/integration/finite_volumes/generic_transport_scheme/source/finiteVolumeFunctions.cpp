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
    const csmp::Index k_key  = model.Database().StorageKey("permeability");
     const csmp::Index pf_key  = model.Database().StorageKey("fluid pressure");

    const csmp::Index pv_key  = model.Database().StorageKey("FV pore volume");
    // const csmp::Index spv_key = model.Database().StorageKey("sector pore volume");
    const csmp::Index ff_key  = model.Database().StorageKey("facet flux");
    const csmp::Index fb_key  = model.Database().StorageKey("flux balance");
   
     Point<dim> vD;
   
    // 0. zeroing sector pore volumes for accumulation in element loop
    // ---------------------------------------------------------------
    gref.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );
   
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref.ElementsBegin(); it!=it_end; ++it )
      {
         const size_t sectors((*it)->Sectors());
          const size_t facets((*it)->Facets());

         FiniteVolumeHelper<dim> helper(*it);

         // element-based total velocity
          if ( initialize_flux ) {
              const double64 K((*it)->Read(k_key));
              vD = -K * helper.GradientOfScalarNodeProperty((*it)->FV()->Barycenter(), pf_key);
          }

         // 1. computing sector pore volumes
         // --------------------------------
         // (scaled by the cell thickness attribute=1 for volumetric elements)
         const double64 phi = (*it)->Read( phi_key ) * (*it)->Read( thi_key );
         for ( size_t i=0U; i<sectors; ++i ) {
              // sector pore volume
              const double64 sector_volume = (*it)->SectorVolume(i);
              // (*it)->Store( i, 0U, spv_key, makeScalar( PLAIN, phi * sector_volume ) );
              // sector volume is added to  pore volume of FV's containing this sector
              double64 pore_volume   = (*it)->N(i)->Read( pv_key );
              // sector volume from FV traits
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( pv_key, makeScalar(PLAIN,pore_volume) );
           }
          if ( initialize_flux ) {
              for ( size_t i=0U; i<facets; ++i ) {
                  Point<dim> facetNormal(helper.NormalOfFacet(i));
                  const double64  facet_flux = dotProduct(vD, facetNormal);
                  (*it)->Store( i, 0U, ff_key, makeScalar( PLAIN, facet_flux ) );
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
                         const double64 facet_flux = sign * eptr->Read( facet, 0U, ff_key );
                         flux_balance += facet_flux;
                      }
                 }
               (*nit)->Store( fb_key, makeScalar((*nit)->Status(fb_key),flux_balance) );
            
               bmin = std::min( bmin, (*nit)->Read( fb_key ) );
               bmax = std::max( bmax, (*nit)->Read( fb_key ) );
            }
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





template<size_t dim>
FiniteVolumeHelper<dim>::FiniteVolumeHelper(Element<dim>* eptr)
    : eptr_(eptr)
{
    if (eptr_->IsLineElement()) {
        element_dim_ = 1;
    }
    else if (eptr_->IsSurfaceElement()) {
        element_dim_ = 2;
    }
    else if (eptr_->IsVolumeElement()) {
        element_dim_ = 3;
    }

    num_nodes_ = eptr->Nodes();
    eptr->CoordinateMatrix();
}



template<size_t dim>
double64
FiniteVolumeHelper<dim>::InterpolateScalarNodeProperty( const Point<dim>& p, const csmp::Index& prop ) const
{
    std::vector<double64> NRST;
  
    CalculateN(p, NRST);

    double64 var(0.0);
    const size_t nodes(eptr_->Nodes());
    for ( size_t i=0; i<nodes; i++ ) {
      var += eptr_->N(i)->Read( prop ) * NRST[i];
    }

    return var;
}




template<size_t dim>
Point<dim>
FiniteVolumeHelper<dim>::GradientOfScalarNodeProperty( const Point<dim>& p, const csmp::Index& prop ) const
{
    std::vector<double64> DN[dim];
    for (unsigned i = 0; i < element_dim_; ++i) {
        DN[i].resize(num_nodes_);
    }
    CalculateDN(p, DN);


    Point<dim> grad(0.);
    for ( size_t i=0U; i<num_nodes_; ++i ) {
        const double64 value_at_node(eptr_->N(i)->Read(prop));
        for (size_t j = 0; j < element_dim_; ++j) {
            grad[j] += DN[j][i] * value_at_node;
        }
    }
    eptr_->FE()->JacobianInverse();
    
    return Point<dim>(eptr_->FE()->JINV * grad.Coordinates());
}


template<size_t dim>
Point<dim>
FiniteVolumeHelper<dim>::NormalOfFacet(size_t iFacet) const
{
    switch (element_dim_) {
        case 1:
        {
            Point<dim> normal;
            const size_t iNrNodes(eptr_->Nodes());
            for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
                const Point<dim> n(eptr_->N(iNode)->Coordinate());
                auto weights = eptr_->FV()->FacetNormalTransformationNodeWeights(iFacet, iNode);
                normal += weights.first * n;
            }
            return normal;
        }

        case 2:
        {
            Point<dim> tangent;
            Point<dim> bitangent;
            const size_t iNrNodes(eptr_->Nodes());
            for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
                const Point<dim> n(eptr_->N(iNode)->Coordinate());
                auto weights = eptr_->FV()->FacetNormalTransformationNodeWeights(iFacet, iNode);
                tangent += weights.first * n;
                bitangent += weights.second * n;
            }
            double64 length = exteriorProductLength(tangent, bitangent);
            tangent.NormalizeLengthTo(1.0);
            Point<dim> normal = bitangent - dotProduct(tangent,bitangent) * tangent;
            normal.NormalizeLengthTo(length);
            return normal;
        }

        case 3:
        {
            Point<dim> v0(0.0);
            Point<dim> v1(0.0);
            const size_t iNrNodes(eptr_->Nodes());
            for (size_t iNode = 0; iNode < iNrNodes; ++iNode) {
                auto xform_weights = eptr_->FV()->FacetNormalTransformationNodeWeights(iFacet, iNode);
                const Point<dim> n(eptr_->N(iNode)->Coordinate());
                v0 += xform_weights.first * n;
                v1 += xform_weights.second * n;
            }
            return crossProduct(v1, v0);
        }
    }
  
  // SKM fix
  return Point<dim>();
}


template<>
void
FiniteVolumeHelper<1u>::CalculateDN(const Point<1u>& p, std::vector<double64>* DN) const
{
    auto fe = eptr_->FE();
    switch (element_dim_) {
        case 1:
            fe->dNr(p[0], DN[0]);
            fe->Jacobian( DN[0] );
            
    }
}


template<>
void
FiniteVolumeHelper<2u>::CalculateDN(const Point<2u>& p, std::vector<double64>* DN) const
{
    auto fe = eptr_->FE();
    switch (element_dim_) {
        case 1:
            fe->dNr(p[0], DN[0]);
            fe->Jacobian( DN[0] );
            break;

        case 2:
            fe->dNr(p[0], p[1], DN[0]);
            fe->dNs(p[0], p[1], DN[1]);
            fe->Jacobian( DN[0], DN[1] );
            break;
    }
}


template<>
void
FiniteVolumeHelper<3u>::CalculateDN(const Point<3u>& p, std::vector<double64>* DN) const
{
    auto fe = eptr_->FE();
    switch (element_dim_) {
        case 3:
            fe->dNr(p[0], p[1], p[2], DN[0]);
            fe->dNs(p[0], p[1], p[2], DN[1]);
            fe->dNt(p[0], p[1], p[2], DN[2]);
            fe->Jacobian( DN[0], DN[1], DN[2] );
            break;
        case 2:
            fe->dNs(p[0], p[1], DN[1]);
            fe->Jacobian( DN[0], DN[1] );
            break;
        case 1:
            fe->dNr(p[0], DN[0]);
            fe->Jacobian( DN[0] );
            break;
    }
}

  
  template<>
  void
  FiniteVolumeHelper<1u>::CalculateN(const Point<1u>& p, std::vector<double64>& N) const
  {
    auto fe = eptr_->FE();
    switch (element_dim_) {
      case 1:
        fe->Nr( p[0], N );
    }
  }
  
  
  template<>
  void
  FiniteVolumeHelper<2u>::CalculateN(const Point<2u>& p, std::vector<double64>& N) const
  {
    auto fe = eptr_->FE();
    switch (element_dim_) {
      case 1:
        fe->Nr( p[0], N );
        break;
        
      case 2:
        fe->Nrs( p[0], p[1], N );
        break;
    }
  }
  
  
  template<>
  void
  FiniteVolumeHelper<3u>::CalculateN(const Point<3u>& p, std::vector<double64>& N) const
  {
    auto fe = eptr_->FE();
    switch (element_dim_) {
      case 3:
        fe->Nrst( p[0], p[1], p[2], N );
        break;
      case 2:
        fe->Nrs( p[0], p[1], N );
        break;
      case 1:
        fe->Nr( p[0], N );
        break;
    }
  }
  
  


template class FiniteVolumeHelper<1u>;
template class FiniteVolumeHelper<2u>;
template class FiniteVolumeHelper<3u>;


} // end csmp
