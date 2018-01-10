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
   
    const csmp::INDEX<SCALAR,ELEMENT> phi_key(model.Database().StorageKey("porosity"));
    const csmp::INDEX<SCALAR,ELEMENT>  thi_key(model.Database().StorageKey("thickness"));
    const csmp::INDEX<VECTOR,ELEMENT>  vt_key(model.Database().StorageKey("velocity"));
    const csmp::INDEX<SCALAR,ELEMENT>  k_key(model.Database().StorageKey("permeability"));
    const csmp::INDEX<SCALAR,NODE>  pf_key(model.Database().StorageKey("fluid pressure"));

    const csmp::INDEX<SCALAR,NODE>  pv_key(model.Database().StorageKey("FV pore volume"));
    // const csmp::Index spv_key = model.Database().StorageKey("sector pore volume");
    const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>  ff_key(model.Database().StorageKey("facet flux"));
    const csmp::INDEX<SCALAR,NODE>  fb_key(model.Database().StorageKey("flux balance"));
   
     Point<dim> vD;
   
    // 0. zeroing sector pore volumes for accumulation in element loop
    // ---------------------------------------------------------------
    gref.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );
   
   FiniteElementHelper<dim> fe;
  
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref.ElementsBegin(); it!=it_end; ++it )
      {
         fe.FiniteElement(*it);
         const size_t sectors((*it)->Sectors());
          const size_t facets((*it)->Facets());

         // element-based total velocity
          if ( initialize_flux ) {
              const double64 K((*it)->Read(k_key));
              vD = -K * fe.ReadGradientAtBarycenter(pf_key);
          }

         // 1. computing sector pore volumes
         // --------------------------------
         // (scaled by the cell thickness attribute=1 for volumetric elements)
         ScalarVariable var_phi, var_thi;
         fe.ReadAtBarycenter( phi_key, var_phi );
         fe.ReadAtBarycenter( thi_key, var_thi );

         const double64 phi = var_phi() * var_thi();
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
                  Point<dim> facetNormal(fe.NormalOfFacet(i));
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


enum InterpolatorType {
    READ_ELMT,
    READ_NODE,
    NODE_TO_ELMT,
    NODE_TO_EIP,
    NODE_TO_FIP,
    NODE_TO_SIP,
    INTERPOLATOR_COUNT
};




  template<InterpolatorType interp>
  struct PropertyInterpolation
  {
    template<size_t dim>
    void Recalculate( const Element<dim>* eptr, size_t element_dim );

    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t idx2, typename VariableTypeTraits<dim,ty>::VariableType& var );
  };

  template<size_t dim>
  void CalculateN(const Element<dim>* eptr, size_t element_dim, const Point<dim>& p, double64* coeff);
  
  template<>
  void CalculateN<1u>(const Element<1u>* eptr, size_t element_dim, const Point<1u>& p, double64* coeff)
  {
    auto fe = eptr->FE();
    switch (element_dim) {
      case 1:
        fe->Nr( p[0], coeff );
    }
  }

  template<>
  void CalculateN<2u>(const Element<2u>* eptr, size_t element_dim, const Point<2u>& p, double64* coeff)
  {
    auto fe = eptr->FE();
    switch (element_dim) {
      case 1:
        fe->Nr( p[0], coeff );
        break;
        
      case 2:
        fe->Nrs( p[0], p[1], coeff );
        break;
    }
  }

  template<>
  void CalculateN<3u>(const Element<3u>* eptr, size_t element_dim, const Point<3u>& p, double64* coeff)
  {
    auto fe = eptr->FE();
    switch (element_dim) {
      case 3:
        fe->Nrst( p[0], p[1], p[2], coeff );
        break;
      case 2:
        fe->Nrs( p[0], p[1], coeff );
        break;
      case 1:
        fe->Nr( p[0], coeff );
        break;
    }
  }
  

  template<>
  struct PropertyInterpolation<READ_ELMT>
  {
    const InterpolatorType type_ = READ_ELMT;
    
    template<size_t dim>
    void Recalculate( const Element<dim>* eptr, size_t element_dim_ )
    {
    }
    
    InterpolatorType Type() const { return type_; }
    
    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t, size_t, typename VariableTypeTraits<dim, ty>::VariableType& var )
    {
      eptr->Read( prop, var );
    }
  };
  

  
  
  template<>
  struct PropertyInterpolation<READ_NODE>
  {
    const InterpolatorType type_ = READ_NODE;
    
    template<size_t dim>
    void Recalculate( const Element<dim>* eptr, size_t element_dim_ )
    {
    }
    
    InterpolatorType Type() const { return type_; }
    
    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t idx1, size_t, typename VariableTypeTraits<dim, ty>::VariableType& var )
    {
      eptr->N(idx1)->Read( prop, var );
    }
  };
  

  template<>
  struct PropertyInterpolation<NODE_TO_ELMT>
  {
    const InterpolatorType type_ = NODE_TO_ELMT;

    std::vector<double64> coeff_;
    
    PropertyInterpolation()
    {
      coeff_.resize(DM_MAX);
    }

    InterpolatorType Type() const { return type_; }

    template<size_t dim>
    void Recalculate( const Element<dim>* eptr, size_t element_dim )
    {
      const size_t num_nodes(eptr->Nodes());
      if (coeff_.size() < num_nodes) {
        coeff_.resize(num_nodes);
      }
      CalculateN(eptr, element_dim, eptr->FV()->Barycenter(), &coeff_[0]);
    }

    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t, size_t, typename VariableTypeTraits<dim, ty>::VariableType& var )
    {
      const size_t num_nodes(eptr->Nodes());
      for (size_t i = 0; i < var.Components(); ++i) {
        var.Component( i, 0. );
      }
      for ( size_t i=0; i<num_nodes; i++ ) {
        var += eptr->N(i)->Read( prop ) * coeff_[i];
      }
    }
  };

  template<>
  struct PropertyInterpolation<NODE_TO_FIP>
  {
    const InterpolatorType type_ = NODE_TO_FIP;
    
    size_t facets_, ips_per_facet_, num_nodes_;
    std::vector<double64> coeff_;

    PropertyInterpolation()
    {
      coeff_.resize(DM_MAX);
    }
    
    InterpolatorType Type() const { return type_; }
    
    template<size_t dim>
    void Recalculate( const Element<dim>* eptr, size_t element_dim )
    {
      auto fv = eptr->FV();
      num_nodes_ = eptr->Nodes();
      facets_ = fv->Facets();
      ips_per_facet_ = fv->IntegrationPointsPerFacet();
      
      const size_t coeff_size = facets_ * ips_per_facet_ * num_nodes_;

      if (coeff_.size() < coeff_size) {
        coeff_.resize(coeff_size);
      }
      size_t offset = 0;
      for (size_t iFacet = 0; iFacet < facets_; ++iFacet) {
        for (size_t iFip = 0; iFip < ips_per_facet_; ++iFip) {
          CalculateN(eptr, element_dim, fv->FacetIntegrationPoint(iFacet, iFip), &coeff_[offset]);
          offset += num_nodes_;
        }
      }
    }
    
    template<size_t dim,VARIABLE_TYPE ty>
    void
    Interpolate( const Index& prop, const Element<dim>* eptr, size_t facet, size_t fip, typename VariableTypeTraits<dim, ty>::VariableType& var )
    {
      auto fv = eptr->FV();
      const size_t num_nodes = eptr->Nodes();
      const size_t offset = (facet * fv->IntegrationPointsPerFacet() + fip) * num_nodes;
      for (size_t i = 0; i < var.Size(); ++i) {
        var.Component( i, 0. );
      }
      for ( size_t i=0; i<num_nodes; i++ ) {
        var += eptr->N(i)->Read( prop ) * coeff_[offset+i];
      }
    }
  };

struct PropertyInterpolators
{
    std::bitset<INTERPOLATOR_COUNT> interp_valid_;
  
    PropertyInterpolation<READ_ELMT> read_elmt_;
    PropertyInterpolation<READ_NODE> read_node_;
    PropertyInterpolation<NODE_TO_ELMT> node_elmt_;
    // PropertyInterpolation<NODE_TO_EIP> node_eip_;
    PropertyInterpolation<NODE_TO_FIP> node_fip_;
    // PropertyInterpolation<NODE_TO_SIP> node_sip_;
};


  template<PLACEMENT from,PLACEMENT to>
  struct InterpolatorDispatch
  {
  };
  
  template<>
  struct InterpolatorDispatch<ELEMENT,ELEMENT>
  {
    typedef PropertyInterpolation<READ_ELMT> Interpolator;
    
    static Interpolator& GetInterpolator(PropertyInterpolators& interps) {
      return interps.read_elmt_;
    }
  };
  
  template<>
  struct InterpolatorDispatch<NODE,NODE>
  {
    typedef PropertyInterpolation<READ_NODE> Interpolator;
    
    static Interpolator& GetInterpolator(PropertyInterpolators& interps) {
      return interps.read_node_;
    }
  };
  
  template<>
  struct InterpolatorDispatch<NODE,ELEMENT>
  {
    typedef PropertyInterpolation<NODE_TO_ELMT> Interpolator;
    
    static Interpolator& GetInterpolator(PropertyInterpolators& interps) {
      return interps.node_elmt_;
    }
  };

  template<>
  struct InterpolatorDispatch<NODE,FACET_INTEGRATION_POINT>
  {
    typedef PropertyInterpolation<NODE_TO_FIP> Interpolator;
    
    static Interpolator& GetInterpolator(PropertyInterpolators& interps) {
      return interps.node_fip_;
    }
  };

template<size_t dim>
struct FiniteElementHelper<dim>::Impl : public PropertyInterpolators
{
  Element<dim>* eptr_;
  size_t element_dim_;
  size_t num_nodes_;
  std::vector<double64> DN_bctr_[dim];
  
  void FiniteElement( Element<dim>* eptr )
  {
    eptr_ = eptr;
    if (eptr_->IsLineElement()) {
      element_dim_ = 1;
    }
    else if (eptr_->IsSurfaceElement()) {
      element_dim_ = 2;
    }
    else if (eptr_->IsVolumeElement()) {
      element_dim_ = 3;
    }

    auto fe = eptr_->FE();
    num_nodes_ = fe->Nodes();

    eptr->CoordinateMatrix();
    interp_valid_.reset();
  }
};





template<size_t dim>
FiniteElementHelper<dim>::FiniteElementHelper()
  : pimpl_(new FiniteElementHelper::Impl())
{
}

  
template<size_t dim>
void FiniteElementHelper<dim>::FiniteElement( Element<dim>* eptr )
  {
      pimpl_->FiniteElement(eptr);
  }


template<size_t dim>
FiniteElementHelper<dim>::~FiniteElementHelper()
{
}

  
template<size_t dim>
Element<dim>* FiniteElementHelper<dim>::FiniteElement( )
  {
      return pimpl_->eptr_;
  }

  
  
template<size_t dim> template<VARIABLE_TYPE ty,PLACEMENT pl>
void FiniteElementHelper<dim>::ReadAtBarycenter( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var )
{
  auto& interpolator = InterpolatorDispatch<pl,ELEMENT>::GetInterpolator(*pimpl_);
  if (!pimpl_->interp_valid_[interpolator.type_]) {
    pimpl_->interp_valid_[interpolator.type_] = true;
    interpolator.Recalculate(pimpl_->eptr_, pimpl_->element_dim_);
  }

  interpolator.template Interpolate<dim,ty>( prop, pimpl_->eptr_, 0, 0, var );
}

template<size_t dim>
template<VARIABLE_TYPE ty,PLACEMENT pl>
void FiniteElementHelper<dim>::ReadAtNode( const csmp::INDEX<ty,pl>& prop, size_t n, typename VariableTypeTraits<dim,ty>::VariableType& var )
{
  auto& interpolator = InterpolatorDispatch<pl,NODE>::GetInterpolator(*pimpl_);
  if (!pimpl_->interp_valid_[interpolator.type_]) {
    pimpl_->interp_valid_[interpolator.type_] = true;
    interpolator.Recalculate(pimpl_->eptr_, pimpl_->element_dim_);
  }
  
  interpolator.template Interpolate<dim,ty>( prop, pimpl_->eptr_, n, 0, var );
}

  
  template<size_t dim>
  template<VARIABLE_TYPE ty,PLACEMENT pl>
  void FiniteElementHelper<dim>::ReadAtFacetIntegrationPoint( const csmp::INDEX<ty,pl>& prop, size_t facet, size_t fip, typename VariableTypeTraits<dim,ty>::VariableType& var )
  {
    auto& interpolator = InterpolatorDispatch<pl,FACET_INTEGRATION_POINT>::GetInterpolator(*pimpl_);
    if (!pimpl_->interp_valid_[interpolator.type_]) {
      pimpl_->interp_valid_[interpolator.type_] = true;
      interpolator.Recalculate(pimpl_->eptr_, pimpl_->element_dim_);
    }
    
    interpolator.template Interpolate<dim,ty>( prop, pimpl_->eptr_, facet, fip, var );
  }
  
  


template<size_t dim>
Point<dim>
FiniteElementHelper<dim>::ReadGradientAtBarycenter( const csmp::INDEX<SCALAR,NODE>& prop )
{
  const size_t num_nodes = pimpl_->num_nodes_;
  const size_t element_dim = pimpl_->element_dim_;
  auto& DN_bctr = pimpl_->DN_bctr_;
  auto eptr = pimpl_->eptr_;
    for (unsigned i = 0; i < element_dim; ++i) {
      if (DN_bctr[i].size() < num_nodes) {
        DN_bctr[i].resize(num_nodes);
      }
    }
    CalculateDN(eptr->FV()->Barycenter(), DN_bctr);

    Point<dim> grad(0.);
    for ( size_t i=0U; i<num_nodes; ++i ) {
        const double64 value_at_node(eptr->N(i)->Read(prop));
        for (size_t j = 0; j < element_dim; ++j) {
            grad[j] += DN_bctr[j][i] * value_at_node;
        }
    }
    eptr->FE()->JacobianInverse();

    return Point<dim>(eptr->FE()->JINV * grad.Coordinates());
}


template<size_t dim>
Point<dim>
FiniteElementHelper<dim>::NormalOfFacet(size_t iFacet)
{
  auto eptr = pimpl_->eptr_;
    switch (pimpl_->element_dim_) {
        case 1:
        {
            Point<dim> normal;
            const size_t iNrNodes(eptr->Nodes());
            for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
                const Point<dim> n(eptr->N(iNode)->Coordinate());
                auto weights = eptr->FV()->FacetNormalTransformationNodeWeights(iFacet, iNode);
                normal += weights.first * n;
            }
            return normal;
        }

        case 2:
        {
            Point<dim> tangent;
            Point<dim> bitangent;
            const size_t iNrNodes(eptr->Nodes());
            for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
                const Point<dim> n(eptr->N(iNode)->Coordinate());
                auto weights = eptr->FV()->FacetNormalTransformationNodeWeights(iFacet, iNode);
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
            const size_t iNrNodes(eptr->Nodes());
            for (size_t iNode = 0; iNode < iNrNodes; ++iNode) {
                auto xform_weights = eptr->FV()->FacetNormalTransformationNodeWeights(iFacet, iNode);
                const Point<dim> n(eptr->N(iNode)->Coordinate());
                v0 += xform_weights.first * n;
                v1 += xform_weights.second * n;
            }
            return crossProduct(v1, v0);
        }
    }
  throw csmp::Exception(ERROR, "FiniteElementHelper::NormalOfFacet", "Element dimension must be 1, 2, or 3");
}


template<>
void
FiniteElementHelper<1u>::CalculateDN(const Point<1u>& p, std::vector<double64>* DN)
{
    auto fe = pimpl_->eptr_->FE();
    switch (pimpl_->element_dim_) {
        case 1:
            fe->dNr(p[0], DN[0]);
            fe->Jacobian( DN[0] );
            
    }
}


template<>
void
FiniteElementHelper<2u>::CalculateDN(const Point<2u>& p, std::vector<double64>* DN)
{
  auto fe = pimpl_->eptr_->FE();
  switch (pimpl_->element_dim_) {
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
FiniteElementHelper<3u>::CalculateDN(const Point<3u>& p, std::vector<double64>* DN)
{
    auto fe = pimpl_->eptr_->FE();
    switch (pimpl_->element_dim_) {
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


#define READ_AT_BARYCENTER(from) \
  template void FiniteElementHelper<1u>::ReadAtBarycenter<SCALAR,from>(INDEX<SCALAR,from> const&, VariableTypeTraits<1u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtBarycenter<SCALAR,from>(INDEX<SCALAR,from> const&, VariableTypeTraits<2u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtBarycenter<SCALAR,from>(INDEX<SCALAR,from> const&, VariableTypeTraits<3u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtBarycenter<VECTOR,from>(INDEX<VECTOR,from> const&, VariableTypeTraits<1u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtBarycenter<VECTOR,from>(INDEX<VECTOR,from> const&, VariableTypeTraits<2u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtBarycenter<VECTOR,from>(INDEX<VECTOR,from> const&, VariableTypeTraits<3u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtBarycenter<TENSOR,from>(INDEX<TENSOR,from> const&, VariableTypeTraits<1u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtBarycenter<TENSOR,from>(INDEX<TENSOR,from> const&, VariableTypeTraits<2u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtBarycenter<TENSOR,from>(INDEX<TENSOR,from> const&, VariableTypeTraits<3u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtBarycenter<ARRAY,from>(INDEX<ARRAY,from> const&, VariableTypeTraits<1u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtBarycenter<ARRAY,from>(INDEX<ARRAY,from> const&, VariableTypeTraits<2u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtBarycenter<ARRAY,from>(INDEX<ARRAY,from> const&, VariableTypeTraits<3u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtBarycenter<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, VariableTypeTraits<1u,FLAGGEDARRAY>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtBarycenter<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, VariableTypeTraits<2u,FLAGGEDARRAY>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtBarycenter<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, VariableTypeTraits<3u,FLAGGEDARRAY>::VariableType&);

#define READ_AT_NODE(from) \
  template void FiniteElementHelper<1u>::ReadAtNode<SCALAR,from>(INDEX<SCALAR,from> const&, size_t, VariableTypeTraits<1u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtNode<SCALAR,from>(INDEX<SCALAR,from> const&, size_t, VariableTypeTraits<2u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtNode<SCALAR,from>(INDEX<SCALAR,from> const&, size_t, VariableTypeTraits<3u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtNode<VECTOR,from>(INDEX<VECTOR,from> const&, size_t, VariableTypeTraits<1u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtNode<VECTOR,from>(INDEX<VECTOR,from> const&, size_t, VariableTypeTraits<2u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtNode<VECTOR,from>(INDEX<VECTOR,from> const&, size_t, VariableTypeTraits<3u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtNode<TENSOR,from>(INDEX<TENSOR,from> const&, size_t, VariableTypeTraits<1u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtNode<TENSOR,from>(INDEX<TENSOR,from> const&, size_t, VariableTypeTraits<2u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtNode<TENSOR,from>(INDEX<TENSOR,from> const&, size_t, VariableTypeTraits<3u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtNode<ARRAY,from>(INDEX<ARRAY,from> const&, size_t, VariableTypeTraits<1u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtNode<ARRAY,from>(INDEX<ARRAY,from> const&, size_t, VariableTypeTraits<2u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtNode<ARRAY,from>(INDEX<ARRAY,from> const&, size_t, VariableTypeTraits<3u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtNode<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, size_t, VariableTypeTraits<1u,FLAGGEDARRAY>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtNode<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, size_t, VariableTypeTraits<2u,FLAGGEDARRAY>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtNode<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, size_t, VariableTypeTraits<3u,FLAGGEDARRAY>::VariableType&);

#define READ_AT_FACET_INTEGRATION_POINT(from) \
  template void FiniteElementHelper<1u>::ReadAtFacetIntegrationPoint<SCALAR,from>(INDEX<SCALAR,from> const&, size_t, size_t, VariableTypeTraits<1u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtFacetIntegrationPoint<SCALAR,from>(INDEX<SCALAR,from> const&, size_t, size_t, VariableTypeTraits<2u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtFacetIntegrationPoint<SCALAR,from>(INDEX<SCALAR,from> const&, size_t, size_t, VariableTypeTraits<3u,SCALAR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtFacetIntegrationPoint<VECTOR,from>(INDEX<VECTOR,from> const&, size_t, size_t, VariableTypeTraits<1u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtFacetIntegrationPoint<VECTOR,from>(INDEX<VECTOR,from> const&, size_t, size_t, VariableTypeTraits<2u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtFacetIntegrationPoint<VECTOR,from>(INDEX<VECTOR,from> const&, size_t, size_t, VariableTypeTraits<3u,VECTOR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtFacetIntegrationPoint<TENSOR,from>(INDEX<TENSOR,from> const&, size_t, size_t, VariableTypeTraits<1u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtFacetIntegrationPoint<TENSOR,from>(INDEX<TENSOR,from> const&, size_t, size_t, VariableTypeTraits<2u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtFacetIntegrationPoint<TENSOR,from>(INDEX<TENSOR,from> const&, size_t, size_t, VariableTypeTraits<3u,TENSOR>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtFacetIntegrationPoint<ARRAY,from>(INDEX<ARRAY,from> const&, size_t, size_t, VariableTypeTraits<1u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtFacetIntegrationPoint<ARRAY,from>(INDEX<ARRAY,from> const&, size_t, size_t, VariableTypeTraits<2u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtFacetIntegrationPoint<ARRAY,from>(INDEX<ARRAY,from> const&, size_t, size_t, VariableTypeTraits<3u,ARRAY>::VariableType&); \
  template void FiniteElementHelper<1u>::ReadAtFacetIntegrationPoint<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, size_t, size_t, VariableTypeTraits<1u,FLAGGEDARRAY>::VariableType&); \
  template void FiniteElementHelper<2u>::ReadAtFacetIntegrationPoint<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, size_t, size_t, VariableTypeTraits<2u,FLAGGEDARRAY>::VariableType&); \
  template void FiniteElementHelper<3u>::ReadAtFacetIntegrationPoint<FLAGGEDARRAY,from>(INDEX<FLAGGEDARRAY,from> const&, size_t, size_t, VariableTypeTraits<3u,FLAGGEDARRAY>::VariableType&);

template class FiniteElementHelper<1u>;
template class FiniteElementHelper<2u>;
template class FiniteElementHelper<3u>;

  READ_AT_BARYCENTER(NODE)
  READ_AT_BARYCENTER(ELEMENT)

  READ_AT_NODE(NODE)

  READ_AT_FACET_INTEGRATION_POINT(NODE)

} // end csmp
