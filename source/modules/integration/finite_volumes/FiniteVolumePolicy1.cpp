#include "FiniteVolumePolicy.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

namespace csmp {

template<template<uint32_t> class CELL>
FiniteVolumePolicy<1U,CELL>::FiniteVolumePolicy( const csmp::FiniteVolumeStencil<1U>* fvptr )
 : fvptr_(fvptr)
 {
 }


template<template<uint32_t> class CELL>
void FiniteVolumePolicy<1U,CELL>::AssignFiniteVolume( const csmp::FiniteVolumeStencil<1U>* fvptr )
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( fvptr != nullptr );
    assert( e     != nullptr );
    assert( fvptr->Geometry() == parseFiniteElementDimension( e->FE_Type() ) );
    assert( fvptr->Sectors()  == e->FE()->CornerNodes() );
    fvptr_ = fvptr;
 }

// FV STENCIL INFO

template<template<uint32_t> class CELL>
const FiniteVolumeStencil<1U>* const FiniteVolumePolicy<1U,CELL>::FV() const
 {
    return fvptr_;
 }

/// Returns number of sectors, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<1U,CELL>::Sectors() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Sectors();
    return 0U;
  }

/// Returns number of facets, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<1U,CELL>::Facets() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Facets();
    return 0;
  }


/// Returns number of integration points per sector, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<1U,CELL>::IntegrationPointsPerSector() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerSector();
    return 0;
  }

/// Returns number of integration points per facet, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<1U,CELL>::IntegrationPointsPerFacet() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerFacet();
    return 0;
  }




template<template<uint32_t> class CELL>
void FiniteVolumePolicy<1U,CELL>::N_AtFacetIntegrationPoint( uint32_t iFacet, uint32_t ip ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( e != nullptr );
    e->FE()->Nr( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ), e->FE()->NRST );
 }


template<template<uint32_t> class CELL>
void FiniteVolumePolicy<1U,CELL>::N_AtSectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( e != nullptr );
    e->FE()->Nr( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ), e->FE()->NRST );
 }


// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS

template<template<uint32_t> class CELL>
void FiniteVolumePolicy<1U,CELL>::Local_dN_At( const Point<1U>& rst ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( e != nullptr );
    e->FE()->dNr( rst[0], e->FE()->DNR );
 }

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::dN_At( const Point<1U>& rst,
                                                 DenseMatrix<DM_MIN>& DN ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( e != nullptr );
    e->FE()->dNr( rst[0], e->FE()->DNR );
    e->FE()->Jacobian( e->FE()->DNR );
    const uint32_t nodes(e->Nodes());
    DN.Resize(1U,nodes);
    for ( auto i{0}; i<nodes; i++ ) DN(0U,i) = e->FE()->DNR[i];
    double detJ = e->FE()->JacobianInverse();
    DN = e->FE()->JINV * DN;
    return detJ;
 }







// PROPERTIES AT FACET AND SECTOR INTEGRATION POINTS

// only for scalars
template<template<uint32_t> class CELL>
double FiniteVolumePolicy<1U,CELL>::PropertyValueAtFacetIntegrationPoint(
                                                 uint32_t iFacet,
                                                 uint32_t ip,
                                                 const csmp::Index& prop_key ) const
{
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interppolating property to integration point
    double  sum(static_cast<double>(0.));

    assert( e != nullptr );
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i=0; i<nodes; i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

} // end PropertyValueAtFacetIntegrationPoint


// only for scalars
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::PropertyValueAtSectorIntegrationPoint(
                                                            uint32_t iSector,
                                                            uint32_t ip,
                                                            const csmp::Index& prop_key ) const
{
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iSector<fvptr_->Sectors());
    assert( ip<fvptr_->IntegrationPointsPerSector());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    assert( e != nullptr );
    N_AtSectorIntegrationPoint( iSector, ip );
    double  sum(0.);
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i=0; i < nodes; i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

 } // end PropertyValueAtVolumeIntegrationPoint


template<template<uint32_t> class CELL>
template<class Var>
void FiniteVolumePolicy<1U,CELL>::PropertyValueAtFacetIntegrationPoint(
                                                 const csmp::Index& prop_key,
                                                 uint32_t iFacet,
                                                 uint32_t ip,
                                                 Var& var ) const
{
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key, var );

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interppolating property to integration point
    assert( e != nullptr );
    Var temp;
    temp.Size( prop_key.dataDepth );
    var.Size( prop_key.dataDepth );
    var = 0;
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i=0; i<nodes; i++ )
      {
          e->N(i)->Read( prop_key, temp );
          var += temp * e->FE()->NRST[i];
      }

} // end PropertyValueAtFacetIntegrationPoint


template<template<uint32_t> class CELL>
template<class Var>
void  FiniteVolumePolicy<1U,CELL>::PropertyValueAtSectorIntegrationPoint(
                                                            const csmp::Index& prop_key,
                                                            uint32_t iSector,
                                                            uint32_t ip,
                                                            Var& var ) const
{
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iSector<fvptr_->Sectors());
    assert( ip<fvptr_->IntegrationPointsPerSector());
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key, var );

    N_AtSectorIntegrationPoint( iSector, ip );
    Var temp;
    temp.Size( prop_key.dataDepth );
    var.Size( prop_key.dataDepth );
    var = 0;
    assert( e != nullptr );
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i=0; i < nodes; i++ )
      {
          e->N(i)->Read( prop_key, temp );
          var += temp * e->FE()->NRST[i];
      }

 } // end PropertyValueAtVolumeIntegrationPoint





// FACET AND SECTOR INTEGRALS

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::FacetIntegral( uint32_t iFacet,
                                                    const csmp::Index& prop_key ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place != ELEMENT ) return e->Read( prop_key );

    return PropertyValueAtFacetIntegrationPoint( iFacet, 0U, prop_key );
 }



template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::SectorIntegral( uint32_t iSector,
                                                     const csmp::Index& prop_key ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    assert( e != nullptr );
    e->CoordinateMatrix();
    double fIntegral(static_cast<double>(0.));

    for ( uint32_t j=0; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
        // 1. take integration point location from FV stencil
        fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );

        // 2. Compute the Jacobian from DNR... (methods of FiniteElement
        e->FE()->JacobianAt( e->FE()->NRST );

        // 3. Compute the determinant of the jacobian
        const double detJ(e->FE()->JacobianDeterminant());

        // 4. Map your local integral to global space by multiplication with detJ
        fIntegral += PropertyValueAtSectorIntegrationPoint( iSector, j, prop_key ) *
                     fvptr_->SectorIntegrationWeight( iSector, j ) * detJ;
     } // Loop over integration points

  return fIntegral;

} // end VolumeIntegral




// PROJECTION ON FACET NORMAL

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::ProjectionOnFacetNormal( uint32_t /*iFacet*/,
                                                              const csmp::Index& prop_key ) const
{
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    VectorVariable<1U>  vc;

    if ( prop_key.place == ELEMENT ) {
         e->Read( prop_key, vc );
         return vc.Length();
      }
    // the facet normal in physical space is 1 and aligned with the element
    // here we just assume that this normal originates from the barycenter
    assert( e != nullptr );
    e->PropertyValueAtBaryCenter( prop_key, vc );
    return vc.Length();

}


template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::ProjectionOnFacetNormal( uint32_t /*iFacet*/,
                                                              const VectorVariable<1U>& vc ) const
 {
    return vc.Length();
 }







// GEOMETRY


// SECTOR VOLUME

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::SectorVolume( uint32_t iSector ) const
 {
    const CELL<1U>* e( static_cast<const CELL<1U>*>(this) );
    assert( iSector < fvptr_->Sectors());

    assert( e != nullptr );
    e->CoordinateMatrix();
    double fVolume(static_cast<double>(0.));
   
    const uint32_t ipoints(fvptr_->IntegrationPointsPerSector());
    for ( uint32_t j=0; j<ipoints; j++ ) {
         fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
         e->FE()->JacobianAt(e->FE()->NRST);
         const double detJ(e->FE()->JacobianDeterminant());
         fVolume += fvptr_->SectorIntegrationWeight(iSector, j)*detJ;
      }
    return fVolume;

} // end SectorVolume


// FACET AREA AND NORMAL

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::FacetArea( uint32_t ) const
{
    return 1.;
}


template<template<uint32_t> class CELL>
Point<1U> FiniteVolumePolicy<1U,CELL>::FacetNormal( uint32_t ) const
{
    //in 1D we only have line elements (Isoparametric Linear Bars)
    return Point<1U>(1.);
}

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::FacetAreaMapped( uint32_t ) const
{
    return 1.;
}


template<template<uint32_t> class CELL>
Point<1U> FiniteVolumePolicy<1U,CELL>::FacetNormalMapped( uint32_t ) const
{
    //in 1D we only have line elements (Isoparametric Linear Bars)
    return Point<1U>(1.);
}
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<1U,CELL>::ParametricFacetArea( uint32_t ) const
 {
    return 1.;
 }

template<template<uint32_t> class CELL>
Point<1U>  FiniteVolumePolicy<1U,CELL>::ParametricFacetNormal( uint32_t ) const
{
   return Point<1U>(1.);
}

template class FiniteVolumePolicy<1U,Element>;
template class FiniteVolumePolicy<1U,Face>;
template class FiniteVolumePolicy<1U,InterFace>;

} // end namespace csmp
