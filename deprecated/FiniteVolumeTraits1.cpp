#include "FiniteVolumeTraits.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

namespace csmp {

template<template<size_t> class SIMPLEX>
FiniteVolumeTraits<1U,SIMPLEX>::FiniteVolumeTraits()
 {
 }


// FV STENCIL INFO

/// Returns number of sectors, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumeTraits<1U,SIMPLEX>::Sectors() const
  {
    const SIMPLEX<1U>* eptr( static_cast<const SIMPLEX<1U>*>(this) );
    assert( eptr != NULL );
    if( eptr->FV_Stencil() )
        return eptr->FV_Stencil()->Sectors();
    return 0;
  }

/// Returns number of facets, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumeTraits<1U,SIMPLEX>::Facets() const
  {
    const SIMPLEX<1U>* eptr( static_cast<const SIMPLEX<1U>*>(this) );
    assert( eptr != NULL );
    if( eptr->FV_Stencil() )
        return eptr->FV_Stencil()->Facets();
    return 0;
  }


/// Returns number of integration points per sector, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumeTraits<1U,SIMPLEX>::IntegrationPointsPerSector() const
  {
    const SIMPLEX<1U>* eptr( static_cast<const SIMPLEX<1U>*>(this) );
    assert( eptr != NULL );
    if( eptr->FV_Stencil() )
        return eptr->FV_Stencil()->IntegrationPointsPerSector();
    return 0;
  }

/// Returns number of integration points per facet, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumeTraits<1U,SIMPLEX>::IntegrationPointsPerFacet() const
  {
    const SIMPLEX<1U>* eptr( static_cast<const SIMPLEX<1U>*>(this) );
    assert( eptr != NULL );
    if( eptr->FV_Stencil() )
        return eptr->FV_Stencil()->IntegrationPointsPerFacet();
    return 0;
  }



// MAPPING BETWEEN LOCAL AND GLOBAL COORDINATES

template<template<size_t> class SIMPLEX>
Point<1U> FiniteVolumeTraits<1U,SIMPLEX>::RstToXYZ( const Point<1U>& rst ) const
{
   const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
   N_At( rst );
   double64  sum(0.);
   for ( size_t i=0; i<e->Nodes(); i++ )
     sum += e->FE()->NRST[i] * e->N(i)->x();

   return Point<1U>(sum);
}




// SHAPE FUNCTIONS AT DIFFERENT POINTS

template<template<size_t> class SIMPLEX>
void FiniteVolumeTraits<1U,SIMPLEX>::N_At( const Point<1U>& rst ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    e->FE()->Nr( rst[0], e->FE()->NRST );
 }


template<template<size_t> class SIMPLEX>
void FiniteVolumeTraits<1U,SIMPLEX>::N_At( const Point<1U>& rst,
                                           std::vector<double64>& IPOL ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    e->FE()->Nr( rst[0], IPOL );
 }

template<template<size_t> class SIMPLEX>
void FiniteVolumeTraits<1U,SIMPLEX>::N_AtFacetIntegrationPoint( size_t iFacet, size_t ip ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    e->FE()->Nr( e->FV_Stencil()->FacetIntegrationPoint( iFacet, ip, 0U ), e->FE()->NRST );
 }


template<template<size_t> class SIMPLEX>
void FiniteVolumeTraits<1U,SIMPLEX>::N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    e->FE()->Nr( e->FV_Stencil()->SectorIntegrationPoint( iSector, ip, 0U ), e->FE()->NRST );
 }


// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS

template<template<size_t> class SIMPLEX>
void FiniteVolumeTraits<1U,SIMPLEX>::Local_dN_At( const Point<1U>& rst ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    e->FE()->dNr( rst[0], e->FE()->DNR );
 }

template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::dN_At( const Point<1U>& rst,
                                             DenseMatrix<DM_MIN>& DN ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    e->FE()->dNr( rst[0], e->FE()->DNR );
    e->FE()->Jacobian( e->FE()->DNR );
    DN.Resize(1U,e->Nodes());
    for ( size_t i=0U; i<e->Nodes(); i++ ) DN(0U,i) = e->FE()->DNR[i];
    double64 detJ = e->FE()->JacobianInverse();
    DN = e->FE()->JINV * DN;
    return detJ;
 }







// PROPERTIES AT FACET AND SECTOR INTEGRATION POINTS

// only for scalars
template<template<size_t> class SIMPLEX>
double64 FiniteVolumeTraits<1U,SIMPLEX>::PropertyValueAtFacetIntegrationPoint(
                                                 size_t iFacet,
                                                 size_t ip,
                                                 const csmp::Index& prop_key ) const
{
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iFacet < e->FV_Stencil()->Facets());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interppolating property to integration point
    double64  sum(static_cast<double64>(0.));

    for ( size_t i=0; i<e->Nodes(); i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

} // end PropertyValueAtFacetIntegrationPoint

// only for scalars
template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::PropertyValueAtSectorIntegrationPoint(
                                                            size_t iSector,
                                                            size_t ip,
                                                            const csmp::Index& prop_key ) const
{
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iSector<e->FV_Stencil()->Sectors());
    assert( ip<e->FV_Stencil()->IntegrationPointsPerSector());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    N_AtSectorIntegrationPoint( iSector, ip );
    double64  sum(static_cast<double64>(0.));
    for ( size_t i=0; i < e->Nodes(); i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

 } // end PropertyValueAtVolumeIntegrationPoint


template<template<size_t> class SIMPLEX>
template<class Var>
void FiniteVolumeTraits<1U,SIMPLEX>::PropertyValueAtFacetIntegrationPoint(
                                                 const csmp::Index& prop_key,
                                                 size_t iFacet,
                                                 size_t ip,
                                                 Var& var ) const
{
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iFacet < e->FV_Stencil()->Facets());
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key, var );

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interppolating property to integration point
    Var temp;
    temp.Size( prop_key.dataDepth );
    var.Size( prop_key.dataDepth );
    var = 0;
    for ( size_t i=0; i<e->Nodes(); i++ )
      {
          e->N(i)->Read( prop_key, temp );
          var += temp * e->FE()->NRST[i];
      }

} // end PropertyValueAtFacetIntegrationPoint

template<template<size_t> class SIMPLEX>
template<class Var>
void  FiniteVolumeTraits<1U,SIMPLEX>::PropertyValueAtSectorIntegrationPoint(
                                                            const csmp::Index& prop_key,
                                                            size_t iSector,
                                                            size_t ip,
                                                            Var& var ) const
{
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iSector<e->FV_Stencil()->Sectors());
    assert( ip<e->FV_Stencil()->IntegrationPointsPerSector());
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key, var );

    N_AtSectorIntegrationPoint( iSector, ip );
    Var temp;
    temp.Size( prop_key.dataDepth );
    var.Size( prop_key.dataDepth );
    var = 0;
    for ( size_t i=0; i < e->Nodes(); i++ )
      {
          e->N(i)->Read( prop_key, temp );
          var += temp * e->FE()->NRST[i];
      }

 } // end PropertyValueAtVolumeIntegrationPoint





// FACET AND SECTOR INTEGRALS

template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::FacetIntegral( size_t iFacet,
                                                         const csmp::Index& prop_key ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iFacet < e->FV_Stencil()->Facets());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place != ELEMENT ) return e->Read( prop_key );

    return PropertyValueAtFacetIntegrationPoint( iFacet, 0U, prop_key );
 }



template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::SectorIntegral( size_t iSector,
                                                          const csmp::Index& prop_key ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iSector < e->FV_Stencil()->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    e->CoordinateMatrix();
    double64 fIntegral(static_cast<double64>(0.));

    for ( size_t j=0; j<e->FV_Stencil()->IntegrationPointsPerSector(); j++ ) {
        // 1. take integration point location from FV stencil
        e->FV_Stencil()->SectorIntegrationPoint( iSector, j, e->FE()->NRST );

        // 2. Compute the Jacobian from DNR... (methods of FiniteElement
        e->FE()->JacobianAt( e->FE()->NRST );

        // 3. Compute the determinant of the jacobian
        const double64 detJ(e->FE()->JacobianDeterminant());

        // 4. Map your local integral to global space by multiplication with detJ
        fIntegral += PropertyValueAtSectorIntegrationPoint( iSector, j, prop_key ) *
                     e->FV_Stencil()->SectorIntegrationWeight( iSector, j ) * detJ;
     } // Loop over integration points

  return fIntegral;

} // end VolumeIntegral




// PROJECTION ON FACET NORMAL

template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::ProjectionOnFacetNormal( size_t /*iFacet*/,
                                                                   const csmp::Index& prop_key ) const
{
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    VectorVariable<1U>  vc;

    if ( prop_key.place == ELEMENT ) {
         e->Read( prop_key, vc );
         return vc.Length();
      }
    // the facet normal in physical space is 1 and aligned with the element
    // here we just assume that this normal originates from the barycenter
    e->PropertyValueAtBaryCenter( prop_key, vc );
    return vc.Length();

}


template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::ProjectionOnFacetNormal( size_t /*iFacet*/,
                                                                   const VectorVariable<1U>& vc ) const
 {
    return vc.Length();
 }







// GEOMETRY


// SECTOR VOLUME

template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::SectorVolume( size_t iSector ) const
 {
    const SIMPLEX<1U>* e( static_cast<const SIMPLEX<1U>*>(this) );
    assert( iSector < e->FV_Stencil()->Sectors());

    e->CoordinateMatrix();
    double64 fVolume(static_cast<double64>(0.));

    for ( size_t j=0; j<e->FV_Stencil()->IntegrationPointsPerSector(); j++ ) {
         e->FV_Stencil()->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
         e->FE()->JacobianAt(e->FE()->NRST);
         const double64 detJ(e->FE()->JacobianDeterminant());
         fVolume += e->FV_Stencil()->SectorIntegrationWeight(iSector, j)*detJ;
      }
    return fVolume;

} // end SectorVolume


// FACET AREA AND NORMAL

template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::FacetArea( size_t ) const
{
    return 1.;
}


template<template<size_t> class SIMPLEX>
Point<1U> FiniteVolumeTraits<1U,SIMPLEX>::FacetNormal( size_t ) const
{
    //in 1D we only have line elements (Isoparametric Linear Bars)
    return Point<1U>(1.);
}

template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::FacetAreaMapped( size_t ) const
{
    return 1.;
}


template<template<size_t> class SIMPLEX>
Point<1U> FiniteVolumeTraits<1U,SIMPLEX>::FacetNormalMapped( size_t ) const
{
    //in 1D we only have line elements (Isoparametric Linear Bars)
    return Point<1U>(1.);
}
template<template<size_t> class SIMPLEX>
double64  FiniteVolumeTraits<1U,SIMPLEX>::ParametricFacetArea( size_t ) const
 {
    return 1.;
 }

template<template<size_t> class SIMPLEX>
Point<1U>  FiniteVolumeTraits<1U,SIMPLEX>::ParametricFacetNormal( size_t ) const
{
   return Point<1U>(1.);
}

template class FiniteVolumeTraits<1U,Element>;
template class FiniteVolumeTraits<1U,Face>;
template class FiniteVolumeTraits<1U,InterFace>;

} // end namespace csmp
