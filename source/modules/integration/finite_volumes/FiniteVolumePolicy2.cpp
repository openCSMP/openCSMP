#include "FiniteVolumePolicy.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

namespace csmp {

template<template<size_t> class SIMPLEX>
FiniteVolumePolicy<2U,SIMPLEX>::FiniteVolumePolicy( const csmp::FiniteVolumeStencil<2U>* fvptr )
 : fvptr_(fvptr)
 {
 }


template<template<size_t> class SIMPLEX>
void FiniteVolumePolicy<2U,SIMPLEX>::AssignFiniteVolume( const csmp::FiniteVolumeStencil<2U>* fvptr )
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( fvptr != nullptr );
    assert( e     != nullptr );
    assert( fvptr->Geometry() == parseFiniteElementDimension( e->FE_Type() ) );
    assert( fvptr->Sectors()  == e->FE()->CornerNodes() );
    fvptr_ = fvptr;
 }

// FV STENCIL INFO

template<template<size_t> class SIMPLEX>
const FiniteVolumeStencil<2U>* const FiniteVolumePolicy<2U,SIMPLEX>::FV_Stencil() const
 {
    return fvptr_;
 }


/// Returns number of sectors, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumePolicy<2U,SIMPLEX>::Sectors() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Sectors();
    return 0;
  }

/// Returns number of facets, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumePolicy<2U,SIMPLEX>::Facets() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Facets();
    return 0;
  }

/// Returns number of integration points per sector, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumePolicy<2U,SIMPLEX>::IntegrationPointsPerSector() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerSector();
    return 0;
  }

/// Returns number of integration points per facet, 0 if no FiniteVolumeStencil assigned
template<template<size_t> class SIMPLEX>
size_t FiniteVolumePolicy<2U,SIMPLEX>::IntegrationPointsPerFacet() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerFacet();
    return 0;
  }


// MAPPING BETWEEN LOCAL AND GLOBAL COORDINATES

template<template<size_t> class SIMPLEX>
Point<2U>  FiniteVolumePolicy<2U,SIMPLEX>::RstToXYZ( const Point<2U>& rst ) const
{
   const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
   N_At( rst );
  
   assert( e  != nullptr );
   double64 sumx(0.), sumy(0.);
   const size_t nodes(e->Nodes());
   for ( size_t i=0U; i<nodes; i++ ) {
       sumx += e->FE()->NRST[i] * e->N(i)->x();
       sumy += e->FE()->NRST[i] * e->N(i)->y();
    }
   // standard RVO
   return Point<2U>(sumx,sumy);
}





// SHAPE FUNCTIONS AT DIFFERENT POINTS

template<template<size_t> class SIMPLEX>
void FiniteVolumePolicy<2U,SIMPLEX>::N_At( const Point<2U>& rst ) const
 {
   const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
   assert( e != nullptr );
    if ( e->IsSurfaceElement() ) {
         e->FE()->Nrs( rst[0], rst[1], e->FE()->NRST );
         return;
      }
    // line element
    e->FE()->Nr( rst[0], e->FE()->NRST );
 }


template<template<size_t> class SIMPLEX>
void FiniteVolumePolicy<2U,SIMPLEX>::N_At( const Point<2U>& rst,
                                           std::vector<double64>& IPOL ) const
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( e != nullptr );
    if ( e->IsSurfaceElement() ) {
         e->FE()->Nrs( rst[0], rst[1], IPOL );
         return;
      }
    // line element
    e->FE()->Nr( rst[0], IPOL );
 }

template<template<size_t> class SIMPLEX>
void FiniteVolumePolicy<2U,SIMPLEX>::N_AtFacetIntegrationPoint( size_t iFacet, size_t ip ) const
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet());

    // surface elements
    assert( e != nullptr );
    if ( e->IsSurfaceElement() ) {
         e->FE()->Nrs( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                       fvptr_->FacetIntegrationPoint( iFacet, ip, 2U ),
                       e->FE()->NRST );
         return;
     }
    // line elements
    e->FE()->Nr( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                 e->FE()->NRST );
 }



template<template<size_t> class SIMPLEX>
void FiniteVolumePolicy<2U,SIMPLEX>::N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( ip < fvptr_->IntegrationPointsPerSector());

    // surface elements
    assert( e != nullptr );
    if ( e->IsSurfaceElement() ) {
          e->FE()->Nrs( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ),
                        fvptr_->SectorIntegrationPoint( iSector, ip, 2U ),
                        e->FE()->NRST );
         return;
    }

    // line elements
    e->FE()->Nr( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ), e->FE()->NRST );
 }



// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS


template<template<size_t> class SIMPLEX>
void FiniteVolumePolicy<2U,SIMPLEX>::Local_dN_At( const Point<2U>& rst ) const
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( e != nullptr );
    if ( e->IsSurfaceElement() ) {
         e->FE()->dNr( rst[0], rst[1], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], e->FE()->DNS );
         return;
      }
    // line element
    e->FE()->dNr( rst[0], e->FE()->DNR );
 }



template<template<size_t> class SIMPLEX>
double64 FiniteVolumePolicy<2U,SIMPLEX>::dN_At( const Point<2U>& rst,
                                                DenseMatrix<DM_MIN>& DN ) const
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( e != nullptr );
    e->CoordinateMatrix();
    if ( e->IsSurfaceElement() ) {
         e->FE()->dNr( rst[0], rst[1], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], e->FE()->DNS );

         e->FE()->Jacobian( e->FE()->DNR, e->FE()->DNS );

         DN.Resize(2U,e->Nodes());
         for ( size_t i=0U; i<e->Nodes(); i++ ) {
              DN(0U,i) = e->FE()->DNR[i];
              DN(1U,i) = e->FE()->DNS[i];
           }
         double64 detJ(e->FE()->JacobianInverse());
         DN = e->FE()->JINV * DN;
         return detJ;
      }

// NOT DONE YET
    // line element
    e->FE()->dNr( rst[0], e->FE()->DNR );
    std::cout <<"\nFiniteVolumePolicy<2U,SIMPLEX>::dN_At: 1D->2D Jacobian is required to get this right."<< std::endl;
    e->FE()->Jacobian( e->FE()->DNR );
    return e->FE()->JacobianInverse();

 } // end dN_At



// PROPERTIES AT FACET AND SECTOR INTEGRATION POINTS

// only for scalars
template<template<size_t> class SIMPLEX>
double64 FiniteVolumePolicy<2U,SIMPLEX>::PropertyValueAtFacetIntegrationPoint(
                                                 size_t iFacet,
                                                 size_t ip,
                                                 const csmp::Index& prop_key ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

  if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

  N_AtFacetIntegrationPoint( iFacet, ip );

  assert( e != nullptr );
  double64  sum(0.);
  const size_t nodes(e->Nodes());
  for ( size_t i=0U; i<nodes; i++ )
    sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

  return sum;

} // end PropertyValueAtFacetIntegrationPoint





// only for scalars
template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::PropertyValueAtSectorIntegrationPoint(
                                                            size_t iSector,
                                                            size_t ip,
                                                            const csmp::Index& prop_key ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( ip < fvptr_->IntegrationPointsPerSector());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    N_AtSectorIntegrationPoint( iSector, ip );

    assert(  e != nullptr );
    double64  sum(static_cast<double64>(0.));
    const size_t nodes(e->Nodes());
    for ( size_t i=0U; i<nodes; i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

 } // end PropertyValueAtVolumeIntegrationPoint

template<template<size_t> class SIMPLEX>
template<class Var>
void FiniteVolumePolicy<2U,SIMPLEX>::PropertyValueAtFacetIntegrationPoint(
                                                 const csmp::Index& prop_key,
                                                 size_t iFacet,
                                                 size_t ip,
                                                 Var& var ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key, var );

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interppolating property to integration point
    Var temp;
    temp.Size( prop_key.dataDepth );
    var.Size( prop_key.dataDepth );
    var = 0;
    assert( e != nullptr );
    const size_t nodes(e->Nodes());
    for ( size_t i=0; i<nodes; i++ ) {
        e->N(i)->Read( prop_key, temp );
        var += temp * e->FE()->NRST[i];
     }

} // end PropertyValueAtFacetIntegrationPoint

template<template<size_t> class SIMPLEX>
template<class Var>
void  FiniteVolumePolicy<2U,SIMPLEX>::PropertyValueAtSectorIntegrationPoint(
                                                            const csmp::Index& prop_key,
                                                            size_t iSector,
                                                            size_t ip,
                                                            Var& var ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
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
    const size_t nodes(e->Nodes());
    for ( size_t i=0; i < nodes; i++ )
    {
        e->N(i)->Read( prop_key, temp );
        var += temp * e->FE()->NRST[i];
    }

 } // end PropertyValueAtVolumeIntegrationPoint




// FACET AND SECTOR INTEGRALS

template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::FacetIntegral( size_t iSector,
                                                         const csmp::Index& prop_key ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    assert( e != nullptr );
    e->CoordinateMatrix();

    double64  integral(0.);

    for ( size_t iFacet=0U; iFacet<fvptr_->FacetsPerSector(iSector); iFacet++ )
      for ( size_t j=0U; j<fvptr_->IntegrationPointsPerFacet(); j++ ) {
           fvptr_->FacetIntegrationPoint( iFacet, j, e->FE()->NRST );
           e->FE()->JacobianAt(e->FE()->NRST);
           const double64 detJ(e->FE()->JacobianDeterminant());
           integral += PropertyValueAtFacetIntegrationPoint( iFacet, j, prop_key ) *
                       fvptr_->FacetIntegrationWeight( iFacet, j ) * detJ;
       }

    return integral;

} // end FacetIntegral




template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::SectorIntegral( size_t iSector,
                                                          const csmp::Index& prop_key ) const
 {
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    assert( e != nullptr );
    e->CoordinateMatrix();

    double64 fIntegral(0.);
    for ( size_t j=0U; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
         fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
         e->FE()->JacobianAt(e->FE()->NRST);
         fIntegral += PropertyValueAtSectorIntegrationPoint( iSector, j, prop_key ) *
                      fvptr_->SectorIntegrationWeight( iSector, j ) *
                      e->FE()->JacobianDeterminant();
      }
   return fIntegral;

} // end SectorIntegral




// PROJECTION ON FACET NORMAL

template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::ProjectionOnFacetNormal(
                                                 size_t iFacet,
                                                 const VectorVariable<2U>& vc
                                                 ) const
{
   assert( iFacet < fvptr_->Facets() );
   return dotProduct(FacetNormal(iFacet),vc.P());

} // end ProjectionOnFacetNormal




template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::ProjectionOnFacetNormal( size_t iFacet,
                                                                   const csmp::Index& prop_key ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( prop_key.type  == VECTOR );
    assert( prop_key.place != INTER_FACE  and prop_key.place != ELEMENT_INTEGRATION_POINT );
    assert( fvptr_->IntegrationPointsPerFacet() == 1U );

    Point<2U>  vecNormal(FacetNormal(iFacet));

    VectorVariable<2U>  vc(PLAIN,PLAIN,0.,0.);

    if ( prop_key.place == ELEMENT ) {
       e->Read( prop_key, vc );
       return vecNormal[0] * vc[0] + vecNormal[1] * vc[1];
    }

    // node properties are nterpolated to facet integration points
    const size_t ip(0U);
    fvptr_->FacetIntegrationPoint( iFacet, ip, e->FE()->NRST );
    N_At( Point<2U>(e->FE()->NRST) );

    assert( e != nullptr );
    double64  sum0(0.);
    double64  sum1(0.);
    const size_t nodes(e->Nodes());
    for ( size_t i=0U; i<nodes; i++ ) {
       e->N(i)->Read( prop_key, vc );
       sum0 += e->FE()->NRST[i] * vc[0];
       sum1 += e->FE()->NRST[i] * vc[1];
    }

    return vecNormal[0] * sum0 + vecNormal[1] * sum1;

} // end ProjectionOnFacetNormal




// GEOMETRY


// SECTOR VOLUME

template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::SectorVolume( size_t iSector ) const
 {
     const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
     assert( iSector < fvptr_->Sectors());

     assert( e != nullptr );
     e->CoordinateMatrix();

     double64 fVolume(static_cast<double64>(0.));
     for ( size_t j=0U; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
          fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
          e->FE()->JacobianAt(e->FE()->NRST);
          fVolume += fvptr_->SectorIntegrationWeight( iSector, j ) *
                     e->FE()->JacobianDeterminant();
       }
     return fVolume;

} // end SectorVolume





// FACET AREA AND NORMAL

template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::FacetArea( size_t iFacet ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());

    assert( e != nullptr );
    if ( e->IsSurfaceElement() )
      return RstToXYZ(fvptr_->FacetPoint(iFacet,0U)).DistanceTo( RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
    
    return 1.; // if it is a line element
} // end FacetArea




template<template<size_t> class SIMPLEX>
Point<2U> FiniteVolumePolicy<2U,SIMPLEX>::FacetNormal( size_t iFacet ) const
{
  const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
  assert( iFacet < fvptr_->Facets());

  if ( e->IsSurfaceElement() ) {
    // the segment midpoint is the origin of the facet
    size_t inside_node, outside_node;
    fvptr_->FacetEdgeNodes( iFacet, inside_node, outside_node );
    Point<2U> org(midPoint(e->N(inside_node)->Coordinate(),
                           e->N(outside_node)->Coordinate()) );

    // the element barycenter is the end-point of the facet
    Point<2U> dest(e->BaryCenter());
    // segment center for rotation
    Point<2U>  m(midPoint(org,dest));
    Point<2U>  v( dest - org );
    Point<2U>  np( v[1], -v[0] );
    // clockwise rotation of segment
    org  = m - (0.5 * np);
    dest = m + (0.5 * np);

    // now the normal length must be normalized to 1
    return (dest - org) / v.Length();
  }

   // the normal of a line element is aligned with it and the line element's length
   // is equivalent to its volume
   return (e->N(1U)->Coordinate() - e->N(0U)->Coordinate()) / e->Volume();

} // end FacetNormal



template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::FacetAreaMapped( size_t iFacet ) const
{
    const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());

    if(e->IsLineElement()) return 1.;

    const Point<2U> fp0( RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
    const Point<2U> fp1( RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );

    return fp0.DistanceTo(fp1);

} // end FacetArea


template<template<size_t> class SIMPLEX>
Point<2U> FiniteVolumePolicy<2U,SIMPLEX>::FacetNormalMapped( size_t iFacet ) const
{
   const SIMPLEX<2U>* e( static_cast<const SIMPLEX<2U>*>(this) );
   assert( iFacet < fvptr_->Facets());

  if(e->IsLineElement()) return Point<2U>(1.,0);

  const Point<2U> fp0( RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
  const Point<2U> fp1( RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );

  const double64 detJ ( 0.5 * sqrt( pow(fp0[0]-fp1[0],2) + pow(fp0[1]-fp1[1],2) ) );

  return Point<2U>( -0.5*(fp0[1]-fp1[1])/detJ, 0.5*(fp0[0]-fp1[0])/detJ );
}


template<template<size_t> class SIMPLEX>
double64  FiniteVolumePolicy<2U,SIMPLEX>::ParametricFacetArea( size_t iFacet ) const
{
  assert(iFacet<fvptr_->Facets());

  double64  area(static_cast<double64>(0.));

  for ( size_t j=0U; j<fvptr_->IntegrationPointsPerFacet(); j++ )
    area += fvptr_->FacetIntegrationWeight( iFacet, j );

  return area;

} // end ParametricFacetArea


template<template<size_t> class SIMPLEX>
Point<2U>  FiniteVolumePolicy<2U,SIMPLEX>::ParametricFacetNormal( size_t iFacet ) const
{
   assert( iFacet < fvptr_->Facets());

   return fvptr_->UnitParametricNormalTo( iFacet );
}

template class FiniteVolumePolicy<2U,Element>;
template class FiniteVolumePolicy<2U,Face>;
template class FiniteVolumePolicy<2U,InterFace>;

} // end namespace csmp
