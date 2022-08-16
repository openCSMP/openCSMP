#include "FiniteVolumePolicy.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

namespace csmp {

template<template<uint32_t> class CELL>
FiniteVolumePolicy<2U,CELL>::FiniteVolumePolicy( const csmp::FiniteVolumeStencil<2U>* fvptr )
 : fvptr_(fvptr)
 {
 }


template<template<uint32_t> class CELL>
void FiniteVolumePolicy<2U,CELL>::AssignFiniteVolume( const csmp::FiniteVolumeStencil<2U>* fvptr )
 {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( fvptr != nullptr );
    assert( e     != nullptr );
    assert( fvptr->Geometry() == parseFiniteElementDimension( e->FE_Type() ) );
    assert( fvptr->Sectors()  == e->FE()->CornerNodes() );
    fvptr_ = fvptr;
 }

// FV STENCIL INFO

template<template<uint32_t> class CELL>
const FiniteVolumeStencil<2U>* const FiniteVolumePolicy<2U,CELL>::FV() const
 {
    return fvptr_;
 }


/// Returns number of sectors, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<2U,CELL>::Sectors() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Sectors();
    return 0;
  }

/// Returns number of facets, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<2U,CELL>::Facets() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Facets();
    return 0;
  }

/// Returns number of integration points per sector, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<2U,CELL>::IntegrationPointsPerSector() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerSector();
    return 0;
  }

/// Returns number of integration points per facet, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<2U,CELL>::IntegrationPointsPerFacet() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerFacet();
    return 0;
  }






// SHAPE FUNCTIONS AT DIFFERENT POINTS

template<template<uint32_t> class CELL>
void FiniteVolumePolicy<2U,CELL>::N_AtFacetIntegrationPoint( uint32_t iFacet, uint32_t ip ) const
 {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet());

    // surface elements
    assert( e != nullptr );
    if ( e->IsSurface() ) {
         e->FE()->Nrs( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                       fvptr_->FacetIntegrationPoint( iFacet, ip, 2U ),
                       e->FE()->NRST );
         return;
     }
    // line elements
    e->FE()->Nr( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                 e->FE()->NRST );
 }


template<template<uint32_t> class CELL>
void FiniteVolumePolicy<2U,CELL>::N_AtFacetIntegrationPoint( uint32_t iFacet, uint32_t ip, std::vector<double>& NRST) const
  {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet());

    // surface elements
    assert( e != nullptr );
    if ( e->IsSurface() ) {
      e->FE()->Nrs( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                    fvptr_->FacetIntegrationPoint( iFacet, ip, 2U ),
                    NRST );
      return;
    }
    // line elements
    e->FE()->Nr( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                 NRST );
  }


template<template<uint32_t> class CELL>
void FiniteVolumePolicy<2U,CELL>::N_AtSectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const
 {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( ip < fvptr_->IntegrationPointsPerSector());

    // surface elements
    assert( e != nullptr );
    if ( e->IsSurface() ) {
          e->FE()->Nrs( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ),
                        fvptr_->SectorIntegrationPoint( iSector, ip, 2U ),
                        e->FE()->NRST );
         return;
    }

    // line elements
    e->FE()->Nr( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ), e->FE()->NRST );
 }



// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS


template<template<uint32_t> class CELL>
void FiniteVolumePolicy<2U,CELL>::Local_dN_At( const Point<2U>& rst ) const
 {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( e != nullptr );
    if ( e->IsSurface() ) {
         e->FE()->dNr( rst[0], rst[1], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], e->FE()->DNS );
         return;
      }
    // line element
    e->FE()->dNr( rst[0], e->FE()->DNR );
 }



template<template<uint32_t> class CELL>
double FiniteVolumePolicy<2U,CELL>::dN_At( const Point<2U>& rst,
                                                DenseMatrix<DM_MIN>& DN ) const
 {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( e != nullptr );
    e->CoordinateMatrix();
    if ( e->IsSurface() ) {
         e->FE()->dNr( rst[0], rst[1], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], e->FE()->DNS );

         e->FE()->Jacobian( e->FE()->DNR, e->FE()->DNS );

         DN.Resize(2U,e->Nodes());
         for ( auto i{0U}; i<e->Nodes(); i++ ) {
              DN(0U,i) = e->FE()->DNR[i];
              DN(1U,i) = e->FE()->DNS[i];
           }
         double detJ(e->FE()->JacobianInverse());
         DN = e->FE()->JINV * DN;
         return detJ;
      }

// NOT DONE YET
    // line element
    e->FE()->dNr( rst[0], e->FE()->DNR );
    std::cout <<"\nFiniteVolumePolicy<2U,CELL>::dN_At: 1D->2D Jacobian is required to get this right."<< std::endl;
    e->FE()->Jacobian( e->FE()->DNR );
    return e->FE()->JacobianInverse();

 } // end dN_At



// PROPERTIES AT FACET AND SECTOR INTEGRATION POINTS

// only for scalars
template<template<uint32_t> class CELL>
double FiniteVolumePolicy<2U,CELL>::PropertyValueAtFacetIntegrationPoint(
                                                 uint32_t iFacet,
                                                 uint32_t ip,
                                                 const csmp::Index& prop_key ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

  if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

  N_AtFacetIntegrationPoint( iFacet, ip );

  assert( e != nullptr );
  double  sum(0.);
  const uint32_t nodes(e->Nodes());
  for ( auto i{0U}; i<nodes; i++ )
    sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

  return sum;

} // end PropertyValueAtFacetIntegrationPoint





// only for scalars
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::PropertyValueAtSectorIntegrationPoint(
                                                            uint32_t iSector,
                                                            uint32_t ip,
                                                            const csmp::Index& prop_key ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( ip < fvptr_->IntegrationPointsPerSector());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    N_AtSectorIntegrationPoint( iSector, ip );

    assert(  e != nullptr );
    double  sum(static_cast<double>(0.));
    const uint32_t nodes(e->Nodes());
    for ( auto i{0U}; i<nodes; i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

 } // end PropertyValueAtVolumeIntegrationPoint




/**
    TODO: test for all relevant cases
*/
template<template<uint32_t> class CELL>
template<class Var>
void FiniteVolumePolicy<2U,CELL>::PropertyValueAtFacetIntegrationPoint(
                                                 const csmp::Index& prop_key,
                                                 uint32_t iFacet,
                                                 uint32_t ip,
                                                 Var& var ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
         e->Read( prop_key, var );
         return;
      }

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interpolating property to integration point
    Var temp;
    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( prop_key.dataDepth );
        var.Resize( prop_key.dataDepth, 0. );
      }
    var = 0;
    assert( e != nullptr );
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i<nodes; i++ ) {
        e->N(i)->Read( prop_key, temp );
        var += temp * e->FE()->NRST[i];
     }

} // end PropertyValueAtFacetIntegrationPoint

template void FiniteVolumePolicy<2U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<2U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<2U>& ) const;
template void FiniteVolumePolicy<2U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<2U>& ) const;

template void FiniteVolumePolicy<2U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<2U>& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<2U>& ) const;

template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<2U>& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<2U>& ) const;

template void FiniteVolumePolicy<2U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;

template void FiniteVolumePolicy<2U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;




template<template<uint32_t> class CELL>
template<class Var>
void  FiniteVolumePolicy<2U,CELL>::PropertyValueAtSectorIntegrationPoint(
                                                            const csmp::Index& prop_key,
                                                            uint32_t iSector,
                                                            uint32_t ip,
                                                            Var& var ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iSector<fvptr_->Sectors());
    assert( ip<fvptr_->IntegrationPointsPerSector());
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
         e->Read( prop_key, var );
         return;
      }

    N_AtSectorIntegrationPoint( iSector, ip );
    Var temp;
    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( prop_key.dataDepth );
        var.Resize( prop_key.dataDepth, 0. );
      }
    var = 0;
  
    assert( e != nullptr );
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i < nodes; i++ )
    {
        e->N(i)->Read( prop_key, temp );
        var += temp * e->FE()->NRST[i];
    }

 } // end PropertyValueAtVolumeIntegrationPoint

template void FiniteVolumePolicy<2U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<2U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<2U>& ) const;
template void FiniteVolumePolicy<2U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<2U>& ) const;

template void FiniteVolumePolicy<2U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<2U>& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<2U>& ) const;

template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<2U>& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<2U>& ) const;

template void FiniteVolumePolicy<2U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;

template void FiniteVolumePolicy<2U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<2U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<2U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;




// FACET AND SECTOR INTEGRALS

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::FacetIntegral( uint32_t iSector,
                                                         const csmp::Index& prop_key ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    assert( e != nullptr );
    e->CoordinateMatrix();

    double  integral(0.);

    for ( uint32_t iFacet=0U; iFacet<fvptr_->FacetsPerSector(iSector); iFacet++ )
      for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerFacet(); j++ ) {
           fvptr_->FacetIntegrationPoint( iFacet, j, e->FE()->NRST );
           e->FE()->JacobianAt(e->FE()->NRST);
           const double detJ(e->FE()->JacobianDeterminant());
           integral += PropertyValueAtFacetIntegrationPoint( iFacet, j, prop_key ) *
                       fvptr_->FacetIntegrationWeight( iFacet, j ) * detJ;
       }

    return integral;

} // end FacetIntegral




template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::SectorIntegral( uint32_t iSector,
                                                          const csmp::Index& prop_key ) const
 {
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    assert( e != nullptr );
    e->CoordinateMatrix();

    double fIntegral(0.);
    for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
         fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
         e->FE()->JacobianAt(e->FE()->NRST);
         fIntegral += PropertyValueAtSectorIntegrationPoint( iSector, j, prop_key ) *
                      fvptr_->SectorIntegrationWeight( iSector, j ) *
                      e->FE()->JacobianDeterminant();
      }
   return fIntegral;

} // end SectorIntegral




// PROJECTION ON FACET NORMAL

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::ProjectionOnFacetNormal(
                                                 uint32_t iFacet,
                                                 const VectorVariable<2U>& vc
                                                 ) const
{
   assert( iFacet < fvptr_->Facets() );
   return dotProduct(FacetNormal(iFacet),vc.P());

} // end ProjectionOnFacetNormal




template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::ProjectionOnFacetNormal( uint32_t iFacet,
                                                                   const csmp::Index& prop_key ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
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
    const uint32_t ip(0U);
    fvptr_->FacetIntegrationPoint( iFacet, ip, e->FE()->NRST );
    e->N_At( Point<2U>(e->FE()->NRST) );

    assert( e != nullptr );
    double  sum0(0.);
    double  sum1(0.);
    const uint32_t nodes(e->Nodes());
    for ( auto i{0U}; i<nodes; i++ ) {
       e->N(i)->Read( prop_key, vc );
       sum0 += e->FE()->NRST[i] * vc[0];
       sum1 += e->FE()->NRST[i] * vc[1];
    }

    return vecNormal[0] * sum0 + vecNormal[1] * sum1;

} // end ProjectionOnFacetNormal




// GEOMETRY


// SECTOR VOLUME

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::SectorVolume( uint32_t iSector ) const
 {
     const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
     assert( iSector < fvptr_->Sectors());

     assert( e != nullptr );
     e->CoordinateMatrix();

     double fVolume(static_cast<double>(0.));
     for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
          fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
          e->FE()->JacobianAt(e->FE()->NRST);
          fVolume += fvptr_->SectorIntegrationWeight( iSector, j ) *
                     e->FE()->JacobianDeterminant();
       }
     return fVolume;

} // end SectorVolume





// FACET AREA AND NORMAL

template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::FacetArea( uint32_t iFacet ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());

    assert( e != nullptr );
    if ( e->IsSurface() )
      return e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)).DistanceTo( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
    
    return 1.; // if it is a line element
} // end FacetArea




template<template<uint32_t> class CELL>
Point<2U> FiniteVolumePolicy<2U,CELL>::FacetNormal( uint32_t iFacet ) const
{
  const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
  assert( iFacet < fvptr_->Facets());

  if ( e->IsSurface() ) {
    // the segment midpoint is the origin of the facet
    uint32_t inside_node, outside_node;
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



template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::FacetAreaMapped( uint32_t iFacet ) const
{
    const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
    assert( iFacet < fvptr_->Facets());

    if(e->IsLine()) return 1.;

    const Point<2U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
    const Point<2U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );

    return fp0.DistanceTo(fp1);

} // end FacetArea


template<template<uint32_t> class CELL>
Point<2U> FiniteVolumePolicy<2U,CELL>::FacetNormalMapped( uint32_t iFacet ) const
{
   const CELL<2U>* e( static_cast<const CELL<2U>*>(this) );
   assert( iFacet < fvptr_->Facets());

  if(e->IsLine()) return Point<2U>(1.,0);

  const Point<2U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
  const Point<2U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );

  // const double detJ ( 0.5 * sqrt( pow(fp0[0]-fp1[0],2) + pow(fp0[1]-fp1[1],2) ) );
    
  const double invDetJ( 1.0 / (0.5 * hypot( fp0[0]-fp1[0], fp0[1]-fp1[1] ) ) );
  return Point<2U>( -0.5*(fp0[1]-fp1[1])*invDetJ, 0.5*(fp0[0]-fp1[0])*invDetJ );
}


template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<2U,CELL>::ParametricFacetArea( uint32_t iFacet ) const
{
  assert(iFacet<fvptr_->Facets());

  double  area(static_cast<double>(0.));

  for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerFacet(); j++ )
    area += fvptr_->FacetIntegrationWeight( iFacet, j );

  return area;

} // end ParametricFacetArea


template<template<uint32_t> class CELL>
Point<2U>  FiniteVolumePolicy<2U,CELL>::ParametricFacetNormal( uint32_t iFacet ) const
{
   assert( iFacet < fvptr_->Facets());

   return fvptr_->UnitParametricNormalTo( iFacet );
}

template class FiniteVolumePolicy<2U,Element>;
template class FiniteVolumePolicy<2U,Face>;
template class FiniteVolumePolicy<2U,InterFace>;

} // end namespace csmp
