#include "FiniteVolumePolicy.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "CSMP_mathUtilities.h"
#include "Exception.h"

// invokes integration over planar- as opposed to warped facets
// because volume conservation cannot be achieved without it
#ifndef PYRAMID_TRIANGULAR_FACETS
#define PYRAMID_TRIANGULAR_FACETS
#endif

namespace csmp {

/** Default constructor of the traits class.

@section application Application

The Traits constructor is called numerous times, it must not allocate memory or do expensive things.
It initialises the private variable of the class, associated with the element and MemoryManager and
updates Coordinate Matrix of the element.

*/
template<template<uint32_t> class CELL>
FiniteVolumePolicy<3U,CELL>::FiniteVolumePolicy( const csmp::FiniteVolumeStencil<3U>* fvptr )
 : fvptr_(fvptr)
 {
 }


/**
    Connects policy to a finite-volume stencil that must match the parent element.
    Check is performed in DEBUG mode that the stencil matches parent element.
*/
template<template<uint32_t> class CELL>
void FiniteVolumePolicy<3U,CELL>::AssignFiniteVolume( const csmp::FiniteVolumeStencil<3U>* fvptr )
 {
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
    assert( fvptr != nullptr );
    assert( e     != nullptr );
    assert( fvptr->Geometry() == parseFiniteElementDimension( e->FE_Type() ) );
    assert( fvptr->Sectors()  == e->FE()->CornerNodes() );
    fvptr_ = fvptr;
 }


// FV STENCIL INFO

template<template<uint32_t> class CELL>
const FiniteVolumeStencil<3U>* const FiniteVolumePolicy<3U,CELL>::FV() const
 {
    return fvptr_;
 }


/// Returns number of sectors, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<3U,CELL>::Sectors() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Sectors();
    return 0;
  }

/// Returns number of facets, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<3U,CELL>::Facets() const
  {
    if( fvptr_ != nullptr ) return fvptr_->Facets();
    return 0;
  }

/// Returns number of integration points per sector, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<3U,CELL>::IntegrationPointsPerSector() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerSector();
    return 0;
  }

/// Returns number of integration points per facet, 0 if no FiniteVolumeStencil assigned
template<template<uint32_t> class CELL>
uint32_t FiniteVolumePolicy<3U,CELL>::IntegrationPointsPerFacet() const
  {
    if( fvptr_ != nullptr ) return fvptr_->IntegrationPointsPerFacet();
    return 0;
  }









/// SKM add-on (no rst vector needed !)
template<template<uint32_t> class CELL>
void FiniteVolumePolicy<3U,CELL>::N_AtFacetIntegrationPoint( uint32_t iFacet, uint32_t ip ) const
 {
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet());

    // volume elements first, since speed matters the most
    assert( e != nullptr );
    if ( e->IsVolumeElement() ) {
         e->FE()->Nrst( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                        fvptr_->FacetIntegrationPoint( iFacet, ip, 1U ),
                        fvptr_->FacetIntegrationPoint( iFacet, ip, 2U ),
                        e->FE()->NRST );
         return; 
      }
      
    // surface elements  
    if ( e->IsSurfaceElement() ) {
      e->FE()->Nrs( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ),
                    fvptr_->FacetIntegrationPoint( iFacet, ip, 1U ),
                    e->FE()->NRST );
         return;
      }
    
    // line elements 
    e->FE()->Nr( fvptr_->FacetIntegrationPoint( iFacet, ip, 0U ), 
                 e->FE()->NRST );
 }
 


 
/// SKM add-on
template<template<uint32_t> class CELL>
void FiniteVolumePolicy<3U,CELL>::N_AtSectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const
 {
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
    assert( iSector < fvptr_->Sectors());
    assert( ip < fvptr_->IntegrationPointsPerSector());

    // volume elements first, since speed matters the most
    assert( e != nullptr );
    if ( e->IsVolumeElement() ) {
         e->FE()->Nrst( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ),
                        fvptr_->SectorIntegrationPoint( iSector, ip, 1U ),
                        fvptr_->SectorIntegrationPoint( iSector, ip, 2U ),
                        e->FE()->NRST );
         return; 
      }
      
    // surface elements  
    if ( e->IsSurfaceElement() ) {
      e->FE()->Nrs( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ),
                    fvptr_->SectorIntegrationPoint( iSector, ip, 1U ),
                    e->FE()->NRST );
         return;
      }
    
    // line elements 
    e->FE()->Nr( fvptr_->SectorIntegrationPoint( iSector, ip, 0U ), 
                 e->FE()->NRST );
 }





// DERIVATIVES OF SHAPE FUNCTIONS AT DIFFERENT POINTS

/** The method gets parametric derivatives of interpolation functions at given
parametric coordinates.

@section application Application

Normally the method is called in connection with all the operations, using
the transformation from parametric to physical space and back. For that
purpose, the user can use public variable of the Element class, containing
derivateves: stl vectors DNR,DNS,DNT. For example, before constructing
the Jacobian matrix of the transformation rst->xyz the user should always
call this method.
*/
template<template<uint32_t> class CELL>
void FiniteVolumePolicy<3U,CELL>::Local_dN_At( const Point<3U>& rst ) const
 {
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
    assert( e != nullptr );
    if ( e->IsVolumeElement() ) {
         e->FE()->dNr( rst[0], rst[1], rst[2], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], rst[2], e->FE()->DNS );
         e->FE()->dNt( rst[0], rst[1], rst[2], e->FE()->DNT );
         return;
      }
    if ( e->IsSurfaceElement() ) {
         e->FE()->dNr( rst[0], rst[1], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], e->FE()->DNS );
         return;
      }
    // line element
    e->FE()->dNr( rst[0], e->FE()->DNR );
 }



/**
     @attention SKM: Caution! - not clear what FEM matrices have to do with FV traits
*/
template<template<uint32_t> class CELL>
double FiniteVolumePolicy<3U,CELL>::dN_At( const Point<3U>& rst,
                                                DenseMatrix<DM_MIN>& DN ) const
 {
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
    const uint32_t nodes(e->Nodes());
    assert( e != nullptr );
    e->CoordinateMatrix();

    if ( e->IsVolumeElement() ) {
         e->FE()->dNr( rst[0], rst[1], rst[2], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], rst[2], e->FE()->DNS );
         e->FE()->dNt( rst[0], rst[1], rst[2], e->FE()->DNT );

         e->FE()->Jacobian( e->FE()->DNR, e->FE()->DNS, e->FE()->DNT );

         DN.Resize(3U,nodes);
         for ( auto i{0U}; i<nodes; i++ ) {
              DN(0U,i) = e->FE()->DNR[i];
              DN(1U,i) = e->FE()->DNS[i];
              DN(2U,i) = e->FE()->DNT[i];
           }
         double detJ(e->FE()->JacobianInverse());
         DN = e->FE()->JINV * DN;
         return detJ;
      }

    if ( e->IsSurfaceElement() ) {
         e->FE()->dNr( rst[0], rst[1], e->FE()->DNR );
         e->FE()->dNs( rst[0], rst[1], e->FE()->DNS );
         // compute Jacobian matrix
         e->FE()->JAC.Resize(2U,3U);
         e->FE()->JAC.Zero();
         for ( auto i{0U}; i<3U; i++ )
           for ( uint32_t j{0U}; j<nodes; j++ ) {
                e->FE()->JAC(0U,i) += e->FE()->DNR[j] *
                                                         e->FE()->XY(j,i);
                e->FE()->JAC(1U,i) += e->FE()->DNS[j] *
                                                         e->FE()->XY(j,i);
             }
         // compute Jacobi J' := "determinant" of the 3x2 Jacobian
         // using equation J' = ( E * g - F^2 )^0.5 see CELL-development-in-CSP.doc equation (15)
         // compute E, F, and g
         double  efg0(0.), efg1(0.), efg2(0.);
         for ( uint32_t i{0U}; i<3U; i++ ) {
              efg0 += e->FE()->JAC(0U,i) * e->FE()->JAC(0U,i);
              efg1 += e->FE()->JAC(0U,i) * e->FE()->JAC(1U,i);
              efg2 += e->FE()->JAC(1U,i) * e->FE()->JAC(1U,i);
           }
         double  detJ(sqrt(efg0 * efg2 - efg1 * efg1));
         double  det_inverse(1. / (detJ * detJ));

         for ( uint32_t i{0U}; i<nodes; i++ ) {
              DN(0U,i)  = det_inverse * e->FE()->JAC(0U,0U) *
                        (efg2 * e->FE()->DNR[i] - efg1 * e->FE()->DNS[i] );
              DN(0U,i) += det_inverse * e->FE()->JAC(1U,0U) *
                        (efg0 * e->FE()->DNS[i] - efg1 * e->FE()->DNR[i] );
              DN(1U,i)  = det_inverse * e->FE()->JAC(0U,1U) *
                        (efg2 * e->FE()->DNR[i] - efg1 * e->FE()->DNS[i] );
              DN(1U,i) += det_inverse * e->FE()->JAC(1U,1U) *
                         (efg0 * e->FE()->DNS[i] - efg1 * e->FE()->DNR[i] );
              DN(2U,i)  = det_inverse * e->FE()->JAC(0U,2U) *
                         (efg2 * e->FE()->DNR[i] - efg1 * e->FE()->DNS[i] );
              DN(2U,i) += det_inverse * e->FE()->JAC(1U,2U) *
                         (efg0 * e->FE()->DNS[i] - efg1 * e->FE()->DNR[i] );
           }
         return detJ;
      }

    // line element, TODO: 2D->3d Jacobian
    e->FE()->dNr( rst[0], e->FE()->DNR );
    std::cerr <<"\nFiniteVolumePolicy<3U,CELL>::dN_At: 2D->3d Jacobian is required to get this right."<< std::endl;
    e->FE()->Jacobian( e->FE()->DNR );
    return e->FE()->JacobianInverse();

 } // end dN_At















// PROPERTIES AT FACET AND SECTOR INTEGRATION POINTS



// only for scalars
/** Returns double value of property at given internal facet integration point.

@param iFacet index (No) of the internal facet isnside the FE
@param ip index (No) of the integration point on the facet
@param prop_key key of the property name

@section application Application

For current implementation index ip is constrained to 0 only, i.e-> the
FVPEM method is working with 1 facet integration point only.
*/
template<template<uint32_t> class CELL>
double FiniteVolumePolicy<3U,CELL>::PropertyValueAtFacetIntegrationPoint(
                                                 uint32_t iFacet,
                                                 uint32_t ip,
                                                 const csmp::Index& prop_key ) const
{
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

    assert( iFacet < fvptr_->Facets());
    assert( ip < fvptr_->IntegrationPointsPerFacet() );
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    // if this is just a property that is constant on the element
    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    // get integration point location and corresponding shape function values
    N_AtFacetIntegrationPoint( iFacet, ip );

    // interpolating property to integration point
    double  sum(static_cast<double>(0.));

    assert( e != nullptr );
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i<nodes; i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

} // end PropertyValueAtFacetIntegrationPoint




/**

Returns double value of property at given volume integration point, inside the
volumetric sector of FE, composing FV.

@param ip index (No) of the integration point on the facet;
@param prop_key key of the property name;

@section application Application

For current implementation index ip is constrained to 0 only, i.e-> the
FVPEM method is working with 1 volume sector integration point only.
*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::PropertyValueAtSectorIntegrationPoint( uint32_t iSector,
                                                                            uint32_t ip,
                                                                            const csmp::Index& prop_key ) const
{
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

    assert( iSector < fvptr_->Sectors());
    assert( ip<fvptr_->IntegrationPointsPerSector());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    // if this is just a property that is constant on the element
    if ( prop_key.place == ELEMENT ) return e->Read( prop_key );

    N_AtSectorIntegrationPoint( iSector, ip );

    assert( e != nullptr );
    // interpolating property to integration point
    double  sum(0.);
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i<nodes; i++ )
      sum += e->FE()->NRST[i] * e->N(i)->Read( prop_key );

    return sum;

 } // end PropertyValueAtVolumeIntegrationPoint





/**
    @attention Does not work for array variables at the moment.
*/
template<template<uint32_t> class CELL>
template<class Var>
void FiniteVolumePolicy<3U,CELL>::PropertyValueAtFacetIntegrationPoint( const csmp::Index& prop_key,
                                                                        uint32_t iFacet,
                                                                        uint32_t ip,
                                                                        Var& var ) const
{
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

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
    var = 0;
    
    // only if we are dealing with an array or flagged array variable
    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( prop_key.dataDepth );
        var.Resize( prop_key.dataDepth, 0. );
      }
      
    assert( e != nullptr );
    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i<nodes; i++ ) {
        e->N(i)->Read( prop_key, temp );
        var += temp * e->FE()->NRST[i];
    }

} // end PropertyValueAtFacetIntegrationPoint

template void FiniteVolumePolicy<3U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<3U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<3U>& ) const;
template void FiniteVolumePolicy<3U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<3U>& ) const;

template void FiniteVolumePolicy<3U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<3U>& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<3U>& ) const;

template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<3U>& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<3U>& ) const;

template void FiniteVolumePolicy<3U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;

template void FiniteVolumePolicy<3U,Element>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtFacetIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;




template<template<uint32_t> class CELL>
template<class Var>
void  FiniteVolumePolicy<3U,CELL>::PropertyValueAtSectorIntegrationPoint( const csmp::Index& prop_key,
                                                                             uint32_t iSector,
                                                                             uint32_t ip,
                                                                             Var& var ) const
{
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
    assert( iSector<fvptr_->Sectors());
    assert( ip<fvptr_->IntegrationPointsPerSector());
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
         e->Read( prop_key, var );
         return;
      }

    N_AtSectorIntegrationPoint( iSector, ip );
    Var temp;
    var = 0;

    if constexpr( TypeMatchesVariableType<Var,ARRAY>::value || TypeMatchesVariableType<Var,FLAGGEDARRAY>::value ) {
        temp.Resize( prop_key.dataDepth );
        var.Resize( prop_key.dataDepth, 0. );
      }

    assert( e != nullptr );

    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i < nodes; i++ ) {
         e->N(i)->Read( prop_key, temp );
         var += temp * e->FE()->NRST[i];
      }

 } // end PropertyValueAtSectorIntegrationPoint

template void FiniteVolumePolicy<3U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<3U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<3U>& ) const;
template void FiniteVolumePolicy<3U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<3U>& ) const;

template void FiniteVolumePolicy<3U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<3U>& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<3U>& ) const;

template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ScalarVariable& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, VectorVariable<3U>& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, TensorVariable<3U>& ) const;

template void FiniteVolumePolicy<3U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, ArrayVariable& ) const;

template void FiniteVolumePolicy<3U,Element>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<3U,Face>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;
template void FiniteVolumePolicy<3U,InterFace>::PropertyValueAtSectorIntegrationPoint( const Index&, uint32_t, uint32_t, FlaggedArrayVariable& ) const;




// FACET AND SECTOR INTEGRALS


/** Method calculates the surface integral for given scalar property and
finite volume sector.

@param iSector index (No) of the target sector of the FE;
@param prop_key key of the integrated property name;

@section application Application

As the method calculates the value of the surface integral for given property and
facet, the value is returned by the method.
*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::FacetIntegral( uint32_t iSector,
                                                         const csmp::Index& prop_key ) const
{
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    assert( e != nullptr );
    e->CoordinateMatrix();

    double  integral(static_cast<double>(0.));
    for ( uint32_t iFacet=0U; iFacet<fvptr_->FacetsPerSector(iSector); iFacet++ )
      for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerFacet(); j++ ) {
             // 1. take integration point location from FV stencil
             fvptr_->FacetIntegrationPoint( iFacet, j, e->FE()->NRST );
             // 2. Compute the Jacobian
             e->FE()->JacobianAt(e->FE()->NRST);

             // 3. Add contribution to integral
             integral += PropertyValueAtFacetIntegrationPoint( iFacet, j, prop_key ) *
                         fvptr_->FacetIntegrationWeight( iFacet, j ) *
                         e->FE()->JacobianDeterminant();
        }
    return integral;

} // end FacetIntegral





/**

Method calculates the value of the volume integral for given scalar property and
volume sector of finite element.

@param iSector index (No) of the internal sector inside the FE;
@param prop_key key of the integrated property name;

@section application Application

As the method calculates the value of the volume integral for the given property and
volume sector, the value is returned by the method.
*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::SectorIntegral( uint32_t iSector,
                                                          const csmp::Index& prop_key ) const
 {
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

    assert( iSector < fvptr_->Sectors());
    assert( prop_key.type  == SCALAR );
    assert( prop_key.place == NODE or prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE );

    assert( e != nullptr );
    e->CoordinateMatrix();

    double fIntegral(0.);
    for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
        fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
        e->FE()->JacobianAt(e->FE()->NRST );
        fIntegral += PropertyValueAtSectorIntegrationPoint( iSector, j, prop_key ) *
                     fvptr_->SectorIntegrationWeight( iSector, j ) *
                     e->FE()->JacobianDeterminant();
      }
    return fIntegral;

} // end VolumeIntegral









// PROJECTION ON FACET NORMAL

/**

 Overloaded method carries out projection of the given vector property on the normal of
 the internal facet. In this implementation of the method, projection is done
 in physical space-> The surface parameter corresponds to local number of surfaces
 per node, 3 or 4  for vertice 5 of the pyramid.

@param iFacet index (No) of the internal facet isnside the FE;
@param vc vector variable (not by database name) for projection;

@return returned value of the determinant of the transformation
of the topologically quadrilateral facet from parametric
space of that surface to physical space;

@section application Application

This overloaded method defines the normal to the facet in physical space, using standard
FE methodology of normal definition, fitting the quadrilateral patch to the
corner points of the facet in physical space-> This method could be rather
computationally expensive, therefore parametric projection on surface normal
is more efficient.

For current implementation index ip is constrained to 0 only, i.e-> the
FVPEM method is working with 1 facet integration point only.

*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::ProjectionOnFacetNormal(
                                                 uint32_t iFacet,
                                                 const VectorVariable<3U>& vc
                                                 ) const
{
    //1. Get physical surface normal based on the definition of the quad surface
    //2.  if this is just a property that is constant on the element, return projection
    return dotProduct(FacetNormal(iFacet), vc.P());

} // end ProjectionOnFacetNormal




/**

 Method carries out projection of the given by name vector property on the normal of
 the internal facet. In this implementation of the method, projection is done
 in physical space-> The surface parameter corresponds to local number of surfaces
 per node, 3 or 4  for vertice 5 of the pyramid.

@param iFacet index (No) of the internal facet inside the FE;
@param prop_key key of the property name;

@return returned value of the determinant of the transformation
of the topologically quadrilateral facet from parametric
space of that surface to physical space;

@section application Application

The method defines the normal to the facet in physical space, using standard
FE methodology of normal definition, fitting the quadrilateral patch to the
corner points of the facet in physical space-> This method could be rather
computationally expensive, therefore parametric projection on surface normal
is more efficient.

For current implementation index ip is constrained to 0 only, i.e-> the
FVPEM method is working with 1 facet integration point only.
*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::ProjectionOnFacetNormal( uint32_t iFacet,
                                                               const csmp::Index& prop_key
                                                             ) const
{
   const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

    assert( iFacet < fvptr_->Facets());
    assert( fvptr_->IntegrationPointsPerFacet() == 1U );
    assert( prop_key.type  == VECTOR );
    assert( prop_key.place != INTER_FACE && prop_key.place != ELEMENT_INTEGRATION_POINT );

    VectorVariable<3U>  vc(PLAIN,PLAIN,PLAIN,0.,0.,0.);

    //2.  if this is just a property that is constant on the element, return projection
    if ( prop_key.place == ELEMENT ) {
         e->Read( prop_key, vc );
         return dotProduct(FacetNormal(iFacet), vc.P());
      }

    //3. If not, interpolate node property to integration point
    const uint32_t ip(0U);
    fvptr_->FacetIntegrationPoint( iFacet, ip, e->FE()->NRST );
    assert( e != nullptr );
    e->N_At( Point<3U>(e->FE()->NRST) );

    double sum0(static_cast<double>(0.));
    double sum1(static_cast<double>(0.));
    double sum2(static_cast<double>(0.));

    const uint32_t nodes(e->Nodes());
    for ( uint32_t i{0U}; i<nodes; i++ ) {
         e->N(i)->Read( prop_key, vc );
         sum0 += e->FE()->NRST[i] * vc[0];
         sum1 += e->FE()->NRST[i] * vc[1];
         sum2 += e->FE()->NRST[i] * vc[2];
      }

    Point<3U>  vecNormal(FacetNormal(iFacet));

    return vecNormal[0] * sum0 + vecNormal[1] * sum1 + vecNormal[2] * sum2;

} // end ProjectionOnFacetNormal




// GEOMETRY


// SECTOR VOLUME

/**

Method calculates the volume for the given volumetric sector.

@param iSector index (No) of the internal facet isnside the FE;

@section application Application

Method automatically takes into account composite nature of the sector at the apex of
the pyramid.
*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::SectorVolume( uint32_t iSector ) const
 {
   const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );

   assert( iSector < fvptr_->Sectors());

   assert( e != nullptr );
   e->CoordinateMatrix();

   double fVolume(0.);
   for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerSector(); j++ ) {
        fvptr_->SectorIntegrationPoint( iSector, j, e->FE()->NRST );
        e->FE()->JacobianAt(e->FE()->NRST);
        fVolume += fvptr_->SectorIntegrationWeight(iSector, j) *
                   e->FE()->JacobianDeterminant();
     }
   return fVolume;
}


// FACET AREA AND NORMAL


/** Method calculates the surface area for the given facet.

@param iFacet index (No) of the internal facet inside the FE;

@section application Application
This method takes the tabulated facet points, converts them to physical space,
and then calculates the area of the given facet.
If the element is a surface element, then only two points describe the facet,
therefore, only the distance between these two points will be calculated.

*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::FacetArea( uint32_t iFacet ) const
{
   const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
   assert( iFacet < fvptr_->Facets());

    assert( e != nullptr );
      if ( e->IsVolumeElement() ) {
#ifdef PYRAMID_TRIANGULAR_FACETS
           if(fvptr_->FacetPoints(iFacet) == 3U )
             return triangleArea(e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)),
                                 e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)),
                                 e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)));
#endif
             return facetArea1(e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)),
                               e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)),
                               e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)),
                               e->RstToXYZ(fvptr_->FacetPoint(iFacet,3U)) );
        }
      else if ( e->IsSurfaceElement() ) {
           return e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)).DistanceTo(
                  e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
        }
      return 1.; //its a line element

} // end FacetArea



/**

Normal is generated in physical space,
using quadrilateral patch, fitted to the facet.
Method returns the determinant of the Jacobian for the surface
at integration point (=facet area).

@param iFacet index (No) of the internal facet is inside the FE;

@return Point containing physical normal coordinates.

@section application Application

This  method defines the normal to the facet in physical space, using standard
FE methodology of normal definition, fitting the quadrilateral patch to the
corner points of the facet in physical space-> This method could be rather
computationally expensive, therefore parametric definition of the surface normal
is more efficient.

For current implementation index ip is constrained to 0 only, i.e-> the
FVPEM method is working with 1 facet integration point only.

*/
template<template<uint32_t> class CELL>
Point<3U>  FiniteVolumePolicy<3U,CELL>::FacetNormal( uint32_t iFacet ) const
{
    const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
#ifndef NDEBUG
    assert( iFacet < fvptr_->Facets());
#endif

  if ( e->IsVolumeElement() ) {
      switch (fvptr_->FacetPoints(iFacet)) {
          case 3:
              return normalOfTriangle ( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)),
                                        e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)),
                                        e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)));
              
          case 4:
              return normalAtFacetCenter( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)),
                                          e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)),
                                          e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)),
                                          e->RstToXYZ(fvptr_->FacetPoint(iFacet,3U)) );

          default:
              throw Exception(ERROR, "FiniteVolumePolicy<3U,CELL>::FacetNormal", "Unexpected facet point count" );

      }
    }

    assert( e != nullptr );
   if ( e->IsSurfaceElement() ) {
       // this version has the best definition
       Point<3U> segment(midPoint(e->N(fvptr_->InsideNode(iFacet))->Coordinate(),
                                  e->N(fvptr_->OutsideNode(iFacet))->Coordinate()) - e->BaryCenter());
       // segment center for rotation
       e->CoordinateMatrix(); // needed by UnitNormal
       e->FE()->UnitNormal( e->FE()->NRST );
       Point<3U>  enrml( e->FE()->NRST );
       Point<3U>  vecNormal(crossProduct( enrml, segment ));
       vecNormal.NormalizeLengthTo(1.);
       return vecNormal;
    }

  // line element: unit normal is parallel to element and point from node i to i+1
  const double elength(e->Volume());
  return Point<3U>( (e->N(1U)->x()-e->N(0U)->x()) / elength,
                    (e->N(1U)->y()-e->N(0U)->y()) / elength,
                    (e->N(1U)->z()-e->N(0U)->z()) / elength );

} // end FacetNormal




template<template<uint32_t> class CELL>
  Point<3U>  FiniteVolumePolicy<3U,CELL>::FacetPoint( uint32_t iFacet, uint32_t iPoint ) const
  {
      const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
      return e->RstToXYZ(fvptr_->FacetPoint(iFacet,iPoint));
  }


template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::FacetAreaMapped( uint32_t iFacet ) const
{
  const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
  if(e->IsLineElement()) return 1.;

  assert( e != nullptr );
  if(e->IsSurfaceElement())
  {
    const Point<3U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
    const Point<3U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );

    const double detJ ( sqrt( pow(fp0[0]-fp1[0],2) + pow(fp0[1]-fp1[1],2) + pow(fp0[2]-fp1[2],2) ));
    return detJ;
  }

  // it is a volume element
#ifdef PYRAMID_TRIANGULAR_FACETS
  if(fvptr_->FacetPoints(iFacet) == 3)
  {
    const Point<3U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
    const Point<3U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
    const Point<3U> fp2( e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)) );

    const double j1( -fp0[0]+fp1[0] ),
                   j2( -fp0[1]+fp1[1] ),
                   j3( -fp0[2]+fp1[2] ),
                   j4( -fp0[0]+fp2[0] ),
                   j5( -fp0[1]+fp2[1] ),
                   j6( -fp0[2]+fp2[2] );

    // 0.5 is the area of the triangle in local space
    const double det( sqrt((j1*j1+j2*j2+j3*j3)*(j4*j4+j5*j5+j6*j6)-pow((j1*j4+j2*j5+j3*j6),2)) );

    const double area = 0.5 * det;

    return area;
  }
#endif

  const Point<3U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
  const Point<3U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
  const Point<3U> fp2( e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)) );
  const Point<3U> fp3( e->RstToXYZ(fvptr_->FacetPoint(iFacet,3U)) );

  const double j1( -0.25* (fp0[0]-fp1[0]-fp2[0]+fp3[0] ) ),
                 j2( -0.25* (fp0[1]-fp1[1]-fp2[1]+fp3[1] ) ),
                 j3( -0.25* (fp0[2]-fp1[2]-fp2[2]+fp3[2] ) ),
                 j4( -0.25* (fp0[0]+fp1[0]-fp2[0]-fp3[0] ) ),
                 j5( -0.25* (fp0[1]+fp1[1]-fp2[1]-fp3[1] ) ),
                 j6( -0.25* (fp0[2]+fp1[2]-fp2[2]-fp3[2] ) );

  // 4 is the area of the quadrilateral in local space
  const double det( sqrt((j1*j1+j2*j2+j3*j3)*(j4*j4+j5*j5+j6*j6)-pow((j1*j4+j2*j5+j3*j6),2)) );

  const double area = 4. * det;

  return area;

}


/**
    Maps the normal from parametric space to physical space using a 1D to 3D Jacobian in a volumetric element
*/
template<template<uint32_t> class CELL>
Point<3U>  FiniteVolumePolicy<3U,CELL>::FacetNormalMapped( uint32_t iFacet ) const
{
  const CELL<3U>* e( static_cast<const CELL<3U>*>(this) );
  assert( e != nullptr );
  if( e->IsLineElement()) return Point<3U>(1.,0.,0.);

  if( e->IsSurfaceElement()) // project to bar element - rotate coplanar to the element by 90 degrees
    {
      const Point<3U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
      const Point<3U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );

      // not needed by RST to XYZ but by UnitNormal
      e->CoordinateMatrix();
      e->FE()->UnitNormal( e->FE()->NRST );
      const double x(e->FE()->NRST[0]), y(e->FE()->NRST[1]), z(e->FE()->NRST[2]);

      const double j1(-0.5*(fp0[0]-fp1[0]));
      const double j2(-0.5*(fp0[1]-fp1[1]));
      const double j3(-0.5*(fp0[2]-fp1[2]));

      Point<3U>  normal(  x*x*j1    + (x*y+z)*j2 + (x*z-y)*j3,
                         (x*y-z)*j1 +  y*y*j2    + (y*z+x)*j3,
                         (z*x+y)*j1 + (z*y-x)*j2 +  z*z*j3 );

      normal.NormalizeLengthTo(1.);

      return normal;
    }

  //e->IsVolumeElement
    double j1,j2,j3,j4,j5,j6;
    switch (fvptr_->FacetPoints(iFacet))
    {
        case 3:
        {
            const Point<3U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
            const Point<3U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
            const Point<3U> fp2( e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)) );
            j1 = -fp0[0]+fp1[0];
            j2 = -fp0[1]+fp1[1];
            j3 = -fp0[2]+fp1[2];
            j4 = -fp0[0]+fp2[0];
            j5 = -fp0[1]+fp2[1];
            j6 = -fp0[2]+fp2[2];
            break;
        }
            
        case 4:
        {
            const Point<3U> fp0( e->RstToXYZ(fvptr_->FacetPoint(iFacet,0U)) );
            const Point<3U> fp1( e->RstToXYZ(fvptr_->FacetPoint(iFacet,1U)) );
            const Point<3U> fp2( e->RstToXYZ(fvptr_->FacetPoint(iFacet,2U)) );
            const Point<3U> fp3( e->RstToXYZ(fvptr_->FacetPoint(iFacet,3U)) );
            
            j1 = -0.25 * (fp0[0]-fp1[0]-fp2[0]+fp3[0] );
            j2 = -0.25 * (fp0[1]-fp1[1]-fp2[1]+fp3[1] );
            j3 = -0.25 * (fp0[2]-fp1[2]-fp2[2]+fp3[2] );
            j4 = -0.25 * (fp0[0]+fp1[0]-fp2[0]-fp3[0] );
            j5 = -0.25 * (fp0[1]+fp1[1]-fp2[1]-fp3[1] );
            j6 = -0.25 * (fp0[2]+fp1[2]-fp2[2]-fp3[2] );
            break;
        }
            
        default:
        {
            throw Exception(ERROR, "FiniteVolumePolicy<3U,CELL>::FacetNormalMapped", "Unexpected facet point count" );
        }
    }

  const double det( fabs((j1*j1+j2*j2+j3*j3)*(j4*j4+j5*j5+j6*j6)-square(j1*j4+j2*j5+j3*j6) ));

  const double jj1 ((j1*j1+j2*j2+j3*j3)/det),
                  jj2 (-(j1*j4+j2*j5+j3*j6)/det),
                  jj3 (jj2),
                  jj4 ((j4*j4+j5*j5+j6*j6)/det);

  Point<3U>  v0( j1*jj2+j4*jj4,
                    j2*jj2+j5*jj4,
                    j3*jj2+j6*jj4 );

  Point<3U>  v1( j1*jj1+j4*jj3,
                    j2*jj1+j5*jj3,
                    j3*jj1+j6*jj3 );

  Point<3U>    normal( crossProduct(v1,v0) );
  normal.NormalizeLengthTo(1.);

  return normal;
}



/** Method calculates the value of the facet area in the parametric space of the given
finite element.

@param iFacet index (No) of the internal facet isnside the FE;

@section application Application

The method is applicable to the algorithm of the parametric calculation of fluxes.
It should be used  together with the method  ParametricProjectionOnFacetNormal().
*/
template<template<uint32_t> class CELL>
double  FiniteVolumePolicy<3U,CELL>::ParametricFacetArea( uint32_t iFacet ) const
{
  assert(iFacet<fvptr_->Facets());

  double  area(static_cast<double>(0.));
  
  // loop that loops only if there are indeed multiple integration points
  for ( uint32_t j{0U}; j<fvptr_->IntegrationPointsPerFacet(); j++ )
    area += fvptr_->FacetIntegrationWeight( iFacet, j );

  return area;

} // end ParametricFacetArea



/**

@return Reference normal is returned.

@param iFacet index (No) of the internal facet isnside the FE;

@return reference to STL vector of the reference normal coordinates.

@section application Application

Application of the method gives the reference normal. This reference normal is used to
find outwards pointing normal to the facet, not taking care of the orientation of the
boundary loop of the facet. Dot product of the reference normal and any normal in physical
or parametric realisation can be used to define outward pointing normal.
*/
template<template<uint32_t> class CELL>
Point<3U>  FiniteVolumePolicy<3U,CELL>::ParametricFacetNormal( uint32_t iFacet ) const
{
  assert(iFacet<fvptr_->Facets());

  // standard RVO
  return fvptr_->UnitParametricNormalTo( iFacet );
  
} // end ParametricFacetNormal



template class FiniteVolumePolicy<3U,Element>;
template class FiniteVolumePolicy<3U,Face>;
template class FiniteVolumePolicy<3U,InterFace>;

} // end namespace csmp
