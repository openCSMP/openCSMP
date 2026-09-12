#ifndef FINITE_VOLUME_POLICY_H
#define FINITE_VOLUME_POLICY_H

#include "Index.h"
#include "Point.h"
#include "FiniteVolumeStencil.h"
#include "Node.h"

#include "QuadrilateralFacet.h"
#include "TriangularFacet.h"

#include "DenseMatrix.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/// Element policy for node-centered finite-volume computations; see full specializations for 1D, 2D, 3D
template<uint32_t, template<uint32_t> class CELL> class FiniteVolumePolicy;

/**

@brief Policy of the Element class which provides it with the functionality needed
for finite-volume stencil formation, i.e. node-centered finite-volume computations.

@author S.K. Matthai
@date 2003

Class is the implementation of the TRAITS design pattern. The traits template 
technique, pioneered by Nathan Myers,is a means of bundling type-dependent 
declarations together.In essence, using traits allows you to mix and match
certain types and and values with contexts that use them in a flexible manner, 
while keeping your code readable and maintainable.  

In our case for example, such code design permits to maintain efficient access 
to the Element class, thus giving the possibility to tie up Finite Volume and 
Finite Element parts of the stencil implementation.  

 
@section functionality Functionality
 
Provides a set of functionality to perform operations, used in the CVCELL 
transport method, implemented in CSMP. This includes calculation of volumes 
of sectors and surfaces of internal division walls - facets - inside the
element, definition of projections of the velocity at the facet integration
point on the normal to the facet, required for flux calculation, etc. 
  
 
@section motivation Motivation

Part of a consistent framework for generation of the Finite Element/Finite 
Volume method in CSP. Uses the concept of the TRAITS design pattern for 
efficient access to the the Element class of CSP. Give the possibility for
elegant coupling and access to Finite Element and Finite Volume frameworks
inside CSMP.
 
 
@section implementation Implementation

Uses the traints desigh pattern. Currently the class is supporting linear
iso-parametric family of Finite Elements with reduced polynomial interpolation
basis. Distinctive feature of the FiniteVolumePolicy class is the 
possibility to access 
 
 
@section examples Application Examples

No specific application examples are required, as the method is multiply  
called the great code CSP, within the implementation of the classes 
ExplicitFiniteVolumeProcessor and others.  

TODO: @todo (3) Put finiteVolumeAuxiliaryFunctions into FiniteVolumePolicy class (A)
TODO: @todo SKM try to do computations in parametric space

 */
template<template<uint32_t> class CELL>
class FiniteVolumePolicy<3U, CELL> {
  public:
    FiniteVolumePolicy( const csmp::FiniteVolumeStencil<3U>* = nullptr );
    FiniteVolumePolicy( const FiniteVolumePolicy& p ) : fvptr_(p.fvptr_) {}
    FiniteVolumePolicy& operator=( const FiniteVolumePolicy& p ) { if ( this != &p ) fvptr_ = p.fvptr_; return *this; }
    /// move semantics
    FiniteVolumePolicy( FiniteVolumePolicy&& p ) noexcept : fvptr_(p.fvptr_) { p.fvptr_ = nullptr; }
    FiniteVolumePolicy& operator=( FiniteVolumePolicy&& p ) noexcept { std::swap(fvptr_,p.fvptr_); p.fvptr_ = nullptr; return *this; }

    /// connect policy to a finite-volume stencil that must match the parent element
    void    AssignFiniteVolume( const csmp::FiniteVolumeStencil<3U>* );
    void    AssignFiniteVolumeNullPtr() { fvptr_ = nullptr; }

    /// access the stencil functionality directly
    const   FiniteVolumeStencil<3U>* const FV() const;
 
    /// the surface patches that constitute the outside walls of the finite volume
    uint32_t Facets()  const;
  
    /// the volumetric partitions of the finite elements that belong to the different node-centered finite volumes
    uint32_t Sectors() const;
  
    /// volume quadrature points: usually one, but there may be multiple
    uint32_t IntegrationPointsPerSector() const;

    /// surface quadrature points: usually one, but there may be multiple
    uint32_t IntegrationPointsPerFacet()  const;

   /// finite-element interpolation function values Ni output to NRST vector (stored by the current finite element) at the numbered facet integration point
    void   N_AtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip ) const;
    void   N_AtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip, std::vector<double>& NRST ) const;
    void   N_AtSectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const;
  
    /// interpolation function derivatives at point specified in local coordinates; result is returned into NRST vector of current finite element
    void   Local_dN_At( const Point<3U>& rst ) const;
  
    /// interpolation function derivatives at point specified in local coordinates returned into NRST vector of current finite element + detJ for integration
    double dN_At( const Point<3U>& rst, DenseMatrix<DM_MIN>& DN )  const;

    /// interpolates value of scalar node variable to facet integration point (XYZ)
    double PropertyValueAtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip, const csmp::Index& ) const;

    /// interpolates value of scalar node variable to sector integration point (XYZ)
    double PropertyValueAtSectorIntegrationPoint( uint32_t iSector, uint32_t ip, const csmp::Index& ) const;
  
    /// interpolate values of any node variable to facet integration point (XYZ)
    template<class Var> requires CsmpVariable<3U, Var>
    void   PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, uint32_t iFacet,  uint32_t ip,  Var& )   const;

    /// interpolate values of any node variable to sector integration point (XYZ)
    template<class Var> requires CsmpVariable<3U, Var>
    void   PropertyValueAtSectorIntegrationPoint( const csmp::Index&, uint32_t iSector, uint32_t ip,  Var& )   const;

    /// integrates value of property over the area of the facet
    double FacetIntegral(  uint32_t iFacet,  const csmp::Index& ) const;

    /// integrates value of property over the sector
    double SectorIntegral( uint32_t iSector, const csmp::Index& ) const;

    /// reads value of vector property and projects it onto facet normal in physical space; returns projected value
    double ProjectionOnFacetNormal( uint32_t iFacet, const csmp::Index& )  const;

    /// projects value of vector property onto facet normal in physical space; returns projection
    double ProjectionOnFacetNormal( uint32_t iFacet, const VectorVariable<3U>& ) const;

    /// returns subvolume of finite element that corresponds to the requested finite volume sector
    double SectorVolume( uint32_t iSector )     const;
  
    /// returns the area of the finite-volume facet in physical space as obtained by construction of facet in physical space
    double FacetArea( uint32_t iFacet )         const;
  
    /// returns the area of the finite volume facet in physical space as obtained by Jacobian transformation of area in parametric space
    double FacetAreaMapped( uint32_t iFacet )   const;

    /// returns the unit normal to the finite-volume facet
    Point<3U> FacetNormal( uint32_t iFacet )       const;
  
    /// transforms the facet normal from parametric to physical space and optionally normalises it to obtain a unit length
    Point<3U> FacetNormalMapped( uint32_t iFacet ) const;

    /// returns the facet area in parametric space
    double ParametricFacetArea( uint32_t iFacet ) const;
  
    /// returns the normal to the facet in parametric space
    Point<3U> ParametricFacetNormal( uint32_t iFacet ) const;
    
    /// returns facet corner point in physical coordinates using RST_to_XYX
    Point<3U> FacetPoint( uint32_t iFacet, uint32_t iPoint ) const;

    friend class Element<3U>;
    friend class Face<3U>;
    friend class InterFace<3U>;

    // get finite-volume facet area, normal and sector volume in parametric space from FiniteVolumeStencil class
private:
    FiniteVolumePolicy( const CELL<3U>& );
    const csmp::FiniteVolumeStencil<3U>*  fvptr_;
};


/// full specialization ( 2D )
template<template<uint32_t> class CELL>
class FiniteVolumePolicy<2U,CELL> {
  public:
    FiniteVolumePolicy( const csmp::FiniteVolumeStencil<2U>* = nullptr );
    FiniteVolumePolicy( const FiniteVolumePolicy& p ) : fvptr_(p.fvptr_) {}
    void       AssignFiniteVolume( const csmp::FiniteVolumeStencil<2U>* );
    void       AssignFiniteVolumeNullPtr() { fvptr_ = nullptr; }
    const FiniteVolumeStencil<2U>* const FV() const;
    uint32_t     Facets()  const;
    uint32_t     Sectors() const;
    uint32_t     IntegrationPointsPerSector() const;
    uint32_t     IntegrationPointsPerFacet()  const;
    void       N_AtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip ) const;
    void       N_AtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip, std::vector<double>& NRST ) const;
    void       N_AtSectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const;
    void       Local_dN_At( const Point<2U>& rst ) const;
    double   dN_At( const Point<2U>& rst, DenseMatrix<DM_MIN>& DN )  const;
    double   PropertyValueAtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip, const csmp::Index& prop_key ) const;
    double   PropertyValueAtSectorIntegrationPoint( uint32_t iSector, uint32_t ip, const csmp::Index& prop_key ) const;
    template<class Var> requires CsmpVariable<2U, Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, uint32_t iFacet, uint32_t ip,  Var& )    const;
    template<class Var> requires CsmpVariable<2U, Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, uint32_t iSector, uint32_t ip, Var& )    const;
    double   FacetIntegral(  uint32_t iFacet,  const csmp::Index& prop_key ) const;
    double   SectorIntegral( uint32_t iSector, const csmp::Index& prop_key ) const;
    double   ProjectionOnFacetNormal( uint32_t iFacet, const csmp::Index& prop_key )  const;
    double   ProjectionOnFacetNormal( uint32_t iFacet, const VectorVariable<2U>& vc ) const;
    double   SectorVolume( uint32_t iSector )     const;
    double   FacetArea( uint32_t iFacet )         const;
    double   FacetAreaMapped( uint32_t iFacet )   const;
    Point<2U>  FacetNormal( uint32_t iFacet )       const;
    Point<2U>  FacetNormalMapped( uint32_t iFacet ) const;
    double   ParametricFacetArea( uint32_t iFacet ) const;
    Point<2U>  ParametricFacetNormal( uint32_t iFacet ) const;

    friend class Element<2U>;
    friend class Face<2U>;
    friend class InterFace<2U>;

private:
    FiniteVolumePolicy( const CELL<2U>& );
    const csmp::FiniteVolumeStencil<2U>*  fvptr_;
};



/// full specialization ( 1D )
template<template<uint32_t> class CELL>
class FiniteVolumePolicy<1U, CELL> {
  public:
    FiniteVolumePolicy( const csmp::FiniteVolumeStencil<1U>* = nullptr );
    FiniteVolumePolicy( const FiniteVolumePolicy& p ) : fvptr_(p.fvptr_) {}
    void     AssignFiniteVolume( const csmp::FiniteVolumeStencil<1U>* );
    void     AssignFiniteVolumeNullPtr() { fvptr_ = nullptr; }
    const FiniteVolumeStencil<1U>* const FV() const;
    uint32_t   Facets()  const;
    uint32_t   Sectors() const;
    uint32_t   IntegrationPointsPerSector() const;
    uint32_t   IntegrationPointsPerFacet()  const;
    void     N_AtFacetIntegrationPoint( uint32_t iFacet, uint32_t ip )   const;
    void     N_AtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip, std::vector<double>& NRST ) const;
    void     N_AtSectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const;
    void     Local_dN_At( const Point<1U>& rst ) const;
    double   dN_At( const Point<1U>& rst, DenseMatrix<DM_MIN>& DN )  const;
    double   PropertyValueAtFacetIntegrationPoint(  uint32_t iFacet,  uint32_t ip, const csmp::Index& prop_key ) const;
    double   PropertyValueAtSectorIntegrationPoint( uint32_t iSector, uint32_t ip, const csmp::Index& prop_key ) const;
    template<class Var> requires CsmpVariable<1U, Var>
    void     PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, uint32_t iFacet, uint32_t ip,  Var& )    const;
    template<class Var> requires CsmpVariable<1U, Var>
    void     PropertyValueAtSectorIntegrationPoint( const csmp::Index&, uint32_t iSector, uint32_t ip,  Var& )   const;
    double   FacetIntegral(  uint32_t iFacet,  const csmp::Index& prop_key ) const;
    double   SectorIntegral( uint32_t iSector, const csmp::Index& prop_key ) const;
    double   ProjectionOnFacetNormal( uint32_t iFacet, const csmp::Index& prop_key )  const;
    double   ProjectionOnFacetNormal( uint32_t iFacet, const VectorVariable<1U>& vc ) const;
    double   SectorVolume( uint32_t iSector )     const;
    double   FacetArea( uint32_t iFacet )         const;
    double   FacetAreaMapped( uint32_t iFacet )   const;
    Point<1U>  FacetNormal( uint32_t iFacet )       const;
    Point<1U>  FacetNormalMapped( uint32_t iFacet ) const;
    double     ParametricFacetArea( uint32_t iFacet ) const;
    Point<1U>  ParametricFacetNormal( uint32_t iFacet ) const;

    friend class Element<1U>;
    friend class Face<1U>;
    friend class InterFace<1U>;

private:
    FiniteVolumePolicy( const CELL<1U>& );
    const csmp::FiniteVolumeStencil<1U>*  fvptr_ = nullptr;
};

} // end namespace 

#endif /* FINITE_VOLUME_POLICY_H */
