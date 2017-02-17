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
template<size_t, template<size_t> class SIMPLEX> class FiniteVolumePolicy;

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
 
Provides a set of functionality to perform operations, used in the CVSIMPLEX 
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
template<template<size_t> class SIMPLEX>
class FiniteVolumePolicy<3U, SIMPLEX> {
  public:
    FiniteVolumePolicy( const csmp::FiniteVolumeStencil<3U>* = nullptr );
    FiniteVolumePolicy( const FiniteVolumePolicy& p ) : fvptr_(p.fvptr_) {}

    /// connect policy to a finite-volume stencil that must match the parent element
    void       AssignFiniteVolume( const csmp::FiniteVolumeStencil<3U>* );

    /// access the stencil functionality directly
    const FiniteVolumeStencil<3U>* const FV_Stencil() const;
 
    /// the surface patches that constitute the outside walls of the finite volume
    size_t     Facets()  const;
  
    /// the volumetric partitions of the finite elements that belong to the different node-centered finite volumes
    size_t     Sectors() const;
  
    /// volume quadrature points: usually one, but there may be multiple
    size_t     IntegrationPointsPerSector() const;

    /// surface quadrature points: usually one, but there may be multiple
    size_t     IntegrationPointsPerFacet()  const;

    /// mapping of integration points from local to global coordinates
    Point<3U>  RstToXYZ( const Point<3U>& rst ) const;

    /// shape function values at sector and facet integration points
    void       N_At( const Point<3U>& rst ) const;
    void       N_At( const Point<3U>& rst, std::vector<double64>& N )  const;
  
    /// finite-element interpolation function values Ni output to NRST vector (stored by the current finite element) at the numbered facet integration point
    void       N_AtFacetIntegrationPoint(  size_t iFacet,  size_t ip ) const;
    void       N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const;
  
    /// interpolation function derivatives at point specified in local coordinates; result is returned into NRST vector of current finite element
    void       Local_dN_At( const Point<3U>& rst ) const;
  
    /// interpolation function derivatives at point specified in local coordinates returned into NRST vector of current finite element + detJ for integration
    double64   dN_At( const Point<3U>& rst, DenseMatrix<DM_MIN>& DN )  const;

    /// interpolates value of scalar node variable to facet integration point (XYZ)
    double64   PropertyValueAtFacetIntegrationPoint(  size_t iFacet,  size_t ip, const csmp::Index& ) const;

    /// interpolates value of scalar node variable to sector integration point (XYZ)
    double64   PropertyValueAtSectorIntegrationPoint( size_t iSector, size_t ip, const csmp::Index& ) const;
  
    /// interpolate values of any node variable to facet integration point (XYZ)
    template<class Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, size_t iFacet,  size_t ip,  Var& )   const;

    /// interpolate values of any node variable to sector integration point (XYZ)
    template<class Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, size_t iSector, size_t ip,  Var& )   const;

    /// integrates value of property over the area of the facet
    double64   FacetIntegral(  size_t iFacet,  const csmp::Index& ) const;

    /// integrates value of property over the sector
    double64   SectorIntegral( size_t iSector, const csmp::Index& ) const;

    /// reads value of vector property and projects it onto facet normal in physical space; returns projected value
    double64   ProjectionOnFacetNormal( size_t iFacet, const csmp::Index& )  const;

    /// projects value of vector property onto facet normal in physical space; returns projection
    double64   ProjectionOnFacetNormal( size_t iFacet, const VectorVariable<3U>& ) const;

    /// returns subvolume of finite element that corresponds to the requested finite volume sector
    double64   SectorVolume( size_t iSector )     const;
  
    /// returns the area of the finite-volume facet in physical space as obtained by construction of facet in physical space
    double64   FacetArea( size_t iFacet )         const;
  
    /// returns the area of the finite volume facet in physical space as obtained by Jacobian transformation of area in parametric space
    double64   FacetAreaMapped( size_t iFacet )   const;

    /// returns the unit normal to the finite-volume facet
    Point<3U>  FacetNormal( size_t iFacet )       const;
  
    /// transforms the facet normal from parametric to physical space and normalises it to obtain a unit length
    Point<3U>  FacetNormalMapped( size_t iFacet ) const;

    /// returns the facet area in parametric space
    double64   ParametricFacetArea( size_t iFacet ) const;
  
    /// returns the normal to the facet in parametric space
    Point<3U>  ParametricFacetNormal( size_t iFacet ) const;

    // get finite-volume facet area, normal and sector volume in parametric space from FiniteVolumeStencil class
private:
    FiniteVolumePolicy( const SIMPLEX<3U>& );
    const csmp::FiniteVolumeStencil<3U>*  fvptr_;
};


/// full specialization ( 2D )
template<template<size_t> class SIMPLEX>
class FiniteVolumePolicy<2U,SIMPLEX> {
  public:
    FiniteVolumePolicy( const csmp::FiniteVolumeStencil<2U>* = nullptr );
    FiniteVolumePolicy( const FiniteVolumePolicy& p ) : fvptr_(p.fvptr_) {}
    void       AssignFiniteVolume( const csmp::FiniteVolumeStencil<2U>* );
    const FiniteVolumeStencil<2U>* const FV_Stencil() const;
    size_t     Facets()  const;
    size_t     Sectors() const;
    size_t     IntegrationPointsPerSector() const;
    size_t     IntegrationPointsPerFacet()  const;
    Point<2U>  RstToXYZ( const Point<2U>& rst ) const;
    void       N_At( const Point<2U>& rst ) const;
    void       N_At( const Point<2U>& rst, std::vector<double64>& N )  const;
    void       N_AtFacetIntegrationPoint(  size_t iFacet,  size_t ip ) const;
    void       N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const;
    void       Local_dN_At( const Point<2U>& rst ) const;
    double64   dN_At( const Point<2U>& rst, DenseMatrix<DM_MIN>& DN )  const;
    double64   PropertyValueAtFacetIntegrationPoint(  size_t iFacet,  size_t ip, const csmp::Index& prop_key ) const;
    double64   PropertyValueAtSectorIntegrationPoint( size_t iSector, size_t ip, const csmp::Index& prop_key ) const;
    template<class Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, size_t iFacet, size_t ip,  Var& )    const;
    template<class Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, size_t iSector, size_t ip, Var& )    const;
    double64   FacetIntegral(  size_t iFacet,  const csmp::Index& prop_key ) const;
    double64   SectorIntegral( size_t iSector, const csmp::Index& prop_key ) const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const csmp::Index& prop_key )  const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const VectorVariable<2U>& vc ) const;
    double64   SectorVolume( size_t iSector )     const;
    double64   FacetArea( size_t iFacet )         const;
    double64   FacetAreaMapped( size_t iFacet )   const;
    Point<2U>  FacetNormal( size_t iFacet )       const;
    Point<2U>  FacetNormalMapped( size_t iFacet ) const;
    double64   ParametricFacetArea( size_t iFacet ) const;
    Point<2U>  ParametricFacetNormal( size_t iFacet ) const;

private:
    FiniteVolumePolicy( const SIMPLEX<2U>& );
    const csmp::FiniteVolumeStencil<2U>*  fvptr_;
};



/// full specialization ( 1D )
template<template<size_t> class SIMPLEX>
class FiniteVolumePolicy<1U, SIMPLEX> {
  public:
    FiniteVolumePolicy( const csmp::FiniteVolumeStencil<1U>* = nullptr );
    FiniteVolumePolicy( const FiniteVolumePolicy& p ) : fvptr_(p.fvptr_) {}
    void       AssignFiniteVolume( const csmp::FiniteVolumeStencil<1U>* );
    const FiniteVolumeStencil<1U>* const FV_Stencil() const;
    size_t     Facets()  const;
    size_t     Sectors() const;
    size_t     IntegrationPointsPerSector() const;
    size_t     IntegrationPointsPerFacet()  const;
    Point<1U>  RstToXYZ( const Point<1U>& rst ) const;
    void       N_At( const Point<1U>& rst ) const;
    void       N_At( const Point<1U>& rst, std::vector<double64>& N )  const;
    void       N_AtFacetIntegrationPoint( size_t iFacet, size_t ip )   const;
    void       N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const;
    void       Local_dN_At( const Point<1U>& rst ) const;
    double64   dN_At( const Point<1U>& rst, DenseMatrix<DM_MIN>& DN )  const;
    double64   PropertyValueAtFacetIntegrationPoint(  size_t iFacet,  size_t ip, const csmp::Index& prop_key ) const;
    double64   PropertyValueAtSectorIntegrationPoint( size_t iSector, size_t ip, const csmp::Index& prop_key ) const;
    template<class Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, size_t iFacet, size_t ip,  Var& )    const;
    template<class Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, size_t iSector, size_t ip,  Var& )   const;
    double64   FacetIntegral(  size_t iFacet,  const csmp::Index& prop_key ) const;
    double64   SectorIntegral( size_t iSector, const csmp::Index& prop_key ) const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const csmp::Index& prop_key )  const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const VectorVariable<1U>& vc ) const;
    double64   SectorVolume( size_t iSector )     const;
    double64   FacetArea( size_t iFacet )         const;
    double64   FacetAreaMapped( size_t iFacet )   const;
    Point<1U>  FacetNormal( size_t iFacet )       const;
    Point<1U>  FacetNormalMapped( size_t iFacet ) const;
    double64   ParametricFacetArea( size_t iFacet ) const;
    Point<1U>  ParametricFacetNormal( size_t iFacet ) const;

private:
    FiniteVolumePolicy( const SIMPLEX<1U>& );
    const csmp::FiniteVolumeStencil<1U>*  fvptr_;
};

} // end namespace 

#endif /* FINITE_VOLUME_POLICY_H */
