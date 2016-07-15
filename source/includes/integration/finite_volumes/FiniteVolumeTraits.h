#ifndef FINITE_VOLUME_TRAITS_H
#define FINITE_VOLUME_TRAITS_H

#include "DenseMatrix.h"
#include "Index.h"

#include "Point.h"
#include "FiniteVolumeStencil.h"

#include "Node.h"

#include "QuadrilateralFacet.h"
#include "TriangularFacet.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/// Element policy for node-centered finite-volume computations; see full specializations for 1D, 2D, 3D
template<size_t, template<size_t> class SIMPLEX> class FiniteVolumeTraits;

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
basis. Distinctive feature of the FiniteVolumeTraits class is the 
possibility to access 
 
 
@section examples Application Examples

No specific application examples are required, as the method is multiply  
called the great code CSP, within the implementation of the classes 
ExplicitFiniteVolumeProcessor and others.  

TODO: @todo (3) Put finiteVolumeAuxiliaryFunctions into FiniteVolumeTraits class (A)
TODO: @todo SKM try to do computations in parametric space

 */
template<template<size_t> class SIMPLEX>
class FiniteVolumeTraits<3U, SIMPLEX> {

  public:

    FiniteVolumeTraits();

    // FV Stencil info
    size_t     Facets()  const;
    size_t     Sectors() const;
    size_t     IntegrationPointsPerSector() const;
    size_t     IntegrationPointsPerFacet()  const;

    // mapping of integration point from local to global coordinates
    Point<3U>  RstToXYZ( const Point<3U>& rst ) const;

    // shape function values and their derivatives at sector and facet integration points
    void       N_At( const Point<3U>& rst ) const;
    void       N_At( const Point<3U>& rst, std::vector<double64>& N )  const;
    // initialize vector NRST (stored by the current finite element)
    void       N_AtFacetIntegrationPoint(  size_t iFacet,  size_t ip ) const;
    void       N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const;
    // shape function derivatives in local coordinates
    void       Local_dN_At( const Point<3U>& rst ) const;
    // shape function derivatives in global coordinates + detJ for integration
    double64   dN_At( const Point<3U>& rst, DenseMatrix<DM_MIN>& DN )  const;

    // interpolate value of scalar node variable to facet/sector integration point (XYZ)
    double64   PropertyValueAtFacetIntegrationPoint(  size_t iFacet,  size_t ip, const csmp::Index& prop_key ) const;
    double64   PropertyValueAtSectorIntegrationPoint( size_t iSector, size_t ip, const csmp::Index& prop_key ) const;
    // interpolate value of any node variable to facet/sector integration point (XYZ)
    template<class Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, size_t iFacet,  size_t ip,  Var& )   const;
    template<class Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, size_t iSector, size_t ip,  Var& )   const;

    // facet and sector integrals
    double64   FacetIntegral(  size_t iFacet,  const csmp::Index& prop_key ) const;
    double64   SectorIntegral( size_t iSector, const csmp::Index& prop_key ) const;

    // projection of vector property onto facet normal in physical space (returns projection and determinant into arg4) (XYZ)
    double64   ProjectionOnFacetNormal( size_t iFacet, const csmp::Index& prop_key )  const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const VectorVariable<3U>& vc ) const;

    // finite volume properties in physical space (XYZ)
    double64   SectorVolume( size_t iSector )     const;
    double64   FacetArea( size_t iFacet )         const;
    Point<3U>  FacetNormal( size_t iFacet )       const;
    double64   FacetAreaMapped( size_t iFacet )   const;
    Point<3U>  FacetNormalMapped( size_t iFacet ) const;

    // finite volume properties in the parametric space of the reference finite element (RST)
    // TODO: refactor so that these methods actually map the properties from parametric space to physical space
    double64   ParametricFacetArea( size_t iFacet ) const;
    const Point<3U>&  ParametricFacetNormal( size_t iFacet ) const;

private:

    FiniteVolumeTraits( const SIMPLEX<3U>& );

};


/// partial specialization ( 2D )

template<template<size_t> class SIMPLEX>
class FiniteVolumeTraits<2U,SIMPLEX> {

  public:

    FiniteVolumeTraits();

    // FV Stencil info
    size_t     Facets()  const;
    size_t     Sectors() const;
    size_t     IntegrationPointsPerSector() const;
    size_t     IntegrationPointsPerFacet()  const;

    // mapping of integration point from local to global coordinates
    Point<2U>  RstToXYZ( const Point<2U>& rst ) const;

    // shape function values and their derivatives at sector and facet integration points
    void       N_At( const Point<2U>& rst ) const;
    void       N_At( const Point<2U>& rst, std::vector<double64>& N )  const;
    // initialize vector NRST (stored by the current finite element)
    void       N_AtFacetIntegrationPoint(  size_t iFacet,  size_t ip ) const;
    void       N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const;
    // shape function derivatives in local coordinates
    void       Local_dN_At( const Point<2U>& rst ) const;
    // shape function derivatives in global coordinates + detJ for integration
    double64   dN_At( const Point<2U>& rst, DenseMatrix<DM_MIN>& DN )  const;


    // interpolate value of scalar node variable to facet/sector integration point (XYZ)
    double64   PropertyValueAtFacetIntegrationPoint(  size_t iFacet,  size_t ip, const csmp::Index& prop_key ) const;
    double64   PropertyValueAtSectorIntegrationPoint( size_t iSector, size_t ip, const csmp::Index& prop_key ) const;
    // interpolate value of any node variable to facet/sector integration point (XYZ)
    template<class Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, size_t iFacet, size_t ip,  Var& )    const;
    template<class Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, size_t iSector, size_t ip, Var& )    const;

    // facte and sector integrals
    double64   FacetIntegral(  size_t iFacet,  const csmp::Index& prop_key ) const;
    double64   SectorIntegral( size_t iSector, const csmp::Index& prop_key ) const;

    // projection of vector property onto facet normal in physical space (returns projection and determinant into arg4) (XYZ)
    double64   ProjectionOnFacetNormal( size_t iFacet, const csmp::Index& prop_key )  const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const VectorVariable<2U>& vc ) const;

    // finite volume properties in physical space (XYZ)
    double64   SectorVolume( size_t iSector )     const;
    double64   FacetArea( size_t iFacet )         const;
    double64   FacetAreaMapped( size_t iFacet )   const;
    Point<2U>  FacetNormal( size_t iFacet )       const;
    Point<2U>  FacetNormalMapped( size_t iFacet ) const;

    // finite volume properties in the parametric space of the reference finite element (RST)
    double64   ParametricFacetArea( size_t iFacet ) const;
    const Point<2U>&  ParametricFacetNormal( size_t iFacet ) const;

private:

    FiniteVolumeTraits( const SIMPLEX<2U>& );

};

/// partial specialization ( 1D )

template<template<size_t> class SIMPLEX>
class FiniteVolumeTraits<1U, SIMPLEX> {

  public:

    FiniteVolumeTraits();

    // FV Stencil info
    size_t     Facets()  const;
    size_t     Sectors() const;
    size_t     IntegrationPointsPerSector() const;
    size_t     IntegrationPointsPerFacet()  const;

    // mapping of integration point from local to global coordinates
    Point<1U>  RstToXYZ( const Point<1U>& rst ) const;

    // shape function values and their derivatives at sector and facet integration points
    void       N_At( const Point<1U>& rst ) const;
    void       N_At( const Point<1U>& rst, std::vector<double64>& N )  const;
    // initialize vector NRST (stored by the current finite element)
    void       N_AtFacetIntegrationPoint( size_t iFacet, size_t ip )   const;
    void       N_AtSectorIntegrationPoint( size_t iSector, size_t ip ) const;
    // shape function derivatives in local coordinates
    void       Local_dN_At( const Point<1U>& rst ) const;
    // shape function derivatives in global coordinates + detJ for integration
    double64   dN_At( const Point<1U>& rst, DenseMatrix<DM_MIN>& DN )  const;

    // interpolate value of scalar node variable to facet/sector integration point (XYZ)
    double64   PropertyValueAtFacetIntegrationPoint(  size_t iFacet,  size_t ip, const csmp::Index& prop_key ) const;
    double64   PropertyValueAtSectorIntegrationPoint( size_t iSector, size_t ip, const csmp::Index& prop_key ) const;
    // interpolate value of any node variable to facet/sector integration point (XYZ)
    template<class Var>
    void       PropertyValueAtFacetIntegrationPoint(  const csmp::Index&, size_t iFacet, size_t ip,  Var& )    const;
    template<class Var>
    void       PropertyValueAtSectorIntegrationPoint( const csmp::Index&, size_t iSector, size_t ip,  Var& )   const;

    // facet and sector integrals
    double64   FacetIntegral(  size_t iFacet,  const csmp::Index& prop_key ) const;
    double64   SectorIntegral( size_t iSector, const csmp::Index& prop_key ) const;

    // projection of vector property onto facet normal in physical space (returns projection and determinant into arg4) (XYZ)
    double64   ProjectionOnFacetNormal( size_t iFacet, const csmp::Index& prop_key )  const;
    double64   ProjectionOnFacetNormal( size_t iFacet, const VectorVariable<1U>& vc ) const;

    // finite volume properties in physical space (XYZ)
    double64   SectorVolume( size_t iSector )     const;
    double64   FacetArea( size_t iFacet )         const;
    double64   FacetAreaMapped( size_t iFacet )   const;
    Point<1U>  FacetNormal( size_t iFacet )       const;
    Point<1U>  FacetNormalMapped( size_t iFacet ) const;

    // finite volume properties in the parametric space of the reference finite element (RST)
    double64   ParametricFacetArea( size_t iFacet )         const;
    const Point<1U>  ParametricFacetNormal( size_t iFacet ) const;

private:

    FiniteVolumeTraits( const SIMPLEX<1U>& );

};

} // end namespace 

#endif
