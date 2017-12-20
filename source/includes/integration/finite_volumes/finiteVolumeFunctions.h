//
//  finiteVolumeFunctions.h
//
//  Created by Stephan Matthai on 9/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#ifndef FINITE_VOLUME_FUNCTIONS_H
#define FINITE_VOLUME_FUNCTIONS_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

#include <bitset>

namespace csmp {

template<size_t> class Model;
template<size_t> class Region;
template<size_t> class Element;
struct Index;

/// sector volume, finite volume, FV pore volume
template<size_t dim> void initializeFiniteVolumeProperties( Model<dim>&, Region<dim>&  );



template<size_t dim, VARIABLE_TYPE vt>
struct VariableTypeTraits
{
};


template<size_t dim>
struct VariableTypeTraits<dim,SCALAR>
{
    typedef ScalarVariable VariableType;
};

template<size_t dim>
struct VariableTypeTraits<dim,VECTOR>
{
    typedef VectorVariable<dim> VariableType;
};

template<size_t dim>
struct VariableTypeTraits<dim,TENSOR>
{
    typedef TensorVariable<dim> VariableType;
};

template<size_t dim>
struct VariableTypeTraits<dim,ARRAY>
{
    typedef ArrayVariable VariableType;
};

template<size_t dim>
struct VariableTypeTraits<dim,FLAGGEDARRAY>
{
    typedef FlaggedArrayVariable VariableType;
};


template<size_t dim>
class FiniteElementHelper
{
public:
    FiniteElementHelper( );
  
    Element<dim>* FiniteElement( );
    void FiniteElement( Element<dim>* eptr );
    ~FiniteElementHelper();

    // Normal of a facet, scaled by facet area
    Point<dim> NormalOfFacet( size_t iFacet );

    // Gradient of a scalar node-based property
    Point<dim> ReadGradientAtBarycenter( const csmp::INDEX<SCALAR,NODE>& prop );

    // Read scalar properties
    template<VARIABLE_TYPE ty, PLACEMENT pl>
    void ReadAtBarycenter( const csmp::INDEX<ty,pl>& prop, typename VariableTypeTraits<dim,ty>::VariableType& var );

    template<VARIABLE_TYPE ty, PLACEMENT pl>
    void ReadAtNode( const csmp::INDEX<ty,pl>& prop, size_t n, typename VariableTypeTraits<dim,ty>::VariableType& var );

    template<VARIABLE_TYPE ty, PLACEMENT pl>
    void ReadAtElementIntegrationPoint( const csmp::INDEX<ty,pl>& prop, size_t ip, typename VariableTypeTraits<dim,ty>::VariableType& var );

    template<VARIABLE_TYPE ty, PLACEMENT pl>
    void ReadAtFacetIntegrationPoint( const csmp::INDEX<ty,pl>& prop, size_t facet, size_t ip, typename VariableTypeTraits<dim,ty>::VariableType& var );

    template<VARIABLE_TYPE ty, PLACEMENT pl>
    void ReadAtSectorIntegrationPoint( const csmp::INDEX<ty,pl>& prop, size_t sector, size_t ip, typename VariableTypeTraits<dim,ty>::VariableType& var );

private:
    void CalculateDN(const Point<dim>& p, std::vector<double64>* DN);

    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // csmp

#endif /* defined(FINITE_VOLUME_UNIVERSAL_FUNCTIONS_H) */
