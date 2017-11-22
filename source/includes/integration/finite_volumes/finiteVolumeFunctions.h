//
//  finiteVolumeFunctions.h
//
//  Created by Stephan Matthai on 9/08/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "Point.h"

#ifndef FINITE_VOLUME_FUNCTIONS_H
#define FINITE_VOLUME_FUNCTIONS_H

namespace csmp {

template<size_t> class Model;
template<size_t> class Region;
template<size_t> class Element;
struct Index;

/// sector volume, finite volume, FV pore volume
template<size_t dim> void initializeFiniteVolumeProperties( Model<dim>&, Region<dim>&  );

template<size_t dim>
class FiniteVolumeHelper
{
public:
    FiniteVolumeHelper( Element<dim>* eptr );

    void SetParametricCoordinate( const Point<dim>& p );

    Point<dim> GradientOfScalarNodeProperty( const csmp::Index& prop ) const;

    Point<dim> NormalOfFacet( size_t iFacet ) const;

private:
    void CalculateDN(const Point<dim>& p);
    
    Element<dim>* eptr_;
    size_t element_dim_;
    size_t num_nodes_;
    std::vector<double64> DN_[dim];
    
};

} // csmp

#endif /* defined(FINITE_VOLUME_UNIVERSAL_FUNCTIONS_H) */
