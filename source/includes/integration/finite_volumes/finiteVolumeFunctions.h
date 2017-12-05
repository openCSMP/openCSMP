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

    // Gradient of a scalar node-based property, at a point in parametric space
    Point<dim> GradientOfScalarNodeProperty( const Point<dim>& p, const csmp::Index& prop ) const;

    // Interpolation of a scalar node-based property, at a point in parametric space
    double64 InterpolateScalarNodeProperty( const Point<dim>& p, const csmp::Index& prop ) const;

    Point<dim> NormalOfFacet( size_t iFacet ) const;

private:
    void CalculateDN(const Point<dim>& p, std::vector<double64>* dn) const;
    void CalculateN(const Point<dim>& p, std::vector<double64>& n) const;

    Element<dim>* eptr_;
    size_t element_dim_;
    size_t num_nodes_;
    
};

} // csmp

#endif /* defined(FINITE_VOLUME_UNIVERSAL_FUNCTIONS_H) */
