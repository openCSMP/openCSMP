//
//  CoordinateTransformer.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 11/17/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef COORDINATE_TRANSFORMER_H
#define COORDINATE_TRANSFORMER_H

#include <list>
#include "Point.h"

namespace csmp {

/// utility to convert between coordinate systems and flip axes
template<size_t dim>
class CoordinateTransformer {
  public:
    CoordinateTransformer();
  
    /// transform point coordinates using all specified operations, @todo make sure that specified sequence is honored
    void Transform( Point<dim>& ) const;
    Point<dim> Transform( const std::vector<double>& ) const;
  
    // Set the transformation operations for later processing
    /// translation along the specified axis
    void Translate( size_t axis, double distance_meters );
  
    /// multiplies coordinates of specified axis by -1
    void FlipAxis( size_t axis );

    /// exchange of axis
    void ExchangeAxes( size_t axis_a, size_t axis_b );
  
    /// returns settings to ones which preserve value of input coordinate
    void Reset();
  
  private:
    double translation_[dim]; ///< in x,y,z direction
    double flip_[dim];        ///< inverts the corresponding coordinate axis
    std::list<std::pair<size_t,size_t> > swapped_axes_;   ///< as the last step
 };

} // end csmp

#endif /* defined(COORDINATE_TRANSFORMER_H) */
