//
//  CoordinateTransformer.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 11/17/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include <algorithm>
#include <cassert>
#include "CoordinateTransformer.h"

using namespace std;

namespace csmp {

/// initialises translations to zero and flip vector to 1
template<uint32_t dim>
CoordinateTransformer<dim>::CoordinateTransformer()
 {
    for ( size_t i=0; i<dim; i++ ) {
         translation_[i] = 0.; ///< in x,y,z direction
         flip_[i]        = 1.; ///< coordinate multiplier
      }
 }



/** 
    processes translations etc. returning the point
*/
template<uint32_t dim>
void CoordinateTransformer<dim>::Transform( Point<dim>& p ) const
 {
    for ( auto i{0}; i<dim; i++ ) {
         p[i] += translation_[i];
         p[i] *= flip_[i];
      }
    if ( !swapped_axes_.empty() )
      for ( list<pair<size_t,size_t> >::const_iterator
            it=swapped_axes_.begin(); it!=swapped_axes_.end(); it++ )
        {
           double swap       = p[ (*it).second ];
           p[ (*it).second ] = p[ (*it).first ];
           p[ (*it).first ]  = swap;
        }
 }
 
 
 
/** 
    processes translations etc. returning the point
*/
template<uint32_t dim>
Point<dim> CoordinateTransformer<dim>::Transform( const vector<double>& vec ) const
 {
    Point<dim>  p(vec);
    for ( auto i{0}; i<dim; i++ ) {
         p[i] += translation_[i];
         p[i] *= flip_[i];
      }
    if ( !swapped_axes_.empty() )
      for ( list<pair<size_t,size_t> >::const_iterator
            it=swapped_axes_.begin(); it!=swapped_axes_.end(); it++ )
        {
           double swap       = p[ (*it).second ];
           p[ (*it).second ] = p[ (*it).first ];
           p[ (*it).first ]  = swap;
        }
   
    return p;
 }
  
  
/*
TODO: add center of rotation into class so that rotations can be processed

/// clockwise rotation in degrees from 0..360
template<uint32_t dim>
void CoordinateTransformer<dim>::Rotate( size_t axis, double angle )
 {
    assert( axis < dim );
    rotation_[axis] = angle;
 }
*/



template<uint32_t dim>
void CoordinateTransformer<dim>::Translate( size_t coord, double distance_meters )
 {
    assert( coord < dim );
    translation_[coord] = distance_meters;
 }




template<uint32_t dim>
void CoordinateTransformer<dim>::FlipAxis( size_t coord )
 {
    assert( coord < dim );
    flip_[coord] *= -1.;
 }
  


/// classical swap
template<uint32_t dim>
void CoordinateTransformer<dim>::ExchangeAxes( size_t axis_a, size_t axis_b )
 {
    // avoiding duplicates
    if ( find( swapped_axes_.begin(), swapped_axes_.end(), make_pair(axis_a,axis_b) ) != swapped_axes_.end() ) return;
    swapped_axes_.push_back( make_pair(axis_a,axis_b) );
 }




/// sets to preserve input coordinate
template<uint32_t dim>
void CoordinateTransformer<dim>::Reset()
 {
    for ( auto i{0}; i<dim; i++ ) {
         translation_[i] = 0.;
         flip_[i] = 1.;
      }
 }
  

template class CoordinateTransformer<1U>;
template class CoordinateTransformer<2U>;
template class CoordinateTransformer<3U>;

} // end csmp
