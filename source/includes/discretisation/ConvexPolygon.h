// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  ConvexPolygon3D.h
//
//  Created by Stephan Matthai on 13/10/2024.
//

#ifndef CONVEX_POLYGON_3D_H
#define CONVEX_POLYGON_3D_H

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "Point.h"

namespace csmp {

template <uint32_t dim>
class ConvexPolygon {
public:
    // Constructors
    ConvexPolygon( const std::vector<Point<dim>>& );
    
    // Constructor: accepts iterators to a range of points
    template <typename Iterator>
    ConvexPolygon( Iterator begin, Iterator end );

    // Member functions
    
    /// construct the normal at the point of interest inside of the polygon (like at a hub vertex thinking of the polygon as a wheel)
    Point<dim> UnitNormal( const Point<dim>& inside_point ) const; 
    bool IsPointInside(const Point<dim>& point) const;
    void Print() const;
    std::vector<Point<dim>> GetCounterClockwiseVertices() const;

    // Const iterators for counterclockwise traversal
    typename std::vector<Point<dim>>::const_iterator CounterClockwiseVerticesBegin() const;
    typename std::vector<Point<dim>>::const_iterator CounterClockwiseVerticesEnd() const;

    // Method to output points to a .vtk file for ParaView visualization
    void Out(const std::string& filename) const;

   private:
     /// computes the unit normal for the first triangle in the polygon as is needed to test for point containment
     Point<dim> UnitNormalFirstTriangle() const;

     std::vector<Point<dim>> vertices;
     
     friend class ConvexPolygon_Test;
};


// Constructor accepting iterators
template <uint32_t dim>
template <typename Iterator>
inline ConvexPolygon<dim>::ConvexPolygon(Iterator begin, Iterator end) {
    std::vector<Point<dim>> points(begin, end);
    if (points.size() < 3) {
        throw std::invalid_argument("A convex polygon requires at least 3 points.");
    }
    vertices = points;
}



/// demonstrate how to use this class
int demo_ConvexPolygon3D();


} // end csmp

#endif /* CONVEX_POLYGON_3D_H */
