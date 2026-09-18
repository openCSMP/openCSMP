// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  ConvexPolygon.cpp
//  Open ACGSS
//
//  Created by Stephan Matthai on 13/10/2024.
//

#include "ConvexPolygon.h"

using namespace std;

namespace csmp {

// ConvexPolygon member function definitions

template <uint32_t dim>
ConvexPolygon<dim>::ConvexPolygon(const std::vector<Point<dim>>& points) {
    if (points.size() < 3) {
        throw std::invalid_argument("A convex polygon requires at least 3 points.");
    }
    vertices = points;
}


template <uint32_t dim>
Point<dim> ConvexPolygon<dim>::UnitNormal( const Point<dim>& point ) const {
//    static_assert(dim == 3, "UnitNormal is only defined for 3D polygons.");

    if (!IsPointInside(point)) {
        throw std::invalid_argument("ConvexPolygon<dim>::UnitNormal: Point must be inside the polygon to calculate normal.");
    }

    Point<dim> averageNormal(0.0); // Initialize the average normal
    size_t numEdges = vertices.size();

    // Loop through each edge to calculate the normal
    for ( uint32_t i = 0u; i < numEdges; ++i) {
        const Point<dim>& v1 = vertices[i];
        const Point<dim>& v2 = vertices[(i + 1) % numEdges];

        // Calculate the edges
        Point<dim> edge1 = v2 - v1; // Edge from v1 to v2
        Point<dim> edge2 = point - v1; // Edge from v1 to the interior point

        // Calculate the normal using the cross product
        Point<dim> normal = crossProduct(edge1, edge2);
        averageNormal += normal; // Accumulate the normals
    }
    averageNormal.NormalizeLengthTo(1.);

    // Normalize the average normal
    return averageNormal;
}




template <uint32_t dim>
Point<dim> ConvexPolygon<dim>::UnitNormalFirstTriangle() const {
//    static_assert(dim == 3, "UnitNormal is only defined for 3D polygons.");
    if (vertices.size() < 3) {
        throw std::runtime_error("ConvexPolygon<dim>::UnitNormalFirstTriangle: A polygon must have at least 3 vertices to define a normal.");
    }

    // Choose the first three vertices to compute the normal
    Point<dim> edge1 = vertices[1] - vertices[0]; // First edge
    Point<dim> edge2 = vertices[2] - vertices[0]; // Second edge

    // Compute the cross product to get the normal vector
    Point<dim> normal = crossProduct(edge1, edge2);

    // Normalize the normal vector
    double length = normal.Length();
    if (length > 0.0) {
        return normal / length; // Return the unit normal
    } else {
        throw std::runtime_error("ConvexPolygon<dim>::UnitNormalFirstTriangle: Computed normal has zero length.");
    }
}



template <uint32_t dim>
bool ConvexPolygon<dim>::IsPointInside( const Point<dim>& point ) const {
    // First, calculate the normal of the polygon using any point in the plane
    Point<dim> normal = UnitNormalFirstTriangle(); // Use the first vertex as a reference

    // Loop through each edge of the polygon
    for (size_t i = 0; i < vertices.size(); ++i) {
        // Get the current vertex and the next vertex (wrap around)
        Point<dim> v1 = vertices[i];
        Point<dim> v2 = vertices[(i + 1) % vertices.size()];

        // Compute the vector from v1 to v2 and from v1 to the point
        Point<dim> edge = v2 - v1;
        Point<dim> toPoint = point - v1;

        // Calculate the cross product to determine the position of the point relative to the edge
        Point<dim> cross_product = crossProduct( edge, toPoint );

        // Check the sign of the cross product with respect to the normal
        // If the dot product is negative, the point is outside the polygon
        if ( dotProduct( normal, cross_product ) < 0. ) {
            return false; // The point is outside the polygon
        }
    }
    return true; // The point is inside the polygon
}



template <uint32_t dim>
void ConvexPolygon<dim>::Print() const {
    for (const auto& v : vertices) {
        std::cout << v << std::endl;
    }
}



template <uint32_t dim>
std::vector<Point<dim>> ConvexPolygon<dim>::GetCounterClockwiseVertices() const {
    // Simplified area calculation on the XY plane (for 2D or projection)
    auto signedArea = [](const Point<dim>& p1, const Point<dim>& p2) {
        return (p1[0] * p2[1] - p2[0] * p1[1]);  // Signed area for 2D points
    };

    double totalArea = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const Point<dim>& current = vertices[i];
        const Point<dim>& next = vertices[(i + 1) % vertices.size()];
        totalArea += signedArea(current, next);
    }

    if (totalArea < 0) {
        std::vector<Point<dim>> ccwVertices = vertices;
        std::reverse(ccwVertices.begin(), ccwVertices.end());
        return ccwVertices;
    }

    return vertices;
}

// Const iterators for traversing the vertices in counterclockwise order
template <uint32_t dim>
typename std::vector<Point<dim>>::const_iterator ConvexPolygon<dim>::CounterClockwiseVerticesBegin() const {
    static std::vector<Point<dim>> ccwVertices = GetCounterClockwiseVertices();
    return ccwVertices.cbegin();
}

template <uint32_t dim>
typename std::vector<Point<dim>>::const_iterator ConvexPolygon<dim>::CounterClockwiseVerticesEnd() const {
    static std::vector<Point<dim>> ccwVertices = GetCounterClockwiseVertices();
    return ccwVertices.cend();
}

// Function to output the points to a .vtk file
template <uint32_t dim>
void ConvexPolygon<dim>::Out(const std::string& filename) const {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Could not open the file for writing.");
    }

    // VTK file header
    outFile << "# vtk DataFile Version 2.0\n";
    outFile << "Convex Polygon\n";
    outFile << "ASCII\n";
    outFile << "DATASET POLYDATA\n";
    
    // Writing vertices
    outFile << "POINTS " << vertices.size() << " float\n";
    for (const auto& vertex : vertices) {
        outFile << vertex[0] << " " << vertex[1] << " ";
        if constexpr (dim == 3) {
            outFile << vertex[2]; // Only add Z coordinate if dim == 3
        } else {
            outFile << "0.0"; // Add 0.0 for Z coordinate if 2D
        }
        outFile << "\n";
    }

    // Writing polygon connectivity
    outFile << "POLYGONS 1 " << vertices.size() + 1 << "\n";
    outFile << vertices.size();  // Number of vertices in the polygon
    for (uint32_t i = 0; i < vertices.size(); ++i) {
        outFile << " " << i;
    }
    outFile << "\n";

    outFile.close();
}

template class ConvexPolygon<1U>;
template class ConvexPolygon<2U>;
template class ConvexPolygon<3U>;



// Main function to demonstrate usage
int demo_ConvexPolygon3D() {

    std::vector<Point<3U>> points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    ConvexPolygon polygon(points);

    std::cout << "Original vertices:" << std::endl;
    polygon.Print();

    // Accessing counterclockwise vertices via iterators
    std::cout << "Counterclockwise vertices:" << std::endl;
    for (auto it = polygon.CounterClockwiseVerticesBegin(); it != polygon.CounterClockwiseVerticesEnd(); ++it) {
        std::cout << *it << std::endl;
    }

    // Output the points to a .vtk file for ParaView
    polygon.Out("convex_polygon.vtk");

    Point<3U> testPoint = {0.5, 0.5, 0};
    std::cout << "Is point inside? " << (polygon.IsPointInside(testPoint) ? "Yes" : "No") << std::endl;

    return 0;
}

} // end csmp



