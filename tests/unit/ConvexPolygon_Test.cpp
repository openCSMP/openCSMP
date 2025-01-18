//
//  ConvexPolygon_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 13/10/2024.
//

#include "ConvexPolygon_Test.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

// Define the member functions outside the class

void ConvexPolygon_Test::run() {
    TestBasicConstruction();
    TestUnitNormalOnVertices();
    TestUnitNormalOnInterior();
    TestIsPointInside();
    TestNormalComputationIn3D();
    TestOutFunction();
}

void ConvexPolygon_Test::TestBasicConstruction() {
    _info("Testing basic construction of ConvexPolygon...");

    // Define points for a simple triangle
    std::vector<Point<3U>> points = {
        Point<3U>(0.0, 0.0, 0.0),
        Point<3U>(1.0, 0.0, 0.0),
        Point<3U>(0.0, 1.0, 0.0),
    };

    ConvexPolygon<3U> polygon(points);
    _test(polygon.vertices.size() == 3 ); //, "Triangle should have 3 vertices.");
}

void ConvexPolygon_Test::TestUnitNormalOnVertices() {
    _info("Testing UnitNormal on vertices...");

    std::vector<Point<3U>> points = {
        Point<3U>(0.0, 0.0, 0.0),
        Point<3U>(1.0, 0.0, 0.0),
        Point<3U>(0.0, 1.0, 0.0),
    };

    ConvexPolygon<3U> polygon(points);
    Point<3U> vertexNormal = polygon.UnitNormal(Point<3U>(0.5, 0.5, 0.0)); // Calculate the average normal

    _info("Calculated normal: " + std::to_string(vertexNormal[0]) + ", " +
         std::to_string(vertexNormal[1]) + ", " +
         std::to_string(vertexNormal[2]));

    _equal( vertexNormal.Length(), 1., 1.0e-14 ); //, "Normal at a vertex should not be zero.");
}

void ConvexPolygon_Test::TestUnitNormalOnInterior() {
    _info("Testing UnitNormal on an interior point...");

    std::vector<Point<3U>> points = {
        Point<3U>(0.0, 0.0, 0.0),
        Point<3U>(1.0, 0.0, 0.0),
        Point<3U>(0.0, 1.0, 0.0),
    };

    ConvexPolygon<3U> polygon(points);
    Point<3U> interiorPoint(0.1, 0.1, 0.0);
    Point<3U> normal = polygon.UnitNormal(interiorPoint);

    _info("Calculated average normal at interior point: " + std::to_string(normal[0]) + ", " +
         std::to_string(normal[1]) + ", " +
         std::to_string(normal[2]));

    _equal( normal.Length(), 1., 1.0e-14 ); // "Average normal at an interior point should not be zero.");
}


void ConvexPolygon_Test::TestIsPointInside() {
    _info("Testing IsPointInside...");

    std::vector<Point<3U>> points = {
        Point<3U>(0.0, 0.0, 0.0),
        Point<3U>(1.0, 0.0, 0.0),
        Point<3U>(0.0, 1.0, 0.0),
    };

    ConvexPolygon<3U> polygon(points);

    _test( polygon.IsPointInside(Point<3U>(0.1, 0.1, 0.0) ) ); // "Point (0.1, 0.1, 0) should be inside.");
    _test(!polygon.IsPointInside(Point<3U>(1.0, 1.0, 0.0) ) ); // "Point (1.0, 1.0, 0) should be outside.");
}



/**
     point data taken from a 3D-space somewhat distorted mesh
     @todo fails since point locations are not entered correctly
 */
void ConvexPolygon_Test::TestNormalComputationIn3D()
 {
    // field data test
    /*
    std::vector<Point<3U>> points = {
        Point<3U>(657284.,-1433.41,-5.73394e6),
        Point<3U>(657284.,-1431.8,-5.73395e6),
        Point<3U>(657281.,-1432.06,-5.73395e6),
        Point<3U>(657277.,-1432.31,-5.73395e6),
        Point<3U>(657277.,-1433.92,-5.73395e6),
        Point<3U>(657277.,-1435.52,-5.73395e6),
        Point<3U>(657281.,-1435.27,-5.73394e6),
        Point<3U>(657284.,-1435.01,-5.73394e6),
    };
    Point<3U> vertex( 657281., -1433.65, -5.73395e6);
    Point<3U> actual_normal( 0.279869, 0.483203, 0.829571 );
    */

    std::vector<Point<3U>> points = {
        Point<3U>( 1., 1., 4. ),
        Point<3U>( 3., 1., 4. ),
        Point<3U>( 3., 4., 4. ),
        Point<3U>( 1., 3.9, 4. ),
    };
    Point<3U> vertex( 2., 2.5, 4. );

    ConvexPolygon<3U> polygon(points);
    
    // node point shared by 4 quadrilaterals
    auto normal = polygon.UnitNormal( vertex );
    
    const double tolerance = numeric_limits<double>::epsilon() * 100.;
    _test( approximatelyEqual(normal.Length(),1.,tolerance) );
    
    Point<3U> actual_normal( 0, 0., 1. );
    
    _test( approximatelyEqual( normal[0], actual_normal[0] ) );
    _test( approximatelyEqual( normal[1], actual_normal[1] ) );
    _test( approximatelyEqual( normal[2], actual_normal[2] ) );
    
     polygon.Out("ConvexPolygon_perimeter.vtk");
    
 } // end



void ConvexPolygon_Test::TestOutFunction() {
    _info("Testing Out function...");

    // This test would require file handling to verify the output.
    _info("Out function test is not implemented.");
}

} // end csmp
