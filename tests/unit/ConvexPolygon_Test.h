//
//  ConvexPolygon_Test.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 13/10/2024.
//

#ifndef CSMP_CONVEX_POLYGON_TEST_H
#define CSMP_CONVEX_POLYGON_TEST_H

#include "Test.h"
#include "ConvexPolygon.h"

namespace csmp {

class ConvexPolygon_Test : public Test {
public:
    virtual void run() override;

private:
    void TestBasicConstruction();
    void TestUnitNormalOnVertices();
    void TestUnitNormalOnInterior();
    void TestIsPointInside();
    void TestOutFunction();
    void TestNormalComputationIn3D();
};

} // csmp

#endif /* CSMP_CONVEX_POLYGON_TEST_H */
