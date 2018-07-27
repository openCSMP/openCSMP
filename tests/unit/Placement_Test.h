#ifndef PLACEMENT_TEST_H
#define PLACEMENT_TEST_H

#include "ANSYS_Model2D.h"

#include "Test.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

namespace csmp {

  class  Placement_Test : public Test
  {

  public:
    Placement_Test();
    ~Placement_Test();
    virtual void run();

  private:
    double64 floatTolerance;
    ANSYS_Model2D* mockModel = nullptr;
    double64 nodeScalar, elementScalar, elementIPScalar, sectorIPScalar, facetIPScalar;
    double64 vx, vy, xx, xy, yx, yy;
    double64 tx, ty, tz;
    double64 arx, ary, arz;

    // target element
    void NodeToElementTest();
    void ElementToElementTest();
    void ElementIPToElementTest();
    void SectorIPToElementTest();
    void FacetIPToElementTest();

    // target elementIP
    void NodeToElementIPTest();
    void ElementToElementIPTest();
    void SectorIPToElementIPTest();
    void FacetIPToElementIPTest();

    // target facetIP
    void NodeToFacetIPTest();
    void ElementToFacetIPTest();
    void SectorIPToFacetIPTest();

    // target sectorIP
    void NodeToSectorIPTest();
    void ElementToSectorIPTest();

    // target node
    void ElementIPToNodeTest(); // extrapolation
    void SectorIPToNodeTest();  // extrapolation
    void FacetIPToNodeTest();   // extrapolation

    // faceFlux
    // normalFacet

      
  };

} // csmp

#endif // PLACEMENT_TEST_H
