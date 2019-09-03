#ifndef PLACEMENT_TEST_2D_H
#define PLACEMENT_TEST_2D_H

#include "ANSYS_Model2D.h"

#include "Test.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

namespace csmp {

  class  Placement_Test_2D : public Test
  {

  public:
    Placement_Test_2D();
    ~Placement_Test_2D();
    virtual void run();

  private:
    double64 floatTolerance;
    ANSYS_Model2D* _model = nullptr;
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
    
    void ElementIPToElementIPTest();
    void FacetIPToFacetIPTest();
    void SectorIPToSectorIPTest();
    void NodeToNodeTest();
    void GradientTest();
    void SectorVolumeTest();
    void FacetNormalTest();
    
    //FiniteVolumePlacement tests
    void FV_ElementToElementTest();
    void FV_FacetAreaAndNormalTest(); 
    void FV_ProjectOntoFacetNormalTest();
    void FV_SectorVolumeTest();
    void FV_InsideOutsideNodeTest();
    void FV_GradientTest();  
  };

} // csmp

#endif // PLACEMENT_TEST_H
