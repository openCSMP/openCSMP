#ifndef FINITE_VOLUME_POLICY_TEST_H
#define FINITE_VOLUME_POLICY_TEST_H

#include "Test.h"
#include "FiniteVolumePolicy.h"
#include "Element.h"

namespace csmp {

/** extension of Element interface for FV computations
    for all elements we need to test:

    - *build the element by applying the following transformation to the orginal reference element:
      - Translate(-3,0,0)
      - Rotate(N0, 40degrees)
      - Shear(N0, N1, 35degrees)
    - *facet and sector weights are as calculated in Rhino
      - *normals are correct
      - *FacetArea() and SectorVolume() give correct results
      - *Sector volumes add to volume of element
      - *projection of vector variable in physical space gives the same result 
        *as the more efficient projection (ParametricProjection)
      - (checked)interpolation functions at sector and facet integration points are correct 
      - (checked)PropertyValues are interpolated to the correct integration points
      - (used in reimpl of FacetArea heavily)rst to xyz conversion is correct 
*/
class FiniteVolumePolicy_Test : public Test {
  public:
    FiniteVolumePolicy_Test();
    ~FiniteVolumePolicy_Test();
    void run(); // runs all the tests for the class (register other methods)
    void IsoparametricLinearLineElement_Test(double64 fTolerance, double64 fToleranceInternal);
    template<size_t dim> void IsoparametricLinearTriangle_Test(double64 fTolerance, double64 fToleranceInternal);
    template<size_t dim> void IsoparametricLinearQuadrilateral_Test(double64 fTolerance, double64 fToleranceInternal);
    void IsoparametricLinearTetrahedron_Test(double64 fTolerance, double64 fToleranceInternal);
    void IsoparametricLinearPyramid_Test(double64 fTolerance, double64 fToleranceInternal);
    void IsoparametricLinearPrism_Test(double64 fTolerance, double64 fToleranceInternal);

    //Revised by J.E.M. 08-09-2010
    void Test_UnitaryIsoparametricLinearHexahedron(double64 fTolerance, double64 fToleranceInternal);
    void Test_IsoparametricLinearHexahedron1(double64 fTolerance, double64 fToleranceInternal);
    void Test_IsoparametricLinearHexahedron2(double64 fTolerance, double64 fToleranceInternal);
    void Test_CreateVSet();
	
private:
	bool m_bTestFacetAreas;
	bool m_bTestParametricFacetArea;
	bool m_bTestFacetNormals;
	bool m_bTestSectorVolumes;
	bool m_bTestParametricFacetNormals;
	bool m_bProjectionOnFacetNormal;
	bool m_bPropertyValueAt;
	bool m_bRSTToXYZ;
};
 
}

#endif /* FINITE_VOLUME_POLICY_TEST_H */
