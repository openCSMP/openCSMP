#ifndef CSP_FINITE_VOLUME_STENCIL_TEST_H
#define CSP_FINITE_VOLUME_STENCIL_TEST_H

#include "CSMP_definitions.h"
#include "FiniteVolumeStencil.h"
#include "FiniteElement.h"
#include "Test.h"

namespace csmp {

class FiniteElement;

class FiniteVolumeStencil_Test : public Test {
  public:
    explicit FiniteVolumeStencil_Test( bool verbose=false );
    ~FiniteVolumeStencil_Test();
    void run(); // runs all the tests for the class (register other methods)
    void displayReferenceCoordinates();
    void facetAndSectorNumbersTest(); // are they right for all element types
    void sectorFacetConnectivityTest(); // facets that delimit a certain sector for all the element types
    void orientationAndLengthOfNormalsTest();
    void normalTransformationTest();
    void locationOfSectorIntegrationPointsTest(); // argument list could deal with multiple integration points
    void weightsAndFacetIntegrationPointsTest(); //test location and weights of facet inetgration points
    void weightsOfSectorIntegrationPointsTest(); // argument list could deal with multiple integration points
    void shapeFunctionsTest();
    void shapeFunctionDerivativesTest();
    // projection tests are done in the EFVT unit test
  
  private:
    std::vector<FiniteVolumeStencil<3U> >  fvs_;
    std::vector<FiniteElement*>            vecFEs_;
    
    template<uint32_t dim>
	  std::vector< Point<dim> > GetPointsOfFacet( const uint32_t& iFacet, const CSMP_FEM_TYPE& elType, 
                                                const FiniteVolumeStencil<dim>& fvs, const FiniteElement & fe );
    const bool verbose_;
};



} //end csmp

#endif
