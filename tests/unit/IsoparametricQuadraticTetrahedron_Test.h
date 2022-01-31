//
//  IsoparametricQuadraticTetrahedron_Test.h
//  
//
//  Created by Stephan Matthai on 7/14/14.
//
//

#ifndef IsoparametricQuadraticTetrahedron_Test_h
#define IsoparametricQuadraticTetrahedron_Test_h

#include "Test.h"
#include "FiniteElement.h"
#include "Element.h"

namespace csmp {

/**
    Testing:
    
    1. Interpolation functions at all nodes = 1; 0 at all other nodes
    2. Placement of integration points is correct
    3. Interpolation functions at up to 1 at all integration points
*/

class IsoparametricQuadraticTetrahedron_Test : public Test {
  public:
    explicit IsoparametricQuadraticTetrahedron_Test( bool verbose );
    ~IsoparametricQuadraticTetrahedron_Test();
  
    virtual void run();
  
    bool TestIntegrationPointPlacement();
    void TestInterpolationFunctionValues( const Element<3U>& );
    void TestSumOfInterpolationFunctionValuesEqualTo1( const Element<3U>& );
    void CheckElementFaceConsistency( const Element<3U>& );
    void CheckInterpolation( const Element<3U>& );
  
    void OutputIntegrationPointsToVTK( const char* file, const Element<3U>& ) const;
    void OutputFaceNormalsToVTK( const char* file, const Element<3U>& ) const;
  
  private:
    void ChangeNodeCoordinatesToTestConfiguration();
    void ChangeNodeCoordinatesToParametric();
  
  private:
    FiniteElement*  ltetra_, *qtetra_;
    Element<3U>     element_, element2_;
    Node<3U>        n0, n1, n2, n3, n4, n5, n6, n7, n8, n9;
    const double  tolerance_factor_;
    bool            verbose_;
};

} // end csmp


#endif /* defined(IsoparametricQuadraticTetrahedron_Test_h) */
