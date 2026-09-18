#ifndef CSMP_VDATA_TEST_H
#define CSMP_VDATA_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class VSet;

// MISSING
// TODO: test read / write of VSet/VData with faces and interfaces
// TODO: test extraction of manifolds

/**
      Tests file I/O and VSet repair / rectification methods, especially the creation of neighbor connectivity that is broken in ANSYS.

    @note Test uses utilities/vsetMakers.h to check the recreation of inter-element connectivity
        
    @author PL Nov 2010 & SKM 23/12022
*/
class VData_Test : public Test {
  public:
    virtual void run();
    
    // test read-write complete VSet, including boundaries and interfaces
    void TestBinaryIO();
    
    // checks whether mesh is still intact after corner elements were split
    // test model FracBox
    bool TestReplacementOfCornerTetrahedra();

    /// for vsetMaker vset 'MeshPatchWithLineElements' tests whether the handcoded neighbor connectivity is reproduced
    void Test_EstablishElementConnectivity2D();

    /**
       checks whether a neighbor-based traversal of line elements is possible after calling this method
       inside of EstablishElementConnectivity2D()
    */
    void Test_CreateConsistentLineElementOrientations2D();
    
    void Test_RecreateConnectivityOfTriangularMesh();

    void Test_RecreateConnectivityOfQuadrilateralMesh();
    
    /// reestablishes neighbor connectivity for 3D model Tetra (6 tets created from one hex)
    void Test_RecreateConnectivityOfTetrahedralMesh();

    /// reestablishes neighbor connectivity for 3D model Rubik cube
    void Test_RecreateConnectivityOfHexahedralMesh();
    
    void Test_RecreateConnectivityOfPrismMesh();

    void Test_InitialiseNodeTopologyIdentifiers();
    
    // mesh refinement and transformation
    void Test_refineSimplexMesh();
    
    void Test_linearToQuadraticMeshConversion();
    
  private:
    const static bool verbose_ = true;  
};


    


} // csmp

#endif // VDATA_TEST_H
