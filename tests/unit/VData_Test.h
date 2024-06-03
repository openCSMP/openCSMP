

#ifndef VDATA_TEST_H
#define VDATA_TEST_H

#include "Test.h"

namespace csmp {

template<uint32_t> class VSet;

// MISSING
// TODO: test all the functionality that recreates neighbors etc.
// TODO: test read / write of VSet/VData with faces and interfaces
// TODO: test extraction of manifolds

/// PL Nov 2010 & SKM 23/12022
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

    /* checks whether a neighbor-based traversal of line elements is possible after calling this method
       inside of EstablishElementConnectivity2D()
    */
    void Test_CreateConsistentLineElementOrientations2D();
    
    /// reestablishes neighbor connectivity for 3D model Tetra (6 tets created from one hex)
    void Test_RecreateConnectivityOfTetrahedralMesh();

    /// reestablishes neighbor connectivity for 3D model Rubik cube
    void Test_RecreateConnectivityOfHexahedralMesh();

    void Test_InitialiseNodeTopologyIdentifiers();
    
  private:
    const static bool verbose_ = true;  
};


    


} // csmp

#endif // VDATA_TEST_H
