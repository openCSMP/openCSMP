// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  MeshManagementUtilities_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 16/1/2026.
//

#ifndef CSMP_MESH_MANAGEMENT_UTILITIES_TEST_H
#define CSMP_MESH_MANAGEMENT_UTILITIES_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp {

/**
   Tests of the free functions declared in meshManagementUtilities.h.

   Test_transformMeshIntoRefinedLinearAndQuadraticMeshes() exercises the functions on the large 'FracBox' model;
   the remaining methods use hand-checkable inputs with exactly known results:

   - Test_StandAloneElements():      single elements built in place (angles between cells, degenerate, inverted,
                                     skewed, interpenetrating and collocated cells)
   - Test_HexahedraMesh():           3x3x3 unit hexahedra, 64 nodes (vsetMakers: create_Hexahedra_VSet)
   - Test_TetrahedraMesh():          cube [-1,1]^3 cut into 6 tetrahedra (create_Tetra_VSet)
   - Test_DisconnectedMesh2D():      2 quadrilaterals and 4 triangles forming two patches (create_Disconnected2D_VSet)
   - Test_MeshWithLineElements2D():  hybrid 2D mesh with line elements (create_MeshPatchWithLineElements_VSet)
   - Test_FracturedBox():            the unrefined 'FracBox' with volume, surface and line elements in one region

   Not covered here: findSplitInterfaceElements (ANSYS_SplitBoundaryMatch_Test), checkNeighborNormalsForConsistentOrientation
   (BoundaryInterface_Test), boundaryAndTopoTypeFlagsToVTU (SplitBoundary_Test), findNodesViaHigherDimensionalNeighbors
   (marked 'NOT USED YET'), and the functions listed under 'open questions' at the top of the source file.

    @author SKM
    @date 16/1/2026
*/
class MeshManagementUtilities_Test : public Test {
  public:
    MeshManagementUtilities_Test() {}
    
    virtual void run();
  
  private:
  
    void Test_transformMeshIntoRefinedLinearAndQuadraticMeshes();
    void Test_StandAloneElements();
    void Test_HexahedraMesh();
    void Test_TetrahedraMesh();
    void Test_DisconnectedMesh2D();
    void Test_MeshWithLineElements2D();
    void Test_FracturedBox();
  
    static const bool verbose_ = true;
  
};

} // end csmp

#endif
