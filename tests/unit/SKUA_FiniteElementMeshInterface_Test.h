//
//  SKUA_FiniteElementMeshInterface_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 2/12/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef SKUA_FINITE_ELEMENT_MESH_INTERFACE_TEST_H
#define SKUA_FINITE_ELEMENT_MESH_INTERFACE_TEST_H

#include "Test.h"
#include "SKUA_FiniteElementMeshInterface.h"
#include "Model.h"

namespace csmp {

class VData;

/** @brief Tests SKUA models from basic to more advanced.

    /// consecutive numbering of nodes and elements
    
    /// mesh consistency (Jacobians, element volume range, non-manifold vertices, triangle boxes
    
    /// pfverts, disambiguated neighbors of lower-dimensional elements
    /// consistent facing directions of surface elements
    /// consecutive numbering of series of line elements
    /// material IDs
    /// test assignment of box boundary flags (sides, edges, corners)
    /// property assignment
    /// use of regions file (is this still needed?)
    // ADVANCED
    /// generation of boundaries
    /// generation of split boundaries
    

*/
class SKUA_FiniteElementMeshInterface_Test : public Test {
  public:
    virtual void run();
    
  private:
    const static bool verbose_ = true;

    /// test that poperty values are assigned in the correct order using variables 'node number' and 'number'
    bool TestNodePropertyAssignment( const Model<3U>& );
    
    bool Test_Boundary_BoxBoundary_FlagConsistency( const Model<3U>& model,
                                                    const std::string& boundary_name );
    
    bool Test_BoxBoundaryFlags();

    bool TestElementNeighborConnectivity( Model<3U>& model );


  // SPLIT-BOUNDARY SPECIFIC TESTS
  
    /// tests whether the number of nodes in InterFace cells is consistent with their parent tetrahedral Element cells
    void Test_NodeNumberingConsistencyBetweenTetrahedraAndInterFaces( const VData& );
  
    /// Manifolds: 1) are they present at all split boundaries, 2) do they have the right number of entries
    void Test_NodesAndTheirManifolds( const Model<3U>& );
    
    /// Tests nodes on opposite sides of InterFaces: 1) are they collocated, 2) do they have the right mapping
    void Test_InterfaceNodes( const Model<3U>& );
    
 // HELPER METHODS
 
    /// for the elements in original order, prints original IDs of neighbor elements or boundary flag values where there is no neighbor
    void PrintOriginalNeighborIDs( const Model<3U>& ) const;
    
    void OutputRegionsToVTK( const Model<3U>&, const char* var_name ) const;
};

} // end csmp

#endif /* SKUA_FINITE_ELEMENT_MESH_INTERFACE_TEST_H */
