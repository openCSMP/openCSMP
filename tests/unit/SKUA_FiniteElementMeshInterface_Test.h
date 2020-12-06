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

namespace csmp {

class SKUA_FiniteElementMeshInterface_Test : public Test {
  public:
    virtual void run();
    
  private:
  /// consecutive numbering of nodes and elements
  /// mesh consistency (Jacobians, element volume range, non-manifold vertices, triangle boxes
  
  /// pfverts, disambiguated neighbors of lower-dimensional elements
  bool TestNeighborConnectivity( Model<3U>& model, bool verbose=true );
  /// consistent facing directions of surface elements
  /// consecutive numbering of series of line elements
  /// material IDs
  /// test assignment of box boundary flags (sides, edges, corners)
  /// property assignment
  /// use of regions file (is this still needed?)
  // ADVANCED
  /// generation of boundaries  
  /// generation of split boundaries
  
  /// for the elements in original order, prints original IDs of neighbor elements or boundary flag values where there is no neighbor
  void PrintOriginalNeighborIDs( const Model<3U>& ) const;
};

} // end csmp

#endif /* SKUA_FINITE_ELEMENT_MESH_INTERFACE_TEST_H */
