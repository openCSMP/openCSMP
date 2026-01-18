//
//  MeshManager_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "MeshManagementUtilities_Test.h"
#include "meshManagementUtilities.h"
#include "CSMP_definitions.h"
#include "ErrorHandler.h"
#include "vsetMakers.h"
#include "VTU_Interface.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
//#include "NodeManifold.h"
#include "compareFloats.h"

#include "VTK_Interface.h"

// #define CSMP_MESH_MANAGER_TEST_DEBUG

using namespace std;

namespace csmp {

void MeshManagementUtilities_Test::run()
 {
    Test_transformMeshIntoRefinedLinearAndQuadraticMeshes();
 }
 
 

void MeshManagementUtilities_Test::Test_transformMeshIntoRefinedLinearAndQuadraticMeshes()
  {
      VSet<3>       mesh;
      ModelTopology topo = create_FracBox( mesh );

      VSet<3>       refined_lin_mesh, quadratic_mesh;
      ModelTopology refined_linear_topo, quadratic_topo;
      
      // removing any variables with a placement other than ELEMENT or FACE since they cannot handled by the tested method yet
      set<string> props_to_delete;
      for ( auto pit=mesh.PropertyValuesBegin(); pit!=mesh.PropertyValuesEnd(); ++pit )
        if ( (*pit).second.Placement() != ELEMENT && (*pit).second.Placement() != FACE )
          props_to_delete.insert( (*pit).first );
      for ( const auto& pit : props_to_delete )
        mesh.RemoveData( pit.c_str() );
  
      transformMeshIntoRefinedLinearAndQuadraticMeshes( mesh, topo, refined_lin_mesh, refined_linear_topo,
                                                        quadratic_mesh, quadratic_topo );
                                                        
      // testing the transformed meshes
      {
        const bool treat_domains_as_regions_and_use_regions_file_if_any{true}; // TODO: change wording
        Model<3> lin_model( refined_linear_topo, refined_lin_mesh, "CSMP-1phase-variables.txt", treat_domains_as_regions_and_use_regions_file_if_any );

        // Testing using MeshManagementUtilities computation
        _test( distance(lin_model.BoundariesBegin(), lin_model.BoundariesEnd()) == 6 );
      }
  }
  
  
} // end csmp
