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
        const Region<3>& mdomain = lin_model.Region("Model");

        _test( distance(lin_model.BoundariesBegin(), lin_model.BoundariesEnd()) == 6 );

        // Testing using MeshManagementUtilities computation
        size_t volume_cells{0}, surface_cells{0}, line_cells{0};
        size_t n_total_cells = currentCellTypes( lin_model.Mesh(), ELEMENT, volume_cells, surface_cells, line_cells );
        size_t n_model_elmts = lin_model.Mesh().Elements();
        _test( n_total_cells == n_model_elmts + lin_model.Mesh().Faces() );
        _test( volume_cells == 10904 ); // after conversion into Face objects
        _test( volume_cells < refined_lin_mesh.Elements() ); // after conversion into Face objects
        _test( line_cells == 0 ); // all removed since not mentioned in regions file
        
        _test( isoparametricElementMesh(lin_model) == true );

        set<size_t> detached_elmts;
        _test( detectElementsWithAllNodesOnBoundary( lin_model.Mesh(), detached_elmts ) == 0 );

        const auto duplicates = detectDuplicateCells<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd(), verbose_ );
        _test( duplicates == 0 );
    
        // Compute parent element barycentre-to-node distances for range of nodes;
        // returns them into vector [e1,e2...e_n,e_sum] with a length of parent elements+1
        vector<vector<double> > distances_and_weight;
        distancesAndWeights<3>( mdomain.NodesBegin(), mdomain.NodesEnd(), distances_and_weight );
        _test( distances_and_weight[25][2] > 0. );

        // trying corners:
        //for supplied edge nodes, find their volumetric parent elements; if find segment ids is on, their local numbers are assigned to Idx of the parent elements
        for ( auto& it : mdomain.CellVector() )
           if ( isCorner( it->AtBoundary(1) /*face=1*/ ) ) {
                 const bool find_segm_ids{ true };
                 vector<Node<3>*>                   edge_nodes;
                 map<Element<3>*,vector<Node<3>*> > segm_parents;
                 size_t n_ele = parentElementsSharingMultipleEdgeNodes( edge_nodes, segm_parents, find_segm_ids );
                 _test( n_ele >= 2 );
                 _test( edge_nodes.size() == 2 );
                 _test( segm_parents.size() == 2 );
             }

        // finds node by point coordinate; returns -1 if not found; @attention tolerance needs to account for single-precision of CAD tools
        auto test_node = mdomain.N(25);
        long idx = findNode( lin_model, test_node->Coordinate(), 1.0e-7, verbose_ );
        _test( idx >= 0 ); // node was found
        _test( test_node->Idx() == static_cast<size_t>(idx) );
        size_t id = findNode( lin_model, test_node->x(), test_node->y(), test_node->z(), 1.0e-7 );
        _test( test_node->Idx() == id );

         set<Node<3>*> nodes;
        _test( detectOrphanNodes( lin_model.Mesh(), nodes ) == 0 );
        _test( findInterconnectedNodeCluster( test_node, nodes ) == lin_model.Mesh().Nodes() ); // all nodes should be connected

        // Relying on the node-to-node connectivity, uses a recursive depth first traversal to discover
        // the interconnected nodes in the supplied range; set of node pairs defining the cel segments must be input
        set<pair<const Node<3>*,const Node<3>*>> validEdges;
        size_t n_inds = findInterconnectedNodes<3>( lin_model.Mesh().Nodes(), validEdges );
        _test( n_inds == lin_model.Mesh().Nodes() );

        double min_dn = minimumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        double max_dn = maximumNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        double avg_dn = averageNodeSpacing<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        _test( min_dn > 0.4 );
        _test( min_dn < max_dn );
        _test( max_dn < 50. );
        _test( avg_dn > min_dn );

        // must be a volumetric element 
        Point<3> bctr = (mdomain.E(1000)->IsEquidimensional()) ? mdomain.E(1000)->BaryCenter() : mdomain.E(2000)->BaryCenter();
        _test( isContainedIn<3>( mdomain, array<double,3>{ bctr[0], bctr[1], bctr[2] } ) == mdomain.E(1000) );

        size_t duplicated_elmts = duplicatesCheck<3,Element>( mdomain.CellsBegin(), mdomain.CellsEnd() );
        _test( duplicated_elmts == 0 );
                          
        size_t connectivity_errors = connectivityCheck<3>( mdomain.CellsBegin(), mdomain.CellsEnd() );
        _test( connectivity_errors == 0 );

        pair<Point<3>,Point<3>> cnrs = boundingBox<3>( mdomain.NodesBegin(), mdomain.NodesEnd() );
        _test( cnrs.first < cnrs.second );
      }
  }
  
} // end csmp
